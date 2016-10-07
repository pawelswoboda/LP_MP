
#ifndef LP_MP_LP_INTERFACE_GUROBI_HXX
#define LP_MP_LP_INTERFACE_GUROBI_HXX

#include "lp_interface.h"

namespace LP_MP {
  
  class LpInterfaceGurobi : public LpInterfaceAdapter {
  public:

    LpInterfaceGurobi() : env_(GRBEnv()),model_(GRBModel(env_)){}
  
    LpInterfaceGurobi(INDEX noVars) : env_(GRBEnv()),model_(GRBModel(env_)),noVars_(noVars),epsilon_(std::numeric_limits<REAL>::infinity()) {
      std::vector<double> obj(noVars_,1);
      std::vector<char> types(noVars_,GRB_INTEGER);
      std::vector<REAL> lb(noVars_,-std::numeric_limits<REAL>::infinity());
      MainVars_ = model_.addVars(&lb[0],NULL,&obj[0],&types[0],NULL,noVars_);
    }
        
    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    LpInterfaceGurobi(FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd, MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd,bool MIP = true)
      : env_(GRBEnv()),model_(GRBModel(env_)),epsilon_(std::numeric_limits<REAL>::infinity()) {

      // Standard Parameter for Gurobi
      model_.getEnv().set(GRB_DoubleParam_TimeLimit,3600);
      model_.getEnv().set(GRB_IntParam_Threads,1);
      model_.getEnv().set(GRB_IntParam_MIPFocus, 1);
      
      noVars_ = 0;
      noAuxVars_ = 0;
      INDEX preAuxVars = 0;
      INDEX noFactors = 0;
      for(auto factorIt = factorBegin; factorIt != factorEnd; ++factorIt) {
        noVars_ += factorIt->size();
        noAuxVars_ += factorIt->GetNumberOfAuxVariables();
        factorIt->SetAuxOffset(preAuxVars);
        preAuxVars = noAuxVars_;
        noFactors++;
      }
      std::vector<double> obj_1(noVars_,1.0);
      std::vector<double> obj_2(noAuxVars_,0.0);
      std::vector<char> types_1;
      if(MIP){
        types_1.resize(noVars_,GRB_INTEGER);
      } else {
        types_1.resize(noVars_,GRB_CONTINUOUS);
      }
      std::vector<char> types_2(noAuxVars_,GRB_CONTINUOUS);
      MainVars_ = model_.addVars(NULL,NULL,&obj_1[0],&types_1[0],NULL,noVars_);
      if( noAuxVars_ > 0){
        std::vector<REAL> lb (noAuxVars_,-std::numeric_limits<REAL>::infinity());
        MainAuxVars_ = model_.addVars(&lb[0],NULL,&obj_2[0],&types_2[0],NULL,noAuxVars_);
      }
      model_.update();

      /* Add Factor Constraints */
      INDEX fixedVars = 0;
      for(auto factorIt = factorBegin; factorIt != factorEnd; ++factorIt) {
        Offset_ = factorIt->GetPrimalOffset();
        size_ = factorIt->size();
        OffsetAux_ = factorIt->GetAuxOffset();
        sizeAux_ = factorIt->GetNumberOfAuxVariables();
        auto pot = factorIt->GetReparametrizedPotential();
        for(INDEX i=0;i<size_;++i) {
          REAL value = pot[i];
          if( std::isfinite(value) ){
            GetVariable(i).set(GRB_DoubleAttr_Obj,value);
          } else {
            GetVariable(i).set(GRB_DoubleAttr_Obj,0);
            GetVariable(i).set(GRB_DoubleAttr_LB,0);
            GetVariable(i).set(GRB_DoubleAttr_UB,0);
            fixedVars++;
          }
        }
        factorIt->CreateConstraints(this);
      }
      //printf("Reparametrization fixed %d variables\n",fixedVars);

      /* Add Message Constraints */
      for(auto messageIt = messageBegin; messageIt != messageEnd; ++messageIt) {
        leftSize_ = messageIt->GetLeftFactor()->size();
        rightSize_ = messageIt->GetRightFactor()->size();
        OffsetLeft_ = messageIt->GetLeftFactor()->GetPrimalOffset();
        OffsetRight_ = messageIt->GetRightFactor()->GetPrimalOffset();
        messageIt->CreateConstraints(this);
      }

      model_.update();
    } 

    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    void ReduceLp(FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd, MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd,REAL epsilon){
      epsilon_ = epsilon;
      for(auto factorIt = factorBegin; factorIt != factorEnd; ++factorIt) {
        Offset_ = factorIt->GetPrimalOffset();
        size_ = factorIt->size();
        OffsetAux_ = factorIt->GetAuxOffset();
        sizeAux_ = factorIt->GetNumberOfAuxVariables();
        factorIt->ReduceLp(this);
      }
      model_.update();
    }
    
    LinExpr CreateLinExpr() { return LinExpr(); }
    
    INDEX GetFactorSize() const { return size_; }
    INDEX GetLeftFactorSize() const { return leftSize_; }
    INDEX GetRightFactorSize() const { return rightSize_; }
        
    LpVariable GetVariable(const INDEX i) const { assert(i < size_);  assert(Offset_ + i < noVars_); return MainVars_[Offset_ + i]; }
    LpVariable GetLeftVariable(const INDEX i) const { assert(i < leftSize_); assert(OffsetLeft_ + i < noVars_); return MainVars_[OffsetLeft_ + i]; }
    LpVariable GetRightVariable(const INDEX i) const { assert(i < rightSize_); assert(OffsetRight_ + i < noVars_); return MainVars_[OffsetRight_ + i]; }

    LpVariable GetAuxVariable(const INDEX i) const { assert(i < sizeAux_); return MainAuxVars_[OffsetAux_ + i]; }

    REAL GetEpsilon() const { return epsilon_; }
    
    REAL GetVariableValue(const INDEX i) const;
    REAL GetObjectiveValue() const;
    REAL GetBestBound() const;
    
    void SetVariableBound(LpVariable v,REAL lb,REAL ub,bool integer = false);
    void SetTimeLimit(REAL t){ model_.getEnv().set(GRB_DoubleParam_TimeLimit,t); }
    void SetNumberOfThreads(INDEX t){ model_.getEnv().set(GRB_IntParam_Threads,t); };
    void SetDisplayLevel(INDEX t){ assert(t <= 1); model_.getEnv().set(GRB_IntParam_OutputFlag,t); };
    
    void addLinearEquality(LinExpr lhs,LinExpr rhs);
    void addLinearInequality(LinExpr lhs,LinExpr rhs);

    template<class factor>
    void addFactor(const factor& f,INDEX offset);
    
    int solve();
    int solve(PrimalSolutionStorage::Element primal);

    void WriteLpModel(std::string name){ model_.write(name); }
    
  private:
    GRBEnv env_;
    GRBModel model_;
    LpVariable* MainVars_;
    LpVariable* MainAuxVars_;
    
    INDEX Offset_,OffsetAux_,OffsetLeft_,OffsetRight_;
    INDEX size_,sizeAux_,leftSize_,rightSize_;
    INDEX noVars_,noAuxVars_;    
    REAL epsilon_;
  };

  template<class factor>
  void LpInterfaceGurobi::addFactor(const factor& f,INDEX offset){
    Offset_ = offset;
    size_ = f.size();
    f.CreateConstraints(this);
  }
  
  int LpInterfaceGurobi::solve(){
    int status = 2;
    try{
      model_.update();

      model_.optimize();
      auto stat = model_.get(GRB_IntAttr_Status);
      //std::cout << "Gurobi Status: " << stat << std::endl;
      if(stat == GRB_OPTIMAL ){
        status = 0;
      }
      else if(model_.get(GRB_IntAttr_SolCount) > 0){
        status = 3;
      }
      else if(stat == GRB_INFEASIBLE){
        status = 1;
      }
      else if(stat == GRB_TIME_LIMIT){
        status = 4;
      }
      
    }
    catch(GRBException e) {
      std::cout << "Error code = " << e.getErrorCode() << std::endl;
      std::cout << e.getMessage() << std::endl;
    }
    return status;
  }

  int LpInterfaceGurobi::solve(PrimalSolutionStorage::Element primal){
     // go through primal and when it is not unknownState, set variables to zero in model
     for(INDEX i=0; i<noVars_; ++i) {
        if(primal[i] == true) {
           GetVariable(i).set(GRB_DoubleAttr_Obj,0);
           GetVariable(i).set(GRB_DoubleAttr_LB,1);
           GetVariable(i).set(GRB_DoubleAttr_UB,1);
        } else if(primal[i] == false) {
           GetVariable(i).set(GRB_DoubleAttr_Obj,0);
           GetVariable(i).set(GRB_DoubleAttr_LB,0);
           GetVariable(i).set(GRB_DoubleAttr_UB,0);
        } else {
           assert(primal[i] == unknownState);
        }
     }
     int status = 2;
     try{
        model_.update();

        model_.optimize();
        status = 0;
     }
     catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
     }
     // now return to original model without fixed variables
     return status;
  }

  REAL LpInterfaceGurobi::GetVariableValue(const INDEX i) const{
    assert(i < noVars_);
    return MainVars_[i].get(GRB_DoubleAttr_X);
  }

  REAL LpInterfaceGurobi::GetObjectiveValue() const{
    return model_.get(GRB_DoubleAttr_ObjVal);
  }
  
  void LpInterfaceGurobi::SetVariableBound(LpVariable v,REAL lb,REAL ub,bool integer){
    v.set(GRB_DoubleAttr_LB,lb);
    v.set(GRB_DoubleAttr_UB,ub);
    if(integer){
      v.set(GRB_CharAttr_VType,'I');
    } else {
      v.set(GRB_CharAttr_VType,'C');
    }
  }
  
  void LpInterfaceGurobi::addLinearEquality(LinExpr lhs,LinExpr rhs){
    model_.addConstr(lhs,GRB_EQUAL,rhs);
  }

    void LpInterfaceGurobi::addLinearInequality(LinExpr lhs,LinExpr rhs){
    model_.addConstr(lhs,GRB_LESS_EQUAL,rhs);
  }
  
  
  REAL LpInterfaceGurobi::GetBestBound() const{
    REAL bound = 0.0;
    if(model_.get(GRB_IntAttr_IsMIP)){
      bound = model_.get(GRB_DoubleAttr_ObjBound);
    }
    return bound;
  }
  
}

#endif // LP_MP_LP_INTERFACE_GUROBI_HXX
