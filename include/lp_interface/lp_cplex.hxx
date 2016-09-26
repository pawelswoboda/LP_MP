#ifndef LP_MP_LP_INTERFACE_CPLEX_HXX
#define LP_MP_LP_INTERFACE_CPLEX_HXX

#include "lp_interface.h"

namespace LP_MP {
  
  class LpInterfaceCplex : public LpInterfaceAdapter {
  public:

    LpInterfaceCplex(INDEX noVars) : env_(IloEnv()),model_(IloModel(env_)),noVars_(noVars) {
      MainVars_ = IloNumVarArray(env_,noVars_,0.0,1.0,LpVariable::Float);
      //cplex_ = IloCplex(model_);
    }
    
    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    LpInterfaceCplex(FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd, MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd,bool MIP = true)
      : env_(IloEnv()),model_(IloModel(env_)) {

      // Standard Parameter for Cplex
      //model_.getEnv().set(GRB_DoubleParam_TimeLimit,3600);
      //model_.getEnv().set(GRB_IntParam_Threads,1);
      
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
      if( MIP ){
        MainVars_ = IloNumVarArray(env_,noVars_,0.0,IloInfinity,LpVariable::Int);
      } else {
        MainVars_ = IloNumVarArray(env_,noVars_,0.0,IloInfinity,LpVariable::Float);
      }
      model_.add(MainVars_);
      if( noAuxVars_ > 0){
        MainAuxVars_ = IloNumVarArray(env_,noAuxVars_,0.0,IloInfinity,LpVariable::Float);
        model_.add(MainAuxVars_);
      }

      /* Add Factor Constraints */
      INDEX fixedVars = 0;
      IloNumArray ObjValues(env_,noVars_); 
      for(auto factorIt = factorBegin; factorIt != factorEnd; ++factorIt) {
        Offset_ = factorIt->GetPrimalOffset();
        size_ = factorIt->size();
        OffsetAux_ = factorIt->GetAuxOffset();
        sizeAux_ = factorIt->GetNumberOfAuxVariables();
        auto pot = factorIt->GetReparametrizedPotential();
        for(INDEX i=0;i<size_;++i) {
          REAL value = pot[i];
          if( std::isfinite(value) ){
            ObjValues[Offset_ + i] = value;
          } else {
            ObjValues[Offset_ + i] = 0.0;
            GetVariable(i).setBounds(0.0,0.0);
            fixedVars++;
          }
        }
        factorIt->CreateConstraints(this);
      }
      printf("Reparametrization fixed %d variables\n",fixedVars);
      IloObjective obj = IloMinimize(env_);
      obj.setLinearCoefs(MainVars_, ObjValues);
      model_.add(obj);
      
      /* Add Message Constraints */
      for(auto messageIt = messageBegin; messageIt != messageEnd; ++messageIt) {
        leftSize_ = messageIt->GetLeftFactor()->size();
        rightSize_ = messageIt->GetRightFactor()->size();
        OffsetLeft_ = messageIt->GetLeftFactor()->GetPrimalOffset();
        OffsetRight_ = messageIt->GetRightFactor()->GetPrimalOffset();
        messageIt->CreateConstraints(this);
      }

      cplex_ = IloCplex(model_);
      cplex_.setParam(IloCplex::TiLim,3600); 
      cplex_.setParam(IloCplex::Threads,1);
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
    }
    
    LinExpr CreateLinExpr() const { return LinExpr(env_); }
    
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
    void SetTimeLimit(REAL t){ cplex_.setParam(IloCplex::TiLim,t); }
    void SetNumberOfThreads(INDEX t){ cplex_.setParam(IloCplex::Threads,t); };
    
    void addLinearEquality(LinExpr lhs,LinExpr rhs);
    void addLinearInequality(LinExpr lhs,LinExpr rhs);

    template<class factor>
    void addFactor(const factor& f,INDEX offset);
    
    int solve();
    int solve(PrimalSolutionStorage::Element primal) { return solve(); }

    void WriteLpModel(std::string name){ cplex_.exportModel(name.c_str()); }
    
  private:
    IloEnv env_;
    IloModel model_;
    IloCplex cplex_;
    IloNumVarArray MainVars_;
    IloNumVarArray MainAuxVars_;
    
    INDEX Offset_,OffsetAux_,OffsetLeft_,OffsetRight_;
    INDEX size_,sizeAux_,leftSize_,rightSize_;
    INDEX noVars_,noAuxVars_;
    REAL epsilon_;
  };

  template<class factor>
  void LpInterfaceCplex::addFactor(const factor& f,INDEX offset){
    Offset_ = offset;
    size_ = f.size();
    f.CreateConstraints(this);
  }
  
  int LpInterfaceCplex::solve(){
    int status = 2;
    try{
      cplex_.solve();
      auto stat = cplex_.getCplexStatus();
      std::cout << "Cplex Status: " << stat << std::endl;
      if(stat == CPX_STAT_OPTIMAL){
        status = 0;
      }
      if(stat == CPX_STAT_INFEASIBLE){
        status = 1;
      }
    }
    catch (IloException& e) {
      std::cerr << "Concert exception caught: " << e << std::endl;
    }
    catch (...) {
      std::cerr << "Unknown exception caught" << std::endl;
    }
    return status;
  }

  REAL LpInterfaceCplex::GetVariableValue(const INDEX i) const{
    assert(i < noVars_);
    return cplex_.getValue(MainVars_[i]);
  }

  REAL LpInterfaceCplex::GetObjectiveValue() const{
    return cplex_.getObjValue();
  }
  
  void LpInterfaceCplex::SetVariableBound(LpVariable v,REAL lb,REAL ub,bool integer){
    v.setBounds(lb,ub);
    if(integer){
      model_.add(IloConversion(env_, v, IloNumVar::Int));
    } else {
      model_.add(IloConversion(env_, v, IloNumVar::Float));
    }
  }
  
  void LpInterfaceCplex::addLinearEquality(LinExpr lhs,LinExpr rhs){
    model_.add(IloRange(env_,0.0,lhs-rhs,0.0));//(lhs,GRB_EQUAL,rhs);
  }

  void LpInterfaceCplex::addLinearInequality(LinExpr lhs,LinExpr rhs){
    model_.add(IloRange(env_,0.0,rhs-lhs,IloInfinity));//addConstr(lhs,GRB_LESS_EQUAL,rhs);
  }
  
  REAL LpInterfaceCplex::GetBestBound() const{
    REAL bound = cplex_.getObjValue();
    if(cplex_.isMIP()){
      bound = cplex_.getBestObjValue();
    }
    return bound;
  }
  
}

#endif // LP_MP_LP_INTERFACE_CPLEX_HXX
