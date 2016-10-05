
#ifndef LP_MP_LP_INTERFACE_LOCALSOLVER_HXX
#define LP_MP_LP_INTERFACE_LOCALSOLVER_HXX

#include "lp_interface.h"


namespace LP_MP {

  class FeasibleSolutionCallback : public localsolver::LSCallback {
  public:
    FeasibleSolutionCallback(){
            
    }

    void callback(localsolver::LocalSolver& solver,localsolver::LSCallbackType type) {
      localsolver::LSSolution sol = solver.getSolution();
      localsolver::LSSolutionStatus stat = sol.getStatus();
      if(stat == localsolver::SS_Feasible){
        solver.stop();
      }
    }
    
  private:

  };
  
  class LpInterfaceLocalSolver : public LpInterfaceAdapter {
  public:

    LpInterfaceLocalSolver(INDEX noVars) : noVars_(noVars),epsilon_(std::numeric_limits<REAL>::infinity()) {
      model_ = localsolver_.getModel();
      MainVars_ = std::vector<LpVariable>(noVars_);
      for(INDEX i=0;i<noVars_;i++){
        MainVars_[i] = model_.boolVar();
      }
    }
    
    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    LpInterfaceLocalSolver(FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd, MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd,bool MIP = true)
      : epsilon_(std::numeric_limits<REAL>::infinity()) {

      try{
        model_ = localsolver_.getModel();
        cb_ = FeasibleSolutionCallback();
        localsolver_.addCallback(localsolver::CT_IterationTicked,&cb_);
        
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
      
        MainVars_ = std::vector<LpVariable>(noVars_);
        if(MIP){
          for(INDEX i=0;i<noVars_;i++){
            MainVars_[i] = model_.boolVar();
          }        
        } else {
          for(INDEX i=0;i<noVars_;i++){
            MainVars_[i] = model_.floatVar(0.0,1.0);
          }   
        }
        if( noAuxVars_ > 0){
          MainAuxVars_ = std::vector<LpVariable>(noAuxVars_);
          for(INDEX i=0;i<noAuxVars_;i++){
            MainAuxVars_[i] = model_.floatVar(-std::numeric_limits<REAL>::max(),std::numeric_limits<REAL>::max());
          }
        }
      
        VariableBounds_.resize(model_.getNbExpressions());
        isBoundSet_.resize(model_.getNbExpressions(),false);
        
        obj_ = CreateLinExpr();
           
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
              obj_ += value*GetVariable(i);
            } else {
              model_.constraint(GetVariable(i) == 0);
              fixedVars++;
            }
          }
          factorIt->CreateConstraints(this);
        }
        model_.minimize(obj_);
      
        /* Add Message Constraints */
        for(auto messageIt = messageBegin; messageIt != messageEnd; ++messageIt) {
          leftSize_ = messageIt->GetLeftFactor()->size();
          rightSize_ = messageIt->GetRightFactor()->size();
          OffsetLeft_ = messageIt->GetLeftFactor()->GetPrimalOffset();
          OffsetRight_ = messageIt->GetRightFactor()->GetPrimalOffset();
          messageIt->CreateConstraints(this);
        }

      } catch (const localsolver::LSException& e) {
        std::cout << "LSException:" << e.getMessage() << std::endl;
        exit(1);
      }
    } 

    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    void ReduceLp(FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd, MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd,REAL epsilon){
      try{
        epsilon_ = epsilon;
        for(auto factorIt = factorBegin; factorIt != factorEnd; ++factorIt) {
          Offset_ = factorIt->GetPrimalOffset();
          size_ = factorIt->size();
          OffsetAux_ = factorIt->GetAuxOffset();
          sizeAux_ = factorIt->GetNumberOfAuxVariables();
          factorIt->ReduceLp(this);
        }
      } catch (const localsolver::LSException& e) {
        std::cout << "LSException:" << e.getMessage() << std::endl;
        //exit(1);
      }
    }
    
    LinExpr CreateLinExpr() { return model_.sum(); }
    
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
    void SetTimeLimit(REAL t){ timelimit_ = t; };
    void SetNumberOfThreads(INDEX t){ noThreads_ = t; }//param_.setNbThreads(t); };
    void SetDisplayLevel(INDEX t){ assert(t <= 1); display_ = t; }  //param_.setVerbosity(t); };
    
    void addLinearEquality(LinExpr lhs,LinExpr rhs);
    void addLinearInequality(LinExpr lhs,LinExpr rhs);

    template<class factor>
    void addFactor(const factor& f,INDEX offset);
    
    int solve();
    int solve(PrimalSolutionStorage::Element primal);

    void WriteLpModel(std::string name){  }
    
  private:
    localsolver::LocalSolver localsolver_;
    localsolver::LSModel model_;
    localsolver::LSPhase phase_; //localsolver_.createPhase();;
    localsolver::LSParam param_ = localsolver_.getParam();
    
    
    std::vector<LpVariable> MainVars_;
    std::vector<LpVariable> MainAuxVars_;
    std::vector<LpVariable> VariableBounds_;
    std::vector<bool> isBoundSet_;

    LinExpr obj_;
    FeasibleSolutionCallback cb_;
    
    INDEX Offset_,OffsetAux_,OffsetLeft_,OffsetRight_;
    INDEX size_,sizeAux_,leftSize_,rightSize_;
    INDEX noVars_,noAuxVars_;    
    REAL epsilon_;
    REAL timelimit_ = 3600;
    REAL noThreads_ = 1;
    INDEX display_ = 0;
  };

  template<class factor>
  void LpInterfaceLocalSolver::addFactor(const factor& f,INDEX offset){
    Offset_ = offset;
    size_ = f.size();
    f.CreateConstraints(this);
  }
  
  int LpInterfaceLocalSolver::solve(){
    int status = 2;
    try{
      for(INDEX i=0;i<VariableBounds_.size();i++){
        if(isBoundSet_[i]){
          model_.constraint(VariableBounds_[i]);
        }
      }
      
      model_.close();
      
      localsolver::LSPhase phase_ = localsolver_.createPhase();//getPhase(0);
      phase_.setTimeLimit(timelimit_);

      param_.setNbThreads(noThreads_);
      param_.setVerbosity(display_);
      param_.setIterationBetweenTicks(2500);
      
      localsolver_.solve();
      localsolver::LSSolution sol = localsolver_.getSolution();
      localsolver::LSSolutionStatus stat = sol.getStatus();
      localsolver::LSStatistics st = localsolver_.getStatistics();
      
      if(stat == localsolver::SS_Optimal ){
        status = 0;
      }
      else if(stat == localsolver::SS_Feasible){
        status = 3;
      }
      else if(st.getRunningTime() >= ((int) timelimit_)){
        status = 4;
      }
      else if(stat == localsolver::SS_Infeasible){
        status = 1;
      }
      
      
    }
    catch (const localsolver::LSException& e) {
      std::cout << "LSException:" << e.getMessage() << std::endl;
    }
    return status;
  }

  int LpInterfaceLocalSolver::solve(PrimalSolutionStorage::Element primal){
    // go through primal and when it is not unknownState, set variables to zero in model
         
    int status = 2;
    try{
      for(INDEX i=0;i<VariableBounds_.size();i++){
        if(isBoundSet_[i]){
          model_.constraint(VariableBounds_[i]);
        }
      }      
      
      for(INDEX i=0; i<noVars_; ++i) {
        if(primal[i] == true) {
          SetVariableBound(MainVars_[i],1,1,false);
        } else if(primal[i] == false) {
          SetVariableBound(MainVars_[i],0,0,false);
        } else {
          assert(primal[i] == unknownState);
        }
      }
      
      model_.close();
      
      localsolver::LSPhase phase_ = localsolver_.createPhase();//getPhase(0);
      phase_.setTimeLimit(timelimit_);

      param_.setNbThreads(noThreads_);
      param_.setVerbosity(display_);
            
      localsolver_.solve();
      localsolver::LSSolution sol = localsolver_.getSolution();
      localsolver::LSSolutionStatus stat = sol.getStatus();
      
      if(stat == localsolver::SS_Optimal ){
        status = 0;
      }
      else if(stat == localsolver::SS_Feasible){
        status = 3;
      }
      else if(stat == localsolver::SS_Infeasible){
        status = 1;
      }
      
    }
    catch (const localsolver::LSException& e) {
      std::cout << "LSException:" << e.getMessage() << std::endl;
    }
    // now return to original model without fixed variables
    return status;
  }

  REAL LpInterfaceLocalSolver::GetVariableValue(const INDEX i) const{
    try{
      assert(i < noVars_);
      if(MainVars_[i].isInt()){
        return MainVars_[i].getIntValue();
      } else {
        return MainVars_[i].getDoubleValue(); 
      }
    } catch (const localsolver::LSException& e) {
      std::cout << "LSException (GetVariableValue):" << e.getMessage() << std::endl;
      exit(1);
    }
  }

  REAL LpInterfaceLocalSolver::GetObjectiveValue() const{
    try{
      localsolver::LSSolution sol = localsolver_.getSolution();
      if( obj_.isInt() ){
        return sol.getIntValue(obj_);
      } else {
        return sol.getDoubleValue(obj_);
      }
    } catch (const localsolver::LSException& e) {
      std::cout << "LSException (GetObjetiveValue):" << e.getMessage() << std::endl;
      exit(1);
    }
  }
  
  void LpInterfaceLocalSolver::SetVariableBound(LpVariable v,REAL lb,REAL ub,bool integer){
    INDEX i = v.getIndex();
    assert(i < VariableBounds_.size());
    if(!isBoundSet_[i]){
      isBoundSet_[i] = true;
    }
    VariableBounds_[i] = model_.and_(lb <= v,v <= ub); 
  }
  
  void LpInterfaceLocalSolver::addLinearEquality(LinExpr lhs,LinExpr rhs){
    model_.constraint(lhs == rhs);
  }

  void LpInterfaceLocalSolver::addLinearInequality(LinExpr lhs,LinExpr rhs){
    model_.constraint(lhs <= rhs);
  }
  
  
  REAL LpInterfaceLocalSolver::GetBestBound() const{
    localsolver::LSSolution sol = localsolver_.getSolution();
    if( obj_.isInt() ){
      return sol.getIntValue(obj_);
    } else {
      return sol.getDoubleValue(obj_);
    }
  }
  
}

#endif // LP_MP_LP_INTERFACE_GUROBI_HXX
