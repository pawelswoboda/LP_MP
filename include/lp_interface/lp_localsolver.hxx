
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

    LpInterfaceLocalSolver() { }

    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    LpInterfaceLocalSolver( FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd, 
                        MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd,
                        REAL epsilon = std::numeric_limits<REAL>::infinity(),
                        bool integer = true)
      : LpInterfaceAdapter(factorBegin,factorEnd,messageBegin,messageEnd,
                          PrimalSolutionStorage(factorBegin,factorEnd),epsilon,integer)
    { 
      Setup(factorBegin,factorEnd,messageBegin,messageEnd);
    } 
    
    
    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    LpInterfaceLocalSolver( FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd, 
                        MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd,
                        PrimalSolutionStorage primal,
                        REAL epsilon = std::numeric_limits<REAL>::infinity(),
                        bool integer = true)
      : LpInterfaceAdapter(factorBegin,factorEnd,messageBegin,messageEnd,primal,epsilon,integer)
    { 
      Setup(factorBegin,factorEnd,messageBegin,messageEnd);
    }
    
    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    void Setup(FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd,
                MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd){
      try{           
        model_ = localsolver_.getModel();
        cb_ = FeasibleSolutionCallback();
        localsolver_.addCallback(localsolver::CT_IterationTicked,&cb_);
        obj_ = CreateLinExpr();
        
        CreateVariables(factorBegin,factorEnd,messageBegin,messageEnd);
        Build(factorBegin,factorEnd,messageBegin,messageEnd);
        
        model_.minimize(obj_);
        
      } catch (const localsolver::LSException& e) {
        std::cout << "LSException:" << e.getMessage() << std::endl;
        exit(1);
      }
       
    }
     
    LinExpr CreateLinExpr() { return model_.sum(); }
    
    void AddObjective(INDEX i,REAL value)
    {
      assert(i < pot_.size());
      if( IsObjective(i) == false){
        obj_ += value*GetVariable(i);
        LpInterfaceAdapter::AddObjective(i,value);
      }
    }    
    
    void CreateMainVariable(REAL lb, REAL ub,bool integer = false){
      if(integer){
        MainVars_.push_back(model_.boolVar());
      } else {
        MainVars_.push_back(model_.floatVar(lb,ub));
      }
    }
    void CreateMainVariables(INDEX n,REAL lb, REAL ub,bool integer = false){
      for(INDEX i=0;i<n;i++){      
        CreateMainVariable(lb,ub,integer);
      }
    }
    
    void CreateAuxVariable(REAL lb, REAL ub,bool integer = false){
     if(integer){
        MainAuxVars_.push_back(model_.boolVar());
      } else {
        MainAuxVars_.push_back(model_.floatVar(lb,ub));
      }
    }
    void CreateAuxVariables(INDEX n,REAL lb, REAL ub,bool integer = false){
       for(INDEX i=0;i<n;i++){      
        CreateAuxVariable(lb,ub,integer);
      }
    } 
 
    REAL GetVariableValue(const INDEX i) const;
    REAL GetObjectiveValue() const;
    REAL GetBestBound() const;
    
    void SetTimeLimit(REAL t){ timelimit_ = t; };
    void SetNumberOfThreads(INDEX t){ noThreads_ = t; }//param_.setNbThreads(t); };
    void SetDisplayLevel(INDEX t){ assert(t <= 1); display_ = t; }  //param_.setVerbosity(t); };
    
    void addLinearEquality(LinExpr lhs,LinExpr rhs);
    void addLinearInequality(LinExpr lhs,LinExpr rhs);
    
    INDEX solve();
    
  private:
    localsolver::LocalSolver localsolver_;
    localsolver::LSModel model_;
    localsolver::LSPhase phase_; //localsolver_.createPhase();;
    localsolver::LSParam param_ = localsolver_.getParam();
        
    LinExpr obj_;
    FeasibleSolutionCallback cb_;
    
    REAL timelimit_ = 3600;
    REAL noThreads_ = 1;
    INDEX display_ = 0;
  };
  
  INDEX LpInterfaceLocalSolver::solve(){
    INDEX status = 2;
    try{      
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
