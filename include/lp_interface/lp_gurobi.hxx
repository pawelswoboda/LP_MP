
#ifndef LP_MP_LP_INTERFACE_GUROBI_HXX
#define LP_MP_LP_INTERFACE_GUROBI_HXX

#include "lp_interface.h"

namespace LP_MP {
  
  class LpInterfaceGurobi : public LpInterfaceAdapter {
  public:

    LpInterfaceGurobi() : env_(GRBEnv()),model_(GRBModel(env_)){ }


    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    LpInterfaceGurobi( FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd, 
                        MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd,
                        REAL epsilon = std::numeric_limits<REAL>::infinity(),
                        bool integer = true)
      : env_(GRBEnv()),model_(GRBModel(env_)),
        LpInterfaceAdapter(factorBegin,factorEnd,messageBegin,messageEnd,
                          PrimalSolutionStorage(factorBegin,factorEnd),epsilon,integer)
    { 
      Setup(factorBegin,factorEnd,messageBegin,messageEnd);
    } 
    
    
    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    LpInterfaceGurobi( FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd, 
                        MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd,
                        PrimalSolutionStorage primal,
                        REAL epsilon = std::numeric_limits<REAL>::infinity(),
                        bool integer = true)
      : env_(GRBEnv()),model_(GRBModel(env_)),
        LpInterfaceAdapter(factorBegin,factorEnd,messageBegin,messageEnd,primal,epsilon,integer)
    { 
      Setup(factorBegin,factorEnd,messageBegin,messageEnd);
    }
    
    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    void Setup(FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd,
                MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd){
                  
      CreateVariables(factorBegin,factorEnd,messageBegin,messageEnd);
      model_.update();
      Build(factorBegin,factorEnd,messageBegin,messageEnd);
      model_.update();
      
      // Default Parameter for Gurobi
      model_.getEnv().set(GRB_DoubleParam_TimeLimit,3600);
      model_.getEnv().set(GRB_IntParam_Threads,1);
      model_.getEnv().set(GRB_IntParam_MIPFocus,1); 
    }
    
    LinExpr CreateLinExpr() { return LinExpr(); }
    
    void AddObjective(INDEX i,REAL value)
    {
      assert(i < pot_.size()); 
      GetVariable(i).set(GRB_DoubleAttr_Obj,value);
      //printf("var (%d) * %.2f \n",Offset_+i,value);
      LpInterfaceAdapter::AddObjective(i,value);      
    }    
    
    void CreateMainVariable(REAL lb, REAL ub,bool integer = false){
      auto TypeFlag = GRB_CONTINUOUS;
      if(integer){ TypeFlag = GRB_INTEGER; }
      MainVars_.push_back(model_.addVar(lb,ub,0.0,TypeFlag));
    }
    void CreateMainVariables(INDEX n,REAL lb, REAL ub,bool integer = false){
      auto TypeFlag = GRB_CONTINUOUS;
      if(integer){ TypeFlag = GRB_INTEGER; }
      for(INDEX i=0;i<n;i++){      
        MainVars_.push_back(model_.addVar(lb,ub,0.0,TypeFlag));
      }
    }
    
    void CreateAuxVariable(REAL lb, REAL ub,bool integer = false){
      auto TypeFlag = GRB_CONTINUOUS;
      if(integer){ TypeFlag = GRB_INTEGER; }
      MainAuxVars_.push_back(model_.addVar(lb,ub,0.0,TypeFlag));
    }
    void CreateAuxVariables(INDEX n,REAL lb, REAL ub,bool integer = false){
      auto TypeFlag = GRB_CONTINUOUS;
      if(integer){ TypeFlag = GRB_INTEGER; }
      for(INDEX i=0;i<n;i++){      
        MainAuxVars_.push_back(model_.addVar(lb,ub,0.0,TypeFlag));
      }
    } 
    
    REAL GetVariableValue(const INDEX i) const;
    REAL GetObjectiveValue() const;
    REAL GetBestBound() const;
    
    void SetTimeLimit(REAL t){ model_.getEnv().set(GRB_DoubleParam_TimeLimit,t); }
    void SetNumberOfThreads(INDEX t){ model_.getEnv().set(GRB_IntParam_Threads,t); };
    void SetDisplayLevel(INDEX t){ assert(t <= 1); model_.getEnv().set(GRB_IntParam_OutputFlag,t); };
    
    void addLinearEquality(LinExpr lhs,LinExpr rhs);
    void addLinearInequality(LinExpr lhs,LinExpr rhs);
    
    INDEX solve();
    
  private:
    GRBEnv env_;
    GRBModel model_;
  };
  
  INDEX LpInterfaceGurobi::solve(){
    INDEX status = 2;
    try{
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
  
  void LpInterfaceGurobi::addLinearEquality(LinExpr lhs,LinExpr rhs){
    model_.addConstr(lhs,GRB_EQUAL,rhs);
  }

    void LpInterfaceGurobi::addLinearInequality(LinExpr lhs,LinExpr rhs){
    model_.addConstr(lhs,GRB_LESS_EQUAL,rhs);
  }

  inline REAL LpInterfaceGurobi::GetVariableValue(const INDEX i) const{
    assert(i < noVars_);
    if(primal_[i] == false){
      return 0.0;
    }
    else{
      return MainVars_[i].get(GRB_DoubleAttr_X);
    }
  }

  REAL LpInterfaceGurobi::GetObjectiveValue() const{
    return model_.get(GRB_DoubleAttr_ObjVal);
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
