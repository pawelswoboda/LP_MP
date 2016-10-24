#ifndef LP_MP_LP_INTERFACE_CPLEX_HXX
#define LP_MP_LP_INTERFACE_CPLEX_HXX

#include "lp_interface.h"

namespace LP_MP {
  
  class LpInterfaceCplex : public LpInterfaceAdapter {
  public:

    LpInterfaceCplex() : env_(IloEnv()),model_(IloModel(env_)),obj_(IloMinimize(env_)){ }

    ~LpInterfaceCplex(){
      env_.end();
    }
    
    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    LpInterfaceCplex( FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd, 
                        MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd,
                        REAL epsilon = std::numeric_limits<REAL>::infinity(),
                        bool integer = true)
      : env_(IloEnv()),model_(IloModel(env_)),obj_(IloMinimize(env_)),
        LpInterfaceAdapter(factorBegin,factorEnd,messageBegin,messageEnd,
                          PrimalSolutionStorage(factorBegin,factorEnd),epsilon,integer)
    { 
      Setup(factorBegin,factorEnd,messageBegin,messageEnd);
    } 
    
    
    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    LpInterfaceCplex( FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd, 
                        MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd,
                        PrimalSolutionStorage primal,
                        REAL epsilon = std::numeric_limits<REAL>::infinity(),
                        bool integer = true)
      : env_(IloEnv()),model_(IloModel(env_)),obj_(IloMinimize(env_)),
        LpInterfaceAdapter(factorBegin,factorEnd,messageBegin,messageEnd,primal,epsilon,integer)
    { 
      Setup(factorBegin,factorEnd,messageBegin,messageEnd);
    }

    
    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    void Setup(FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd,
                MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd){

      ObjVars_ = IloNumVarArray(env_);

      CreateVariables(factorBegin,factorEnd,messageBegin,messageEnd);

      ObjValues_ = IloNumArray(env_,MainVars_.size());
      
      Build(factorBegin,factorEnd,messageBegin,messageEnd);
            
      obj_.setLinearCoefs(ObjVars_,ObjValues_);
      model_.add(obj_);
      
      // Default Parameter for CPLEX
      cplex_ = IloCplex(model_);
      cplex_.setParam(IloCplex::TiLim,3600); 
      cplex_.setParam(IloCplex::Threads,1);
      cplex_.setParam(IloCplex::MIPEmphasis, CPX_MIPEMPHASIS_HIDDENFEAS); // CPX_MIPEMPHASIS_FEASIBILITY
      cplex_.setParam(IloCplex::Param::MIP::Strategy::Probe, 0); // probing
    }
    
    LinExpr CreateLinExpr() { return LinExpr(env_); }

    void AddObjective(INDEX i,REAL value)
    {
      assert(i < pot_.size());
      ObjValues_[Offset_+i]=value;
      LpInterfaceAdapter::AddObjective(i,value);      
    }  
    
    void CreateMainVariable(REAL lb, REAL ub,bool integer = false){
      auto TypeFlag = LpVariable::Float;
      if(integer){ TypeFlag = LpVariable::Int; }
      MainVars_.push_back(IloNumVar(env_,lb,ub,TypeFlag));
      ObjVars_.add(MainVars_.back());
    }
    void CreateMainVariables(INDEX n,REAL lb, REAL ub,bool integer = false){
      auto TypeFlag = LpVariable::Float;
      if(integer){ TypeFlag = LpVariable::Int; }
      for(INDEX i=0;i<n;i++){      
        MainVars_.push_back(IloNumVar(env_,lb,ub,TypeFlag));
      }
    }
    
    void CreateAuxVariable(REAL lb, REAL ub,bool integer = false){
      auto TypeFlag = LpVariable::Float;
      if(integer){ TypeFlag = LpVariable::Int; }
      MainAuxVars_.push_back(IloNumVar(env_,lb,ub,TypeFlag));
    }
    void CreateAuxVariables(INDEX n,REAL lb, REAL ub,bool integer = false){
      auto TypeFlag = LpVariable::Float;
      if(integer){ TypeFlag = LpVariable::Int; }
      for(INDEX i=0;i<n;i++){      
        MainAuxVars_.push_back(IloNumVar(env_,lb,ub,TypeFlag));
      }
    } 
    
    REAL GetVariableValue(const INDEX i) const;
    REAL GetObjectiveValue() const;
    REAL GetBestBound() const;
    
    
    void SetTimeLimit(REAL t){ cplex_.setParam(IloCplex::TiLim,t); }
    void SetNumberOfThreads(INDEX t){ cplex_.setParam(IloCplex::Threads,t); };
    void SetDisplayLevel(INDEX t){ 
      assert(t <= 1); 
      if( t == 0){ cplex_.setOut(env_.getNullStream()); }
    };
    
    void addLinearEquality(LinExpr lhs,LinExpr rhs);
    void addLinearInequality(LinExpr lhs,LinExpr rhs);
    
    INDEX solve();
    
  private:
    IloEnv env_;
    IloModel model_;
    IloCplex cplex_;
    IloObjective obj_;
    
    IloNumVarArray ObjVars_;
    IloNumArray ObjValues_;
  };
  
  INDEX LpInterfaceCplex::solve(){
    INDEX status = 2;
    try{
      cplex_.solve();
      auto stat = cplex_.getCplexStatus();
      //std::cout << "Cplex Status: " << stat << std::endl;
      if(stat == CPX_STAT_OPTIMAL){
        status = 0;
      }
      else if(cplex_.getSolnPoolNsolns() > 0){
        status = 3;
      }
      else if(stat == CPX_STAT_INFEASIBLE){
        status = 1;
      }
      else if(stat == CPX_STAT_ABORT_TIME_LIM){
        status = 4;
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
    if(primal_[i] == false){
      return 0.0;
    } else {
      return cplex_.getValue(ObjVars_[i]);//MainVars_[i]);;
    }
  }

  REAL LpInterfaceCplex::GetObjectiveValue() const{
    return cplex_.getObjValue();
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
