#ifndef LP_MP_LP_INTERFACE_HXX
#define LP_MP_LP_INTERFACE_HXX

#include <gurobi_c++.h>

namespace LP_MP {

  class LpInterfaceGurobi {
  public:

    typedef GRBVar LpVariable;
    typedef GRBLinExpr LinExpr;
    
    
    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    LpInterfaceGurobi(FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd, MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd)
    {
      env_ = GRBEnv();
      model_ = GRBModel(env);
      
      /* TODO: How to get number of all vars  */
      noVars_ = 10; // <--- TODO
      std::vector<double> obj(noVars_,1);
      std::vector<char> types(noVars_,GRB_INTEGER);
      MainVars_ = model.addVars(NULL,NULL,&obj[0],&types[0],NULL,noVars_);
      model_.update();

      /* TODO: Factor Constraints */
      for(auto factorIt = factorBegin; factorIt != factorEnd; ++factorIt) {
        Offset_ = factorIt->GetPrimalOffset();

        factorIt->CreateConstraints(this);
      }

      /* TODO: Message Constraints */
      for(auto messageIt = messageBegin; messageIt != messageEnd; ++messageIt) {
        Offset_ = messageIt->GetPrimalOffset();

        messageIt->CreateConstraints(this);
      }
      
    } 

    LpVariable CreateAuxiliaryVariable(REAL lb,REAL ub,bool integer);
    LpVariable GetVariable(const INDEX i){ return MainVars_[Offset_ + i]; }
    void addLinearEquation(LinExpr lhs,LinExpr rhs);
    void addLinearInequality(LinExpr lhs,LinExpr rhs);
    
  private:
    GRBEnv env_;
    GRBModel model_;
    LpVariable* MainVars_;
    
    INDEX Offset_;
    INDEX noVars_;    
  };

  void LpInterfaceGurobi::addLinearEquation(LinExpr lhs,LinExpr rhs){
    model_.addConstr(lhs,GRB_EQUAL,rhs);
    //model_.update(); // We may just update it right before optimization
  }

    void LpInterfaceGurobi::addLinearInequality(LinExpr lhs,LinExpr rhs){
    model_.addConstr(lhs,GRB_LESS_EQUAL,rhs);
    //model_.update(); // We may just update it right before optimization
  }
  
  typename LpInterfaceGurobi::LpVariable
  LpInterfaceGurobi::CreateAuxiliaryVariable(REAL lb,REAL ub,bool integer = false){
    LpVariable v;
    if( integer ){
      v = model_.addVar(lb,ub,0.0,GRB_INTEGER);
    } else {
      v = model_.addVar(lb,ub,0.0,GRB_CONTINUOUS);
    }
    model_.update();
    return v;
  }
  
}

#endif // LP_MP_LP_INTERFACE_HXX
