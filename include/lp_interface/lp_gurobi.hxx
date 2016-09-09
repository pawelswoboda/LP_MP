
#ifndef LP_MP_LP_INTERFACE_GUROBI_HXX
#define LP_MP_LP_INTERFACE_GUROBI_HXX

#include "lp_interface.h"

namespace LP_MP {
  
  class LpInterfaceGurobi : public LpInterfaceAdapter {
  public:

    LpInterfaceGurobi() : env_(GRBEnv()),model_(GRBModel(env_)) { }
    
    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    LpInterfaceGurobi(FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd, MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd)
      : env_(GRBEnv()),model_(GRBModel(env_)) {
      
      noVars_ = 0; 
      for(auto factorIt = factorBegin; factorIt != factorEnd; ++factorIt) {
        noVars_ += factorIt->size();
      }
      std::vector<double> obj(noVars_,1);
      std::vector<char> types(noVars_,GRB_INTEGER);
      MainVars_ = model_.addVars(NULL,NULL,&obj[0],&types[0],NULL,noVars_);
      model_.update();

      /* Add Factor Constraints */
      for(auto factorIt = factorBegin; factorIt != factorEnd; ++factorIt) {
        Offset_ = factorIt->GetPrimalOffset();
        printf("offset -> %d \n",Offset_);
        size_ = factorIt->size();
        for(INDEX i=0;i<size_;++i) {
          REAL value = factorIt->GetReparametrizedPotential()[i];
          if( std::isfinite(value) ){
            GetVariable(i).set(GRB_DoubleAttr_Obj,value);
          } else {
            GetVariable(i).set(GRB_DoubleAttr_Obj,0);
            GetVariable(i).set(GRB_DoubleAttr_LB,0);
            GetVariable(i).set(GRB_DoubleAttr_UB,0);
          }
        }
        factorIt->CreateConstraints(this);
      }

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

    LpVariable CreateAuxiliaryVariable(REAL lb,REAL ub,bool integer = false);

    INDEX GetFactorSize() const { return size_; }
    INDEX GetLeftFactorSize() const { return leftSize_; }
    INDEX GetRightFactorSize() const { return rightSize_; }
        
    LpVariable GetVariable(const INDEX i) const { assert(i < size_);  assert(Offset_ + i < noVars_); return MainVars_[Offset_ + i]; }
    LpVariable GetLeftVariable(const INDEX i) const { assert(i < leftSize_); assert(OffsetLeft_ + i < noVars_); return MainVars_[OffsetLeft_ + i]; }
    LpVariable GetRightVariable(const INDEX i) const { assert(i < rightSize_); assert(OffsetRight_ + i < noVars_); return MainVars_[OffsetRight_ + i]; }
    
    void SetVariableBound(LpVariable v,REAL lb,REAL ub,bool integer = false);
    
    void addLinearEquality(LinExpr lhs,LinExpr rhs);
    void addLinearInequality(LinExpr lhs,LinExpr rhs);
    
    void solve();
    void solve(PrimalSolutionStorage::Element primal) { solve(); }

    void WriteLpModel(std::string name){ model_.write(name); }
    
  private:
    GRBEnv env_;
    GRBModel model_;
    LpVariable* MainVars_;
    
    INDEX Offset_,OffsetLeft_,OffsetRight_;
    INDEX size_,leftSize_,rightSize_;
    INDEX noVars_;    
  };

  void LpInterfaceGurobi::solve(){
    
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

  LpVariable
  LpInterfaceGurobi::CreateAuxiliaryVariable(REAL lb,REAL ub,bool integer){
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

#endif // LP_MP_LP_INTERFACE_GUROBI_HXX
