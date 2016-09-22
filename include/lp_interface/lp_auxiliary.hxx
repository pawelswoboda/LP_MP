#ifndef LP_MP_LP_INTERFACE_AUX_HXX
#define LP_MP_LP_INTERFACE_AUX_HXX

#include "config.hxx"

namespace LP_MP {

  class LpVariable {
  public:
    LpVariable() : index(0),coeff(0) {}
    LpVariable(INDEX idx,REAL obj) : index(idx),coeff(obj) {}

    LpVariable operator* (REAL value){
      coeff = value;
      return *this;
    }
    
    INDEX index;
    REAL coeff;
  };
    
  class LinExpr {
  public:
    LinExpr(){};
    
    LinExpr& operator+=(const LpVariable& var){
      vars_.push_back(var);
      return *this;
    }
    
    LinExpr& operator+=(const REAL value){
      constant_ += value;
      return *this;
    }

    LinExpr& operator-=(const LpVariable& var){
      assert(false);
      vars_.push_back(var);
      return *this;
    }
    
    LinExpr& operator-=(const REAL value){
      constant_ -= value;
      return *this;
    }
    
  private:
    std::vector<LpVariable> vars_;
    REAL constant_;
  };
  
}

#endif // LP_MP_LP_INTERFACE_AUX_HXX
