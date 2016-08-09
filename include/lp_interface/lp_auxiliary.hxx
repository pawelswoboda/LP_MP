#ifndef LP_MP_LP_INTERFACE_AUX_HXX
#define LP_MP_LP_INTERFACE_AUX_HXX

#include "config.hxx"

namespace LP_MP {

  class LpVariable {
  public:
    LpVariable() : index(0),coeff(0),IsVar(false) {}
    INDEX index;
    REAL coeff;
    bool IsVar;
  };

  class LinExpr {
  public:
    LinExpr(){};
  private:
    std::vector<LpVariable> vars_;
  };
  
}

#endif // LP_MP_LP_INTERFACE_AUX_HXX
