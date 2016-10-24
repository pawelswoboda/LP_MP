#ifndef LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX
#define LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX

#include "LP_MP.h"
#include "minConv.hxx"
#include <math.h> 

namespace LP_MP{

  using MinConv = discrete_tomo::MinConv<REAL,INDEX>;
  
  class DiscreteTomographyFactorCounting{

  public:
 
    enum class NODE {left,right,up,reg};
    
    DiscreteTomographyFactorCounting(INDEX numberOfLabels,INDEX numberOfVarsLeft,INDEX numberOfVarsRight,INDEX SumBound);
    
    //template<typename REPAM_ARRAY>
    //INDEX ComputeOptimalLabeling(const REPAM_ARRAY& repam) const { }; //--required

    template<typename REPAM_ARRAY>
    REAL LowerBound(const REPAM_ARRAY& repam) const; //--required

    //template<typename REPAM_ARRAY>
    //void MaximizePotentialAndComputePrimal(const REPAM_ARRAY& repam, PrimalSolutionStorage::Element primal) const;
    
    //template<typename REPAM_ARRAY>
    //void NarrowPrimal(const REPAM_ARRAY& repam, const REAL epsilon, typename PrimalSolutionStorage::Element primal);
    
    void PropagatePrimal(PrimalSolutionStorage::Element primal) const;
    
    template<typename REPAM_ARRAY>
    REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const; //--required

    INDEX size() const;
    
    //void WritePrimal(const PrimalSolutionStorage::Element, std::ofstream& fs) const; //

    // Lp Interface related
    void CreateConstraints(LpInterfaceAdapter* lp) const;
    void CreateAuxVariables(LpInterfaceAdapter* lp) const {
      lp->CreateAuxVariables(leftSize_*rightSize_,0.0,1.0);
    }
    
    INDEX getSize(NODE) const;

    template<typename REPAM_ARRAY>
    REAL eval(INDEX,INDEX,INDEX,const REPAM_ARRAY& repam);
    
  private:

    const INDEX numberOfLabels_,numberOfVarsLeft_,numberOfVarsRight_,SumBound_;
    INDEX upSize_,leftSize_,rightSize_,regSize_;
  };
  
  template<typename REPAM_ARRAY>
  REAL DiscreteTomographyFactorCounting::eval(INDEX up,INDEX left,INDEX right,const REPAM_ARRAY& repam){
    assert(up < upSize_);
    assert(left < leftSize_);
    assert(right < rightSize_);
    assert(repam.size() == (upSize_ + leftSize_ + rightSize_ + regSize_));

    auto xa = [&](INDEX idx){ return idx % numberOfLabels_;  };
    auto xb = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return xa(idx); };
    auto z = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return (idx - xa(idx))/numberOfLabels_; };
    
    if( xa(up) == xa(left) &&
        xb(up) == xb(right)&&
        z(left) + z(right) == z(up) )
      { return
          repam[up] +
          repam[upSize_ + left] +
          repam[upSize_ + leftSize_ + right] +
          repam[upSize_ + leftSize_ + rightSize_ + xb(left) + xa(right)*pow(numberOfLabels_,2)]; }
    else
      { return std::numeric_limits<REAL>::infinity(); }
    
  }

  inline INDEX DiscreteTomographyFactorCounting::size() const {
    return getSize(DiscreteTomographyFactorCounting::NODE::up) 
      + getSize(DiscreteTomographyFactorCounting::NODE::left) 
      + getSize(DiscreteTomographyFactorCounting::NODE::right) 
      + getSize(DiscreteTomographyFactorCounting::NODE::reg); 
  }

  
  inline INDEX DiscreteTomographyFactorCounting::getSize(NODE n) const {
    if( n == NODE::left ){ return leftSize_;  }
    if( n == NODE::right ){ return rightSize_;  }
    if( n == NODE::up ){ return upSize_;  }
    if( n == NODE::reg ){ return regSize_;  }
    return 0;
  }
   
  inline DiscreteTomographyFactorCounting::DiscreteTomographyFactorCounting(INDEX numberOfLabels,INDEX numberOfVarsLeft,INDEX numberOfVarsRight,INDEX SumBound)
    : numberOfLabels_(numberOfLabels),numberOfVarsLeft_(numberOfVarsLeft),numberOfVarsRight_(numberOfVarsRight),SumBound_(SumBound) {
    assert(numberOfLabels_ > 1);
    assert(numberOfVarsLeft_ > 0);
    assert(numberOfVarsRight_ > 0);
    assert(SumBound_ > 0);
    
    upSize_ = pow(numberOfLabels_,2)*std::min(((numberOfVarsLeft_+numberOfVarsRight_)*(numberOfLabels_-1)+1),SumBound);
    leftSize_ = pow(numberOfLabels_,2)*std::min((numberOfVarsLeft_*(numberOfLabels_-1)+1),SumBound);
    rightSize_ = pow(numberOfLabels_,2)*std::min((numberOfVarsRight_*(numberOfLabels_-1)+1),SumBound);
    
    regSize_ = pow(numberOfLabels_,2);
  }

  inline void DiscreteTomographyFactorCounting::CreateConstraints(LpInterfaceAdapter* lp) const {
    REAL inf = std::numeric_limits<REAL>::infinity();
    
    auto xa = [&](INDEX idx){ return idx % numberOfLabels_;  };
    auto xb = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return xa(idx); };
    auto xz = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return (idx - xa(idx))/numberOfLabels_; };
   
    LinExpr lhs_all = lp->CreateLinExpr();
    lhs_all += 0;
    std::vector<LinExpr> lhs_up(upSize_);
    std::vector<LinExpr> lhs_left(leftSize_);
    std::vector<LinExpr> lhs_right(rightSize_);
    std::vector<LinExpr> lhs_reg(regSize_);
    
    auto InitVector = [&](std::vector<LinExpr>& v){
      for(INDEX i=0;i<v.size();i++){
        v[i] = lp->CreateLinExpr();
        v[i] += 0.0;
      }
    };
    InitVector(lhs_up);
    InitVector(lhs_left);
    InitVector(lhs_right);
    InitVector(lhs_reg);

    INDEX u_prim = upSize_;
    INDEX l_prim = leftSize_;
    INDEX r_prim = rightSize_;
    INDEX p_prim = regSize_;

    auto GetPrims = [&](INDEX size,INDEX offset,INDEX &prim){
      INDEX count = 0;
      for(INDEX i=offset;i<offset+size;i++){
        if(lp->IsPrimal(i) == true){
          prim = i;
          count++;
        }
      }
      assert(count <= 1);
    };
    GetPrims(upSize_,0,u_prim);
    GetPrims(leftSize_,upSize_,l_prim);
    GetPrims(rightSize_,upSize_+leftSize_,r_prim);
    GetPrims(regSize_,upSize_+leftSize_+rightSize_,p_prim);
    
    INDEX z_max = upSize_/pow(numberOfLabels_,2);
     
    REAL epsi = lp->GetEpsilon();
    auto repam = lp->GetRepam();
    REAL lb = LowerBound(repam);
    
    INDEX size = upSize_ + leftSize_ + rightSize_ + regSize_;
    
    INDEX noVars = 0;
    for(INDEX i=0;i<leftSize_;i++){
      for(INDEX j=0;j<rightSize_;j++){
        LpVariable var;
        {
          // up variable constraint
          // sum_{ i(z) + j(z) = k(z) && i(a) = k(a) && j(b) = k(b) } eta(i,j) = mu_u(k)
          INDEX lz = xz(i);
          INDEX rz = xz(j);
          INDEX z = lz + rz;
          if( z < z_max ){
            
            INDEX a = xa(i);
            INDEX b = xb(i);
            INDEX c = xa(j);
            INDEX d = xb(j);

            INDEX rz = xz(j);
            INDEX uz = lz+rz; 

            INDEX l = a + b*numberOfLabels_ + lz*pow(numberOfLabels_,2);
            INDEX r = c + d*numberOfLabels_ + rz*pow(numberOfLabels_,2);
            INDEX u = a + d*numberOfLabels_ + uz*pow(numberOfLabels_,2);
            INDEX p = b + c*numberOfLabels_;
            
            assert(l < leftSize_);
            assert(r < rightSize_);
            assert(u < upSize_);
            assert(p < pow(numberOfLabels_,2));
            
            REAL value = repam[u] + repam[upSize_+l] + repam[upSize_+leftSize_+r] + repam[upSize_+leftSize_+rightSize_+p];
             
            if( value >= lb + epsi ){
              continue;
            }
            if( u_prim < upSize_ && u != u_prim ){
              continue;
            }
            if( l_prim < leftSize_ && l != l_prim ){
              continue;
            }
            if( r_prim < rightSize_ && r != r_prim ){
              continue;
            }
            if( p_prim < regSize_ && p != p_prim ){
              continue;
            }
            
            var = lp->GetAuxVariable(i + j*leftSize_);
            lhs_up[u] += var;
            
            lp->AddObjective(u,repam[u]);
            lp->AddObjective(upSize_+l,repam[upSize_+l]);
            lp->AddObjective(upSize_+leftSize_+r,repam[upSize_+leftSize_+r]);
            lp->AddObjective(upSize_+leftSize_+rightSize_+p,repam[upSize_+leftSize_+rightSize_+p]);
            
          } else {
            continue;
          }
        }
        {
          // sum_j eta(i,j) = mu_l(i) 
          lhs_left[i] += var;
          
          // sum_i eta(i,j) = mu_r(j) 
          lhs_right[j] += var;
        }
        {
          // sum_{i,j} eta(i,j) = 1
          lhs_all += var;
          noVars++;
        }
        {
          // pairwise potential
          // sum_{i(b) = k(a) && j(a) = k(b)} eta(i,j) = mu_p(k)
          INDEX a = xb(i);
          INDEX b = xa(j);
          
          lhs_reg[a+b*numberOfLabels_] += var;
        }
      }
    }
    
    auto AddAllConstraints = [&]( std::vector<LinExpr> lhs,
                                  INDEX offset){
      for(INDEX i=0;i<lhs.size();i++){
        LinExpr rhs = lp->CreateLinExpr();
        rhs += lp->GetVariable(offset + i);
        lp->addLinearEquality(lhs[i],rhs);
      }
    };
  
    LinExpr rhs_all = lp->CreateLinExpr();
    rhs_all += 1;
    lp->addLinearEquality(lhs_all,rhs_all);
    
    AddAllConstraints(lhs_up,0);
    AddAllConstraints(lhs_left,upSize_);
    AddAllConstraints(lhs_right,upSize_+leftSize_);
    AddAllConstraints(lhs_reg,upSize_+leftSize_+rightSize_);
    
  }
  
  template<typename REPAM_ARRAY>
  REAL DiscreteTomographyFactorCounting::LowerBound(const REPAM_ARRAY& repam) const{
    assert(repam.size() == (upSize_ + leftSize_ + rightSize_ + regSize_));

    const INDEX noLabels = numberOfLabels_;
    const INDEX up_size = upSize_/pow(noLabels,2);
    const INDEX left_size = leftSize_/pow(noLabels,2);
    const INDEX right_size = rightSize_/pow(noLabels,2);
    
    REAL bound = std::numeric_limits<REAL>::infinity();
    
    if( 5000 < upSize_ ){
      auto op = [&](INDEX i,INDEX j){ return (i+j < up_size) ? i+j : up_size;  }; // 0 <= i+j < up_size
      
      for(INDEX i=0;i<pow(noLabels,4);i++){
        INDEX idx = i;
        INDEX a = idx % noLabels;
        idx = ( idx - a )/noLabels;
        INDEX b = idx % noLabels;
        idx = ( idx - b )/noLabels;
        INDEX c = idx % noLabels;
        idx = ( idx - c )/noLabels;
        INDEX d = idx % noLabels;

        auto z_up = [&](INDEX k){
          return repam[a + d*noLabels + k*pow(noLabels,2)];  };
        auto z_left = [&](INDEX k){
          return repam[up_size*pow(noLabels,2) + a + b*noLabels + k*pow(noLabels,2)];  };
        auto z_right = [&](INDEX k){
          return repam[up_size*pow(noLabels,2) + left_size*pow(noLabels,2) + c + d*noLabels + k*pow(noLabels,2)];  };

        REAL reg = repam[up_size*pow(noLabels,2) + left_size*pow(noLabels,2) + right_size*pow(noLabels,2) + b + c*noLabels];
        assert(reg > -std::numeric_limits<REAL>::max());

        MinConv mc(z_left,z_right,left_size,right_size,up_size);
        mc.CalcConv(op,z_left,z_right);

        for(INDEX k=0;k<up_size;k++){
          assert(k == op(mc.getIdxA(k),mc.getIdxB(k)));
          assert(k == (mc.getIdxA(k) + mc.getIdxB(k)));
          assert(!std::isnan(reg));

          REAL val = mc.getConv(k) + z_up(k)  + reg;
          
          assert(!std::isnan(val));
          assert(val > -eps);
          bound = std::min(bound,val);
          assert(bound > -eps);
        }
      }
    } else{
      for(INDEX lz=0;lz<left_size;lz++){
        for(INDEX rz=0;rz<right_size && rz+lz<up_size;rz++){
          for(INDEX i=0;i<pow(noLabels,4);i++){
            INDEX idx = i;
            INDEX a = idx % noLabels;
            idx = (idx - a)/noLabels;
            INDEX b = idx % noLabels;
            idx = (idx - b)/noLabels;
            INDEX c = idx % noLabels;
            INDEX d = (idx - c)/noLabels;

            assert(a + d*noLabels + (lz+rz)*pow(noLabels,2) < up_size*pow(noLabels,2));
            assert(a + b*noLabels + lz*pow(noLabels,2) < left_size*pow(noLabels,2));
            assert(c + d*noLabels + rz*pow(noLabels,2) < right_size*pow(noLabels,2));
            assert(b + c*noLabels < pow(noLabels,2));
            
            REAL value = repam[a + d*noLabels + (lz+rz)*pow(noLabels,2)];
            value += repam[up_size*pow(noLabels,2) + a + b*noLabels + lz*pow(noLabels,2)];
            value += repam[up_size*pow(noLabels,2) + left_size*pow(noLabels,2) + c + d*noLabels + rz*pow(noLabels,2)];
            value += repam[up_size*pow(noLabels,2) + left_size*pow(noLabels,2) + right_size*pow(noLabels,2) + b + c*noLabels];

            assert(value > -eps);
            bound = std::min(bound,value);
          }
        }
      }
    }
    return bound;
  }

  /*
  template<typename REPAM_ARRAY>
  inline void DiscreteTomographyFactorCounting::NarrowPrimal(const REPAM_ARRAY& repam, const REAL epsilon, typename PrimalSolutionStorage::Element primal)
  {
     const REAL min_val = LowerBound(repam);

     // go over all labelings and set to false all entries that lead to labelings with value greater than min_val + epsilon
     // do it as follows: given some entry in repam, go over all combinations such that it is used, and see whether it always leads to labeling costs larger than min_val + epsilon. If so, set primal in the corresponding place to false

     INDEX sum_up = std::accumulate(repam.begin(),repam.begin()+upSize_,0);
     INDEX sum_left = std::accumulate(repam.begin()+upSize_,repam.begin()+upSize_+leftSize_,0);
     INDEX sum_right = std::accumulate(repam.begin()+upSize_+leftSize_,repam.begin()+upSize_+leftSize_+rightSize_,0);
     INDEX sum_reg = std::accumulate(repam.begin()+upSize_+leftSize_+rightSize_,repam.end(),0);

     // (i) up
     for(INDEX i=0; i<upSize_; ++i) {

     }
     // (ii) left
     // (iii) right
     // (iv) reg
     INDEX z_up_size = upSize_/pow(numberOfLabels_,2);
     INDEX z_left_size = leftSize_/pow(numberOfLabels_,2);
     INDEX z_right_size = rightSize_/pow(numberOfLabels_,2);

     auto op = [&](INDEX i,INDEX j){ return (i+j < z_up_size) ? i+j : z_up_size;  };
     for( INDEX i=0;i<pow(numberOfLabels_,4);i++ ){
        INDEX idx = i;
        INDEX a = idx % numberOfLabels_;
        idx = ( idx - a )/numberOfLabels_;
        INDEX b = idx % numberOfLabels_;
        idx = ( idx - b )/numberOfLabels_;
        INDEX c = idx % numberOfLabels_;
        idx = ( idx - c )/numberOfLabels_;
        INDEX d = idx % numberOfLabels_;

        auto z_up = [&](INDEX k){ return repam[a + d*numberOfLabels_ + k*pow(numberOfLabels_,2)]; };
        auto z_left = [&](INDEX k){ return repam[upSize_ + a + b*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };
        auto z_right = [&](INDEX k){ return repam[upSize_ + leftSize_ + c + d*numberOfLabels_ + k*pow(numberOfLabels_,2)];  };

        REAL reg = repam[upSize_ + leftSize_ + rightSize_ + b + c*numberOfLabels_];
        assert(!std::isnan(reg));
        assert(reg > -std::numeric_limits<REAL>::max() );

        MinConv mc(z_left,z_right,z_left_size,z_right_size,z_up_size);
        mc.CalcConv(op,z_left,z_right);

        REAL m_new = 0;
        for( INDEX j=0;j<z_up_size;j++ ){
           assert(j == op(mc.getIdxA(j),mc.getIdxB(j)));
           assert(z_up(j) > -std::numeric_limits<REAL>::max() );
           assert(!std::isnan(z_up(j)));
	
        m_new = mc.getConv(j)+z_up(j)+reg;
        if(m_new > min_val + epsilon) {
           //primal[] = false;
        }
      }
    }
  }
  */

  inline void DiscreteTomographyFactorCounting::PropagatePrimal(PrimalSolutionStorage::Element primal) const{
    
    struct IdxLbl {
      IdxLbl(INDEX x,INDEX y,INDEX l) : idx(x),c(y),noLabels(l) {
        INDEX i=idx;
        a = i % noLabels;
        i = (i-a) / noLabels;
        b = i % noLabels;
        z = (i-b) / noLabels;
      }
      
      const INDEX idx,c,noLabels;
      INDEX a = 0;
      INDEX b = 0;
      INDEX z = 0;
    };

    auto getOptLabel = [&](INDEX s,INDEX t){
      INDEX noTrue = 0;
      INDEX noUnkwn = 0;
      INDEX opt = 0;
      for(INDEX i=0;i<s;i++){
        if( primal[t + i] == true ){ noTrue++; opt=i; }
        if( primal[t + i] == unknownState ){ noUnkwn++; opt=i; }
      }
      assert(noTrue <=1);
      assert(noTrue != 1 || noUnkwn == 0);
      assert(noTrue == 0 || noUnkwn == 0);

      if( noTrue == 1 ){ return IdxLbl(opt,1,numberOfLabels_); }
      else if( noUnkwn == 1 && noTrue == 0  ){
        primal[t+opt]=true;
        return IdxLbl(opt,1,numberOfLabels_);
      }
      else if( noUnkwn > 1 && noTrue == 0 ){
        return IdxLbl(opt,0,numberOfLabels_);
      }
      else{ 
        assert(false);
        return IdxLbl(opt,0,numberOfLabels_);
      } // not possible!
    };
    
    /* Get primal label */
    auto up    = getOptLabel(upSize_   ,0);
    auto left  = getOptLabel(leftSize_ ,upSize_);
    auto right = getOptLabel(rightSize_,upSize_ + leftSize_);
    auto reg   = getOptLabel(regSize_  ,upSize_ + leftSize_ + rightSize_);

    auto CalcIdx = [&](INDEX a,INDEX b,INDEX z){ return a + b*numberOfLabels_ + z*((INDEX)pow(numberOfLabels_,2)); };
    auto Set2False = [&](INDEX s,INDEX t){ for(INDEX i=0;i<s;i++){primal[t+i]=false;} };
    INDEX count = up.c + left.c + right.c + reg.c;
    
    /* If exact one label is missing, we can calculate it */
    if( count == 3 ){
      if(up.c == 0 && left.b == reg.a && right.a == reg.b
         && CalcIdx(left.a,right.b,left.z+right.z) < upSize_ ){   // up + consistency
        Set2False(upSize_,0);
        primal[CalcIdx(left.a,right.b,left.z+right.z)]=true; }
      if(left.c == 0 && up.b == right.b && right.a == reg.b
         && up.z >= right.z
         && CalcIdx(up.a,reg.a,up.z-right.z) < leftSize_){ // left + consistency
        Set2False(leftSize_,upSize_);
        primal[upSize_ + CalcIdx(up.a,reg.a,up.z-right.z)]=true; }
      if(right.c == 0 && up.a == left.a && left.b == reg.a   // right + consistency
         && CalcIdx(up.b,reg.b,up.z-left.z)
         && up.z >= right.z ){
        Set2False(rightSize_,upSize_ + leftSize_);
        primal[upSize_ + leftSize_ + CalcIdx(up.b,reg.b,up.z-left.z)]=true; }
      if(reg.c == 0 && up.a == left.a && up.b == right.b 
         && up.z == (left.z + right.z) ){                    // pairwise + consistency
        Set2False(regSize_,upSize_ + leftSize_ + rightSize_);
        primal[upSize_ + leftSize_ + rightSize_ + CalcIdx(left.b,right.a,0)]=true; }
    } 
    /* Just for debugging */
    if( count == 4){
      assert(left.a == up.a && up.b == right.b);
      assert(left.b == reg.a && reg.b == right.a);
      assert(up.z == (left.z + right.z));
    }
    
  }

  template<typename REPAM_ARRAY>
  REAL DiscreteTomographyFactorCounting::EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const{
    assert(repam.size() == (upSize_ + leftSize_ + rightSize_ + regSize_));
    REAL val = 0;
    INDEX count = 0;
    INDEX up = 0;
    INDEX left = 0;
    INDEX right = 0;
    INDEX reg = 0;

    auto xa = [&](INDEX idx){ return idx % numberOfLabels_;  };
    auto xb = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return xa(idx); };
    auto z = [&](INDEX idx){ idx = (idx - xa(idx))/numberOfLabels_; return (idx - xa(idx))/numberOfLabels_; };
    
    auto updateVal = [&](INDEX s,INDEX t,INDEX& idx){
      for(INDEX i=0;i<s;i++){
        if(primal[t+i] == true){
          val += repam[t+i];
          idx = i;
          count++;
        }
      }
    };

    INDEX size = 0;      updateVal(upSize_,size,up);
    size += upSize_;     updateVal(leftSize_,size,left);
    size += leftSize_;   updateVal(rightSize_,size,right);
    size += rightSize_;  updateVal(regSize_,size,reg);
    
    assert(count <= 4);

    if( xa(up) == xa(left) &&
        xb(up) == xb(right)&&
        z(left) + z(right) == z(up) &&
        xa(reg) == xb(left) &&
        xb(reg) == xa(right) &&
        count == 4 )
      {  return val; }
    else{
      return std::numeric_limits<REAL>::infinity();
    }
      
  }
    
} // end namespace LP_MP

#endif // LPMP_DTOMOGRAPHY_FACTOR_COUNTING_HXX

