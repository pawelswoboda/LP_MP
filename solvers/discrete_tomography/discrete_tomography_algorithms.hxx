#ifndef LP_MP_DT_ALGORITHMS_HXX
#define LP_MP_DT_ALGORITHMS_HXX

#include "LP_MP.h"
#include "minConv.hxx"
#include <math.h> 

namespace LP_MP {
  namespace DiscreteTomo {

    static INDEX AlgorithmThreshold = 5000;
    
    // MinConv
    template<typename FACTOR, typename REPAM_ARRAY, typename MSG>
    void MessageCalculation_MinConv_Up(FACTOR* const f, const REPAM_ARRAY& repam, MSG& msg,INDEX noLabels);
    
    template<typename FACTOR, typename REPAM_ARRAY, typename MSG>
    void MessageCalculation_MinConv_Left(FACTOR* const f, const REPAM_ARRAY& repam, MSG& msg,INDEX noLabels);

    template<typename FACTOR, typename REPAM_ARRAY, typename MSG>
    void MessageCalculation_MinConv_Right(FACTOR* const f, const REPAM_ARRAY& repam, MSG& msg,INDEX noLabels);

    template<typename FACTOR, typename REPAM_ARRAY, typename MSG>
    void MessageCalculation_MinConv_Reg(FACTOR* const f, const REPAM_ARRAY& repam, MSG& msg,INDEX noLabels);
    //---

    // Naive
    
    template<typename FACTOR, typename REPAM_ARRAY, typename MSG>
    void MessageCalculation_Naive_Up(FACTOR* const f, const REPAM_ARRAY& repam, MSG& msg,const INDEX noLabels);
    
    template<typename FACTOR, typename REPAM_ARRAY, typename MSG>
    void MessageCalculation_Naive_Left(FACTOR* const f, const REPAM_ARRAY& repam, MSG& msg,INDEX noLabels);

    template<typename FACTOR, typename REPAM_ARRAY, typename MSG>
    void MessageCalculation_Naive_Right(FACTOR* const f, const REPAM_ARRAY& repam, MSG& msg,INDEX noLabels);

    template<typename FACTOR, typename REPAM_ARRAY, typename MSG>
    void MessageCalculation_Naive_Reg(FACTOR* const f, const REPAM_ARRAY& repam, MSG& msg,INDEX noLabels);

    //---
    
    
    
    /* ---- Impl ----  */

    template<typename FACTOR, typename REPAM_ARRAY, typename MSG>
    void MessageCalculation_Naive_Up(FACTOR* const f, const REPAM_ARRAY& repam, MSG& msg,const INDEX noLabels){
      const INDEX up_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::up)/pow(noLabels,2);
      const INDEX right_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::right)/pow(noLabels,2);
      const INDEX left_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::left)/pow(noLabels,2);
      const INDEX reg_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::reg);

      assert((up_size+right_size+left_size)*pow(noLabels,2)+reg_size == repam.size());
       
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

            INDEX u = a + d*noLabels + (lz+rz)*pow(noLabels,2);
            assert(u<msg.size());
            msg[u] = std::min(msg[u],value);
          }
        }
      }
    }

    template<typename FACTOR, typename REPAM_ARRAY, typename MSG>
    void MessageCalculation_Naive_Left(FACTOR* const f, const REPAM_ARRAY& repam, MSG& msg,INDEX noLabels){
      const INDEX up_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::up)/pow(noLabels,2);
      const INDEX right_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::right)/pow(noLabels,2);
      const INDEX left_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::left)/pow(noLabels,2);
      const INDEX reg_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::reg);

      assert((up_size+right_size+left_size)*pow(noLabels,2)+reg_size == repam.size());
 
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

            INDEX l =a + b*noLabels + lz*pow(noLabels,2);
            assert(l < msg.size());
            msg[l] = std::min(msg[l],value);
          }
        }
      }
      
    }

    template<typename FACTOR, typename REPAM_ARRAY, typename MSG>
    void MessageCalculation_Naive_Right(FACTOR* const f, const REPAM_ARRAY& repam, MSG& msg,INDEX noLabels){
      const INDEX up_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::up)/pow(noLabels,2);
      const INDEX right_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::right)/pow(noLabels,2);
      const INDEX left_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::left)/pow(noLabels,2);
      const INDEX reg_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::reg);

      assert((up_size+right_size+left_size)*pow(noLabels,2)+reg_size == repam.size());
      
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

            INDEX r = c + d*noLabels + rz*pow(noLabels,2);
            assert(r < msg.size());
            msg[r] = std::min(msg[r],value);
          }
        }
      }
      
    }

    template<typename FACTOR, typename REPAM_ARRAY, typename MSG>
    void MessageCalculation_Naive_Reg(FACTOR* const f, const REPAM_ARRAY& repam, MSG& msg,INDEX noLabels){
      const INDEX up_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::up)/pow(noLabels,2);
      const INDEX right_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::right)/pow(noLabels,2);
      const INDEX left_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::left)/pow(noLabels,2);
      const INDEX reg_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::reg);

      assert((up_size+right_size+left_size)*pow(noLabels,2)+reg_size == repam.size());

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

            INDEX reg = b + c*noLabels;
            assert(reg < msg.size());
            msg[reg] = std::min(msg[reg],value);
          }
        }
      }
    }
    
    template<typename FACTOR, typename REPAM_ARRAY, typename MSG>
    void MessageCalculation_MinConv_Up(FACTOR* const f, const REPAM_ARRAY& repam, MSG& msg,INDEX noLabels){
      INDEX up_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::up)/pow(noLabels,2);
      INDEX right_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::right)/pow(noLabels,2);
      INDEX left_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::left)/pow(noLabels,2);

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
          assert(k == (mc.getIdxA(k) + mc.getIdxB(k)));
          assert(!std::isnan(reg));

          REAL val = mc.getConv(k) + z_up(k)  + reg;
          INDEX kidx = a + noLabels*d + k*pow(noLabels,2);

          assert(kidx < (up_size*pow(noLabels,2)));
          assert(!std::isnan(val));
          msg[kidx] = std::min(msg[kidx],val);
        }
      }
    }
    
    template<typename FACTOR, typename REPAM_ARRAY, typename MSG>
    void MessageCalculation_MinConv_Left(FACTOR* const f, const REPAM_ARRAY& repam, MSG& msg,INDEX noLabels){

      INDEX left_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::left)/pow(noLabels,2);
      INDEX right_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::right)/pow(noLabels,2);
      INDEX up_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::up)/pow(noLabels,2);
      
      auto op = [&](INDEX i,INDEX j){ // 0 <= i-j < left_size
        if( i < j ){ // i-j < 0
          return left_size;
        }
        else{
          return (i-j < left_size) ? i-j : left_size;
        }
      };
      
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

        MinConv mc(z_up,z_right,up_size,right_size,left_size);
        mc.CalcConv(op,z_up,z_right);

        for(INDEX k=0;k<left_size;k++){
          assert(k == (mc.getIdxA(k) - mc.getIdxB(k)));
          assert(!std::isnan(reg));

          REAL val = mc.getConv(k) + z_left(k) + reg;
          INDEX kidx = a + noLabels*b + k*pow(noLabels,2);

          assert(kidx < (left_size*pow(noLabels,2)));
          assert(!std::isnan(val));

          msg[kidx] = std::min(msg[kidx],val);
        }
      }
    }

    template<typename FACTOR, typename REPAM_ARRAY, typename MSG>
    void MessageCalculation_MinConv_Right(FACTOR* const f, const REPAM_ARRAY& repam, MSG& msg,INDEX noLabels){
      INDEX left_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::left)/pow(noLabels,2);
      INDEX right_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::right)/pow(noLabels,2);
      INDEX up_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::up)/pow(noLabels,2);
     
      auto op = [&](INDEX i,INDEX j){ // 0 <= i-j <= right_size
        if( i < j ){ 
          return right_size;
        }
        else{
          return (i-j < right_size) ? i-j : right_size;
        }
      };

        
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

        MinConv mc(z_up,z_left,up_size,left_size,right_size);
        mc.CalcConv(op,z_up,z_left);

        for(INDEX k=0;k<right_size;k++){
          assert(k == op(mc.getIdxA(k),mc.getIdxB(k)));//(mc.getIdxA(k) - mc.getIdxB(k)));
          assert(!std::isnan(reg));

          REAL val = mc.getConv(k) + z_right(k) + reg;
          INDEX kidx = c + noLabels*d + k*pow(noLabels,2);

          assert(kidx < (right_size*pow(noLabels,2)));
          assert(!std::isnan(val));

          msg[kidx] = std::min(msg[kidx],val);
        }
      }
    }

    template<typename FACTOR, typename REPAM_ARRAY, typename MSG>
    void MessageCalculation_MinConv_Reg(FACTOR* const f, const REPAM_ARRAY& repam, MSG& msg,INDEX noLabels){
      INDEX left_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::left)/pow(noLabels,2);
      INDEX right_size =(*f).getSize(DiscreteTomographyFactorCounting::NODE::right)/pow(noLabels,2);
      INDEX up_size = (*f).getSize(DiscreteTomographyFactorCounting::NODE::up)/pow(noLabels,2);

      auto op = [&](INDEX i,INDEX j){ return (i+j < up_size) ? i+j : up_size; };

      for(INDEX i=0;i<pow(noLabels,2);i++){
        REAL m = std::numeric_limits<REAL>::infinity();
        INDEX b = i % noLabels;
        INDEX c = ((i-b)/noLabels) % noLabels;

        REAL reg = repam[up_size*pow(noLabels,2) + left_size*pow(noLabels,2) + right_size*pow(noLabels,2) + b + c*noLabels];

        for(INDEX j=0;j<pow(noLabels,2);j++){
          INDEX a = j % noLabels;
          INDEX d = ((j-a)/noLabels) % noLabels;

          auto z_up = [&](INDEX k){ return repam[a + d*noLabels + k*pow(noLabels,2)];  };
          auto z_left = [&](INDEX k){ return repam[up_size*pow(noLabels,2) + a + b*noLabels + k*pow(noLabels,2)];  };
          auto z_right = [&](INDEX k){ return repam[up_size*pow(noLabels,2) + left_size*pow(noLabels,2) + c + d*noLabels + k*pow(noLabels,2)];  };

          MinConv mc(z_left,z_right,left_size,right_size,up_size);
          mc.CalcConv(op,z_left,z_right);

          for(INDEX k=0;k<up_size;k++){
            assert( k == (mc.getIdxA(k) + mc.getIdxB(k)));

            REAL val = mc.getConv(k) + z_up(k);
            m = std::min(m,val);
          }
        }
        assert( (m + reg) > -eps);
        msg[i] = m + reg;
      }
    }
    
  }

}

#endif // LP_MP_DT_ALGORITHMS_HXX
