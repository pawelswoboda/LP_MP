#ifndef LP_MP_DTOMOGRAPHY_COUNTING_NAIVE_HXX
#define LP_MP_DTOMOGRAPHY_COUNTING_NAIVE_HXX

#include "LP_MP.h"

namespace LP_MP{
class DiscreteTomographyFactorCountingNaive { // this factor does not perform any message passing, it is just for building the LP-model
public:
   DiscreteTomographyFactorCountingNaive(const INDEX no_labels, const INDEX no_vars, const INDEX sum)
      : no_labels_(no_labels),
      no_vars_(no_vars),
      sum_(sum)
   {}

   template<typename REPAM_ARRAY>
   static REAL LowerBound(const REPAM_ARRAY& repamPot) {
      return 0.0;
   }

   template<typename REPAM_ARRAY>
   REAL EvaluatePrimal(const REPAM_ARRAY& repam, const PrimalSolutionStorage::Element primal) const
   {
      INDEX sum=0;
      for(INDEX i=0; i<no_vars_; ++i) {
         for(INDEX s=1; s<no_labels_; ++s) {
            sum += s*primal[i*no_labels_ + s];
         }
      }
      if(sum_ == sum) {
         return 0.0;
      } else {
         return std::numeric_limits<REAL>::infinity();
      }
   }
   void CreateConstraints(LpInterfaceAdapter* lp) const;

   INDEX GetNumberOfLabels() const { return no_labels_; }

   INDEX size() const { return no_labels_*no_vars_; }

   struct Primal {
      SIGNED_INDEX residual; // what is the residual of currently labelled unaries regarding the projection constraint
      INDEX no_unaries;
      INDEX no_labelled_unaries;
   };
private:
   const INDEX no_labels_;
   const INDEX no_vars_;
   const INDEX sum_;
};

inline void
DiscreteTomographyFactorCountingNaive::CreateConstraints(LpInterfaceAdapter* lp) const
{
   LinExpr lhs = lp->CreateLinExpr();
   for(INDEX i=0; i<no_vars_; ++i) {
      for(INDEX s=1; s<no_labels_; ++s) {
         lhs += REAL(s)*lp->GetVariable(i*no_labels_ + s);
      }
   }
   LinExpr rhs = lp->CreateLinExpr();
   rhs += (REAL) sum_;
   lp->addLinearEquality(lhs,rhs);
}

class DiscreteTomographyUnaryToFactorCountingNaiveMessage {
public:
   DiscreteTomographyUnaryToFactorCountingNaiveMessage(const INDEX idx) : idx_(idx) {}
   template<class LEFT_FACTOR_TYPE,class RIGHT_FACTOR_TYPE>
   void CreateConstraints(LpInterfaceAdapter* lp,LEFT_FACTOR_TYPE* LeftFactor,RIGHT_FACTOR_TYPE* RightFactor) const;

   template<typename A1, typename A2>
   void RepamLeft(A1& l, const A2& msgs)
   {
   }

   // message is used only for primal rounding
   template<typename RIGHT_FACTOR, typename MSG>
   void ReceiveRestrictedMessageFromRight(const RIGHT_FACTOR& r, PrimalSolutionStorage::Element primal)
   {

   }

private:
   const INDEX idx_; // the position in the factor counting
};

template<class LEFT_FACTOR_TYPE,class RIGHT_FACTOR_TYPE>
void DiscreteTomographyUnaryToFactorCountingNaiveMessage::CreateConstraints(LpInterfaceAdapter* lp,LEFT_FACTOR_TYPE* LeftFactor,RIGHT_FACTOR_TYPE* RightFactor) const
{
   const INDEX no_labels = RightFactor->GetNumberOfLabels();
   for(INDEX i=0; i<no_labels; i++){
      LinExpr lhs = lp->CreateLinExpr();
      LinExpr rhs = lp->CreateLinExpr();
      lhs += lp->GetLeftVariable(i);
      rhs += lp->GetRightVariable(idx_*no_labels + i);
      lp->addLinearEquality(lhs,rhs);
   }
}

template<class FMC,
           INDEX MRF_PROBLEM_CONSTRUCTOR_NO,
           INDEX UNARY_FACTOR_NO,
           INDEX COUNTING_FACTOR_NO,
           INDEX MESSAGE_NO>
class DiscreteTomographyNaiveConstructor {
public:
   using MrfConstructorType =
      typename meta::at_c<typename FMC::ProblemDecompositionList,MRF_PROBLEM_CONSTRUCTOR_NO>;
   using UnaryFactorContainer =
      typename meta::at_c<typename FMC::FactorList,UNARY_FACTOR_NO>;
   using DiscreteTomographyCountingFactorContainer =
      typename meta::at_c<typename FMC::FactorList,COUNTING_FACTOR_NO>;
   using DiscreteTomographyCountingMessage =
      typename meta::at_c<typename FMC::MessageList, MESSAGE_NO>::MessageContainerType;

   DiscreteTomographyNaiveConstructor(Solver<FMC>& s) : s_(s) {}

   void SetNumberOfLabels(const INDEX noLabels) { noLabels_ = noLabels; }

   void AddProjection(const std::vector<INDEX>& projectionVar, const std::vector<REAL>& projectionCost)
   {
      for(INDEX i=0; i<projectionCost.size()-1; ++i) {
         assert(projectionCost[i] == std::numeric_limits<REAL>::infinity());
      }
      assert(projectionCost.back() == 0.0);
      AddProjection(projectionVar,projectionCost.size()-1);
   }

   void AddProjection(const std::vector<INDEX>& projectionVar, const INDEX sum)
   {
      // create new counting factor
      auto f = DiscreteTomographyFactorCountingNaive(noLabels_, projectionVar.size(), sum);
      auto* fc = new DiscreteTomographyCountingFactorContainer(f, std::vector<REAL>(f.size(),0.0));
      s_.GetLP().AddFactor(fc);
      // connect unaries with counting factor
      auto& mrf = s_.template GetProblemConstructor<0>();
      for(INDEX i=0; i<projectionVar.size(); ++i) {
         auto* u = mrf.GetUnaryFactor(projectionVar[i]);
         auto *m = new DiscreteTomographyCountingMessage(DiscreteTomographyUnaryToFactorCountingNaiveMessage(i), u, fc, u->size());
         s_.GetLP().AddMessage(m);
      }
   }
private:
   Solver<FMC>& s_;
   INDEX noLabels_;
};

} // end namespace LP_MP


#endif // LP_MP_DTOMOGRAPHY_COUNTING_NAIVE_HXX
