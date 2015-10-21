#ifndef POTTS_PROBLEM_CONSTRUCTION_HXX
#define POTTS_PROBLEM_CONSTRUCTION_HXX

#include "problem_decomposition.hxx"

#include <string>
#include <sstream>

namespace LP_MP {

template<class FMC, INDEX UNARY_FACTOR_NO, INDEX PAIRWISE_FACTOR_NO, INDEX LEFT_MESSAGE_NO, INDEX RIGHT_MESSAGE_NO>
class PottsProblemConstructor {
   using UnaryFactor = meta::at_c<typename FMC::factor_list, UNARY_FACTOR_NO>;
   using UnaryFactorOp = typename UnaryFactor::FactorType;
   using PairwiseFactor = meta::at_c<typename FMC::factor_list, PAIRWISE_FACTOR_NO>;
   using LeftMessage = meta::at_c<meta::at_c<typename FMC::msg_list, LEFT_MESSAGE_NO>, 0>;
   using RightMessage = meta::at_c<meta::at_c<typename FMC::msg_list, RIGHT_MESSAGE_NO>, 0>;
public:
   PottsProblemConstructor() {}
   void ReadLine(const std::string& line) {
      std::string potts;
      INDEX i1, i2;
      REAL cost;
      std::stringstream s{line};
      s >> potts >> i1 >> i2 >> cost;
      if(s.fail()) throw std::runtime_error("incompatible line");
      if(potts != "potts") throw std::runtime_error("line must start with word 'potts'");
      SetNumberUnaryFactors(std::max(i1,i2));
      AddPairwiseFactor(std::vector<REAL>{0.0,cost,cost,0.0});
      LinkUnaryPairwiseFactor(pairwiseFactor_.back(), unaryFactor_[i1], unaryFactor_[i2]);
   }

   PottsProblemConstructor& SetNumberUnaryFactors(const INDEX i)
   {
      for(INDEX j=unaryFactor_.size(); j<i; ++j) {
         AddUnaryFactor(std::vector<REAL>(2,0));
      }
      return *this;
   }

   PottsProblemConstructor& AddUnaryFactor(const std::vector<REAL>& cost) { 
      // make this cleaner in multiplex_factor: provide default constructor for simplex type
      static_assert(1<0,"");
      unaryFactor_.push_back( new UnaryFactor(UnaryFactorOp(cost,std::vector<INDEX>{},1), cost) ); 
      return *this;
   }
   PottsProblemConstructor& AddPairwiseFactor(const std::vector<REAL>& cost) { 
      pairwiseFactor_.push_back( new PairwiseFactor(cost) );
      return *this;
   }
   PottsProblemConstructor& LinkUnaryPairwiseFactor(PairwiseFactor* const p, UnaryFactor* const left, UnaryFactor* right) {
      leftMessage_.push_back( new LeftMessage( left, p ) );
      rightMessage_.push_back( new RightMessage( right, p ) );
      return *this;
   }

   const std::vector<UnaryFactor*>& GetUnaryFactors() { return unaryFactor_; }
   const std::vector<PairwiseFactor*>& GetPairwiseFactors() { return pairwiseFactor_; }

   void Construct(ProblemDecomposition<FMC>& pd) {}
   void Construct(LP* lp) {
      for(auto unaryIt : unaryFactor_) {
         lp->AddFactor(unaryIt);
      }
      for(auto pairwiseIt: pairwiseFactor_) {
         lp->AddFactor(pairwiseIt);
      }
      for(auto messageIt : leftMessage_) {
         lp->AddMessage(messageIt);
      }
      for(auto messageIt : rightMessage_) {
         lp->AddMessage(messageIt);
      }
   }

private:
   std::vector<UnaryFactor*> unaryFactor_;
   std::vector<PairwiseFactor*> pairwiseFactor_;

   std::vector<LeftMessage*> leftMessage_;
   std::vector<RightMessage*> rightMessage_;
};

} // end namespace LP_MP

#endif // POTTS_PROBLEM_CONSTRUCTION_HXX

