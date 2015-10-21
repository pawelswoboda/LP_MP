#ifndef MARGINAL_SUMMATION_REPLICATOR_HXX
#define MARGINAL_SUMMATION_REPLICATOR_HXX

#include "problem_decomposition.hxx"

namespace LP_MP {

// left problem is mrf problem with two labels
// right problem is assignment problem. 


template<class FMC, INDEX REPLICATOR_FACTOR_NO, INDEX LEFT_MESSAGE_NO, INDEX RIGHT_MESSAGE_NO, INDEX PROBLEM_CONSTRUCTOR_1_NO, INDEX PROBLEM_CONSTRUCTOR_2_NO>
class MarginalSummationReplicator {
   using ReplicatorFactor = meta::at_c<typename FMC::factor_list, REPLICATOR_FACTOR_NO>;
   
   using LeftMargMessage = meta::at_c<meta::at_c<typename FMC::msg_list, LEFT_MESSAGE_NO>, 0>;
   using LeftFactor = typename LeftMargMessage::LEFT_FACTOR_TYPE;
   using LeftMargMessageOp = typename LeftMargMessage::MessageType;

   using RightMargMessage = meta::at_c<meta::at_c<typename FMC::msg_list, RIGHT_MESSAGE_NO>, 0>;
   using RightFactor = typename RightMargMessage::RIGHT_FACTOR_TYPE;
   using RightMargMessageOp = typename RightMargMessage::MessageType;

   using ProblemConstructor1 = meta::at_c<typename FMC::problem_decomposition, PROBLEM_CONSTRUCTOR_1_NO>;
   using ProblemConstructor2 = meta::at_c<typename FMC::problem_decomposition, PROBLEM_CONSTRUCTOR_2_NO>;

public:
   MarginalSummationReplicator() {}
   void ReadLine(const std::string line) {
      string marg_sum;
      INDEX i1, i2;
      std::stringstream s{line};
      s >> marg_sum >> l >> r;
      if(s.fail()) throw std::runtime_error("incompatible line");
      if(marg_sum != "marg") throw std::runtime_error("line must start with word 'marg'");
      AddMarginalSummation(i1,i2);
   }

   void Construct(ProblemDecomposition<FMC>& pd) {}

   // left_no is left factor number,
   // right_no is right factor number
   // do zrobienia: templatize this for Get{Left|Right}AssignmentFactors
   void AddMarginalSummation(ProblemDecomposition<FMC>* pd, const INDEX left_no, const INDEX right_no) {
      /*
      ReplicatorFactor* replicator_factor = rf_.find[right_no];
      ProblemConstructor1& pc1 = pd->template GetProblemConstructor<ProblemConstructor1>();
      RightFactor* assignment_factor = pc1.GetLeftAssignmentFactors().operator[](right_no);
      if(replicator_factor == rf_.end()) {
         replicator_factor = new ReplicatorFactor(2);
         rf_[right_no] = replicator_factor;
         rmm_.push_back( new RightMargMessage(RightMargMessage(), replicator_factor, assignment_factor,2) ); // do zrobienia: do the same for GetRightAssignment...
      }

      ProblemConstructor2& pc2 = pd->template GetProblemConstructor<ProblemConstructor2>();
      LeftFactor* left_factor = pc2.GetUnaryFactors().operator[](left_no);
      lmm_.push_back( new LeftMargMessage(LeftMargMessage(), left_factor, replicator_factor, 2) );
      */
   }

   void RegisterMarginalSummation(const INDEX i, const INDEX j) 
   {
      if(i >= marg_summation_indices_.size()) {
         marg_summation_indices_.resize(i+1, std::vector<INDEX>{});
      }
      marg_summation_indices_[i].push_back(j);
   } 



   void AddMarginalSummationLeft(LeftFactor* l, RightFactor* r) {}

private:
   std::map<INDEX, ReplicatorFactor*> rf_; // key refers to factor in assignment problem
   std::vector<LeftMargMessage*> lmm_; // left marginal message
   std::vector<RightMargMessage*> rmm_; // right marginal message

   std::vector<std::vector<INDEX> > marg_summation_indices; // if left index i sums with right index j, then there exists c with marg_summation_indices[i][c] = j
};

}

#endif // MARGINAL_SUMMATION_REPLICATOR_HXX
