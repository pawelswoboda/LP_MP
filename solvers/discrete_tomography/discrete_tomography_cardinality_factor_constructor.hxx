#ifndef LP_MP_DISCRETE_TOMOGRAPHY_CARDINALITY_FACTOR_CONSTRUCTOR_HXX
#define LP_MP_DISCRETE_TOMOGRAPHY_CARDINALITY_FACTOR_CONSTRUCTOR_HXX

#include "LP_MP.h"
#include "tree_decomposition.hxx"
#include <deque>

namespace LP_MP {

template<typename FMC,
   INDEX MRF_PROBLEM_CONSTRUCTOR_NO,
   typename CARDINALITY_FACTOR,
   typename UNARY_CARDINALITY_MESSAGE,
   typename CARDINALITY_MESSAGE>
class discrete_tomography_cardinality_factor_constructor {
public:
   using MrfConstructorType =
      typename meta::at_c<typename FMC::ProblemDecompositionList,MRF_PROBLEM_CONSTRUCTOR_NO>;
   using cardinality_factor_type = CARDINALITY_FACTOR;
   using unary_cardinality_message_type = UNARY_CARDINALITY_MESSAGE;
   using cardinality_message_type = CARDINALITY_MESSAGE;

    template<typename SOLVER>
    discrete_tomography_cardinality_factor_constructor(SOLVER& s) :
       mrf_constructor_(s.template GetProblemConstructor<MRF_PROBLEM_CONSTRUCTOR_NO>()),
       lp_(&s.GetLP()) 
   {}

    void SetNumberOfLabels(const INDEX noLabels) {}

    template<typename ITERATOR>
    cardinality_factor_type* AddProjection(ITERATOR projection_var_begin, ITERATOR projection_var_end, const INDEX sum, LP_tree* tree = nullptr)
    {
       const INDEX no_variables = std::distance(projection_var_begin, projection_var_end);
       assert(no_variables > 1);

       // build cardinality factors recursively: first add bottom layer and connect it with the unaries. Then recursively join two factors until no factor is left.
       std::deque<cardinality_factor_type*> queue;
       for(INDEX i=0; i+1<no_variables; i+=2) {
          const INDEX var = *(projection_var_begin+i);
          const INDEX next_var = *(projection_var_begin+i+1);
          auto* u1 = mrf_constructor_.GetUnaryFactor(var);
          auto* u2 = mrf_constructor_.GetUnaryFactor(next_var);
          auto* f = new cardinality_factor_type(mrf_constructor_.GetNumberOfLabels(var), mrf_constructor_.GetNumberOfLabels(next_var));
          lp_->AddFactor(f);
          queue.push_back(f);
          auto* m1 = new unary_cardinality_message_type(u1, f, Chirality::left);
          lp_->AddMessage(m1);
          auto* m2 = new unary_cardinality_message_type(u2, f, Chirality::right);
          lp_->AddMessage(m2);

          if(tree != nullptr) {
             tree->AddMessage(m1, Chirality::right);
             tree->AddMessage(m2, Chirality::right);
          }
       }

       // if last variable was not civered (#variables is odd) connect last variable and last cardinality factor to each other
       if(no_variables%2 == 1) {
          const INDEX var = *(projection_var_begin + no_variables-1);
          auto* u2 = mrf_constructor_.GetUnaryFactor(var);

          auto* f_left = queue.back();
          queue.pop_back();

          auto* f = new cardinality_factor_type(f_left->GetFactor()->min_conv.size(), mrf_constructor_.GetNumberOfLabels(var));
          f_left->GetFactor()->set_reference_to_min_conv(f->GetFactor()->left);
          lp_->AddFactor(f);
          queue.push_back(f);

          auto* m1 = new cardinality_message_type(f_left,f,Chirality::left);
          lp_->AddMessage(m1);
          auto* m2 = new unary_cardinality_message_type(u2, f, Chirality::right);
          lp_->AddMessage(m2);

          if(tree != nullptr) {
             tree->AddMessage(m1, Chirality::right);
             tree->AddMessage(m2, Chirality::right);
          }
       }

       while(queue.size() > 1) {
          const INDEX q_size = queue.size();
          auto* f_front = queue.front();
          if(q_size%2 == 1) { 
             queue.pop_front(); 
             assert(queue.size() == q_size-1);
          }
          assert(queue.size()%2 == 0);
          for(INDEX i=q_size%2; i<q_size; i+=2) {
             auto* f_left = queue.front();
             queue.pop_front();
             auto* f_right = queue.front();
             queue.pop_front();

             auto* f = new cardinality_factor_type(f_left->GetFactor()->min_conv.size(), f_right->GetFactor()->min_conv.size());
             f_left->GetFactor()->set_reference_to_min_conv(f->GetFactor()->left);
             f_right->GetFactor()->set_reference_to_min_conv(f->GetFactor()->right);
             lp_->AddFactor(f);
             queue.push_back(f);

             auto* m1 = new cardinality_message_type(f_left,f,Chirality::left);
             lp_->AddMessage(m1);
             auto* m2 = new cardinality_message_type(f_right,f,Chirality::right);
             lp_->AddMessage(m2);

             if(tree != nullptr) {
                tree->AddMessage(m1, Chirality::right);
                tree->AddMessage(m2, Chirality::right);
             } 
          }
          if(q_size%2 == 1) { 
             queue.push_front(f_front);
          } 
       }
       assert(queue.size() == 1);
       auto* f_top = queue.back();
       queue.pop_back();

       // possibly do something more effective at top level: only min sum needs to be computed.
       auto& cost_vec = f_top->GetFactor()->min_conv;
       const INDEX sum_top_size = f_top->GetFactor()->left.size() + f_top->GetFactor()->right.size() -1;
       assert(sum+1 <= sum_top_size);
       cost_vec = vector<REAL>(sum_top_size,-std::numeric_limits<REAL>::infinity());
       cost_vec[sum] = 0.0;

       if(tree != nullptr) {
          tree->init();
       }

       return f_top;
    }

  private:
    LP* lp_;
    MrfConstructorType& mrf_constructor_;
};

} // end namespace LP_MP

#endif // LP_MP_DISCRETE_TOMOGRAPHY_CARDINALITY_FACTOR_CONSTRUCTOR_HXX

