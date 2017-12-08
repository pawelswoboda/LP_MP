#ifndef LP_MP_combiLP_HXX
#define LP_MP_combiLP_HXX

#include "LP_MP.h"
#include "DD_ILP.hxx"
#include <vector>
#include <unordered_set>

namespace LP_MP {

  template<typename EXTERNAL_SOLVER, typename BASE_LP_SOLVER>
  class combiLP : public BASE_LP_SOLVER {
    public:
      using BASE_LP_SOLVER::BASE_LP_SOLVER; 

    void End()
    {
      BASE_LP_SOLVER::End();

      bool consistent = true;
      //std::vector<bool> factor_consistent(this->GetNumberOfFactors()), true; // possibly not needed
      std::vector<FactorTypeAdapter*> ILP_factors;
      std::unordered_map<FactorTypeAdapter*,INDEX> ILP_factor_address_to_index;
      DD_ILP::external_solver_interface<EXTERNAL_SOLVER> s;
      std::vector<typename DD_ILP::variable_counters> external_variable_counter;
      std::unordered_set<MessageTypeAdapter*> ILP_messages;

      auto factor_in_ILP = [&](auto* f) {
        return ILP_factor_address_to_index.find(f) != ILP_factor_address_to_index.end(); 
      };

      auto add_factor_to_ILP = [&](auto* f) {
        consistent=false;
        const INDEX factor_index = this->factor_address_to_index_[f];
        if(!factor_in_ILP(f)) {
          ILP_factor_address_to_index.insert(std::make_pair(f, ILP_factors.size()));
          ILP_factors.push_back(f);
          external_variable_counter.push_back(s.get_variable_counters());
          f->construct_constraints(s); 
        }
      };

      auto add_message_to_ILP = [&](auto* m) {
        if(ILP_messages.find(m) == ILP_messages.end()) {
          ILP_messages.insert(m);
          auto* l = m->GetLeftFactor();
          auto* r = m->GetRightFactor();
          assert(factor_in_ILP(l) && factor_in_ILP(r));
          auto left_factor_index = ILP_factor_address_to_index[l].second;
          auto right_factor_index = ILP_factor_address_to_index[r].second;
          m->construct_constraints(s, external_variable_counter[left_factor_index], external_variable_counter[right_factor_index]);
        } 
      };

      // check wheter factor is locally optimal
      for(INDEX i=0; i<this->GetNumberOfFactors(); ++i) {
        auto* f = this->f_[i];
        f->MaximizePotentialAndComputePrimal(); // should this be called here or shall we assume that this has been done already during computation?
        if(f->LowerBound() <= f->EvaluatePrimal() - eps) {
          add_factor_to_ILP(f);
        }
      }

      // check whether factors agree with each other
      for(INDEX i=0; i<this->GetNumberOfMessages(); ++i) {
        auto* m = this->m_[i];
        if(!m->CheckPrimalConsistency()) {
          add_factor_to_ILP(m->GetLeftFactor());
          add_factor_to_ILP(m->GetRightFactor());
        }
      }

      while(!consistent) {
        // solve subproblem on factors that are not consistent
        s.init_variable_loading();
        for(auto* f : ILP_factors) {
          f->load_costs(s);
        }
        for(INDEX i=0; i<this->GetNumberOfMessages(); ++i) {
          auto* m = this->m_[i];
          auto* l = m->GetLeftFactor();
          auto* r = m->GetRightFactor();
          if(factor_in_ILP(l) && factor_in_ILP(r)) {
            add_message_to_ILP(m); 
          }
        }
        s.solve();
        s.init_variable_loading();
        for(auto* f : ILP_factors) {
          f->convert_primal(s);
        } 

        // check whether solutions agree between ILP and LP part
        for(INDEX i=0; i<this->GetNumberOfMessages(); ++i) {
          auto* m = this->m_[i];
          auto* l = m->GetLeftFactor();
          auto* r = m->GetRightFactor();
          if(factor_in_ILP(l) != factor_in_ILP(r)) {
            if(!m->CheckPrimalConsistency()) {
              add_factor_to_ILP(l);
              add_factor_to_ILP(r);
            }
          }
        }
      }
    }
  };

} // namespace LP_MP

#endif // LP_MP_combiLP_HXX
