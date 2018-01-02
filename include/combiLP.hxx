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
          // to do: set initial solution
        }
      };

      auto add_message_to_ILP = [&](auto* m) {
        if(ILP_messages.find(m) == ILP_messages.end()) {
          consistent = false;
          ILP_messages.insert(m);
          auto* l = m->GetLeftFactor();
          auto* r = m->GetRightFactor();
          assert(factor_in_ILP(l) && factor_in_ILP(r));
          auto left_factor_index = ILP_factor_address_to_index[l];
          auto right_factor_index = ILP_factor_address_to_index[r];
          m->construct_constraints(s, external_variable_counter[left_factor_index], external_variable_counter[right_factor_index]);
        } 
      };

      // check whether factor is locally optimal
      for(INDEX i=0; i<this->GetNumberOfFactors(); ++i) {
        auto* f = this->f_[i];
        //f->MaximizePotentialAndComputePrimal(); // should this be called here or shall we assume that this has been done already during computation?
        assert(f->LowerBound() <= f->EvaluatePrimal() + eps);
        if(f->LowerBound() < f->EvaluatePrimal() - eps) {
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

      INDEX combiLP_iteration = 0;
      //std::ofstream ofs("solution_iteration.txt", std::ofstream::out);
      while(!consistent) {
        //for(INDEX i=0; i<100; ++i) {
          //if(factor_in_ILP(this->f_[i])) ofs << "+";
          //else ofs << "-";
          //ofs << static_cast<FMC_SRMP::UnaryFactor*>(this->f_[i])->GetFactor()->primal() << " ";
        //}
        //ofs << std::endl;

        consistent = true;
        std::cout << "solve ILP subproblem with " << ILP_factors.size() << " factors (out of " << this->f_.size() << ")\n";
        //for(auto* f : ILP_factors) {
        //  std::cout << this->factor_address_to_index_[f] << " ";
        //}
        //std::cout << "\n";
        // solve subproblem on factors that are not consistent
        s.init_variable_loading();
        for(auto* f : ILP_factors) {
          f->load_costs(s);
        }
        for(auto* m : this->m_) {
          auto* l = m->GetLeftFactor();
          auto* r = m->GetRightFactor();
          if(factor_in_ILP(l) && factor_in_ILP(r)) {
            add_message_to_ILP(m);
          } else {
            assert(m->CheckPrimalConsistency() == true);
        }
        s.write_to_file("combiLP_ILP_part_iteration" + std::to_string(combiLP_iteration) + ".lp");
        const bool solved = s.solve();
        assert(solved);
        s.init_variable_loading();
        for(auto* f : ILP_factors) {
          f->convert_primal(s);
        } 
        for(auto* m : ILP_messages) {
          assert(m->CheckPrimalConsistency() == true);
        }

        // check whether solutions agree between ILP and LP part
        for(auto* m : this->m_) {
          auto* l = m->GetLeftFactor();
          auto* r = m->GetRightFactor();
          if(factor_in_ILP(l) != factor_in_ILP(r)) {
            if(!m->CheckPrimalConsistency()) {
              add_factor_to_ILP(l);
              add_factor_to_ILP(r);
            }
          }
        }

        for(auto* f : ILP_factors) {
          f->propagate_primal_through_messages(); 
        } 
        std::cout << "primal cost obtained by combiLP = " << this->EvaluatePrimal() << "\n";
        ++combiLP_iteration;
      }
    }
  };
} // namespace LP_MP

#endif // LP_MP_combiLP_HXX
