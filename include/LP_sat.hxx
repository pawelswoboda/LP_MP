#ifndef LP_SAT_HXX
#define LP_SAT_HXX

#include "LP_MP.h"
#include "sat_solver.hxx"

namespace LP_MP {

template<typename BASE_LP_CLASS>
class LP_sat : public BASE_LP_CLASS
{
private:
  struct sat_th {
    sat_th() 
      : th(5.0)
    {}

    void adjust_th(const bool feasible)
    {
      if(feasible) {
        th *= 0.8;
      } else {
        th *= 2.0;
      } 
    }

    REAL th;
  };


public:
  using LP_type = LP_sat<BASE_LP_CLASS>;

  LP_sat(TCLAP::CmdLine& cmd) 
    : BASE_LP_CLASS(cmd), 
    sat_reduction_mode_arg_("","satReductionMode","how to reduce sat problem",false,"interleaved","{interleaved|static}",cmd)
  {
     //sat_.set_no_simplify(); // seems to make solver much faster
  }

  ~LP_sat()
  {
     // if sat solver is still running, wait until it ends
     //assert(!sat_computation_running());
     if(sat_computation_running()) {
        sat_handle_.wait(); 
        collect_sat_result();
     }
  }

  
  virtual INDEX AddFactor(FactorTypeAdapter* f) 
  {
     // we must wait for the sat solver to finish
     if(sat_computation_running()) {
        sat_handle_.wait();
        collect_sat_result();
     }

     sat_var_.push_back( sat_.size()+1 );
     //sat_var_.push_back(sat_.nVars());
     //std::cout << "number of variables in sat_ = " << sat_var_[sat_var_.size()-1] << "\n";
     f->construct_sat_clauses(sat_);
     INDEX n = BASE_LP_CLASS::AddFactor(f);
     f->init_primal();

     sat_dirty_ = true;

     return n;
  }

   virtual INDEX AddMessage(MessageTypeAdapter* m)
   {
      // we must wait for the sat solver to finish
      if(sat_computation_running()) {
         sat_handle_.wait();
         collect_sat_result();
      }

      assert(this->factor_address_to_index_.find(m->GetLeftFactor()) != this->factor_address_to_index_.end());
      assert(this->factor_address_to_index_.find(m->GetRightFactor()) != this->factor_address_to_index_.end());
      const INDEX left_factor_number = this->factor_address_to_index_[m->GetLeftFactor()];
      const INDEX right_factor_number = this->factor_address_to_index_[m->GetRightFactor()];
      m->construct_sat_clauses(sat_, sat_var_[left_factor_number], sat_var_[right_factor_number]);

      sat_dirty_ = true;

      return BASE_LP_CLASS::AddMessage(m);
   }

   void solve_sat_problem_async(const sat_vec& assumptions, sat_th* th)
   {
      assert(!sat_computation_running());
      sat_handle_ = std::async(std::launch::async, solve_sat_problem_async_impl, &sat_, assumptions, th);
   }

   static bool solve_sat_problem_async_impl(sat_solver* sat_ptr, sat_vec assumptions, sat_th* th)
   {
      const bool feasible = sat_ptr->solve(assumptions.begin(), assumptions.end());

      if(debug()) { std::cout << "sat with threshold " << th->th << ": " << (feasible ? "feasible" : "infeasible") << "\n"; }
      th->adjust_th(feasible);
      if(debug()) { std::cout << "new threshold = " << th->th << "\n"; }

      return feasible;
   }

   bool sat_computation_running() const
   {
      if(!sat_handle_.valid()) { // do not collect sat result in first round
         return false;
      }
      const auto sat_state = sat_handle_.wait_for(std::chrono::seconds(0));
      assert(sat_state != std::future_status::deferred); // this should not happen as we launch primal computation immediately.
      return sat_state != std::future_status::ready;
   }

   void collect_sat_result()
   {
      assert(!sat_computation_running());
      const bool feasible = sat_handle_.get();
      if(!feasible) { return; }

      if(debug()) { 
         std::cout << "collect sat result with threshold = ";
         if(cur_sat_reduction_direction_ == Direction::forward) { 
            std::cout << forward_sat_th_.th;
         } else {
            std::cout << backward_sat_th_.th;
         } 
         std::cout << "\n"; // = " << th.th << "\n";
      }

      if(!sat_dirty_) {

         for(sat_var i=0; i<sat_.size(); ++i) {
            //for(sat_var i=0; i<sat_.nVars(); ++i) {
            //assert(sat_.get_model()[i] == CMSat::l_True || sat_.get_model()[1] == CMSat::l_False);
            //}
         }
         // convert sat solution to original solution format and compute primal cost
         for(INDEX i=0; i<this->f_.size(); ++i) {
            assert(this->factor_address_to_index_[this->f_[i]] == i);
            this->f_[i]->convert_primal(sat_, sat_var_[i]);
         }
         if(verbosity >= 2) {
            const REAL primal_cost = this->EvaluatePrimal();
            std::cout << "sat solution cost = " << primal_cost << "\n"; 
         }
      }
   }



   void ComputeForwardPassAndPrimal(const INDEX iteration)
   {
      const auto omega = this->get_omega();
      if(cur_sat_reduction_direction_ == Direction::forward && !sat_computation_running()) {
         compute_pass_reduce_sat(this->forwardUpdateOrdering_.begin(), this->forwardUpdateOrdering_.end(), omega.forward.begin(), forward_sat_th_);
         cur_sat_reduction_direction_ = Direction::backward; 
      } else {
         this->ComputeForwardPass();
      }
   }
   void ComputeBackwardPassAndPrimal(const INDEX iteration)
   {
      const auto omega = this->get_omega();
      if(cur_sat_reduction_direction_ == Direction::backward && !sat_computation_running()) {
         compute_pass_reduce_sat(this->backwardUpdateOrdering_.begin(), this->backwardUpdateOrdering_.end(), omega.backward.begin(), backward_sat_th_);
         cur_sat_reduction_direction_ = Direction::forward; 
      } else {
         for(auto it = this->backwardUpdateOrdering_.begin(); it != this->backwardUpdateOrdering_.end(); ++it) {
            assert((*it)->no_send_messages() == omega.backward[ std::distance(this->backwardUpdateOrdering_.begin(), it) ].size());
         }
         this->ComputeBackwardPass();
      }
   }
   void ComputePassAndPrimal(const INDEX iteration)
   {
      assert(false); // is currently not used
      ComputeForwardPassAndPrimal(2*iteration+1);
      ComputeBackwardPassAndPrimal(2*iteration+2);
   }


   template<typename FACTOR_ITERATOR, typename WEIGHT_ITERATOR>
   void compute_pass_reduce_sat(FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end, WEIGHT_ITERATOR omega_begin, sat_th& th)
   {
      assert(!sat_computation_running());

      if(sat_handle_.valid()) { // do not collect sat result in first round
        collect_sat_result();
      }

      sat_vec assumptions;
      if(sat_reduction_mode_arg_.getValue() == "interleaved") {
        assumptions = reduce_sat_interleaved(factor_begin, factor_end, omega_begin, th.th);
      } else if(sat_reduction_mode_arg_.getValue() == "static") {
#ifdef LP_MP_PARALLEL
        std::cout << "not implemented yet!\n";
        assert(false);
        //this->ComputePassSynchronized(factor_begin, factor_end, omega_begin);
#else
        this->ComputePass(factor_begin, factor_end, omega_begin);
#endif
        assumptions = reduce_sat_static(th.th);
      }

      // run sat solver on reduced problem asynchronously
      if(!sat_handle_.valid()) { 
         if(verbosity >= 2) { std::cout << "start sat calculation\n"; }
         solve_sat_problem_async(assumptions, &th);
         sat_dirty_ = false;
      } else { 
        if(verbosity >= 2) { std::cout << "restart sat calculation\n"; }
        solve_sat_problem_async(assumptions, &th);
        sat_dirty_ = false;
      }
   }

   template<typename FACTOR_ITERATOR, typename WEIGHT_ITERATOR>
   std::vector<sat_literal> reduce_sat_interleaved(FACTOR_ITERATOR factor_begin, FACTOR_ITERATOR factor_end, WEIGHT_ITERATOR omega_begin, const REAL th)
   {
      sat_vec assumptions;
      for(auto it=factor_begin; it!=factor_end; ++it, ++omega_begin) {
         const INDEX factor_number = this->factor_address_to_index_[*it];
         (*it)->UpdateFactorSAT(*omega_begin, th, sat_var_[factor_number], assumptions);
      }
      return std::move(assumptions); 
   }

   std::vector<sat_literal> reduce_sat_static(const REAL th)
   {
      sat_vec assumptions;
      for(INDEX i=0; i<this->f_.size(); ++i) {
        this->f_[i]->reduce_sat(th, sat_var_[i], assumptions);
      }
      return std::move(assumptions); 
   }
private:

   std::vector<sat_var> sat_var_;
   sat_solver sat_;
   decltype(std::async(std::launch::async, solve_sat_problem_async_impl, nullptr, sat_vec{}, nullptr)) sat_handle_;

   sat_th forward_sat_th_, backward_sat_th_;

   TCLAP::ValueArg<std::string> sat_reduction_mode_arg_;

   Direction cur_sat_reduction_direction_ = Direction::forward;
   // possibly not needed anymore
   bool sat_dirty_ = false; // sat solver is run asynchronously. When factor graph changes, then sat solution cannot be read in anymore and has to be discarded. This flag signifies this case
};


} // end namespace LP_MP

#endif // LP_SAT_HXX
