#ifndef LP_MP_LP_EXTERNAL_INTERFACE_HXX
#define LP_MP_LP_EXTERNAL_INTERFACE_HXX

#include "DD_ILP.hxx"
#include "LP_MP.h"
#include "external_solver_interface.hxx"

// interface to DD_ILP object which builds up the given LP_MP problem for various other solvers.
// there are two functions a factor must provide, so that the export can take place:
//   template<typename EXTERNAL_SOLVER>
//   void construct_constraints(EXTERNAL_SOLVER& s, TYPES...) const
// and
//   template<typename EXTERNAL_SOLVER>
//   void convert_to_primal(EXTERNAL_SOLVER& s, TYPES...)
// where TYPES is a list of parameters of type EXTERNAL_SOLVER::{variable|vector|matrix|tensor} whose number and type matches the respective arguments given to export_variables() 
// likewise, a message must also provide a function
// template<typename EXTERNAL_SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
// void construct_constraints(EXTERNAL_SOLVER&s, LEFT_FACTOR& left, LEFT_TYPES..., RIGHT_FACTOR& right, RIGHT_TYPES...)

namespace LP_MP {

template<typename BASE_EXTERNAL_SOLVER, typename BASE_LP_SOLVER>
class LP_external_solver : public BASE_LP_SOLVER {
public:
  using external_solver = DD_ILP::external_solver_interface<BASE_EXTERNAL_SOLVER>;

  using BASE_LP_SOLVER::BASE_LP_SOLVER;

  virtual INDEX AddFactor(FactorTypeAdapter* f)
  {
    const INDEX factor_no = BASE_LP_SOLVER::AddFactor(f);
    external_variable_counter_.push_back(s_.get_variable_counters());
    f->construct_constraints(s_);
    return factor_no;
  }

  virtual INDEX AddMessage(MessageTypeAdapter* m)
  {
    const INDEX message_no = BASE_LP_SOLVER::AddMessage(m); 

    const INDEX left_factor_no = this->factor_address_to_index_[m->GetLeftFactorTypeAdapter()];
    assert(left_factor_no < this->GetNumberOfFactors() && left_factor_no < external_variable_counter_.size());

    const INDEX right_factor_no = this->factor_address_to_index_[m->GetRightFactorTypeAdapter()];
    assert(right_factor_no < this->GetNumberOfFactors() && right_factor_no < external_variable_counter_.size());

    m->construct_constraints(s_, external_variable_counter_[left_factor_no], external_variable_counter_[right_factor_no]);

    return message_no;
  }

  const external_solver& get_external_solver() const { return s_; }

  void solve()
  {
    // go over all factors and write costs to external problem
    s_.init_variable_loading();
    for(INDEX i=0; i<this->GetNumberOfFactors(); ++i) {
      this->f_[i]->load_costs(s_);
    } 

    s_.solve();
  }

  void write_to_file(const std::string& filename) { 
    s_.init_variable_loading();
    for(INDEX i=0; i<this->GetNumberOfFactors(); ++i) {
      this->f_[i]->load_costs(s_);
    } 
    s_.write_to_file(filename); 
  }

private:
   external_solver s_; 

   std::vector<typename DD_ILP::variable_counters> external_variable_counter_;
};


} // end namespace LP_MP

#endif // LP_MP_LP_EXTERNAL_INTERFACE_HXX
