#ifndef LP_MP_partial_external_interface_HXX
#define LP_MP_partial_external_interface_HXX

#include "DD_ILP.hxx"
#include "LP_MP.h"
#include "external_solver_interface.hxx"

namespace LP_MP {

// This class mimics an `LP_MP::LP` but does not inherit from it. This allows
// reusing the very same factors and messages and computing their primal values
// with an external solver.
template<typename EXTERNAL_SOLVER>
class partial_external_solver {
public:
  template<typename FACTOR_CONTAINER_TYPE>
  void add_factor(FACTOR_CONTAINER_TYPE* f) {
    if (!has_factor(f)) {
      dirty_ = true;
      factor_address_to_index_.insert(std::make_pair(f, f_.size()));
      f_.push_back(f);
      assert(factor_address_to_index_.size() == f_.size());

      external_variable_counter_.push_back(s_.get_variable_counters());
      f->construct_constraints(s_);
    }
  }

  template<typename MESSAGE_CONTAINER_TYPE>
  void add_message(MESSAGE_CONTAINER_TYPE* m) {
    if (!has_message(m)) {
      dirty_ = true;
      auto* l = m->GetLeftFactor();
      auto* r = m->GetRightFactor();
      assert(has_factor(l) && has_factor(r));
      auto li = factor_address_to_index_[l];
      auto ri = factor_address_to_index_[r];
      m->construct_constraints(s_, external_variable_counter_[li], external_variable_counter_[ri]);
    }
  }

  template<class LP_TYPE>
  void add_messages(const LP_TYPE &LP) {
    LP.for_each_message([this](auto* m) {
      if (has_factor(m->GetLeftFactor()) && has_factor(m->GetRightFactor()))
        add_message(m);
    });
  }

  bool has_factor(FactorTypeAdapter* f) {
    return factor_address_to_index_.find(f) != factor_address_to_index_.end();
  }

  bool has_message(AbstractMessageContainer* m) {
    return m_.find(m) != m_.end();
  }

  INDEX GetNumberOfFactors() const { return f_.size(); }
  INDEX GetNumberOfMessages() const { return f_.size(); }

  bool solve() {
    bool result = true;

    if (dirty_) {
      s_.init_variable_loading();
      for (auto* f : f_)
        f->load_costs(s_);
      result = s_.solve();

      s_.init_variable_loading();
      for (auto* f : f_)
        f->convert_primal(s_);

      dirty_ = false;
    }

    return result;
  }

  void write_to_file(const std::string& filename) {
    s_.init_variable_loading();
    for (auto* f : f_)
      f->load_costs(s_);
    s_.write_to_file(filename);
  }

  bool dirty () const { return dirty_; }

private:
  DD_ILP::external_solver_interface<EXTERNAL_SOLVER> s_;
  std::vector<FactorTypeAdapter*> f_;
  std::unordered_set<AbstractMessageContainer*> m_;
  std::unordered_map<FactorTypeAdapter*, INDEX> factor_address_to_index_;
  std::vector<typename DD_ILP::variable_counters> external_variable_counter_;
  bool dirty_;
};

} // end namespace LP_MP

#endif // LP_MP_partial_external_interface_HXX

/* vim: set ts=2 sts=2 sw=2 et: */
