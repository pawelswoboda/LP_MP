#ifndef LP_MP_factor_archive_HXX
#define LP_MP_factor_archive_HXX

#include <unordered_map>

namespace LP_MP {

namespace serialization_functor {

struct dual {
  template<typename ARCHIVE>
  void operator()(FactorTypeAdapter* f, ARCHIVE &a) {
    f->serialize_dual(a);
  }
};

struct primal {
  template<typename ARCHIVE>
  void operator()(FactorTypeAdapter* f, ARCHIVE &a) {
    f->serialize_primal(a);
  }
};

} // namespace serialization_functor

template<typename SERIALIZATON_FUNCTOR>
class factor_archive {
public:
  using factor_archive_type = factor_archive<SERIALIZATON_FUNCTOR>;

  factor_archive() { }

  template<typename FACTOR_ITERATOR>
  factor_archive(FACTOR_ITERATOR begin, FACTOR_ITERATOR end)
  {
    allocate_archive aa;
    for (auto it = begin; it != end; ++it) {
      factor_to_index_.insert(std::make_pair(*it, aa.size()));
      functor(*it, aa);
    }
    archive_.aquire_memory(aa.size());

    save_archive sa(archive_);
    for (auto it = begin; it != end; ++it) {
      functor(*it, sa);
    }
  }

  factor_archive(LP &lp)
  : factor_archive(lp.begin(), lp.end()) { }

  void load_factor(FactorTypeAdapter* f) {
    access<load_archive>(f);
  }

  void save_factor(FactorTypeAdapter* f) {
    access<save_archive>(f);
  }

  bool operator==(const factor_archive_type& rhs) const {
    return archive_ == rhs.archive_;
  }

  static bool check_factor_equality(factor_archive_type& fa1, factor_archive_type& fa2, FactorTypeAdapter *f) {
    if (fa1.factor_to_index_.find(f) == fa1.factor_to_index_.end() ||
        fa2.factor_to_index_.find(f) == fa2.factor_to_index_.end())
    {
      return false;
    }

    SERIALIZATON_FUNCTOR functor;
    allocate_archive aa;
    functor(f, aa);

    fa1.archive_.reset_cur();
    fa1.archive_.advance(fa1.factor_to_index_[f]);

    fa2.archive_.reset_cur();
    fa2.archive_.advance(fa1.factor_to_index_[f]);

    return std::memcmp(fa1.archive_.cur_address(), fa2.archive_.cur_address(), aa.size()) == 0;
  }

private:
  SERIALIZATON_FUNCTOR functor;
  serialization_archive archive_;
  std::unordered_map<FactorTypeAdapter*, INDEX> factor_to_index_;

  template<typename ARCHIVE>
  void access(FactorTypeAdapter* f) {
    assert(factor_to_index_.find(f) != factor_to_index_.end());
    ARCHIVE a(archive_);
    archive_.advance(factor_to_index_[f]);
    functor(f, a);
  }
};

} // namespace LP_MP

#endif // LP_MP_factor_archive_HXX

/* vim: set ts=2 sts=2 sw=2 et: */
