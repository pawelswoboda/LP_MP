#ifndef LP_MP_LP_BCFW_BUNDLE_HXX
#define LP_MP_LP_BCFW_BUNDLE_HXX

#include "tree_decomposition.hxx"
#include "SVM.h"

namespace LP_MP {

class LP_tree_BCFW_Bundle : public LP_tree {
public:
   // for the Frank Wolfe implementation
   // to do: change the SVM implementation and make these methods virtual instead of static.
   // to do: Lagrangean factors are shared between trees. Hence, one must in each tree use these factors with costs divided by number of appearances in different trees.
   // _y is the primal labeling to be computed
   // wi is the Lagrangean variables
   static double max_fn(double* wi, BCFW_Bundle::YPtr _y, BCFW_Bundle::TermData term_data)
   {
      LP_tree_BCFW_Bundle* t = (LP_tree_BCFW_Bundle*) term_data;

      // first add weights to problem
      // we only need to add Lagrange variables to Lagrangean_factors_ (others are not shared)
      t->add_weights(wi, -1.0);

      // compute optimal labeling
      t->compute_subgradient();

      // remove weights again
      t->add_weights(wi, +1.0);

      // store primal solution in archive
      void* y = (void*) _y;
      serialization_archive ar(y, t->primal_size_in_bytes());
      save_archive s_ar(ar);
      for(auto f : t->factors_) {
         f->serialize_primal(s_ar);
      }
      ar.release_memory(); // so that memory will not be automatically deallocated

      return -t->primal_cost(); // SVN considers maximization problem!
   }

   static bool compare_fn(BCFW_Bundle::YPtr _y1, BCFW_Bundle::YPtr _y2, BCFW_Bundle::TermData term_data)
   {
      // the primal is a binary archive constructed by serializing the primal solutions
      LP_tree_BCFW_Bundle* t = (LP_tree_BCFW_Bundle*) term_data;
      const INDEX size = t->primal_size_in_bytes();
      return std::memcmp((void*) _y1, (void*) _y2, size) == 0;
   }

   // copy values provided by subgradient from all factors into ai
   static void copy_fn(double* ai, BCFW_Bundle::YPtr _y, BCFW_Bundle::TermData term_data)
   {
      LP_tree_BCFW_Bundle* t = (LP_tree_BCFW_Bundle*) term_data;

      std::fill(ai, ai+t->dual_size(), double(0.0));

      // read in primal solution from which to compute subgradient
      char* mem = (char*) _y;
      serialization_archive ar(_y, t->primal_size_in_bytes());
      load_archive l_ar(ar);
      for(auto f : t->factors_) { // theoretically, we would only need to load primals of Lagrangean factors!
         f->serialize_primal(l_ar);
      } 
      ar.release_memory();

      for(auto L : t->Lagrangean_factors_) {
         L.copy_fn(ai);
      }
   }

   static double dot_product_fn(double* wi, BCFW_Bundle::YPtr _y, BCFW_Bundle::TermData term_data)
   {
      LP_tree_BCFW_Bundle* t = (LP_tree_BCFW_Bundle*) term_data;

      // read in primal solution
      // to do: only primal solution associated with Lagrangean factors needs to be read in
      char* mem = (char*) _y;
      serialization_archive ar(_y, t->primal_size_in_bytes());
      load_archive l_ar(ar);
      for(auto f : t->factors_) {
         f->serialize_primal(l_ar);
      }

      double v = 0.0;
      for(auto L : t->Lagrangean_factors_) {
         v += L.dot_product_fn(wi);
      }

      ar.release_memory();

      return v;
   }
};

// solve problem with proximal bundle with proximal steps implemented with a multi-plane block coordinate Frank-Wolfe method
class LP_BCFW_Bundle : public LP_with_trees {
private:
   typename BCFW_Bundle::SVM* build_up_solver(const REAL lambda = 0.1)
   {
      auto* svm = new BCFW_Bundle::SVM(this->no_Lagrangean_vars(), trees_.size(), LP_tree_BCFW_Bundle::max_fn, LP_tree_BCFW_Bundle::copy_fn, LP_tree_BCFW_Bundle::compare_fn, LP_tree_BCFW_Bundle::dot_product_fn, nullptr, false);//int d, int n, MaxFn max_fn, CopyFn copy_fn, CompareFn compare_fn, DotProductFn dot_product_fn, DotProductKernelFn dot_product_kernel_fn, bool zero_lower_bound);
      //svm->SetParams(lambda, 1.0, 1.0); // lambda, mu, kappa

      for(INDEX i=0; i<trees_.size(); ++i) {
         auto& t = trees_[i];
         const INDEX primal_size_in_bytes = t.primal_size_in_bytes();

         svm->SetTerm(i, &t, t.mapping().size()-1, t.primal_size_in_bytes(), &t.mapping()[0]); // although mapping is of length di + 1 (the last entry being di itself, its length must be given as di!
      }

      svm->options.gap_threshold = 0.0001;
      svm->options.iter_max = 100000;

      return svm;
   }


public:
   using LP_with_trees::LP_with_trees;

   void ComputePass(const INDEX iteration)
   {
      std::cout << "compute pass fw\n";
      // compute descent with quadratic term centered at current reparametrization.
      // unfortunately, the SVM solver has to be built up from scratch in every iteration

      auto* svm = build_up_solver();
      const REAL lb = this->LowerBound();
      //SVM_FW_visitor visitor({0.0,lb});
      //auto visitor_func = std::bind(&SVM_FW_visitor::visit, &visitor, std::placeholders::_1);
      //svm->options.callback_fn = visitor_func;
      double* const w = svm->Solve();
      add_weights(w, -1.0);
      delete svm;
      std::cout << "after lower bound = " << this->LowerBound() << "\n";
   }
};


} // namespace LP_MP

#endif // LP_MP_LP_BCFW_BUNDLE_HXX

