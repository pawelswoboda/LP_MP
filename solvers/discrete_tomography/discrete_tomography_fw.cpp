#include "discrete_tomography.h"
#include "visitors/standard_visitor.hxx"
#include "SVM-NEW/SVM.h"

#include <chrono>

using namespace LP_MP;

static LP_tree mrf_chain;
static LP_tree dt_chain1;
static LP_tree dt_chain2;
static INDEX no_labels;

int kappa_SCALE = 1;
int mu_SCALE = 1;

double lambda = 1.0;
double mu = 1.0;
double kappa = 1.0;



// YPtr is an array of ints holding the labeling of unaries
// TermData holds an LP_tree with discrete tomography subproblems, the unary indices
struct sub_problem {
   template<typename VEC>
   sub_problem(LP_tree t, VEC& u) : dt_chain(t), no_unaries(u.size() + 1) 
   {
      assert(no_unaries > 1);
      unaries = new INDEX[no_unaries]; 
      pairwise_factors = new FMC_DT::dt_sequential_pairwise_factor*[no_unaries-1];
      for(INDEX i=0; i<no_unaries-1; ++i) {
         pairwise_factors[i] = u[i];
      }
   }
   LP_tree dt_chain;
   INDEX no_unaries;
   INDEX* unaries;
   FMC_DT::dt_sequential_pairwise_factor** pairwise_factors; // one less than number of unaries
   INDEX chain_no;
};

double max_fn(double* wi, YPtr _y, double kappa, TermData term_data) // maximization oracle. Must copy argmax_y <a^{iy},[PAD(wi) kappa]> to y, and return the free term a^{iy}[d].
{
   sub_problem* sp = (sub_problem*) term_data;
   const REAL sign = sp->chain_no == 0 ? 1.0 : -1.0;
   std::cout << "\ncompute subgradient on " << sp << "\n";
   for(INDEX i=0; i<9; ++i) {
      std::cout << sign*wi[i] << ", ";
   }
   std::cout << "\n";
   // first move the unary costs to the pairwise factors
   INDEX c=0;
   for(INDEX i=0; i<sp->no_unaries-1; ++i) {
      auto* f = sp->pairwise_factors[i]->GetFactor();
      for(INDEX x1=0; x1<no_labels; ++x1) {
         for(INDEX x2=0; x2<no_labels; ++x2) {
            assert(!std::isnan(wi[c]));
            (*f).reg(x1,x2) += sign*wi[c];
         } 
         ++c;
      }
   }
   // last pairwise is reparametrized twice
   {
      auto* f = sp->pairwise_factors[ sp->no_unaries-2 ]->GetFactor();
      for(INDEX x2=0; x2<no_labels; ++x2) {
         for(INDEX x1=0; x1<no_labels; ++x1) {
            assert(!std::isnan(wi[c]));
            (*f).reg(x1,x2) += sign*wi[c];
         } 
         ++c;
      } 
   }
   assert(c == 9);

   sp->dt_chain.compute_subgradient();
   std::cout << "lb = " << sp->dt_chain.lower_bound() << ", kappa = " << kappa << "\n";
   INDEX* y = (INDEX*) _y;
   for(INDEX i=0; i<sp->no_unaries-1; ++i) {
      y[i] = sp->pairwise_factors[i]->GetFactor()->state_[0]; 
   }
   y[sp->no_unaries-1] = sp->pairwise_factors[sp->no_unaries-2]->GetFactor()->state_[1];

   // now remove the Lagrangean variables from unary factors again
   c=0;
   for(INDEX i=0; i<sp->no_unaries-1; ++i) {
      auto* f = sp->pairwise_factors[i]->GetFactor();
      for(INDEX x1=0; x1<no_labels; ++x1) {
         for(INDEX x2=0; x2<no_labels; ++x2) {
            assert(!std::isnan(wi[c]));
            (*f).reg(x1,x2) -= sign*wi[c];
         } 
         ++c;
      }
   }
   // last pairwise is reparametrized twice
   {
      auto* f = sp->pairwise_factors[ sp->no_unaries-2 ]->GetFactor();
      for(INDEX x2=0; x2<no_labels; ++x2) {
         for(INDEX x1=0; x1<no_labels; ++x1) {
            assert(!std::isnan(wi[c]));
            (*f).reg(x1,x2) -= sign*wi[c];
         } 
         ++c;
      } 
   }
   // free term is excluding wi
   const REAL cost = - sp->dt_chain.primal_cost()*kappa; // we minimize, but SVM expects maximization
   std::cout << "cost without weights = " << -cost << ", ";
   std::cout << "sol = ";
   for(INDEX i=0; i<3; ++i) {
      std::cout << y[i] << ", ";
   }
   std::cout << "\n";
   return cost;
}

static bool compare_fn(YPtr _y1, YPtr _y2, TermData term_data)
{
   sub_problem* sp = (sub_problem*) term_data;
   INDEX* y1 = (INDEX*) _y1;
   INDEX* y2 = (INDEX*) _y2;
   for(INDEX i=0; i<sp->no_unaries; ++i) {
      if(y1[i] != y2[i]) {
         return false;
      } 
   }
   return true;
}

static void copy_fn(double* ai, YPtr _y, TermData term_data)
{
   sub_problem* sp = (sub_problem*) term_data;
   const REAL sign = sp->chain_no == 0 ? 1.0 : -1.0;
   std::cout << "\ncopy fn on " << sp << "\n";
   std::fill(ai, ai+9, 0.0);
   INDEX* y = (INDEX*) _y;
   //for(INDEX i=0; i<sp->no_unaries-1; ++i) {
   //   sp->pairwise_factors[i]->GetFactor()->state_[0] = y[i];
   //   sp->pairwise_factors[i]->GetFactor()->state_[1] = y[i+1];
   //}
   for(INDEX i=0; i<sp->no_unaries; ++i) {
      INDEX label = y[i];
      std::cout << label << ", ";
      ai[i*no_labels + label] = sign*1.0;
   }
   std::cout << "\n";

   // free term
   //ai[9] = -sp->dt_chain.primal_cost(); 
}

static double dot_product_fn(double* wi, YPtr _y, TermData term_data)
{
   sub_problem* s = (sub_problem*) term_data;
   const REAL sign = s->chain_no == 0 ? 1.0 : -1.0;
   //std::cout << "dot product on " << s << "\n";
   INDEX* y = (INDEX*) _y;
   double v = 0;
   INDEX c=0;
   for(INDEX i=0; i<s->no_unaries; ++i) {
      v += sign*wi[i*no_labels + y[i]];
   }

   //for(INDEX i=0; i<s->no_unaries-1; ++i) {
   //   s->pairwise_factors[i]->GetFactor()->state_[0] = y[i];
   //   s->pairwise_factors[i]->GetFactor()->state_[1] = y[i+1];
   //}
   //v += -s->dt_chain.primal_cost();
   return v;
}


int main()
{


   Solver<FMC_DT,LP,StandardVisitor> solver;
   auto& mrf_constructor = solver.GetProblemConstructor<0>();
   mrf_constructor.AddUnaryFactor({0,0,0});
   mrf_constructor.AddUnaryFactor({0,0,0});
   mrf_constructor.AddUnaryFactor({0,0,0});
   mrf_constructor.AddEmptyPairwiseFactor(0,1);
   mrf_constructor.AddEmptyPairwiseFactor(1,2);
   auto* f_01_c = mrf_constructor.GetPairwiseFactor(0,1);
   auto* f_12_c = mrf_constructor.GetPairwiseFactor(1,2);
   auto* f_01 = f_01_c->GetFactor();
   auto* f_12 = f_12_c->GetFactor();
   for(INDEX i=0; i<2; ++i) {
      for(INDEX j=0; j<2; ++j) {
         if(i != j) {
            (*f_01).cost(i,j) = 1.0;
            (*f_12).cost(i,j) = 1.0;
         } else {
            assert((*f_01)(i,j) == 0.0);
            assert((*f_12)(i,j) == 0.0);
         }
      }
   }
   mrf_chain = mrf_constructor.add_tree({mrf_constructor.GetUnaryFactor(0), mrf_constructor.GetUnaryFactor(1), mrf_constructor.GetUnaryFactor(2)}, {f_01_c, f_12_c});
   
   auto& dt_constructor = solver.GetProblemConstructor<1>();
   dt_constructor.SetNumberOfLabels(3);
   no_labels = 3;
   dt_constructor.AddProjection({0,1,2}, {5.0,5.0,5.0,0.5}, &dt_chain1);
   // get pairwise factors out of dt_chain 
   auto pairwise1 = dt_chain1.get_factors<FMC_DT::dt_sequential_pairwise_factor>();
   //std::reverse(pairwise1.begin(), pairwise1.end());
   dt_constructor.AddProjection({0,1,2}, {10.0,2.0,3.0,4.0}, &dt_chain2);
   auto pairwise2 = dt_chain2.get_factors<FMC_DT::dt_sequential_pairwise_factor>();
   //std::reverse(pairwise2.begin(), pairwise2.end());

   

	SVM s (3*no_labels, 2, max_fn, copy_fn, compare_fn, dot_product_fn, nullptr, false);//int d, int n, MaxFn max_fn, CopyFn copy_fn, CompareFn compare_fn, DotProductFn dot_product_fn, DotProductKernelFn dot_product_kernel_fn, bool zero_lower_bound);
	s.SetParams(lambda, mu / mu_SCALE, kappa / kappa_SCALE);

   sub_problem sp1(dt_chain1, pairwise1);
   sp1.chain_no = 0;
   sub_problem sp2(dt_chain2, pairwise2);
   sp2.chain_no = 1;

   // to do: push all information from pairwise potentials into discrete tomography chains.

   s.SetTerm(0, &sp1, 3*no_labels, 3*sizeof(INDEX), nullptr );
   s.SetTerm(1, &sp2, 3*no_labels, 3*sizeof(INDEX), nullptr );

	s.options.gap_threshold = 0.001;
	s.options.iter_max = 100;

   s.Solve();

   return 0;
}

