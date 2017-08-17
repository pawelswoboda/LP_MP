#include "tree_decomposition.hxx"
#include <iostream>
#include <array>
#include "assert.h"
#include "SVM-NEW/SVM.h"

#include <chrono>


int kappa_SCALE = 1;
int mu_SCALE = 1;

double lambda = 1.0;
double mu = 1.0;
double kappa = 1.0;

/*

struct sub_problem_test {
   std::array<double,2> cost; // simplex with two entries.
   double sign;
};

double max_fn_test(double* wi, YPtr _y, double kappa, TermData term_data) // maximization oracle. Must copy argmax_y <a^{iy},[PAD(wi) kappa]> to y, and return the free term a^{iy}[d].
{
   //assert(kappa == 1.0);
   sub_problem_test* sp = (sub_problem_test*) term_data;
   const double sign = sp->sign;
   assert(sign == 1.0 || sign == -1.0);
   size_t* y = (size_t*) _y;

   if(sign*wi[0] + kappa*sp->cost[0] > sign*wi[1] + kappa*sp->cost[1]) {
      *y = 0;
   } else {
      *y = 1;
   }
   std::cout << "compute solution on " << sp << " = " << *y << " with w = (" << wi[0] << "," << wi[1] << ")\n";
   return kappa*sp->cost[*y]; 
}

static bool compare_fn_test(YPtr _y1, YPtr _y2, TermData term_data)
{
   sub_problem_test* sp = (sub_problem_test*) term_data;
   size_t* y1 = (size_t*) _y1;
   assert(*y1 == 0 || *y1 == 1);
   size_t* y2 = (size_t*) _y2;
   assert(*y2 == 0 || *y2 == 1);
   return *y1 != *y2;
}

static void copy_fn_test(double* ai, YPtr _y, TermData term_data)
{
   sub_problem_test* sp = (sub_problem_test*) term_data;
   const double sign = sp->sign;
   assert(sign == 1.0 || sign == -1.0);
   size_t* y = (size_t*) _y;
   assert(*y == 0 || *y == 1);
   ai[1-*y] = 0.0;
   ai[*y] = sign*1.0;
}

static double dot_product_fn_test(double* wi, YPtr _y, TermData term_data)
{
   sub_problem_test* sp = (sub_problem_test*) term_data;
   const double sign = sp->sign;
   assert(sign == 1.0 || sign == -1.0);
   size_t* y = (size_t*) _y;
   assert(*y == 0 || *y == 1);
   return sign*wi[*y];
}

int main()
{ 
	SVM s (2, 2, max_fn_test, copy_fn_test, compare_fn_test, dot_product_fn_test, nullptr, true);//int d, int n, MaxFn max_fn, CopyFn copy_fn, CompareFn compare_fn, DotProductFn dot_product_fn, DotProductKernelFn dot_product_kernel_fn, bool zero_lower_bound);
	s.SetParams(1.0, 1.0, 1.0);

   sub_problem_test sp1;
   sp1.sign = 1.0;
   sp1.cost = {9.0,10.0};
   sub_problem_test sp2;
   sp2.sign = -1.0;
   sp2.cost = {20.0,0.0};

   s.SetTerm(0, &sp1, 2, sizeof(size_t), nullptr );
   s.SetTerm(1, &sp2, 2, sizeof(size_t), nullptr );

	s.options.gap_threshold = 0.001;
	s.options.iter_max = 10000;

   double* w = s.Solve();

   std::cout << "\n\nsolution = ";
   for(int i=0; i<2; ++i) {
      std::cout << w[i] << ",";
   }
   std::cout << "\n";

   return 0; 
}
*/

/*
static LP_tree mrf_chain;
static LP_tree dt_chain1;
static LP_tree dt_chain2;
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
            (*f_01).cost(i,j) = 0.0;
            (*f_12).cost(i,j) = 0.0;
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
   const REAL scaling = 0.01;
   dt_constructor.AddProjection({0,1,2}, {scaling*5.0, scaling*5.0, scaling*5.0, scaling*0.5}, &dt_chain1);
   // get pairwise factors out of dt_chain 
   auto pairwise1 = dt_chain1.get_factors<FMC_DT::dt_sequential_pairwise_factor>();
   //std::reverse(pairwise1.begin(), pairwise1.end());
   dt_constructor.AddProjection({0,1,2}, {scaling*10.0, scaling*2.0, scaling*3.0, scaling*4.0}, &dt_chain2);
   auto pairwise2 = dt_chain2.get_factors<FMC_DT::dt_sequential_pairwise_factor>();
   //std::reverse(pairwise2.begin(), pairwise2.end());

   

	SVM s (3*no_labels, 2, max_fn, copy_fn, compare_fn, dot_product_fn, nullptr, false);//int d, int n, MaxFn max_fn, CopyFn copy_fn, CompareFn compare_fn, DotProductFn dot_product_fn, DotProductKernelFn dot_product_kernel_fn, bool zero_lower_bound);
	s.SetParams(lambda, mu / mu_SCALE, kappa / kappa_SCALE);
   assert(mu/mu_SCALE == 1.0);

   sub_problem sp1(dt_chain1, pairwise1);
   sp1.sign = 1.0;
   sub_problem sp2(dt_chain2, pairwise2);
   sp2.sign = -1.0;

   sp1.dt_chain.compute_subgradient();
   sp2.dt_chain.compute_subgradient();
   const double initial_lb = sp1.dt_chain.lower_bound() + sp2.dt_chain.lower_bound();

   s.SetTerm(0, &sp1, 3*no_labels, 3*sizeof(INDEX), nullptr );
   s.SetTerm(1, &sp2, 3*no_labels, 3*sizeof(INDEX), nullptr );

	s.options.gap_threshold = 0.000001;
	s.options.iter_max = 100;

   double* w = s.Solve();

   std::cout << "\n\nsolution = ";
   for(INDEX i=0; i<9; ++i) {
      std::cout << w[i] << ",";
   }
   std::cout << "\n\n";

   //solver.Solve();
   //std::cout << "message passing cost = " << solver.lower_bound() << "\n";

   INDEX y[3];
   const double lb = max_fn(w, y, 1.0, &sp1) + dot_product_fn(w,y,&sp1) + max_fn(w, y, 1.0, &sp2) + dot_product_fn(w,y,&sp2);
   std::cout << "\n\nlower bound for original problem = " << -lb << "\n";
   std::cout << "initial lower bound for original problem = " << initial_lb << "\n";
   solver.Solve();
   std::cout << "optimal lower bound for original problem = " << solver.lower_bound() << "\n";
   return 0; 
}
*/



#include "discrete_tomography.h"
#include "visitors/standard_visitor.hxx"




using namespace LP_MP;

static INDEX no_labels;
static INDEX no_vars;

// YPtr is an array of ints holding the labeling of unaries
// TermData holds an LP_tree with discrete tomography subproblems, the unary indices
struct sub_problem {
   template<typename VEC>
   sub_problem(LP_tree t, VEC& u) : dt_chain(t), no_unaries(u.size() + 1) 
   {
      assert(no_unaries > 1);
      pairwise_factors = new FMC_DT::dt_sequential_pairwise_factor*[no_unaries-1];
      for(INDEX i=0; i<no_unaries-1; ++i) {
         pairwise_factors[i] = u[i];
      }
   }
   ~sub_problem()
   {
      delete[] pairwise_factors;
   }


   void add_weights(double* const wi, const REAL multiplier)
   {
      assert(this->no_unaries > 2);
      const REAL mult = multiplier*this->sign; 
      //std::cout << add_sign << "," << mult << "\n";
      // first move the unary costs to the pairwise factors
      for(INDEX i=0; i<this->no_unaries-1; ++i) {
         auto* f = this->pairwise_factors[i]->GetFactor();
         for(INDEX x1=0; x1<no_labels; ++x1) {
            for(INDEX x2=0; x2<no_labels; ++x2) {
               assert(!std::isnan(wi[this->unaries[i]*no_labels + x1]));
               (*f).reg(x1,x2) += mult*wi[this->unaries[i]*no_labels + x1];
            } 
         }
      }
      // last pairwise is reparametrized twice
      {
         auto* f = this->pairwise_factors[ this->no_unaries-2 ]->GetFactor();
         for(INDEX x1=0; x1<no_labels; ++x1) {
            for(INDEX x2=0; x2<no_labels; ++x2) {
               assert(!std::isnan(wi[this->unaries[this->no_unaries-1]*no_labels + x2]));
               (*f).reg(x1,x2) += mult*wi[this->unaries[this->no_unaries-1]*no_labels + x2];
            } 
         } 
      } 
   }

   LP_tree dt_chain;
   INDEX no_unaries;
   std::vector<INDEX> unaries;
   FMC_DT::dt_sequential_pairwise_factor** pairwise_factors; // one less than number of unaries
   REAL sign;
};

static double dot_product_fn(double* wi, YPtr _y, TermData term_data);
double max_fn(double* wi, YPtr _y, double kappa, TermData term_data) // maximization oracle. Must copy argmax_y <a^{iy},[PAD(wi) kappa]> to y, and return the free term a^{iy}[d].
{
   //assert(kappa == 1.0);
   assert(kappa > 0.0);
   //if(kappa != 1.0) { std::cout << kappa << "\n"; }
   assert(kappa > std::numeric_limits<double>::min());
   assert(kappa < std::numeric_limits<double>::max());
   sub_problem* sp = (sub_problem*) term_data;

   sp->add_weights(wi, -1.0/kappa);

   sp->dt_chain.compute_subgradient();
   INDEX* y = (INDEX*) _y;
   for(INDEX i=0; i<sp->no_unaries-1; ++i) {
      y[i] = sp->pairwise_factors[i]->GetFactor()->state_[0]; 
   }
   y[sp->no_unaries-1] = sp->pairwise_factors[sp->no_unaries-2]->GetFactor()->state_[1];
   for(INDEX i=0; i<sp->no_unaries; ++i) {
      assert(y[i] < no_labels);
   }

   const REAL cost_test = -sp->dt_chain.lower_bound();
   // now remove the Lagrangean variables from unary factors again
   sp->add_weights(wi, 1.0/kappa);

   assert(std::abs(cost_test - (-sp->dt_chain.primal_cost() + 1.0/kappa*dot_product_fn(wi, y, term_data))) <= eps);

   // free term excluding wi
   const REAL cost = -sp->dt_chain.primal_cost()*kappa; // we minimize, but SVM expects maximization
   return cost;
}

static bool compare_fn(YPtr _y1, YPtr _y2, TermData term_data)
{
   sub_problem* sp = (sub_problem*) term_data;
   //return (!memcmp(_y1,_y2, sizeof(INDEX)*sp->no_unaries));
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
   const REAL sign = sp->sign;
   assert(sign == 1.0 || sign == -1.0);
   //std::cout << "\ncopy fn on " << sp << "\n";
   std::fill(ai, ai+no_labels*no_vars, 0.0);
   INDEX* y = (INDEX*) _y;
   for(INDEX i=0; i<sp->no_unaries; ++i) {
      const INDEX label = y[i];
      assert(label < no_labels);
      const INDEX var = sp->unaries[i];
      assert(var < no_vars);
      //std::cout << label << ", ";
      ai[var*no_labels + label] = sign*1.0;
   }
   //std::cout << "\n";
}

static double dot_product_fn(double* wi, YPtr _y, TermData term_data)
{
   sub_problem* sp = (sub_problem*) term_data;
   const REAL sign = sp->sign;
   assert(sign == 1.0 || sign == -1.0);
   //std::cout << "dot product on " << s << "\n";
   INDEX* y = (INDEX*) _y;
   double v = 0;
   //std::cout << sp->no_unaries << ": ";
   for(INDEX i=0; i<sp->no_unaries; ++i) {
      const INDEX label = y[i];
      assert(label < no_labels);
      const INDEX var = sp->unaries[i];
      assert(var < no_vars);
      v += sign*wi[var*no_labels + label];
      //std::cout << "(" << var << "," << label << "," << wi[var*no_labels + label] << ") ";
   }
   //std::cout << "; " << v << "\n";
   return v;
}

int main(int argc, char**argv)
{
   std::string filename(argv[1]);

   double scaling = 1.0;
   if(argc > 2) {
      scaling = std::stod(argv[2]);
   }
   std::cout << "scaling = " << scaling << "\n";


   //MpRoundingSolver<Solver<FMC_DT,LP_sat<LP>,StandardVisitor>> solver;
   Solver<FMC_DT,LP_FW_DS,StandardVisitor> solver;

   pegtl::file_parser problem(filename);

   UaiMrfInput::MrfInput mrfInput;
   bool ret = problem.parse< UaiMrfInput::grammar, UaiMrfInput::action>(mrfInput);
   if(ret != true) {
      throw std::runtime_error("could not read mrf problem in uai format for discrete tomography");
   }
   UaiMrfInput::build_mrf(solver.template GetProblemConstructor<0>(), mrfInput);
   auto& mrf = solver.template GetProblemConstructor<0>();

   solver.template GetProblemConstructor<1>().SetNumberOfLabels(mrfInput.cardinality_[0]);
   no_labels = mrfInput.cardinality_[0];
   no_vars = mrf.GetNumberOfVariables();
   std::cout << no_vars << ", " << no_labels << "\n";
   std::cout << "number of MRF factors = " << solver.GetLP().GetNumberOfFactors() << "\n";

   LP_MP::DiscreteTomographyTextInput::Projections p;
   ret = problem.parse< LP_MP::DiscreteTomographyTextInput::grammar, LP_MP::DiscreteTomographyTextInput::action>(p);
   if(ret != true) {
      throw std::runtime_error("could not read projection constraints for discrete tomography");
   }

   SVM* s = new SVM(mrf.GetNumberOfVariables()*no_labels, p.projectionVar.size(), max_fn, copy_fn, compare_fn, dot_product_fn, nullptr, false);//int d, int n, MaxFn max_fn, CopyFn copy_fn, CompareFn compare_fn, DotProductFn dot_product_fn, DotProductKernelFn dot_product_kernel_fn, bool zero_lower_bound);
	s->SetParams(lambda, mu / mu_SCALE, kappa / kappa_SCALE);
   assert(mu/mu_SCALE == 1.0);

   assert(p.projectionVar.size() == p.projectionCost.size());
   std::vector<sub_problem*> sp_vec;
   for(INDEX i=0; i<p.projectionVar.size(); ++i) {
      LP_tree t;
      solver.template GetProblemConstructor<1>().AddProjection(p.projectionVar[i], p.projectionCost[i], &t);
      solver.GetLP().add_tree(t);

      auto pairwise = t.get_factors<FMC_DT::dt_sequential_pairwise_factor>();
      sub_problem* sp = new sub_problem(t, pairwise);
      t.compute_subgradient();
      if(p.projectionVar[i][0] + 1 == p.projectionVar[i][1]) { // horizontal projection
         sp->sign = 1.0; 
      } else {
         sp->sign = -1.0; 
      }
      //std::cout << "sign = " << sp->sign << "\n";

      sp->unaries = p.projectionVar[i];

      s->SetTerm(i, sp, mrf.GetNumberOfVariables()*no_labels, p.projectionVar[i].size()*sizeof(INDEX), nullptr );
      sp_vec.push_back(sp); 
   }

   // decompose mrf into trees automatically. Do this after adding projection, since they can add pairwise factors as well.
   auto trees = mrf.compute_forest_cover();
   for(auto& tree : trees) {
      solver.GetLP().add_tree(tree);
   }


   std::cout << "number of factors = " << solver.GetLP().GetNumberOfFactors() << "\n";

   solver.Solve();
   return 0;

   s->options.gap_threshold = 0.000001;
	s->options.iter_max = 1000;

   // push pairwise into dt factors
   for(INDEX i=0; i<mrf.GetNumberOfPairwiseFactors(); ++i) {
      auto* f = mrf.GetPairwiseFactor(i);
      for(INDEX x1=0; x1<no_labels; ++x1) {
         for(INDEX x2=0; x2<no_labels; ++x2) {
            f->GetFactor()->cost(x1,x2) *= scaling;
         }
      }
      for(INDEX m=0; m<f->GetNoMessages(); ++m) {
         auto * msg = f->GetMessage(m);
         if(auto* msg_t = dynamic_cast<typename FMC_DT::dt_pairwise_pairwise_message*> (msg)) {
            msg_t->ReceiveMessageFromLeftContainer();
         }
      } 
   }

   //solver.Solve();

   double* w = s->Solve();

   //std::cout << "final Lagrange mulp = ";
   //for(INDEX i=0; i<no_labels*mrf.GetNumberOfVariables(); ++i) {
   //   std::cout << w[i] << ",";
   //}
   //std::cout << "\n";
   //REAL mult_norm = 0.0;
   //for(INDEX i=0; i<no_labels*mrf.GetNumberOfVariables(); ++i) {
   //   mult_norm += w[i]*w[i];
   //}
   //std::cout << "Lagrange mult. norm = " << mult_norm << "\n";

   for(INDEX iter=0; iter<400; ++iter) {
      // put weights into trees
      for(INDEX i=0; i<sp_vec.size(); ++i) {
         sp_vec[i]->add_weights(w, -1.0);
      }

      delete s;

      /*
      solver.Solve();
      for(INDEX i=0; i<mrf.GetNumberOfPairwiseFactors(); ++i) {
         auto* f = mrf.GetPairwiseFactor(i);
         for(INDEX m=0; m<f->GetNoMessages(); ++m) {
            auto * msg = f->GetMessage(m);
            if(auto* msg_t = dynamic_cast<typename FMC_DT::dt_pairwise_pairwise_message*> (msg)) {
               msg_t->ReceiveMessageFromLeftContainer();
            }
         } 
      }
      */ 

      s = new SVM(mrf.GetNumberOfVariables()*no_labels, p.projectionVar.size(), max_fn, copy_fn, compare_fn, dot_product_fn, nullptr, false);
      s->SetParams(lambda, mu / mu_SCALE, kappa / kappa_SCALE);
      //s->SetParams(lambda/REAL(iter+2), mu / mu_SCALE, kappa / kappa_SCALE);
      for(INDEX i=0; i<sp_vec.size(); ++i) {
         s->SetTerm(i, sp_vec[i], mrf.GetNumberOfVariables()*no_labels, p.projectionVar[i].size()*sizeof(INDEX), nullptr );
      }
      s->options.gap_threshold = 0.000001;
      s->options.iter_max = 100;
      w = s->Solve();

   }



   //INDEX y[mrf.GetNumberOfVariables()];
   //double lb = 0.0;
   //for(INDEX i=0; i<sp_vec.size(); ++i) {
   //   lb += max_fn(w, y, 1.0, sp_vec[i]) + dot_product_fn(w,y,sp_vec[i]);
   //}
   //std::cout << "\n\nlower bound for original problem = " << -lb << "\n";


   solver.Solve();
   //std::cout << "optimal lower bound for original problem = " << solver.lower_bound() << "\n";
   delete s;
   return 0;
}


