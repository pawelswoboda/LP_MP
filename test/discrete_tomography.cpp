#include "catch.hpp"
#include <vector>
#include "visitors/standard_visitor.hxx"
#include "solvers/discrete_tomography/discrete_tomography.h"

using namespace LP_MP;
using namespace std;

// do zrobienia: this will not be needed anymore, once we include reparametrization into factor directly.
template<typename FACTOR>
struct FactorContainerMock : public std::vector<REAL> {
   FactorContainerMock(FACTOR& f, std::vector<REAL>& cost) 
      : f_(f),
      std::vector<REAL>(cost)
   {}
   FACTOR* GetFactor() { return &f_; }
   FACTOR& f_;
};


TEST_CASE( "discrete tomography", "[dt]" ) {
   // other test ideas:
   // -if gurobi is present, solve tree model and other and check whether primal computed by ILP has correct solution and message passing gave correct lower/upper bounds
   // -test whether individual operations of factors/messages are correct. some operations of counting factors/messages already implemented.

   const INDEX noLabels = 3;
   std::vector<REAL> PottsCost = {0.0,1.0,1.0,
                                  1.0,0.0,1.0,
                                  1.0,1.0,0.0};

   std::vector<INDEX> projectionVar {0,1,2,3};
   std::vector<REAL> projectionCost {10,10,0,100};


   
   SECTION( "counting factors/messages" ) { // a small triangle
      // below we check non-negativity: if all factors have only non-negative potentials, reparametrized ones should stay non-negative as well.
      DiscreteTomographyFactorCounting left_fac(noLabels,1,1,projectionCost.size());
      DiscreteTomographyFactorCounting right_fac(noLabels,1,1,projectionCost.size());
      DiscreteTomographyFactorCounting top_fac(noLabels,2,2,projectionCost.size());
      DiscreteTomographyMessageCounting<DIRECTION::left> left_msg(noLabels,2,2,projectionCost.size());
      DiscreteTomographyMessageCounting<DIRECTION::right> right_msg(noLabels,2,2,projectionCost.size());

      std::vector<REAL> left_msg_var(left_msg.size(),0.0);
      std::vector<REAL> right_msg_var(right_msg.size(),0.0);

      // factors holds costs as follows:
      // size(up) + size(left) + size(right) + noLabels^2
      // upper label*sum + left label*sum + right label*sum + regularizer
      std::vector<REAL> left_pot(left_fac.size());
      std::vector<REAL> right_pot(right_fac.size());
      std::vector<REAL> top_pot(top_fac.size());
      for(INDEX i=0; i<PottsCost.size(); ++i) {
         left_pot[left_pot.size() - pow(noLabels,2) + i] = PottsCost[i];
         right_pot[left_pot.size() - pow(noLabels,2) + i] = PottsCost[i];
         top_pot[left_pot.size() - pow(noLabels,2) + i] = PottsCost[i];
      }

      for(INDEX i=0;i<projectionCost.size();i++){
         for(INDEX j=0;j<pow(noLabels,2);j++){
            top_pot[j+i*pow(noLabels,2)] = projectionCost[i];
         }
      }


      FactorContainerMock<DiscreteTomographyFactorCounting> left_fac_cont(left_fac, left_pot);
      FactorContainerMock<DiscreteTomographyFactorCounting> right_fac_cont(right_fac, right_pot);
      FactorContainerMock<DiscreteTomographyFactorCounting> top_fac_cont(top_fac, top_pot);

      left_msg.MakeLeftFactorUniform(&left_fac, left_pot, left_msg_var, 1.0);
      REQUIRE(*std::max_element(left_msg_var.begin(), left_msg_var.end()) <= eps);
      for(INDEX i=0; i<left_msg_var.size(); ++i) {
         left_msg.RepamLeft(left_fac_cont, left_msg_var[i], i);
         REQUIRE(*std::min_element(left_pot.begin(), left_pot.end()) >= -eps);
         left_msg.RepamRight(top_fac_cont, -left_msg_var[i], i);
         REQUIRE(*std::min_element(top_pot.begin(), top_pot.end()) >= -eps);
      }

      right_msg.MakeLeftFactorUniform(&right_fac, right_pot, right_msg_var, 1.0);
      REQUIRE(*std::max_element(right_msg_var.begin(), right_msg_var.end()) <= eps);
      for(INDEX i=0; i<right_msg_var.size(); ++i) {
         right_msg.RepamLeft(right_fac_cont, right_msg_var[i], i);
         REQUIRE(*std::min_element(right_pot.begin(), right_pot.end()) >= -eps);
         left_msg.RepamRight(top_fac_cont, -left_msg_var[i], i);
         REQUIRE(*std::min_element(top_pot.begin(), top_pot.end()) >= -eps);
      }

      auto top_pot_copy = top_pot;
      std::fill(left_msg_var.begin(), left_msg_var.end(), 0.0);
      left_msg.MakeRightFactorUniform(&top_fac, top_pot_copy, left_msg_var, 0.5);
      for(INDEX i=0; i<left_msg_var.size(); ++i) {
         left_msg.RepamLeft(left_fac_cont, -left_msg_var[i], i);
         REQUIRE(*std::min_element(left_pot.begin(), left_pot.end()) >= -eps);
         left_msg.RepamRight(top_fac_cont, left_msg_var[i], i);
         REQUIRE(*std::min_element(top_pot.begin(), top_pot.end()) >= -eps);
      }
      right_msg.MakeRightFactorUniform(&top_fac, top_pot_copy, right_msg_var, 0.5);
      std::fill(right_msg_var.begin(), right_msg_var.end(), 0.0);
      for(INDEX i=0; i<right_msg_var.size(); ++i) {
         right_msg.RepamLeft(right_fac_cont, -right_msg_var[i], i);
         REQUIRE(*std::min_element(right_pot.begin(), right_pot.end()) >= -eps);
         left_msg.RepamRight(top_fac_cont, left_msg_var[i], i);
         REQUIRE(*std::min_element(top_pot.begin(), top_pot.end()) >= -eps);
      }

   }

   std::vector<std::string> options = {
      {"discrete tomography test"},
      {"--inputFile"}, {""}, // note: we dot not have an input file, but argument is mandatory
      {"--maxIter"}, {"10"},
      {"--lowerBoundComputationInterval"}, {"1"}
   };
   
   VisitorSolver<Solver<FMC_DT>,StandardVisitor> s(options);
   auto& mrf = s.GetProblemConstructor<0>();
   auto& dt = s.GetProblemConstructor<1>();
   dt.SetNumberOfLabels(noLabels);


   SECTION( "tree" ) { // build a tree for a single projection and optimize. Check, whether global optimum could be achieved.
      // first build chain with four nodes
      mrf.AddUnaryFactor(0,std::vector<REAL>{0.0,1.0,0.5});
      mrf.AddUnaryFactor(1,std::vector<REAL>{0.0,1.0,0.5});
      mrf.AddUnaryFactor(2,std::vector<REAL>{0.0,1.0,0.5});
      mrf.AddUnaryFactor(3,std::vector<REAL>{0.0,1.0,0.5});

      // now link with pairwise Potts factors
      mrf.AddPairwiseFactor(0,1,PottsCost);
      mrf.AddPairwiseFactor(1,2,PottsCost);
      mrf.AddPairwiseFactor(2,3,PottsCost);

      dt.AddProjection(projectionVar,projectionCost);
      //s.Solve();
   }

}
