#include "catch.hpp"
#include <vector>
#include "visitors/standard_visitor.hxx"
#include "lp_interface/lp_gurobi.hxx"
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
   REAL mp_lb = 0.0;

   
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
      std::vector<REAL> left_pot(left_fac.size(),0.0);
      std::vector<REAL> right_pot(right_fac.size(),0.0);
      std::vector<REAL> top_pot(top_fac.size(),0.0);
      
      // Leaf Consistency
      INDEX eq = 1 + noLabels + pow(noLabels,2);
      REAL inf = std::numeric_limits<REAL>::infinity();
      REQUIRE(left_fac.getSize(DiscreteTomographyFactorCounting::NODE::left) == left_fac.getSize(DiscreteTomographyFactorCounting::NODE::right));
      REQUIRE(left_fac.getSize(DiscreteTomographyFactorCounting::NODE::right) == right_fac.getSize(DiscreteTomographyFactorCounting::NODE::left));
      REQUIRE(right_fac.getSize(DiscreteTomographyFactorCounting::NODE::left) == right_fac.getSize(DiscreteTomographyFactorCounting::NODE::right));
      for(INDEX i=0;i<left_fac.getSize(DiscreteTomographyFactorCounting::NODE::left);i++){
        if( i % eq != 0){
          left_pot[left_fac.getSize(DiscreteTomographyFactorCounting::NODE::up) + i] = inf;
          left_pot[left_fac.getSize(DiscreteTomographyFactorCounting::NODE::up) + left_fac.getSize(DiscreteTomographyFactorCounting::NODE::left) + i] = inf;
          right_pot[left_fac.getSize(DiscreteTomographyFactorCounting::NODE::up) + i] = inf;
          right_pot[left_fac.getSize(DiscreteTomographyFactorCounting::NODE::up) + left_fac.getSize(DiscreteTomographyFactorCounting::NODE::left) + i] = inf;
        }
      }
            
      // Regularizer
      for(INDEX i=0; i<PottsCost.size(); ++i) {
         left_pot[left_pot.size() - pow(noLabels,2) + i] = PottsCost[i];
         right_pot[right_pot.size() - pow(noLabels,2) + i] = PottsCost[i];
         top_pot[top_pot.size() - pow(noLabels,2) + i] = PottsCost[i];
      }

      // Counting Constraint
      for(INDEX i=0;i<projectionCost.size();i++){
         for(INDEX j=0;j<pow(noLabels,2);j++){
            top_pot[j+i*pow(noLabels,2)] = projectionCost[i];
         }
      }

      FactorContainerMock<DiscreteTomographyFactorCounting> left_fac_cont(left_fac, left_pot);
      FactorContainerMock<DiscreteTomographyFactorCounting> right_fac_cont(right_fac, right_pot);
      FactorContainerMock<DiscreteTomographyFactorCounting> top_fac_cont(top_fac, top_pot);

      auto PrintVec = [](std::vector<REAL> f){
        for(INDEX i=0;i<f.size();i++){
          printf("%03d -> %10.2f \n",i,f[i]);
        }
      };
      
      auto upSize = top_fac.getSize(DiscreteTomographyFactorCounting::NODE::up);
      auto leftSize = top_fac.getSize(DiscreteTomographyFactorCounting::NODE::left);
      auto rightSize = top_fac.getSize(DiscreteTomographyFactorCounting::NODE::right);
      REQUIRE(leftSize == rightSize);
      REQUIRE(top_fac_cont.size() == upSize + leftSize + rightSize + pow(noLabels,2));
      for(INDEX i=0;i<leftSize;i++){
        REQUIRE(top_fac_cont[upSize + i] == top_fac_cont[upSize + leftSize + i]);
      }
      
      auto SendUpMessages = [&](){
        std::fill(left_msg_var.begin(), left_msg_var.end(), 0.0);
        left_msg.MakeLeftFactorUniform(&left_fac, left_fac_cont, left_msg_var, 1.0);
        //printf("LeftUpMessage\n");
        //PrintVec(left_msg_var);
        REQUIRE(*std::max_element(left_msg_var.begin(), left_msg_var.end()) <= eps);
        for(INDEX i=0; i<left_msg_var.size(); ++i) {
          left_msg.RepamLeft(left_fac_cont, left_msg_var[i], i);
          REQUIRE(left_fac.LowerBound(left_fac_cont) >= -eps);
          left_msg.RepamRight(top_fac_cont, -left_msg_var[i], i);
          REQUIRE(top_fac.LowerBound(top_fac_cont) >= -eps);
        }

        std::fill(right_msg_var.begin(), right_msg_var.end(), 0.0);
        right_msg.MakeLeftFactorUniform(&right_fac, right_fac_cont, right_msg_var, 1.0);
        //printf("RightUpMessage\n");
        //PrintVec(right_msg_var);
        REQUIRE(*std::max_element(right_msg_var.begin(), right_msg_var.end()) <= eps);
        for(INDEX i=0; i<right_msg_var.size(); ++i) {
          right_msg.RepamLeft(right_fac_cont, right_msg_var[i], i);
          REQUIRE(right_fac.LowerBound(right_fac_cont) >= -eps);
          right_msg.RepamRight(top_fac_cont, -right_msg_var[i], i);
          REQUIRE(top_fac.LowerBound(top_fac_cont) >= -eps);
        }
        
        REQUIRE(left_msg_var.size() == right_msg_var.size());
        for(INDEX i=0;i<left_msg_var.size();i++){
          REQUIRE(left_msg_var[i] == right_msg_var[i]);
        }
        for(INDEX i=0;i<left_fac_cont.size();i++){
          REQUIRE(left_fac_cont[i] == right_fac_cont[i]);
        }
        
        REQUIRE(leftSize == rightSize);
        for(INDEX i=0;i<leftSize;i++){
          REQUIRE(top_fac_cont[upSize + i] == top_fac_cont[upSize + leftSize + i]);
        }
        
      };
      
      auto SendDownMessages = [&](){
        auto top_fac_cont_copy = top_fac_cont;
        
        std::fill(left_msg_var.begin(), left_msg_var.end(), 0.0);
        left_msg.MakeRightFactorUniform(&top_fac, top_fac_cont_copy, left_msg_var, 0.5);
        //printf("LeftDownMessage\n");
        //PrintVec(left_msg_var);
        for(INDEX i=0; i<left_msg_var.size(); ++i) {
          left_msg.RepamLeft(left_fac_cont, -left_msg_var[i], i);
          REQUIRE(left_fac.LowerBound(left_fac_cont) >= -eps);
          left_msg.RepamRight(top_fac_cont, left_msg_var[i], i);
          REQUIRE(top_fac.LowerBound(top_fac_cont) >= -eps);
        }
        
        std::fill(right_msg_var.begin(), right_msg_var.end(), 0.0);
        right_msg.MakeRightFactorUniform(&top_fac, top_fac_cont_copy, right_msg_var, 0.5);
        //printf("RightDownMessage\n");
        //PrintVec(right_msg_var);
        for(INDEX i=0; i<right_msg_var.size(); ++i) {
          right_msg.RepamLeft(right_fac_cont, -right_msg_var[i], i);
          REQUIRE(right_fac.LowerBound(right_fac_cont) >= -eps);
          right_msg.RepamRight(top_fac_cont, right_msg_var[i], i);
          REQUIRE(top_fac.LowerBound(top_fac_cont) >= -eps);
        }
      };
      

      
      SendUpMessages();
      SendDownMessages(); 
      
      //printf("Calc Lowerbound\n");
      mp_lb = top_fac.LowerBound(top_fac_cont) + left_fac.LowerBound(left_fac_cont) + right_fac.LowerBound(right_fac_cont);
      REQUIRE(mp_lb == 1.0);
      //printf("lb = %10.4f \n",lb);
      
      
   }

   std::vector<std::string> options = {
     {"discrete tomography test"},
      //{"--lpFile"},{"model.lp"},
      {"-i"}, {""}, // note: we dot not have an input file, but argument is mandatory
      //{"--maxIter"}, {"10"},
      //{"--lowerBoundComputationInterval"}, {"1"}
   };
   
   LpRoundingSolver<FMC_DT,LpInterfaceGurobi> s(options);
   auto& mrf = s.GetProblemConstructor<0>();
   auto& dt = s.GetProblemConstructor<1>();
   dt.SetNumberOfLabels(noLabels);


   SECTION( "LP Interface" ) { // LP Interface
      // first build chain with four nodes
      mrf.AddUnaryFactor(0,std::vector<REAL>{0.0,0.0,0.0});
      mrf.AddUnaryFactor(1,std::vector<REAL>{0.0,0.0,0.0});
      mrf.AddUnaryFactor(2,std::vector<REAL>{0.0,0.0,0.0});
      mrf.AddUnaryFactor(3,std::vector<REAL>{0.0,0.0,0.0});

      // now link with pairwise Potts factors
      mrf.AddPairwiseFactor(0,1,PottsCost);
      mrf.AddPairwiseFactor(1,2,PottsCost);
      mrf.AddPairwiseFactor(2,3,PottsCost);

      dt.AddProjection(projectionVar,projectionCost);
      s.Solve();

      auto LpSolver = s.GetLpSolver();

      std::vector<REAL> vars(4,0.0);
      REAL counting = 0.0;
      
      for(INDEX i=0;i<4;i++){
        for(INDEX j=0;j<noLabels;j++){
          vars[i] += j*LpSolver->GetVariableValue(i*noLabels + j);
        }
        counting += vars[i];
      }

      REQUIRE(std::fabs(counting-2) < eps);
      
   }

   SECTION(" Test CreateConstraint  "){
     DiscreteTomographyFactorCounting bigFactor(3,16,16,30);
     //std::vector<REAL> bigFactor_pot(bigFactor.size(),0.0);
     //FactorContainerMock<DiscreteTomographyFactorCounting> bigFactor_cont(bigFactor, bigFactor_pot);
     
     LpInterfaceGurobi LpBuilder(100*bigFactor.size());
     for(INDEX i=0;i<1;i++){
       printf("factor (%5d)\n",(int) i);
       LpBuilder.addFactor(bigFactor,i*bigFactor.size());
     }
     
   }
   
   
   
}

