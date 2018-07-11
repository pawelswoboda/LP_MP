#include "catch.hpp"
#include "factors_messages.hxx"
#include "solvers/graphical_model/graphical_model.h"
#include "visitors/standard_visitor.hxx"

using namespace LP_MP;

TEST_CASE("factor and message interface", "[]") {
      using FMC = FMC_SRMP;
      using VisitorType = StandardVisitor;
      using SolverType = MpRoundingSolver<Solver<FMC,LP,VisitorType>>;

      std::vector<std::string> solver_options = {
        {"factor messages test"},
        {"--inputFile"},{""}
      };

      SolverType s(solver_options);
      auto& mrf = s.template GetProblemConstructor<0>();

      std::vector<REAL> unary(2);
      matrix<REAL> pairwise(2,2);

      mrf.AddUnaryFactor(unary);
      mrf.AddUnaryFactor(unary);
      mrf.AddUnaryFactor(unary);
      mrf.AddUnaryFactor(unary);
      mrf.AddUnaryFactor(unary);

      mrf.AddPairwiseFactor(0,1,pairwise);
      mrf.AddPairwiseFactor(1,2,pairwise);
      mrf.AddPairwiseFactor(2,3,pairwise);
      mrf.AddPairwiseFactor(0,3,pairwise);

      mrf.AddPairwiseFactor(0,4,pairwise);
      mrf.AddPairwiseFactor(1,4,pairwise);
      mrf.AddPairwiseFactor(2,4,pairwise);
      mrf.AddPairwiseFactor(3,4,pairwise);

      SECTION("factor connections") {
        auto is_f0_neighbor = [&](auto* f) {
          return mrf.GetPairwiseFactor(0,1) == f || mrf.GetPairwiseFactor(0,3) == f || mrf.GetPairwiseFactor(0,4) == f;
        };

        auto is_f4_neighbor = [&](auto* f) {
          return mrf.GetPairwiseFactor(0,4) == f || mrf.GetPairwiseFactor(1,4) == f || mrf.GetPairwiseFactor(2,4) == f || mrf.GetPairwiseFactor(3,4) == f;
        };

        auto* f0 = mrf.GetUnaryFactor(0);
        REQUIRE(f0->GetNoMessages() == 3);
        REQUIRE(f0->no_send_messages() == 3);
        auto m0_it = f0->begin();
        REQUIRE(is_f0_neighbor(m0_it.GetConnectedFactor()));
        ++m0_it;
        REQUIRE(is_f0_neighbor(m0_it.GetConnectedFactor()));
        REQUIRE(m0_it != f0->begin());
        REQUIRE(m0_it != f0->end());
        ++m0_it;
        REQUIRE(is_f0_neighbor(m0_it.GetConnectedFactor()));
        REQUIRE(m0_it != f0->begin());
        REQUIRE(m0_it != f0->end());
        ++m0_it;
        REQUIRE(m0_it == f0->end());


        auto* f4 = mrf.GetUnaryFactor(4);
        REQUIRE(f4->GetNoMessages() == 4);
        REQUIRE(f4->no_send_messages() == 4);
        auto m4_it = f4->begin();
        REQUIRE(is_f4_neighbor(m4_it.GetConnectedFactor()));
        ++m4_it;
        REQUIRE(is_f4_neighbor(m4_it.GetConnectedFactor()));
        REQUIRE(m4_it != f4->begin());
        REQUIRE(m4_it != f4->end());
        ++m4_it;
        REQUIRE(is_f4_neighbor(m4_it.GetConnectedFactor()));
        REQUIRE(m4_it != f4->begin());
        REQUIRE(m4_it != f4->end());
        ++m4_it;
        REQUIRE(is_f4_neighbor(m4_it.GetConnectedFactor()));
        REQUIRE(m4_it != f4->begin());
        REQUIRE(m4_it != f4->end());
        ++m4_it;
        REQUIRE(m4_it == f4->end());

      }
}

