#include "catch.hpp"
#include <string>
#include <fstream>
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"
#include "solvers/graph_matching/graph_matching.h"

// do zrobienia: make those problems harder. Currently, the assignment constraints are not actually needed, as solution is a matching already.

TEST_CASE("graph matching", "[solving a graph matching test problem]") {

// test input in the format used by Torresani, Kolmogorov and Rother
const std::string torresani_test_input = 
R"(c example of graph matching problem
p 3 2 4 2
a 0   0 0    1.0
a 1   1 0    -2.0
a 2   1 1    1.5
a 3   2 1    4.0
e     0 2    -0.5
e     1 3    0.5
n0 0 1
n0 1 2
n1 0 1
)";

// test input in the modified uai format
const std::string uai_test_input =
R"(MARKOV
3
2 2 3
3
1 0
2 0 1
2 1 2

2
 0.436 0.564

4
 1e13 0.872
 0.920 0.080

6
 1e13 0.333 0.457
 0.811 0.000 0.189 

matching
0 0 slack 
1 0 slack
2 0 1 slack
)";

// to do: devise harder test instances which need more iterations to converge

std::vector<std::string> solver_options = {
   {"graph matching test"},
   {"--maxIter"}, {"5"},
   {"--timeout"}, {"60"}, // one minute
   {"--lowerBoundComputationInterval"}, {"1"},
   {"--standardReparametrization"}, {"anisotropic"},
   {"--roundingReparametrization"}, {"anisotropic"}
};

// algorithms considered
using FMC_GM_LEFT = FMC_GM<PairwiseConstruction::Left>;
using FMC_GM_RIGHT = FMC_GM<PairwiseConstruction::Right>;
using FMC_MCF_LEFT = FMC_MCF<PairwiseConstruction::Left>;
using FMC_MCF_RIGHT = FMC_MCF<PairwiseConstruction::Right>;
using FMC_MCF_BOTH_SIDES = FMC_MCF<PairwiseConstruction::BothSides>;
using FMC_MP_LEFT = FMC_MP<PairwiseConstruction::Left>;
using FMC_MP_RIGHT = FMC_MP<PairwiseConstruction::Right>;
using FMC_MP_BOTH_SIDES = FMC_MP<PairwiseConstruction::BothSides>;
using FMC_HUNGARIAN_BP_LEFT = FMC_HUNGARIAN_BP<PairwiseConstruction::Left>;
using FMC_HUNGARIAN_BP_RIGHT = FMC_HUNGARIAN_BP<PairwiseConstruction::Right>;
using FMC_HUNGARIAN_BP_BOTH_SIDES = FMC_HUNGARIAN_BP<PairwiseConstruction::BothSides>;

using VisitorType = StandardVisitor;

   SECTION("input in format used by Torresani et al") {

      std::string tmp_file_name = std::tmpnam(nullptr);
      std::cout << "temporary file name: " << tmp_file_name << '\n';

      std::ofstream tmp_file;
      tmp_file.open (tmp_file_name);
      tmp_file << torresani_test_input;
      tmp_file.close();

      auto torresani_options = solver_options;
      torresani_options.push_back("--inputFile");
      torresani_options.push_back(tmp_file_name);

      const REAL opt_val = 2.0;

      SECTION( "mp left" ) {
         VisitorSolver<MpRoundingSolver<FMC_MP_LEFT>,VisitorType> s(torresani_options);
         s.ReadProblem(TorresaniEtAlInput::ParseProblemMP<FMC_MP_LEFT>);
         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - -opt_val) < LP_MP::eps);
      }
      SECTION( "mp right" ) {
         VisitorSolver<MpRoundingSolver<FMC_MP_RIGHT>,VisitorType> s(torresani_options);
         s.ReadProblem(TorresaniEtAlInput::ParseProblemMP<FMC_MP_RIGHT>);
         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - -opt_val) < LP_MP::eps);
      }
      SECTION( "mp both sides" ) {
         VisitorSolver<MpRoundingSolver<FMC_MP_BOTH_SIDES>,VisitorType> s(torresani_options);
         s.ReadProblem(TorresaniEtAlInput::ParseProblemMP<FMC_MP_BOTH_SIDES>);
         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - -opt_val) < LP_MP::eps);
      }
      SECTION( "mcf left" ) {
         VisitorSolver<MpRoundingSolver<FMC_MCF_LEFT>,VisitorType> s(torresani_options);
         s.ReadProblem(TorresaniEtAlInput::ParseProblemMCF<FMC_MCF_LEFT>);
         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - -opt_val) < LP_MP::eps);
      }
      SECTION( "mcf right" ) {
         VisitorSolver<MpRoundingSolver<FMC_MCF_RIGHT>,VisitorType> s(torresani_options);
         s.ReadProblem(TorresaniEtAlInput::ParseProblemMCF<FMC_MCF_RIGHT>);
         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - -opt_val) < LP_MP::eps);
      }
      SECTION( "mcf both sides" ) {
         VisitorSolver<MpRoundingSolver<FMC_MCF_BOTH_SIDES>,VisitorType> s(torresani_options);
         s.ReadProblem(TorresaniEtAlInput::ParseProblemMCF<FMC_MCF_BOTH_SIDES>);
         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - -opt_val) < LP_MP::eps);
      }
      SECTION( "gm left" ) {
         VisitorSolver<MpRoundingSolver<FMC_GM_LEFT>,VisitorType> s(torresani_options);
         s.ReadProblem(TorresaniEtAlInput::ParseProblemGM<FMC_GM_LEFT>);
         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - -opt_val) < LP_MP::eps);
      }
      SECTION( "gm right" ) {
         VisitorSolver<MpRoundingSolver<FMC_GM_RIGHT>,VisitorType> s(torresani_options);
         s.ReadProblem(TorresaniEtAlInput::ParseProblemGM<FMC_GM_RIGHT>);
         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - -opt_val) < LP_MP::eps);
      }
      SECTION( "hungarian bp left" ) {
         VisitorSolver<MpRoundingSolver<FMC_HUNGARIAN_BP_LEFT>,VisitorType> s(torresani_options);
         s.ReadProblem(TorresaniEtAlInput::ParseProblemMCF<FMC_HUNGARIAN_BP_LEFT>);
         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - -opt_val) < LP_MP::eps);
      }
      SECTION( "hungarian bp right" ) {
         VisitorSolver<MpRoundingSolver<FMC_HUNGARIAN_BP_RIGHT>,VisitorType> s(torresani_options);
         s.ReadProblem(TorresaniEtAlInput::ParseProblemMCF<FMC_HUNGARIAN_BP_RIGHT>);
         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - -opt_val) < LP_MP::eps);
      }
      SECTION( "hungarian bp both sides" ) {
         VisitorSolver<MpRoundingSolver<FMC_HUNGARIAN_BP_BOTH_SIDES>,VisitorType> s(torresani_options);
         s.ReadProblem(TorresaniEtAlInput::ParseProblemMCF<FMC_HUNGARIAN_BP_BOTH_SIDES>);
         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - -opt_val) < LP_MP::eps);
      }

      std::remove(tmp_file_name.c_str());
   }

   SECTION("input in uai format") {

      std::string tmp_file_name = std::tmpnam(nullptr);
      std::cout << "temporary file name: " << tmp_file_name << '\n';

      std::ofstream tmp_file;
      tmp_file.open (tmp_file_name);
      tmp_file << uai_test_input;
      tmp_file.close();

      auto uai_options = solver_options;
      uai_options.push_back("--inputFile");
      uai_options.push_back(tmp_file_name);

      const REAL opt_val = 0.833;

      SECTION( "mp left" ) {
         VisitorSolver<MpRoundingSolver<FMC_MP_LEFT>,VisitorType> s(uai_options);
         s.ReadProblem(UaiGraphMatchingInput::ParseProblemMP<FMC_MP_LEFT>);
         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - opt_val) < LP_MP::eps);
      }
      SECTION( "mcf left" ) {
         VisitorSolver<MpRoundingSolver<FMC_MCF_LEFT>,VisitorType> s(uai_options);
         s.ReadProblem(UaiGraphMatchingInput::ParseProblemMCF<FMC_MCF_LEFT>);
         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - opt_val) < LP_MP::eps);
      }
      SECTION( "gm left" ) {
         VisitorSolver<MpRoundingSolver<FMC_GM_LEFT>,VisitorType> s(uai_options);
         s.ReadProblem(UaiGraphMatchingInput::ParseProblemGM<FMC_GM_LEFT>);
         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - opt_val) < LP_MP::eps);
      }
      SECTION( "hungarian bp left" ) {
         VisitorSolver<MpRoundingSolver<FMC_HUNGARIAN_BP_LEFT>,VisitorType> s(uai_options);
         s.ReadProblem(UaiGraphMatchingInput::ParseProblemMCF<FMC_HUNGARIAN_BP_LEFT>);
         s.Solve();
         REQUIRE(std::abs(s.lower_bound() - opt_val) < LP_MP::eps);
      }
      std::remove(tmp_file_name.c_str());
   }
}
