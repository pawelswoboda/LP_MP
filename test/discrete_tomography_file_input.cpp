#include "catch.hpp"
#include <vector>
#include <string>
#include <cstdio>
#include <fstream>
#include "visitors/standard_visitor.hxx"
#include "solvers/discrete_tomography/discrete_tomography.h"

using namespace LP_MP;
using namespace std;

TEST_CASE( "discrete tomography file input", "[dt file input]" ) {

   // write test_problem to a temporary fiel under /tmp, read from it and test whether it could be solved.

   const std::string test_problem = 
R"(MARKOV
9
3 3 3 3 3 3 3 3 3 
21
2 0 1
2 1 2
2 3 4
2 4 5
2 6 7
2 7 8
2 0 3
2 1 4
2 2 5
2 3 6
2 4 7
2 5 8
1 0
1 1
1 2
1 3
1 4
1 5
1 6
1 7
1 8

9
0.0000 1.0000 2.0000 1.0000 0.0000 1.0000 2.0000 1.0000 0.0000 
9
0.0000 1.0000 2.0000 1.0000 0.0000 1.0000 2.0000 1.0000 0.0000 
9
0.0000 1.0000 2.0000 1.0000 0.0000 1.0000 2.0000 1.0000 0.0000 
9
0.0000 1.0000 2.0000 1.0000 0.0000 1.0000 2.0000 1.0000 0.0000 
9
0.0000 1.0000 2.0000 1.0000 0.0000 1.0000 2.0000 1.0000 0.0000 
9
0.0000 1.0000 2.0000 1.0000 0.0000 1.0000 2.0000 1.0000 0.0000 
9
0.0000 1.0000 2.0000 1.0000 0.0000 1.0000 2.0000 1.0000 0.0000 
9
0.0000 1.0000 2.0000 1.0000 0.0000 1.0000 2.0000 1.0000 0.0000 
9
0.0000 1.0000 2.0000 1.0000 0.0000 1.0000 2.0000 1.0000 0.0000 
9
0.0000 1.0000 2.0000 1.0000 0.0000 1.0000 2.0000 1.0000 0.0000 
9
0.0000 1.0000 2.0000 1.0000 0.0000 1.0000 2.0000 1.0000 0.0000 
9
0.0000 1.0000 2.0000 1.0000 0.0000 1.0000 2.0000 1.0000 0.0000 
3
0.0000 0.0000 0.0000 
3
0.0000 0.0000 0.0000 
3
0.0000 0.0000 0.0000 
3
0.0000 0.0000 0.0000 
3
0.0000 0.0000 0.0000 
3
0.0000 0.0000 0.0000 
3
0.0000 0.0000 0.0000 
3
0.0000 0.0000 0.0000 
3
0.0000 0.0000 0.0000 

PROJECTIONS
6 + 7 + 8 = (Inf,Inf,Inf,0)
3 + 4 + 5 = (Inf,0)
0 + 1 + 2 = (Inf,Inf,Inf,0)
5 + 7 + 8 = (Inf,0)
2 + 4 + 6 = (Inf,Inf,Inf,Inf,Inf,0)
0 + 1 + 3 = (Inf,0)
2 + 5 + 8 = (Inf,Inf,Inf,0)
1 + 4 + 7 = (Inf,0)
0 + 3 + 6 = (Inf,Inf,Inf,0))";

   std::string tmp_file_name = std::tmpnam(nullptr);
   std::cout << "temporary file name: " << tmp_file_name << '\n';

   std::ofstream tmp_file;
   tmp_file.open (tmp_file_name);
   tmp_file << test_problem;
   tmp_file.close();


   std::vector<std::string> options = {
      {"discrete tomography test"},
      {"-i"}, {tmp_file_name}, // note: we dot not have an input file, but argument is mandatory
      {"--maxIter"}, {"10"}
      //{"--lowerBoundComputationInterval"}, {"1"}
   };

   VisitorSolver<Solver<FMC_DT>,StandardVisitor> s(options);
   s.ReadProblem(DiscreteTomographyTextInput::ParseProblem<FMC_DT>);
   s.Solve();
   REQUIRE(s.lower_bound() >= 0);

   std::remove(tmp_file_name.c_str());
}
