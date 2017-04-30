#include "cell_tracking.h"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"
using FMC_INST = LP_MP::FMC_CELL_TRACKING_MOTHER_MACHINE;
//LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_SAT(FMC_INST, LP_MP::cell_tracking_parser::ParseProblemMotherMachine, StandardVisitor);
LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_SAT_CONCURRENT(FMC_INST, LP_MP::cell_tracking_parser::ParseProblemMotherMachine, StandardVisitor);
/*
#include <cryptominisat5/cryptominisat.h>
#include <assert.h>
#include <vector>
using std::vector;
using namespace CMSat;

int main()
{
  SATSolver solver;
  vector<Lit> clause;

  //Let's use 4 threads
  solver.set_num_threads(4);

  //We need 3 variables
  solver.new_vars(3);

  //adds "1 0"
  clause.push_back(Lit(0, false));
  solver.add_clause(clause);

  //adds "-2 0"
  clause.clear();
  clause.push_back(Lit(1, true));
  solver.add_clause(clause);

  //adds "-1 2 3 0"
  clause.clear();
  clause.push_back(Lit(0, true));
  clause.push_back(Lit(1, false));
  clause.push_back(Lit(2, false));
  solver.add_clause(clause);

  lbool ret = solver.solve();
  assert(ret == l_True);
  assert(solver.get_model()[0] == l_True);
  assert(solver.get_model()[1] == l_False);
  assert(solver.get_model()[2] == l_True);
  std::cout
    << "Solution is: "
    << solver.get_model()[0]
    << ", " << solver.get_model()[1]
    << ", " << solver.get_model()[2]
    << std::endl;

  solver.solve();

  return 0;
}

*/
