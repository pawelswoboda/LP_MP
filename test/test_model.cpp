#include "config.hxx"
#include "factors_messages.hxx"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"
#include "test.h"
#include "LP_external_interface.hxx"
#include "test_model.hxx"
#include <random>

using namespace LP_MP; 

// build simple model and check properties of factors/messages
int main()
{

  //std::vector<std::string> options({{"test solver"}, {"--maxIter"}, {"5"}});
  //MpRoundingSolver<Solver<test_FMC, LP, StandardVisitor>> s(options);
  Solver<test_FMC, LP_external_solver<DD_ILP::problem_export, LP>, StandardVisitor> s;
  auto& lp = s.GetLP();
  auto* f1 = new typename test_FMC::factor(0,1);
  auto* f2 = new typename test_FMC::factor(1,0);
  auto* f3 = new typename test_FMC::factor(0,0);
  lp.AddFactor(f1);
  lp.AddFactor(f2);
  lp.AddFactor(f3);

  auto* m12 = new typename test_FMC::message(f1,f2);
  auto* m13 = new typename test_FMC::message(f1,f3);
  lp.AddMessage(m12);
  lp.AddMessage(m13);

  test(lp.GetNumberOfFactors() == 3);
  test(lp.GetNumberOfMessages() == 2);
  test(lp.GetFactor(0) == f1);
  test(lp.GetFactor(1) == f2);
  test(lp.GetFactor(2) == f3);

  test(f1->no_messages() == 2);
  test(f1->no_send_messages() == 2);

  test(f2->no_messages() == 1);
  test(f2->no_send_messages() == 0);

  test(f3->no_messages() == 1);
  test(f3->no_send_messages() == 0);

  s.GetLP().solve();
  s.GetLP().get_external_solver().write_to_file("test_problem.lp");
  //s.Solve();
}
