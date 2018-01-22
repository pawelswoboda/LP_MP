#include "test_model.hxx"
#include "LP_conic_bundle.hxx"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"

using namespace LP_MP;

int main(int argc, char** argv)
{
  Solver<test_FMC, LP_conic_bundle, StandardVisitor> s;
  auto& lp = s.GetLP();

  build_test_model(lp);

  s.Solve();
}



