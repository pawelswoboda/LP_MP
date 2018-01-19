#include "test_model.hxx"
#include "LP_FWMAP.hxx"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"

using namespace LP_MP;

int main(int argc, char** argv)
{
  Solver<test_FMC, LP_tree_FWMAP, StandardVisitor> s;
  auto& lp = s.GetLP();


  auto* f1 = new typename test_FMC::factor(0,1);
  lp.AddFactor(f1);
  {
    factor_tree t1;
    auto* f2 = new typename test_FMC::factor(1,0);
    auto* f3 = new typename test_FMC::factor(0,0);
    lp.AddFactor(f2);
    lp.AddFactor(f3);
    auto* m12 = new typename test_FMC::message(f1,f2);
    auto* m13 = new typename test_FMC::message(f1,f3);
    lp.AddMessage(m12);
    lp.AddMessage(m13);
    t1.AddMessage(m12, Chirality::left);
    t1.AddMessage(m13, Chirality::left);
    t1.init();
    lp.add_tree(t1);
  }

  {
    factor_tree t2;
    auto* f2 = new typename test_FMC::factor(1,0);
    auto* f3 = new typename test_FMC::factor(0,0);
    lp.AddFactor(f2);
    lp.AddFactor(f3);
    auto* m12 = new typename test_FMC::message(f1,f2);
    auto* m23 = new typename test_FMC::message(f2,f3);
    lp.AddMessage(m12);
    lp.AddMessage(m23);
    t2.AddMessage(m12, Chirality::right);
    t2.AddMessage(m23, Chirality::left);
    t2.init();
    lp.add_tree(t2);
  }

  s.Solve();
}


