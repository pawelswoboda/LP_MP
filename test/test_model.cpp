#include "config.hxx"
#include "factors_messages.hxx"
#include "solver.hxx"
#include "visitors/standard_visitor.hxx"
#include "test.h"
#include "LP_external_interface.hxx"
#include <random>

using namespace LP_MP;

struct test_factor {
  test_factor(const REAL x, const REAL y)
    : cost(2)
  {
    cost[0] = x;
    cost[1] = y;
  }
  REAL LowerBound() const { return cost.min(); }
  REAL EvaluatePrimal() const { return cost[primal]; }

  void MaximizePotentialAndComputePrimal()
  {
    primal = std::min_element(cost.begin(),cost.end()) - cost.begin();
  }

  void init_primal() { primal = 0; }

  template<typename ARCHIVE> void serialize_dual(ARCHIVE& ar) {};
  template<typename ARCHIVE> void serialize_primal(ARCHIVE& ar) {}; 

  auto export_variables() { return std::tie(cost); } 

  template<typename SOLVER>
  void construct_constraints(SOLVER& s, typename SOLVER::vector v)
  {
    s.add_simplex_constraint(v.begin(), v.end()); 
  }

  template<typename SOLVER>
  void convert_primal(SOLVER& s, typename SOLVER::vector v)
  {
    if(s.solution(v[0])) { primal = 0; }
    else { primal = 1; }
  }

  vector<REAL> cost;
  INDEX primal;
};

struct test_message {

  template<typename LEFT_FACTOR, typename MSG>
  void RepamLeft(LEFT_FACTOR& l, MSG msg)
  {
    assert(msg.size() == 2);
    l.cost[0] += msg[0];
    l.cost[1] += msg[1]; 
  }

  template<typename RIGHT_FACTOR, typename MSG>
  void RepamRight(RIGHT_FACTOR& r, MSG msg)
  {
    assert(msg.size() == 2);
    r.cost[0] += msg[0];
    r.cost[1] += msg[1]; 
  }

  template<typename LEFT_FACTOR, typename MSG>
  void send_message_to_right(const LEFT_FACTOR& l, MSG msg, const REAL omega)
  {
    msg -= omega*l.cost; 
  }

  template<typename RIGHT_FACTOR, typename MSG>
  void send_message_to_left(const RIGHT_FACTOR& r, MSG msg, const REAL omega)
  {
    msg -= omega*r.cost; 
  }

  template<typename SOLVER, typename LEFT_FACTOR, typename RIGHT_FACTOR>
  void construct_constraints(SOLVER& s, LEFT_FACTOR& l, typename SOLVER::vector v_left, RIGHT_FACTOR& r, typename SOLVER::vector v_right)
  {
    s.make_equal(v_left[0], v_right[0]);
    s.make_equal(v_left[1], v_right[1]);
    //s.make_equal(v_left.begin(), v_left.end(), v_right.begin(), v_right.end()); 
  }
};

struct FMC { // factor message connection
  constexpr static const char* name = "test model";
  using factor = FactorContainer<test_factor, FMC, 0>;
  using message = MessageContainer<test_message, 0, 0, message_passing_schedule::left, variableMessageNumber, atMostOneMessage, FMC, 0>;
  using FactorList = meta::list<factor>;
  using MessageList = meta::list<message>;
  using ProblemDecompositionList = meta::list<>;
};

// build simple model and check properties of factors/messages
int main()
{

  //std::vector<std::string> options({{"test solver"}, {"--maxIter"}, {"5"}});
  //MpRoundingSolver<Solver<FMC, LP, StandardVisitor>> s(options);
  Solver<FMC, LP_external_solver<DD_ILP::problem_export, LP>, StandardVisitor> s;
  auto& lp = s.GetLP();
  auto* f1 = new typename FMC::factor(0,1);
  auto* f2 = new typename FMC::factor(1,0);
  auto* f3 = new typename FMC::factor(0,0);
  lp.AddFactor(f1);
  lp.AddFactor(f2);
  lp.AddFactor(f3);

  auto* m12 = new typename FMC::message(f1,f2);
  auto* m13 = new typename FMC::message(f1,f3);
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
