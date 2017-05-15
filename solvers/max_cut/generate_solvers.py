#!/usr/bin/python

preamble = """
#include "max_cut.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
int main(int argc, char* argv[])
"""

FMC = [ 
      'FMC_CYCLE_MAX_CUT',
      'FMC_ODD_BICYCLE_WHEEL_MAX_CUT',
      'FMC_CYCLE_MAX_CUT',
      'FMC_ODD_BICYCLE_WHEEL_MAX_CUT'
      ]

parse_fun = [
      'max_cut_simple_text_format::ParseProblemMaxCut',
      'max_cut_simple_text_format::ParseProblemMaxCut',
      'max_cut_simple_text_format::ParseProblemQUBO',
      'max_cut_simple_text_format::ParseProblemQUBO'
      ]

file_name = [
      'max_cut_cycle.cpp',
      'max_cut_cycle_odd_bicycle_wheel.cpp',
      'QUBO_cycle.cpp',
      'QUBO_cycle_odd_bicycle_wheel.cpp' 
      ]


main_body = ["ProblemConstructorRoundingSolver<Solver<" + e[0] + ",LP,StandardTighteningVisitor>> solver(argc,argv);\nsolver.ReadProblem(" + e[1] + "<Solver<" + e[0] + ",LP,StandardTighteningVisitor>>);\nreturn solver.Solve();\n" for e in zip(FMC,parse_fun)]
#main_body = ["LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_MP_ROUNDING(" + e[0] + ",TorresaniEtAlInput::" + e[1] + "<" + e[0] + ">,StandardTighteningVisitor);" for e in zip(FMC,parse_fun)]

for e in zip(main_body,file_name):
   f = open(e[1],'w')
   f.write(preamble)
   f.write("\n{\n");
   f.write(e[0])
   f.write("\n}\n");
   f.close()
