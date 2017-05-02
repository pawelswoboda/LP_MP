#!/usr/bin/python

preamble = """
#include "multicut.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
int main(int argc, char* argv[])
"""

FMC = [ 
      'FMC_MULTICUT<MessageSendingType::SRMP>',
      'FMC_MULTICUT<MessageSendingType::SRMP>',
      'FMC_ODD_WHEEL_MULTICUT<MessageSendingType::SRMP>',
      'FMC_LIFTED_MULTICUT',
      'FMC_LIFTED_MULTICUT',
      'FMC_LIFTED_MULTICUT',
      'FMC_MULTIWAY_CUT'
      ]

parse_fun = [
      'MulticutTextInput::ParseProblem',
      'MulticutOpenGmInput::ParseProblem',
      'MulticutOpenGmInput::ParseProblem',
      'MulticutTextInput::ParseLiftedProblem',
      'MulticutH5Input::ParseLiftedProblem',
      'MulticutH5Input::ParseLiftedProblem',
      'MulticutOpenGmInput::ParsePottsProblem'
      ]

parse_fun_param = [
      '',
      '',
      '',
      '',
      ',false',
      ',true',
      ''
      ]

file_name = [
      'multicut_srmp.cpp',
      'multicut_opengm_srmp_cycle.cpp',
      'multicut_opengm_srmp_cycle_odd_wheel.cpp',
      'lifted_multicut_text.cpp',
      'lifted_multicut_h5.cpp',
      'lifted_multicut_h5_grid.cpp',
      'potts_via_multiway_cut.cpp'
      ]


main_body = ["CombinedMPProblemConstructorRoundingSolver<Solver<" + e[0] + ",LP,StandardTighteningVisitor>> solver(argc,argv);\nsolver.ReadProblem(" + e[1] + "<Solver<" + e[0] + ",LP,StandardTighteningVisitor>" + e[2] + ">);\nreturn solver.Solve();\n" for e in zip(FMC,parse_fun, parse_fun_param)]
#main_body = ["LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_MP_ROUNDING(" + e[0] + ",TorresaniEtAlInput::" + e[1] + "<" + e[0] + ">,StandardTighteningVisitor);" for e in zip(FMC,parse_fun)]

for e in zip(main_body,file_name):
   f = open(e[1],'w')
   f.write(preamble)
   f.write("\n{\n");
   f.write(e[0])
   f.write("\n}\n");
   f.close()
