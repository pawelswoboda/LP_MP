#!/usr/bin/python

preamble = """
#include "cell_tracking.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
int main(int argc, char* argv[])
"""

FMC = [ 
      'FMC_CELL_TRACKING_MOTHER_MACHINE',
      'FMC_CONSERVATION_TRACKING'
      ]

parse_fun = [
      'cell_tracking_parser::ParseProblemMotherMachine',
      'conservation_tracking_parser::ParseProblem'
      ]

file_name = [
      'cell_tracking.cpp',
      'conservation_tracking.cpp'
      ]


main_body = ["MpRoundingSolver<Solver<" + e[0] + ",LP_sat<LP_concurrent<LP>>,StandardVisitor>> solver(argc,argv);\nsolver.ReadProblem(" + e[1] + "<Solver<" + e[0] + ",LP_sat<LP_concurrent<LP>>,StandardVisitor>>);\nreturn solver.Solve();\n" for e in zip(FMC,parse_fun)]
#main_body = ["LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_MP_ROUNDING(" + e[0] + ",TorresaniEtAlInput::" + e[1] + "<" + e[0] + ">,StandardTighteningVisitor);" for e in zip(FMC,parse_fun)]

for e in zip(main_body,file_name):
   f = open(e[1],'w')
   f.write(preamble)
   f.write("\n{\n");
   f.write(e[0])
   f.write("\n}\n");
   f.close()
