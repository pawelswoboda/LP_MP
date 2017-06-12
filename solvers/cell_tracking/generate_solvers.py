#!/usr/bin/python

preamble = """
#include "cell_tracking.h"
#include "visitors/standard_visitor.hxx"
using namespace LP_MP;
int main(int argc, char* argv[])
"""

FMC = [ 
      'FMC_CELL_TRACKING_MOTHER_MACHINE',
      'FMC_CONSERVATION_TRACKING',
      'FMC_CELL_TRACKING'
      ]

parse_fun = [
      'cell_tracking_parser_mother_machine::ParseProblemMotherMachine',
      'conservation_tracking_parser::ParseProblem',
      'cell_tracking_parser_2d::ParseProblem'

      ]

file_name = [
      'cell_tracking_mother_machine.cpp',
      'conservation_tracking.cpp',
      'cell_tracking.cpp'
      ]


#main_body = ["MpRoundingSolver<Solver<" + e[0] + ",LP_sat<LP>,StandardTighteningVisitor>> solver(argc,argv);\nsolver.ReadProblem(" + e[1] + "<Solver<" + e[0] + ",LP_sat<LP>,StandardTighteningVisitor>>);\nreturn solver.Solve();\n" for e in zip(FMC,parse_fun)]
main_body = ["Solver<" + e[0] + ",LP_sat<LP>,StandardTighteningVisitor> solver(argc,argv);\nsolver.ReadProblem(" + e[1] + "<Solver<" + e[0] + ",LP_sat<LP>,StandardTighteningVisitor>>);\nreturn solver.Solve();\n" for e in zip(FMC,parse_fun)]

for e in zip(main_body,file_name):
   f = open(e[1],'w')
   f.write(preamble)
   f.write("\n{\n");
   f.write(e[0])
   f.write("\n}\n");
   f.close()
