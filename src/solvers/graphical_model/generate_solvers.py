#!/usr/bin/python

preamble = """
#include "graphical_model.h"
#include "visitors/standard_visitor.hxx"
int main(int argc, char* argv[])
"""

FMC = [ 
      'FMC_SRMP',
      'FMC_SRMP_T',
      'FMC_MPLP',
      'FMC_SRMP',
      'FMC_MPLP'
      ]

parse_fun = [
      'HDF5Input::ParseProblem',
      'HDF5Input::ParseProblem',
      'HDF5Input::ParseProblem',
      'UaiMrfInput::ParseProblem',
      'UaiMrfInput::ParseProblem'
      ]

file_name = [
      'srmp_opengm.cpp',
      'srmp_tightening_opengm.cpp',
      'mplp_opengm.cpp',
      'srmp_uai.cpp',
      'mplp_uai.cpp',
      ]


main_body = ["MpRoundingSolver<Solver<" + e[0] + ",LP,StandardTighteningVisitor>> solver(argc,argv);\nsolver.ReadProblem(" + e[1] + "<Solver<" + e[0] + ",LP,StandardTighteningVisitor>>);\nreturn solver.Solve();\n" for e in zip(FMC,parse_fun)]
#main_body = ["LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_MP_ROUNDING(" + e[0] + ",TorresaniEtAlInput::" + e[1] + "<" + e[0] + ">,StandardTighteningVisitor);" for e in zip(FMC,parse_fun)]

for e in zip(main_body,file_name):
   f = open(e[1],'w')
   f.write(preamble)
   f.write("\n{\n");
   f.write(e[0])
   f.write("\n}\n");
   f.close()
