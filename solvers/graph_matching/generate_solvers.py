#!/usr/bin/python

preamble = """
#include "graph_matching.h"
#include "visitors/standard_visitor.hxx"
int main(int argc, char* argv[])
"""

FMC = [ 
      'FMC_MCF<PairwiseConstruction::Left>',
      'FMC_MCF<PairwiseConstruction::Right>',
      'FMC_MCF<PairwiseConstruction::BothSides>',
      'FMC_MP<PairwiseConstruction::Left>',
      'FMC_MP<PairwiseConstruction::Right>',
      'FMC_MP<PairwiseConstruction::Left>',
      'FMC_GM<PairwiseConstruction::Left>',
      'FMC_GM<PairwiseConstruction::Right>',
      'FMC_HUNGARIAN_BP<PairwiseConstruction::Left>',
      'FMC_HUNGARIAN_BP<PairwiseConstruction::Right>',
      'FMC_HUNGARIAN_BP<PairwiseConstruction::BothSides>',
      'FMC_MCF_T<PairwiseConstruction::Left>',
      'FMC_MCF_T<PairwiseConstruction::Right>',
      'FMC_MCF_T<PairwiseConstruction::BothSides>',
      'FMC_MP_T<PairwiseConstruction::Left>',
      'FMC_MP_T<PairwiseConstruction::Right>',
      'FMC_MP_T<PairwiseConstruction::Left>',
      'FMC_GM_T<PairwiseConstruction::Left>',
      'FMC_GM_T<PairwiseConstruction::Right>',
      'FMC_HUNGARIAN_BP_T<PairwiseConstruction::Left>',
      'FMC_HUNGARIAN_BP_T<PairwiseConstruction::Right>',
      'FMC_HUNGARIAN_BP_T<PairwiseConstruction::BothSides>'
      ]

parse_fun = [
      'ParseProblemMCF',
      'ParseProblemMCF',
      'ParseProblemMCF',
      'ParseProblemMP',
      'ParseProblemMP',
      'ParseProblemMP',
      'ParseProblemGM',
      'ParseProblemGM',
      'ParseProblemHungarian',
      'ParseProblemHungarian',
      'ParseProblemHungarian',
      'ParseProblemMCF',
      'ParseProblemMCF',
      'ParseProblemMCF',
      'ParseProblemMP',
      'ParseProblemMP',
      'ParseProblemMP',
      'ParseProblemGM',
      'ParseProblemGM',
      'ParseProblemHungarian',
      'ParseProblemHungarian',
      'ParseProblemHungarian'
      ]

file_name = [
      "graph_matching_via_mcf_left.cpp",
      "graph_matching_via_mcf_right.cpp",
      "graph_matching_via_mcf_both_sides.cpp",
      "graph_matching_via_mp_left.cpp",
      "graph_matching_via_mp_right.cpp",
      "graph_matching_via_mp_both_sides.cpp",
      "graph_matching_via_gm_left.cpp",
      "graph_matching_via_gm_right.cpp",
      "hungarian_bp_left.cpp",
      "hungarian_bp_right.cpp",
      "hungarian_bp_both_sides.cpp",
      "graph_matching_via_mcf_left_tightening.cpp",
      "graph_matching_via_mcf_right_tightening.cpp",
      "graph_matching_via_mcf_both_sides_tightening.cpp",
      "graph_matching_via_mp_left_tightening.cpp",
      "graph_matching_via_mp_right_tightening.cpp",
      "graph_matching_via_mp_both_sides_tightening.cpp",
      "graph_matching_via_gm_left_tightening.cpp",
      "graph_matching_via_gm_right_tightening.cpp",
      "hungarian_bp_left_tightening.cpp",
      "hungarian_bp_right_tightening.cpp",
      "hungarian_bp_both_sides_tightening.cpp"
      ]


main_body = ["MpRoundingSolver<Solver<" + e[0] + ",LP,StandardTighteningVisitor>> solver(argc,argv);\nsolver.ReadProblem(TorresaniEtAlInput::" + e[1] + "<Solver<" + e[0] + ",LP,StandardTighteningVisitor>>);\nreturn solver.Solve();\n" for e in zip(FMC,parse_fun)]
#main_body = ["LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_MP_ROUNDING(" + e[0] + ",TorresaniEtAlInput::" + e[1] + "<" + e[0] + ">,StandardTighteningVisitor);" for e in zip(FMC,parse_fun)]

for e in zip(main_body,file_name):
   f = open(e[1],'w')
   f.write(preamble)
   f.write("\n{\n");
   f.write(e[0])
   f.write("\n}\n");
   f.close()
