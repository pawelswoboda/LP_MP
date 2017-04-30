#include "evaluate.hxx"
#include "visitors/sqlite_visitor.hxx"
#include "solvers/graph_matching/graph_matching.h"

using VisitorType = SqliteVisitor<StandardTighteningVisitor>;

// algorithms considered
using FMC_GM_LEFT = FMC_GM<PairwiseConstruction::Left>;
using FMC_GM_RIGHT = FMC_GM<PairwiseConstruction::Right>;
using FMC_MCF_LEFT = FMC_MCF<PairwiseConstruction::Left>;
using FMC_MCF_RIGHT = FMC_MCF<PairwiseConstruction::Right>;
using FMC_MCF_BOTH_SIDES = FMC_MCF<PairwiseConstruction::BothSides>;
using FMC_MP_LEFT = FMC_MP<PairwiseConstruction::Left>;
using FMC_MP_RIGHT = FMC_MP<PairwiseConstruction::Right>;
using FMC_MP_BOTH_SIDES = FMC_MP<PairwiseConstruction::BothSides>;
using FMC_HUNGARIAN_BP_LEFT = FMC_HUNGARIAN_BP<PairwiseConstruction::Left>;
using FMC_HUNGARIAN_BP_RIGHT = FMC_HUNGARIAN_BP<PairwiseConstruction::Right>;
using FMC_HUNGARIAN_BP_BOTH_SIDES = FMC_HUNGARIAN_BP<PairwiseConstruction::BothSides>;

// algorithms with tightening
using FMC_GM_LEFT_T = FMC_GM_T<PairwiseConstruction::Left>;
using FMC_GM_RIGHT_T = FMC_GM_T<PairwiseConstruction::Right>;
using FMC_MCF_LEFT_T = FMC_MCF_T<PairwiseConstruction::Left>;
using FMC_MCF_RIGHT_T = FMC_MCF_T<PairwiseConstruction::Right>;
using FMC_MCF_BOTH_SIDES_T = FMC_MCF_T<PairwiseConstruction::BothSides>;
using FMC_MP_LEFT_T = FMC_MP_T<PairwiseConstruction::Left>;
using FMC_MP_RIGHT_T = FMC_MP_T<PairwiseConstruction::Right>;
using FMC_MP_BOTH_SIDES_T = FMC_MP_T<PairwiseConstruction::BothSides>;
using FMC_HUNGARIAN_BP_LEFT_T = FMC_HUNGARIAN_BP_T<PairwiseConstruction::Left>;
using FMC_HUNGARIAN_BP_RIGHT_T = FMC_HUNGARIAN_BP_T<PairwiseConstruction::Right>;
using FMC_HUNGARIAN_BP_BOTH_SIDES_T = FMC_HUNGARIAN_BP_T<PairwiseConstruction::BothSides>;

// Torresani et al input
static auto TORRESANI_GM_LEFT_INPUT = TorresaniEtAlInput::ParseProblemGM<Solver<FMC_GM_LEFT,LP,VisitorType>>;
static auto TORRESANI_GM_RIGHT_INPUT = TorresaniEtAlInput::ParseProblemGM<Solver<FMC_GM_RIGHT,LP,VisitorType>>;
static auto TORRESANI_MCF_LEFT_INPUT = TorresaniEtAlInput::ParseProblemMCF<Solver<FMC_MCF_LEFT,LP,VisitorType>>;
static auto TORRESANI_MCF_RIGHT_INPUT = TorresaniEtAlInput::ParseProblemMCF<Solver<FMC_MCF_RIGHT,LP,VisitorType>>;
static auto TORRESANI_MCF_BOTH_SIDES_INPUT = TorresaniEtAlInput::ParseProblemMCF<Solver<FMC_MCF_BOTH_SIDES,LP,VisitorType>>;
static auto TORRESANI_MP_LEFT_INPUT = TorresaniEtAlInput::ParseProblemMP<Solver<FMC_MP_LEFT,LP,VisitorType>>;
static auto TORRESANI_MP_RIGHT_INPUT = TorresaniEtAlInput::ParseProblemMP<Solver<FMC_MP_RIGHT,LP,VisitorType>>;
static auto TORRESANI_MP_BOTH_SIDES_INPUT = TorresaniEtAlInput::ParseProblemMP<Solver<FMC_MP_BOTH_SIDES,LP,VisitorType>>;
static auto TORRESANI_HUNGARIAN_BP_LEFT_INPUT = TorresaniEtAlInput::ParseProblemHungarian<Solver<FMC_HUNGARIAN_BP_LEFT,LP,VisitorType>>;
static auto TORRESANI_HUNGARIAN_BP_RIGHT_INPUT = TorresaniEtAlInput::ParseProblemHungarian<Solver<FMC_HUNGARIAN_BP_RIGHT,LP,VisitorType>>;
static auto TORRESANI_HUNGARIAN_BP_BOTH_SIDES_INPUT = TorresaniEtAlInput::ParseProblemHungarian<Solver<FMC_HUNGARIAN_BP_BOTH_SIDES,LP,VisitorType>>;

static auto TORRESANI_GM_LEFT_INPUT_T = TorresaniEtAlInput::ParseProblemGM<Solver<FMC_GM_LEFT_T,LP,VisitorType>>;
static auto TORRESANI_GM_RIGHT_INPUT_T = TorresaniEtAlInput::ParseProblemGM<Solver<FMC_GM_RIGHT_T,LP,VisitorType>>;
static auto TORRESANI_MCF_LEFT_INPUT_T = TorresaniEtAlInput::ParseProblemMCF<Solver<FMC_MCF_LEFT_T,LP,VisitorType>>;
static auto TORRESANI_MCF_RIGHT_INPUT_T = TorresaniEtAlInput::ParseProblemMCF<Solver<FMC_MCF_RIGHT_T,LP,VisitorType>>;
static auto TORRESANI_MCF_BOTH_SIDES_INPUT_T = TorresaniEtAlInput::ParseProblemMCF<Solver<FMC_MCF_BOTH_SIDES_T,LP,VisitorType>>;
static auto TORRESANI_MP_LEFT_INPUT_T = TorresaniEtAlInput::ParseProblemMP<Solver<FMC_MP_LEFT_T,LP,VisitorType>>;
static auto TORRESANI_MP_RIGHT_INPUT_T = TorresaniEtAlInput::ParseProblemMP<Solver<FMC_MP_RIGHT_T,LP,VisitorType>>;
static auto TORRESANI_MP_BOTH_SIDES_INPUT_T = TorresaniEtAlInput::ParseProblemMP<Solver<FMC_MP_BOTH_SIDES_T,LP,VisitorType>>;
static auto TORRESANI_HUNGARIAN_BP_LEFT_INPUT_T = TorresaniEtAlInput::ParseProblemHungarian<Solver<FMC_HUNGARIAN_BP_LEFT_T,LP,VisitorType>>;
static auto TORRESANI_HUNGARIAN_BP_RIGHT_INPUT_T = TorresaniEtAlInput::ParseProblemHungarian<Solver<FMC_HUNGARIAN_BP_RIGHT_T,LP,VisitorType>>;
static auto TORRESANI_HUNGARIAN_BP_BOTH_SIDES_INPUT_T = TorresaniEtAlInput::ParseProblemHungarian<Solver<FMC_HUNGARIAN_BP_BOTH_SIDES_T,LP,VisitorType>>;

// UAI input
//static auto UAI_GM_INPUT = UaiGraphMatchingInput::ParseProblemGM<FMC_GM<>>;
//static auto UAI_MCF_INPUT = UaiGraphMatchingInput::old_format::ParseProblem<FMC_MCF<>>;
//static auto UAI_MP_INPUT = UaiGraphMatchingInput::ParseProblem<FMC_MP<>>;

//std::vector<std::string> graphMatchingTestDataset = {
//   {"../../solvers/graph_matching/Large_QAP.txt"}
//};

std::string house_prefix = "graph_matching_datasets/tkr_pami13_data/house/";
std::vector<std::string> graphMatchingHouseDatasets = {
   {house_prefix + "energy_house_frame10frame100.txt"},
   {house_prefix + "energy_house_frame10frame95.txt"},
   {house_prefix + "energy_house_frame10frame96.txt"},
   {house_prefix + "energy_house_frame10frame97.txt"},
   {house_prefix + "energy_house_frame10frame98.txt"},
   {house_prefix + "energy_house_frame10frame99.txt"},
   {house_prefix + "energy_house_frame11frame100.txt"},
   {house_prefix + "energy_house_frame11frame101.txt"},
   {house_prefix + "energy_house_frame11frame96.txt"},
   {house_prefix + "energy_house_frame11frame97.txt"},
   {house_prefix + "energy_house_frame11frame98.txt"},
   {house_prefix + "energy_house_frame11frame99.txt"},
   {house_prefix + "energy_house_frame12frame100.txt"},
   {house_prefix + "energy_house_frame12frame101.txt"},
   {house_prefix + "energy_house_frame12frame102.txt"},
   {house_prefix + "energy_house_frame12frame97.txt"},
   {house_prefix + "energy_house_frame12frame98.txt"},
   {house_prefix + "energy_house_frame12frame99.txt"},
   {house_prefix + "energy_house_frame13frame100.txt"},
   {house_prefix + "energy_house_frame13frame101.txt"},
   {house_prefix + "energy_house_frame13frame102.txt"},
   {house_prefix + "energy_house_frame13frame103.txt"},
   {house_prefix + "energy_house_frame13frame98.txt"},
   {house_prefix + "energy_house_frame13frame99.txt"},
   {house_prefix + "energy_house_frame14frame100.txt"},
   {house_prefix + "energy_house_frame14frame101.txt"},
   {house_prefix + "energy_house_frame14frame102.txt"},
   {house_prefix + "energy_house_frame14frame103.txt"},
   {house_prefix + "energy_house_frame14frame104.txt"},
   {house_prefix + "energy_house_frame14frame99.txt"},
   {house_prefix + "energy_house_frame15frame100.txt"},
   {house_prefix + "energy_house_frame15frame101.txt"},
   {house_prefix + "energy_house_frame15frame102.txt"},
   {house_prefix + "energy_house_frame15frame103.txt"},
   {house_prefix + "energy_house_frame15frame104.txt"},
   {house_prefix + "energy_house_frame15frame105.txt"},
   {house_prefix + "energy_house_frame16frame101.txt"},
   {house_prefix + "energy_house_frame16frame102.txt"},
   {house_prefix + "energy_house_frame16frame103.txt"},
   {house_prefix + "energy_house_frame16frame104.txt"},
   {house_prefix + "energy_house_frame16frame105.txt"},
   {house_prefix + "energy_house_frame17frame102.txt"},
   {house_prefix + "energy_house_frame17frame103.txt"},
   {house_prefix + "energy_house_frame17frame104.txt"},
   {house_prefix + "energy_house_frame17frame105.txt"},
   {house_prefix + "energy_house_frame18frame103.txt"},
   {house_prefix + "energy_house_frame18frame104.txt"},
   {house_prefix + "energy_house_frame18frame105.txt"},
   {house_prefix + "energy_house_frame19frame104.txt"},
   {house_prefix + "energy_house_frame19frame105.txt"},
   {house_prefix + "energy_house_frame1frame86.txt"},
   {house_prefix + "energy_house_frame1frame87.txt"},
   {house_prefix + "energy_house_frame1frame88.txt"},
   {house_prefix + "energy_house_frame1frame89.txt"},
   {house_prefix + "energy_house_frame1frame90.txt"},
   {house_prefix + "energy_house_frame1frame91.txt"},
   {house_prefix + "energy_house_frame20frame105.txt"},
   {house_prefix + "energy_house_frame2frame87.txt"},
   {house_prefix + "energy_house_frame2frame88.txt"},
   {house_prefix + "energy_house_frame2frame89.txt"},
   {house_prefix + "energy_house_frame2frame90.txt"},
   {house_prefix + "energy_house_frame2frame91.txt"},
   {house_prefix + "energy_house_frame2frame92.txt"},
   {house_prefix + "energy_house_frame3frame88.txt"},
   {house_prefix + "energy_house_frame3frame89.txt"},
   {house_prefix + "energy_house_frame3frame90.txt"},
   {house_prefix + "energy_house_frame3frame91.txt"},
   {house_prefix + "energy_house_frame3frame92.txt"},
   {house_prefix + "energy_house_frame3frame93.txt"},
   {house_prefix + "energy_house_frame4frame89.txt"},
   {house_prefix + "energy_house_frame4frame90.txt"},
   {house_prefix + "energy_house_frame4frame91.txt"},
   {house_prefix + "energy_house_frame4frame92.txt"},
   {house_prefix + "energy_house_frame4frame93.txt"},
   {house_prefix + "energy_house_frame4frame94.txt"},
   {house_prefix + "energy_house_frame5frame90.txt"},
   {house_prefix + "energy_house_frame5frame91.txt"},
   {house_prefix + "energy_house_frame5frame92.txt"},
   {house_prefix + "energy_house_frame5frame93.txt"},
   {house_prefix + "energy_house_frame5frame94.txt"},
   {house_prefix + "energy_house_frame5frame95.txt"},
   {house_prefix + "energy_house_frame6frame91.txt"},
   {house_prefix + "energy_house_frame6frame92.txt"},
   {house_prefix + "energy_house_frame6frame93.txt"},
   {house_prefix + "energy_house_frame6frame94.txt"},
   {house_prefix + "energy_house_frame6frame95.txt"},
   {house_prefix + "energy_house_frame6frame96.txt"},
   {house_prefix + "energy_house_frame7frame92.txt"},
   {house_prefix + "energy_house_frame7frame93.txt"},
   {house_prefix + "energy_house_frame7frame94.txt"},
   {house_prefix + "energy_house_frame7frame95.txt"},
   {house_prefix + "energy_house_frame7frame96.txt"},
   {house_prefix + "energy_house_frame7frame97.txt"},
   {house_prefix + "energy_house_frame8frame93.txt"},
   {house_prefix + "energy_house_frame8frame94.txt"},
   {house_prefix + "energy_house_frame8frame95.txt"},
   {house_prefix + "energy_house_frame8frame96.txt"},
   {house_prefix + "energy_house_frame8frame97.txt"},
   {house_prefix + "energy_house_frame8frame98.txt"},
   {house_prefix + "energy_house_frame9frame94.txt"},
   {house_prefix + "energy_house_frame9frame95.txt"},
   {house_prefix + "energy_house_frame9frame96.txt"},
   {house_prefix + "energy_house_frame9frame97.txt"},
   {house_prefix + "energy_house_frame9frame98.txt"},
   {house_prefix + "energy_house_frame9frame99.txt"}
};

std::string hotel_prefix = "graph_matching_datasets/tkr_pami13_data/hotel/";
std::vector<std::string> graphMatchingHotelDatasets = {
   {hotel_prefix + "energy_hotel_frame15frame22.txt"},
   {hotel_prefix + "energy_hotel_frame15frame29.txt"},
   {hotel_prefix + "energy_hotel_frame15frame36.txt"},
   {hotel_prefix + "energy_hotel_frame15frame43.txt"},
   {hotel_prefix + "energy_hotel_frame15frame50.txt"},
   {hotel_prefix + "energy_hotel_frame15frame57.txt"},
   {hotel_prefix + "energy_hotel_frame15frame64.txt"},
   {hotel_prefix + "energy_hotel_frame15frame71.txt"},
   {hotel_prefix + "energy_hotel_frame15frame78.txt"},
   {hotel_prefix + "energy_hotel_frame15frame85.txt"},
   {hotel_prefix + "energy_hotel_frame15frame92.txt"},
   {hotel_prefix + "energy_hotel_frame15frame99.txt"},
   {hotel_prefix + "energy_hotel_frame1frame15.txt"},
   {hotel_prefix + "energy_hotel_frame1frame22.txt"},
   {hotel_prefix + "energy_hotel_frame1frame29.txt"},
   {hotel_prefix + "energy_hotel_frame1frame36.txt"},
   {hotel_prefix + "energy_hotel_frame1frame43.txt"},
   {hotel_prefix + "energy_hotel_frame1frame50.txt"},
   {hotel_prefix + "energy_hotel_frame1frame57.txt"},
   {hotel_prefix + "energy_hotel_frame1frame64.txt"},
   {hotel_prefix + "energy_hotel_frame1frame71.txt"},
   {hotel_prefix + "energy_hotel_frame1frame78.txt"},
   {hotel_prefix + "energy_hotel_frame1frame85.txt"},
   {hotel_prefix + "energy_hotel_frame1frame8.txt"},
   {hotel_prefix + "energy_hotel_frame1frame92.txt"},
   {hotel_prefix + "energy_hotel_frame1frame99.txt"},
   {hotel_prefix + "energy_hotel_frame22frame29.txt"},
   {hotel_prefix + "energy_hotel_frame22frame36.txt"},
   {hotel_prefix + "energy_hotel_frame22frame43.txt"},
   {hotel_prefix + "energy_hotel_frame22frame50.txt"},
   {hotel_prefix + "energy_hotel_frame22frame57.txt"},
   {hotel_prefix + "energy_hotel_frame22frame64.txt"},
   {hotel_prefix + "energy_hotel_frame22frame71.txt"},
   {hotel_prefix + "energy_hotel_frame22frame78.txt"},
   {hotel_prefix + "energy_hotel_frame22frame85.txt"},
   {hotel_prefix + "energy_hotel_frame22frame92.txt"},
   {hotel_prefix + "energy_hotel_frame22frame99.txt"},
   {hotel_prefix + "energy_hotel_frame29frame36.txt"},
   {hotel_prefix + "energy_hotel_frame29frame43.txt"},
   {hotel_prefix + "energy_hotel_frame29frame50.txt"},
   {hotel_prefix + "energy_hotel_frame29frame57.txt"},
   {hotel_prefix + "energy_hotel_frame29frame64.txt"},
   {hotel_prefix + "energy_hotel_frame29frame71.txt"},
   {hotel_prefix + "energy_hotel_frame29frame78.txt"},
   {hotel_prefix + "energy_hotel_frame29frame85.txt"},
   {hotel_prefix + "energy_hotel_frame29frame92.txt"},
   {hotel_prefix + "energy_hotel_frame29frame99.txt"},
   {hotel_prefix + "energy_hotel_frame36frame43.txt"},
   {hotel_prefix + "energy_hotel_frame36frame50.txt"},
   {hotel_prefix + "energy_hotel_frame36frame57.txt"},
   {hotel_prefix + "energy_hotel_frame36frame64.txt"},
   {hotel_prefix + "energy_hotel_frame36frame71.txt"},
   {hotel_prefix + "energy_hotel_frame36frame78.txt"},
   {hotel_prefix + "energy_hotel_frame36frame85.txt"},
   {hotel_prefix + "energy_hotel_frame36frame92.txt"},
   {hotel_prefix + "energy_hotel_frame36frame99.txt"},
   {hotel_prefix + "energy_hotel_frame43frame50.txt"},
   {hotel_prefix + "energy_hotel_frame43frame57.txt"},
   {hotel_prefix + "energy_hotel_frame43frame64.txt"},
   {hotel_prefix + "energy_hotel_frame43frame71.txt"},
   {hotel_prefix + "energy_hotel_frame43frame78.txt"},
   {hotel_prefix + "energy_hotel_frame43frame85.txt"},
   {hotel_prefix + "energy_hotel_frame43frame92.txt"},
   {hotel_prefix + "energy_hotel_frame43frame99.txt"},
   {hotel_prefix + "energy_hotel_frame50frame57.txt"},
   {hotel_prefix + "energy_hotel_frame50frame64.txt"},
   {hotel_prefix + "energy_hotel_frame50frame71.txt"},
   {hotel_prefix + "energy_hotel_frame50frame78.txt"},
   {hotel_prefix + "energy_hotel_frame50frame85.txt"},
   {hotel_prefix + "energy_hotel_frame50frame92.txt"},
   {hotel_prefix + "energy_hotel_frame50frame99.txt"},
   {hotel_prefix + "energy_hotel_frame57frame64.txt"},
   {hotel_prefix + "energy_hotel_frame57frame71.txt"},
   {hotel_prefix + "energy_hotel_frame57frame78.txt"},
   {hotel_prefix + "energy_hotel_frame57frame85.txt"},
   {hotel_prefix + "energy_hotel_frame57frame92.txt"},
   {hotel_prefix + "energy_hotel_frame57frame99.txt"},
   {hotel_prefix + "energy_hotel_frame64frame71.txt"},
   {hotel_prefix + "energy_hotel_frame64frame78.txt"},
   {hotel_prefix + "energy_hotel_frame64frame85.txt"},
   {hotel_prefix + "energy_hotel_frame64frame92.txt"},
   {hotel_prefix + "energy_hotel_frame64frame99.txt"},
   {hotel_prefix + "energy_hotel_frame71frame78.txt"},
   {hotel_prefix + "energy_hotel_frame71frame85.txt"},
   {hotel_prefix + "energy_hotel_frame71frame92.txt"},
   {hotel_prefix + "energy_hotel_frame71frame99.txt"},
   {hotel_prefix + "energy_hotel_frame78frame85.txt"},
   {hotel_prefix + "energy_hotel_frame78frame92.txt"},
   {hotel_prefix + "energy_hotel_frame78frame99.txt"},
   {hotel_prefix + "energy_hotel_frame85frame92.txt"},
   {hotel_prefix + "energy_hotel_frame85frame99.txt"},
   {hotel_prefix + "energy_hotel_frame8frame15.txt"},
   {hotel_prefix + "energy_hotel_frame8frame22.txt"},
   {hotel_prefix + "energy_hotel_frame8frame29.txt"},
   {hotel_prefix + "energy_hotel_frame8frame36.txt"},
   {hotel_prefix + "energy_hotel_frame8frame43.txt"},
   {hotel_prefix + "energy_hotel_frame8frame50.txt"},
   {hotel_prefix + "energy_hotel_frame8frame57.txt"},
   {hotel_prefix + "energy_hotel_frame8frame64.txt"},
   {hotel_prefix + "energy_hotel_frame8frame71.txt"},
   {hotel_prefix + "energy_hotel_frame8frame78.txt"},
   {hotel_prefix + "energy_hotel_frame8frame85.txt"},
   {hotel_prefix + "energy_hotel_frame8frame92.txt"},
   {hotel_prefix + "energy_hotel_frame8frame99.txt"},
   {hotel_prefix + "energy_hotel_frame92frame99.txt"}
};

std::string hassan_prefix = "../../../solvers/graph_matching/Hassan/";
std::vector<std::string> graphMatchingHassanDatasets = {
   {hassan_prefix + "board_torresani.txt"},
   {hassan_prefix + "books_torresani.txt"},
   {hassan_prefix + "hammer_torresani.txt"},
   {hassan_prefix + "party_torresani.txt"},
   {hassan_prefix + "table_torresani.txt"},
   //{hassan_prefix + "tea_torresani.txt"}, 
   {hassan_prefix + "walking_torresani.txt"}
};

std::string worms_prefix = "../../../solvers/graph_matching/graph_matching_datasets/allWorms-16-03-11-1745-dd/";
std::vector<std::string> graphMatchingWormsDatasets = {
   { worms_prefix + "C18G1_2L1_1-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "cnd1threeL1_1213061-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "cnd1threeL1_1228061-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "cnd1threeL1_1229061-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "cnd1threeL1_1229062-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "cnd1threeL1_1229063-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "eft3RW10035L1_0125071-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "eft3RW10035L1_0125072-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "eft3RW10035L1_0125073-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "egl5L1_0606074-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "elt3L1_0503071-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "elt3L1_0503072-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "elt3L1_0504073-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "hlh1fourL1_0417071-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "hlh1fourL1_0417075-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "hlh1fourL1_0417076-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "hlh1fourL1_0417077-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "hlh1fourL1_0417078-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "mir61L1_1228061-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "mir61L1_1228062-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "mir61L1_1229062-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "pha4A7L1_1213061-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "pha4A7L1_1213062-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "pha4A7L1_1213064-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "pha4B2L1_0125072-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "pha4I2L_0408071-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "pha4I2L_0408072-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "pha4I2L_0408073-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "unc54L1_0123071-lowThresh-more-hyp.surf-16-03-11-1745.dd"},
   { worms_prefix + "unc54L1_0123072-lowThresh-more-hyp.surf-16-03-11-1745.dd"}
};

/*
std::string worms_prefix = "../../../solvers/graph_matching/graph_matching_datasets/allWorms-03-04-1750-uai/";
std::vector<std::string> graphMatchingWormsDatasets = {
   { worms_prefix + "C18G1_2L1_1-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "cnd1threeL1_1213061-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "cnd1threeL1_1228061-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "cnd1threeL1_1229061-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "cnd1threeL1_1229062-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "cnd1threeL1_1229063-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "eft3RW10035L1_0125071-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "eft3RW10035L1_0125072-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "eft3RW10035L1_0125073-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "egl5L1_0606074-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "elt3L1_0503071-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "elt3L1_0503072-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "elt3L1_0504073-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "hlh1fourL1_0417071-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "hlh1fourL1_0417075-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "hlh1fourL1_0417076-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "hlh1fourL1_0417077-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "hlh1fourL1_0417078-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "mir61L1_1228061-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "mir61L1_1228062-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "mir61L1_1229062-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "pha4A7L1_1213061-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "pha4A7L1_1213062-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "pha4A7L1_1213064-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "pha4B2L1_0125072-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "pha4I2L_0408071-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "pha4I2L_0408072-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "pha4I2L_0408073-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "unc54L1_0123071-lowThresh-more-hyp.surf-16-03-04-1750.uai"},
   { worms_prefix + "unc54L1_0123072-lowThresh-more-hyp.surf-16-03-04-1750.uai"}
};
*/

std::string car_prefix = "../../../solvers/graph_matching/graph_matching_datasets/car/";
std::vector<std::string> graphMatchingCarDatasets = {
   { car_prefix + "car1.txt" },
   { car_prefix + "car2.txt" },
   { car_prefix + "car3.txt" },
   { car_prefix + "car4.txt" },
   { car_prefix + "car5.txt" },
   { car_prefix + "car6.txt" },
   { car_prefix + "car7.txt" },
   { car_prefix + "car8.txt" },
   { car_prefix + "car9.txt" },
   { car_prefix + "car10.txt" },
   { car_prefix + "car11.txt" },
   { car_prefix + "car12.txt" },
   { car_prefix + "car13.txt" },
   { car_prefix + "car14.txt" },
   { car_prefix + "car15.txt" },
   { car_prefix + "car16.txt" },
   { car_prefix + "car17.txt" },
   { car_prefix + "car18.txt" },
   { car_prefix + "car19.txt" },
   { car_prefix + "car20.txt" },
   { car_prefix + "car21.txt" },
   { car_prefix + "car22.txt" },
   { car_prefix + "car23.txt" },
   { car_prefix + "car24.txt" },
   { car_prefix + "car25.txt" },
   { car_prefix + "car26.txt" },
   { car_prefix + "car27.txt" },
   { car_prefix + "car28.txt" },
   { car_prefix + "car29.txt" },
   { car_prefix + "car30.txt" }
};

std::string motor_prefix = "../../../solvers/graph_matching/graph_matching_datasets/motor/";
std::vector<std::string> graphMatchingMotorDatasets = {
   { motor_prefix + "motor1.txt" },
   { motor_prefix + "motor2.txt" },
   { motor_prefix + "motor3.txt" },
   { motor_prefix + "motor4.txt" },
   { motor_prefix + "motor5.txt" },
   { motor_prefix + "motor6.txt" },
   { motor_prefix + "motor7.txt" },
   { motor_prefix + "motor8.txt" },
   { motor_prefix + "motor9.txt" },
   { motor_prefix + "motor10.txt" },
   { motor_prefix + "motor11.txt" },
   { motor_prefix + "motor12.txt" },
   { motor_prefix + "motor13.txt" },
   { motor_prefix + "motor14.txt" },
   { motor_prefix + "motor15.txt" },
   { motor_prefix + "motor16.txt" },
   { motor_prefix + "motor17.txt" },
   { motor_prefix + "motor18.txt" },
   { motor_prefix + "motor19.txt" },
   { motor_prefix + "motor20.txt" }
};

int main()
{
   // emulate command line options
   std::vector<std::string> options = {
      {"--maxIter"}, {"1000"},
      {"--timeout"}, {"3600"}, // one hour
      {"--minDualImprovement"}, {"0.001"},
      {"--minDualImprovementInterval"}, {"20"},
      {"--lowerBoundComputationInterval"}, {"10"},
      {"--primalComputationInterval"}, {"10"},
      {"--overwriteDbRecord"}, // do zrobienia: possibly deactivate this. Then we do not overwrite
      {"--databaseFile"}, {"graph_matching.db"}
   };

   std::vector<std::string> tightening_options = options;
   tightening_options.push_back("--tighten");
   tightening_options.push_back("--tightenIteration");
   tightening_options.push_back("700");
   tightening_options.push_back("--tightenInterval");
   tightening_options.push_back("20");
   tightening_options.push_back("--tightenConstraintsPercentage");
   tightening_options.push_back("0.1");
   tightening_options.push_back("--tightenReparametrization");
   tightening_options.push_back("damped_uniform");
   tightening_options.push_back("--tightenMinDualImprovement");
   tightening_options.push_back("0.02");
   tightening_options.push_back("--tightenMinDualImprovementInterval");
   tightening_options.push_back("20");

   std::vector<std::string> uniform_options = tightening_options;
   uniform_options.push_back("--standardReparametrization");
   uniform_options.push_back("uniform");
   uniform_options.push_back("--roundingReparametrization");
   uniform_options.push_back("uniform");

   std::vector<std::string> anisotropic_options = tightening_options;
   anisotropic_options.push_back("--standardReparametrization");
   anisotropic_options.push_back("anisotropic");
   anisotropic_options.push_back("--roundingReparametrization");
   anisotropic_options.push_back("anisotropic");



   std::vector<std::string> mcf_options = anisotropic_options;

   std::vector<std::string> mp_options = anisotropic_options;

   std::vector<std::string> gm_options = anisotropic_options;

   std::vector<std::string> hbp_options = uniform_options;


   std::vector<std::string> graphMatchingDatasets;
   //for(auto& d : graphMatchingHotelDatasets) {graphMatchingDatasets.push_back(d);}
   //for(auto& d : graphMatchingHouseDatasets) {graphMatchingDatasets.push_back(d);}
   //graphMatchingDatasets.push_back(graphMatchingHotelDatasets[54]);


   // small evaluation for overview paper
   {
      // attention: better results can be obtained by choosing right side for house and hotel, but then primal rounding does not work currently for mcf
      //std::vector<std::string> hotel_dataset = {{hotel_prefix + "energy_hotel_frame1frame64.txt"}};
      //RunSolver<FMC_MCF_LEFT,VisitorSolver<MpRoundingSolver<Solver<FMC_MCF_LEFT,LP>,VisitorType>>(TORRESANI_MCF_LEFT_INPUT, {graphMatchingCarDatasets[0]}, anisotropic_options,"car","AMCF-O");
      //RunSolver<FMC_MP_LEFT,VisitorSolver<MpRoundingSolver<Solver<FMC_MP_LEFT>,VisitorType>>(TORRESANI_MP_LEFT_INPUT, {graphMatchingCarDatasets[0]}, anisotropic_options,"car","AMP-O");
      //RunSolver<FMC_HUNGARIAN_BP_LEFT,VisitorSolver<MpRoundingSolver<Solver<FMC_HUNGARIAN_BP_LEFT>,VisitorType>>(TORRESANI_HUNGARIAN_BP_LEFT_INPUT, graphMatchingCarDatasets, uniform_options,"car","HUNGARIAN-BP-O");
      //RunSolver<FMC_GM_LEFT,VisitorSolver<MpRoundingSolver<Solver<FMC_GM_LEFT>,VisitorType>>(TORRESANI_GM_LEFT_INPUT, graphMatchingCarDatasets, anisotropic_options,"car","GM-O");

      //std::vector<std::string> car_dataset = {{car_prefix + "car6.txt"}};
      //RunSolver<FMC_MCF_LEFT,VisitorSolver<MpRoundingSolver<Solver<FMC_MCF_LEFT>,VisitorType>>(TORRESANI_MCF_LEFT_INPUT, graphMatchingMotorDatasets, anisotropic_options,"motor","AMCF-O");
      //RunSolver<FMC_MP_LEFT,VisitorSolver<MpRoundingSolver<Solver<FMC_MP_LEFT>,VisitorType>>(TORRESANI_MP_LEFT_INPUT, graphMatchingMotorDatasets, anisotropic_options,"motor","AMP-O");
      //RunSolver<FMC_HUNGARIAN_BP_LEFT,VisitorSolver<MpRoundingSolver<Solver<FMC_HUNGARIAN_BP_LEFT>,VisitorType>>(TORRESANI_HUNGARIAN_BP_LEFT_INPUT, graphMatchingMotorDatasets, uniform_options,"motor","HUNGARIAN-BP-O");
      //RunSolver<FMC_GM_LEFT,VisitorSolver<MpRoundingSolver<Solver<FMC_GM_LEFT>,VisitorType>>(TORRESANI_GM_LEFT_INPUT, graphMatchingMotorDatasets, anisotropic_options,"motor","GM-O");

      //std::vector<std::string> worms_dataset = {{worms_prefix + "C18G1_2L1_1-lowThresh-more-hyp.surf-16-03-11-1745.dd"}};
      //RunSolver<FMC_MCF_LEFT,VisitorSolver<MpRoundingSolver<Solver<FMC_MCF_LEFT>,VisitorType>>(TORRESANI_MCF_LEFT_INPUT, graphMatchingWormsDatasets, anisotropic_options,"worms","AMCF-O");
      //RunSolver<FMC_MP_LEFT,VisitorSolver<MpRoundingSolver<Solver<FMC_MP_LEFT>,VisitorType>>(TORRESANI_MP_LEFT_INPUT, graphMatchingWormsDatasets, anisotropic_options,"worms","AMP-O");
      //RunSolver<FMC_HUNGARIAN_BP_LEFT,VisitorSolver<MpRoundingSolver<Solver<FMC_HUNGARIAN_BP_LEFT>,VisitorType>>(TORRESANI_HUNGARIAN_BP_LEFT_INPUT, graphMatchingWormsDatasets, uniform_options,"worms","HUNGARIAN-BP-O");
      //RunSolver<FMC_GM_LEFT,VisitorSolver<MpRoundingSolver<Solver<FMC_GM_LEFT>,VisitorType>>(TORRESANI_GM_LEFT_INPUT, graphMatchingWormsDatasets, anisotropic_options,"worms","GM-O");
   }

   {
   RunSolver<FMC_MP_BOTH_SIDES_T, MpRoundingSolver<Solver<FMC_MP_BOTH_SIDES_T,LP,VisitorType>>>(TORRESANI_MP_BOTH_SIDES_INPUT_T,graphMatchingHotelDatasets,mp_options,"hotel","AMP-B");
   RunSolver<FMC_MP_BOTH_SIDES_T, MpRoundingSolver<Solver<FMC_MP_BOTH_SIDES_T,LP,VisitorType>>>(TORRESANI_MP_BOTH_SIDES_INPUT_T,graphMatchingHouseDatasets,mp_options,"house","AMP-B");
   //RunSolver<FMC_MP_BOTH_SIDES_T, MpRoundingSolver<Solver<FMC_MP_BOTH_SIDES_T,LP,VisitorType>>(TORRESANI_MP_BOTH_SIDES_INPUT_T,graphMatchingHassanDatasets,mp_options,"Hassan","AMP-B");

   //RunSolver<FMC_MP_LEFT_T, MpRoundingSolver<Solver<FMC_MP_LEFT_T,LP,VisitorType>>(TORRESANI_MP_LEFT_INPUT_T,graphMatchingHotelDatasets,mp_options,"hotel","AMP-O");
   //RunSolver<FMC_MP_LEFT_T, MpRoundingSolver<Solver<FMC_MP_LEFT_T,LP,VisitorType>>(TORRESANI_MP_LEFT_INPUT_T,graphMatchingHouseDatasets,mp_options,"house","AMP-O");
//   RunSolver<FMC_MP_LEFT_T, MpRoundingSolver<Solver<FMC_MP_LEFT_T,LP,VisitorType>>(TORRESANI_MP_LEFT_INPUT_T,graphMatchingHassanDatasets,mp_options,"Hassan","AMP-O");

   //RunSolver<FMC_MP_RIGHT_T, MpRoundingSolver<Solver<FMC_MP_RIGHT_T,LP,VisitorType>>(TORRESANI_MP_RIGHT_INPUT_T,graphMatchingHotelDatasets,mp_options,"hotel","AMP-I");
   //RunSolver<FMC_MP_RIGHT_T, MpRoundingSolver<Solver<FMC_MP_RIGHT_T,LP,VisitorType>>(TORRESANI_MP_RIGHT_INPUT_T,graphMatchingHouseDatasets,mp_options,"house","AMP-I");
   //RunSolver<FMC_MP_RIGHT_T, MpRoundingSolver<Solver<FMC_MP_RIGHT_T,LP,VisitorType>>(TORRESANI_MP_RIGHT_INPUT_T,graphMatchingHassanDatasets,mp_options,"Hassan","AMP-I");

   RunSolver<FMC_MCF_BOTH_SIDES_T, MpRoundingSolver<Solver<FMC_MCF_BOTH_SIDES_T,LP,VisitorType>>>(TORRESANI_MCF_BOTH_SIDES_INPUT_T,graphMatchingHotelDatasets,mcf_options,"hotel","AMCF-B");
   RunSolver<FMC_MCF_BOTH_SIDES_T, MpRoundingSolver<Solver<FMC_MCF_BOTH_SIDES_T,LP,VisitorType>>>(TORRESANI_MCF_BOTH_SIDES_INPUT_T,graphMatchingHouseDatasets,mcf_options,"house","AMCF-B");
   //RunSolver<FMC_MCF_BOTH_SIDES_T, MpRoundingSolver<Solver<FMC_MCF_BOTH_SIDES_T,LP,VisitorType>>(TORRESANI_MCF_BOTH_SIDES_INPUT_T,graphMatchingHassanDatasets,mcf_options,"Hassan","AMCF-B");

   //RunSolver<FMC_MCF_LEFT_T, MpRoundingSolver<Solver<FMC_MCF_LEFT_T,LP,VisitorType>>(TORRESANI_MCF_LEFT_INPUT_T,graphMatchingHotelDatasets,mcf_options,"hotel","AMCF-O");
   //RunSolver<FMC_MCF_LEFT_T, MpRoundingSolver<Solver<FMC_MCF_LEFT_T,LP,VisitorType>>(TORRESANI_MCF_LEFT_INPUT_T,graphMatchingHouseDatasets,mcf_options,"house","AMCF-O");
//   RunSolver<FMC_MCF_LEFT_T, MpRoundingSolver<Solver<FMC_MCF_LEFT_T,LP,VisitorType>>(TORRESANI_MCF_LEFT_INPUT_T,graphMatchingHassanDatasets,mcf_options,"Hassan","AMCF-O");

   //RunSolver<FMC_MCF_RIGHT_T, MpRoundingSolver<Solver<FMC_MCF_RIGHT_T,LP,VisitorType>>(TORRESANI_MCF_RIGHT_INPUT_T,graphMatchingHotelDatasets,mcf_options,"hotel","AMCF-I");
   //RunSolver<FMC_MCF_RIGHT_T, MpRoundingSolver<Solver<FMC_MCF_RIGHT_T,LP,VisitorType>>(TORRESANI_MCF_RIGHT_INPUT_T,graphMatchingHouseDatasets,mcf_options,"house","AMCF-I");
   //RunSolver<FMC_MCF_RIGHT_T, MpRoundingSolver<Solver<FMC_MCF_RIGHT_T,LP,VisitorType>>(TORRESANI_MCF_RIGHT_INPUT_T,graphMatchingHassanDatasets,mcf_options,"Hassan","AMCF-I");
   
   RunSolver<FMC_GM_LEFT_T, MpRoundingSolver<Solver<FMC_GM_LEFT_T,LP,VisitorType>>>(TORRESANI_GM_LEFT_INPUT_T,graphMatchingHotelDatasets,gm_options,"hotel","GM-O");
   RunSolver<FMC_GM_LEFT_T, MpRoundingSolver<Solver<FMC_GM_LEFT_T,LP,VisitorType>>>(TORRESANI_GM_LEFT_INPUT_T,graphMatchingHouseDatasets,gm_options,"house","GM-O");
//   RunSolver<FMC_GM_LEFT_T, MpRoundingSolver<Solver<FMC_GM_LEFT_T,LP,VisitorType>>(TORRESANI_GM_LEFT_INPUT_T,graphMatchingHassanDatasets,gm_options,"Hassan","GM-O");
   
   //RunSolver<FMC_GM_RIGHT_T, MpRoundingSolver<Solver<FMC_GM_RIGHT_T,LP,VisitorType>>(TORRESANI_GM_RIGHT_INPUT_T,graphMatchingHotelDatasets,gm_options,"hotel","GM-I");
   //RunSolver<FMC_GM_RIGHT_T, MpRoundingSolver<Solver<FMC_GM_RIGHT_T,LP,VisitorType>>(TORRESANI_GM_RIGHT_INPUT_T,graphMatchingHouseDatasets,gm_options,"house","GM-I");
   //RunSolver<FMC_GM_RIGHT_T, MpRoundingSolver<Solver<FMC_GM_RIGHT_T,LP,VisitorType>>(TORRESANI_GM_RIGHT_INPUT_T,graphMatchingHassanDatasets,gm_options,"Hassan","GM-I");

   //RunSolver<FMC_HUNGARIAN_BP_RIGHT_T, MpRoundingSolver<Solver<FMC_HUNGARIAN_BP_RIGHT_T,LP,VisitorType>>(TORRESANI_HUNGARIAN_BP_RIGHT_INPUT_T,graphMatchingHotelDatasets,hbp_options,"hotel","HUNGARIAN_BP-I");
   //RunSolver<FMC_HUNGARIAN_BP_RIGHT_T, MpRoundingSolver<Solver<FMC_HUNGARIAN_BP_RIGHT_T,LP,VisitorType>>(TORRESANI_HUNGARIAN_BP_RIGHT_INPUT_T,graphMatchingHouseDatasets,hbp_options,"house","HUNGARIAN_BP-I");
   //RunSolver<FMC_HUNGARIAN_BP_RIGHT_T, MpRoundingSolver<Solver<FMC_HUNGARIAN_BP_RIGHT_T,LP,VisitorType>>(TORRESANI_HUNGARIAN_BP_RIGHT_INPUT_T,graphMatchingHassanDatasets,hbp_options,"Hassan","HUNGARIAN_BP-I");

   RunSolver<FMC_HUNGARIAN_BP_BOTH_SIDES_T, MpRoundingSolver<Solver<FMC_HUNGARIAN_BP_BOTH_SIDES_T,LP,VisitorType>>>(TORRESANI_HUNGARIAN_BP_BOTH_SIDES_INPUT_T,graphMatchingHotelDatasets,hbp_options,"hotel","HUNGARIAN_BP-B");
   RunSolver<FMC_HUNGARIAN_BP_BOTH_SIDES_T, MpRoundingSolver<Solver<FMC_HUNGARIAN_BP_BOTH_SIDES_T,LP,VisitorType>>>(TORRESANI_HUNGARIAN_BP_BOTH_SIDES_INPUT_T,graphMatchingHouseDatasets,hbp_options,"house","HUNGARIAN_BP-B");
   //RunSolver<FMC_HUNGARIAN_BP_BOTH_SIDES_T, MpRoundingSolver<Solver<FMC_HUNGARIAN_BP_BOTH_SIDES_T,LP,VisitorType>>(TORRESANI_HUNGARIAN_BP_BOTH_SIDES_INPUT_T,graphMatchingHassanDatasets,hbp_options,"Hassan","HUNGARIAN_BP-B");

   //RunSolver<FMC_HUNGARIAN_BP_LEFT_T, MpRoundingSolver<Solver<FMC_HUNGARIAN_BP_LEFT_T,LP,VisitorType>>(TORRESANI_HUNGARIAN_BP_LEFT_INPUT_T,graphMatchingHotelDatasets,hbp_options,"hotel","HUNGARIAN_BP-O");
   //RunSolver<FMC_HUNGARIAN_BP_LEFT_T, MpRoundingSolver<Solver<FMC_HUNGARIAN_BP_LEFT_T,LP,VisitorType>>(TORRESANI_HUNGARIAN_BP_LEFT_INPUT_T,graphMatchingHouseDatasets,hbp_options,"house","HUNGARIAN_BP-O");
//   RunSolver<FMC_HUNGARIAN_BP_LEFT_T, MpRoundingSolver<Solver<FMC_HUNGARIAN_BP_LEFT_T,LP,VisitorType>>(TORRESANI_HUNGARIAN_BP_LEFT_INPUT_T,graphMatchingHassanDatasets,hbp_options,"Hassan","HUNGARIAN_BP-O");

   //RunSolver<FMC_MP_LEFT_T, MpRoundingSolver<Solver<FMC_MP_LEFT_T,LP,VisitorType>>(TORRESANI_MP_LEFT_INPUT_T,graphMatchingCarDatasets,mp_options,"car","AMP-O");
   //RunSolver<FMC_MCF_LEFT_T, MpRoundingSolver<Solver<FMC_MCF_LEFT_T,LP,VisitorType>>(TORRESANI_MCF_LEFT_INPUT_T,graphMatchingCarDatasets,mcf_options,"car","AMCF-O");
   //RunSolver<FMC_GM_LEFT_T, MpRoundingSolver<Solver<FMC_GM_LEFT_T,LP,VisitorType>>(TORRESANI_GM_LEFT_INPUT_T,graphMatchingCarDatasets,gm_options,"car","GM-O");
   //RunSolver<FMC_HUNGARIAN_BP_LEFT_T, MpRoundingSolver<Solver<FMC_HUNGARIAN_BP_LEFT_T,LP,VisitorType>>(TORRESANI_HUNGARIAN_BP_LEFT_INPUT_T,graphMatchingCarDatasets,hbp_options,"car","HUNGARIAN_BP-O");

   //RunSolver<FMC_MP_LEFT_T, MpRoundingSolver<Solver<FMC_MP_LEFT_T,LP,VisitorType>>(TORRESANI_MP_LEFT_INPUT_T,graphMatchingMotorDatasets,mp_options,"motor","AMP-O");
   //RunSolver<FMC_MCF_LEFT_T, MpRoundingSolver<Solver<FMC_MCF_LEFT_T,LP,VisitorType>>(TORRESANI_MCF_LEFT_INPUT_T,graphMatchingMotorDatasets,mcf_options,"motor","AMCF-O");
   //RunSolver<FMC_GM_LEFT_T, MpRoundingSolver<Solver<FMC_GM_LEFT_T,LP,VisitorType>>(TORRESANI_GM_LEFT_INPUT_T,graphMatchingMotorDatasets,gm_options,"motor","GM-O");
   //RunSolver<FMC_HUNGARIAN_BP_LEFT_T, MpRoundingSolver<Solver<FMC_HUNGARIAN_BP_LEFT_T,LP,VisitorType>>(TORRESANI_HUNGARIAN_BP_LEFT_INPUT_T,graphMatchingMotorDatasets,hbp_options,"motor","HUNGARIAN_BP-O");

   RunSolver<FMC_MP_BOTH_SIDES_T, MpRoundingSolver<Solver<FMC_MP_BOTH_SIDES_T,LP,VisitorType>>>(TORRESANI_MP_BOTH_SIDES_INPUT_T,graphMatchingCarDatasets,mp_options,"car","AMP-B");
   RunSolver<FMC_MCF_BOTH_SIDES_T, MpRoundingSolver<Solver<FMC_MCF_BOTH_SIDES_T,LP,VisitorType>>>(TORRESANI_MCF_BOTH_SIDES_INPUT_T,graphMatchingCarDatasets,mcf_options,"car","AMCF-B");
   RunSolver<FMC_GM_LEFT_T, MpRoundingSolver<Solver<FMC_GM_LEFT_T,LP,VisitorType>>>(TORRESANI_GM_LEFT_INPUT_T,graphMatchingCarDatasets,gm_options,"car","GM-B");
   RunSolver<FMC_HUNGARIAN_BP_BOTH_SIDES_T, MpRoundingSolver<Solver<FMC_HUNGARIAN_BP_BOTH_SIDES_T,LP,VisitorType>>>(TORRESANI_HUNGARIAN_BP_BOTH_SIDES_INPUT_T,graphMatchingCarDatasets,hbp_options,"car","HUNGARIAN_BP-B");

   RunSolver<FMC_MP_BOTH_SIDES_T, MpRoundingSolver<Solver<FMC_MP_BOTH_SIDES_T,LP,VisitorType>>>(TORRESANI_MP_BOTH_SIDES_INPUT_T,graphMatchingMotorDatasets,mp_options,"motor","AMP-B");
   RunSolver<FMC_MCF_BOTH_SIDES_T, MpRoundingSolver<Solver<FMC_MCF_BOTH_SIDES_T,LP,VisitorType>>>(TORRESANI_MCF_BOTH_SIDES_INPUT_T,graphMatchingMotorDatasets,mcf_options,"motor","AMCF-B");
   RunSolver<FMC_GM_LEFT_T, MpRoundingSolver<Solver<FMC_GM_LEFT_T,LP,VisitorType>>>(TORRESANI_GM_LEFT_INPUT_T,graphMatchingMotorDatasets,gm_options,"motor","GM-B");
   RunSolver<FMC_HUNGARIAN_BP_BOTH_SIDES_T, MpRoundingSolver<Solver<FMC_HUNGARIAN_BP_BOTH_SIDES_T,LP,VisitorType>>>(TORRESANI_HUNGARIAN_BP_BOTH_SIDES_INPUT_T,graphMatchingMotorDatasets,hbp_options,"motor","HUNGARIAN_BP-B");

   //RunSolver<FMC_MP_LEFT_T, MpRoundingSolver<Solver<FMC_MP_LEFT_T,LP,VisitorType>>(TORRESANI_MP_LEFT_INPUT_T,graphMatchingWormsDatasets,mp_options,"worms","AMP-O");
   //RunSolver<FMC_MCF_LEFT_T, MpRoundingSolver<Solver<FMC_MCF_LEFT_T,LP,VisitorType>>(TORRESANI_MCF_LEFT_INPUT_T,graphMatchingWormsDatasets,mcf_options,"worms","AMCF-O");
   //RunSolver<FMC_GM_LEFT_T, MpRoundingSolver<Solver<FMC_GM_LEFT_T,LP,VisitorType>>(TORRESANI_GM_LEFT_INPUT_T,graphMatchingWormsDatasets,gm_options,"worms","GM-O");
   //RunSolver<FMC_HUNGARIAN_BP_LEFT_T, MpRoundingSolver<Solver<FMC_HUNGARIAN_BP_LEFT_T,LP,VisitorType>>(TORRESANI_HUNGARIAN_BP_LEFT_INPUT_T,graphMatchingWormsDatasets,hbp_options,"worms","HUNGARIAN_BP-O");
   }
   return 0;
}
