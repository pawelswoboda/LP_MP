#include "evaluate.hxx"
#include "visitors/sqlite_visitor.hxx"
#include "solvers/graph_matching/graph_matching.h"

// algorithms considered
using FMC_GM_LEFT = FMC_GM<PairwiseConstruction::Left>;
using FMC_GM_RIGHT = FMC_GM<PairwiseConstruction::Right>;
using FMC_MCF_LEFT = FMC_MCF<PairwiseConstruction::Left>;
using FMC_MCF_RIGHT = FMC_MCF<PairwiseConstruction::Right>;
using FMC_MCF_BOTH_SIDES = FMC_MCF<PairwiseConstruction::BothSides>;
using FMC_MP_LEFT = FMC_MP<PairwiseConstruction::Left>;
using FMC_MP_RIGHT = FMC_MP<PairwiseConstruction::Right>;
using FMC_MP_BOTH_SIDES = FMC_MP<PairwiseConstruction::BothSides>;

// Torresani et al input
static auto TORRESANI_GM_LEFT_INPUT = TorresaniEtAlInput::ParseProblem<FMC_GM_LEFT>;
static auto TORRESANI_GM_RIGHT_INPUT = TorresaniEtAlInput::ParseProblem<FMC_GM_RIGHT>;
static auto TORRESANI_MCF_LEFT_INPUT = TorresaniEtAlInput::ParseProblem<FMC_MCF_LEFT>;
static auto TORRESANI_MCF_RIGHT_INPUT = TorresaniEtAlInput::ParseProblem<FMC_MCF_RIGHT>;
static auto TORRESANI_MCF_BOTH_SIDES_INPUT = TorresaniEtAlInput::ParseProblem<FMC_MCF_BOTH_SIDES>;
static auto TORRESANI_MP_LEFT_INPUT = TorresaniEtAlInput::ParseProblem<FMC_MP_LEFT>;
static auto TORRESANI_MP_RIGHT_INPUT = TorresaniEtAlInput::ParseProblem<FMC_MP_RIGHT>;
static auto TORRESANI_MP_BOTH_SIDES_INPUT = TorresaniEtAlInput::ParseProblem<FMC_MP_BOTH_SIDES>;

// UAI input
static auto UAI_GM_INPUT = UAIInput::ParseProblem<FMC_GM<>>;
static auto UAI_MCF_INPUT = UAIInput::ParseProblem<FMC_MCF<>>;
static auto UAI_MP_INPUT = UAIInput::ParseProblem<FMC_MP<>>;

//std::vector<std::string> graphMatchingTestDataset = {
//   {"../../solvers/graph_matching/Large_QAP.txt"}
//};

std::string house_prefix = "../../../solvers/graph_matching/tkr_pami13_data/house/";
std::vector<std::string> graphMatchingHouseDatasets = {
   {house_prefix + "energy_house_frame10frame100.txt"},
   {house_prefix + "energy_house_frame10frame95.txt"},
   /*
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
   */
   {house_prefix + "energy_house_frame9frame99.txt"}
};

std::string hotel_prefix = "../../../solvers/graph_matching/tkr_pami13_data/hotel/";
std::vector<std::string> graphMatchingHotelDatasets = {
   {hotel_prefix + "energy_hotel_frame15frame22.txt"},
   {hotel_prefix + "energy_hotel_frame15frame29.txt"},
   /*
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
   */
   {hotel_prefix + "energy_hotel_frame92frame99.txt"}
};

std::string hassan_prefix = "../../../solvers/graph_matching/Abu_Hassan/torresani/";
std::vector<std::string> graphMatchingHassanDatasets = {
   {hassan_prefix + "board_torresani.txt"},
   {hassan_prefix + "books_torresani.txt"},
   {hassan_prefix + "hammer_torresani.txt"},
   {hassan_prefix + "party_torresani.txt"},
   {hassan_prefix + "table_torresani.txt"},
   {hassan_prefix + "tea_torresani.txt"},
   {hassan_prefix + "walking_torresani.txt"}
};

std::string worms_prefix = "../../../solvers/graph_matching/allWorms-03-04-1750-uai/";
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


// emulate command line options
static std::vector<std::string> options = {
   {"--maxIter"}, {"1000"},
   {"--minDualImprovement"}, {"0.00001"},
   {"--lowerBoundComputationInterval"}, {"1"},
   {"--primalComputationInterval"}, {"5"},
   {"--standardReparametrization"}, {"anisotropic"},
   {"--roundingReparametrization"}, {"uniform"},
   {"--databaseFile"}, {"graph_matching.db"},
   {"--overwriteDbRecord"} // do zrobienia: possibly deactivate this. Then we do not overwrite
};

int main()
{
   std::vector<std::string> graphMatchingDatasets;
   //for(auto& d : graphMatchingHotelDatasets) {graphMatchingDatasets.push_back(d);}
   //for(auto& d : graphMatchingHouseDatasets) {graphMatchingDatasets.push_back(d);}
   //graphMatchingDatasets.push_back(graphMatchingHotelDatasets[54]);

   //InitTikzFiles(graphMatchingDatasets);
   // do zrobienia: renamce dataset of Hassan instances

   RunSolver<FMC_MP_BOTH_SIDES,SqliteVisitor>(TORRESANI_MP_BOTH_SIDES_INPUT,graphMatchingHotelDatasets,options,"hotel","AMP-B");
   RunSolver<FMC_MP_BOTH_SIDES,SqliteVisitor>(TORRESANI_MP_BOTH_SIDES_INPUT,graphMatchingHouseDatasets,options,"house","AMP-B");
   //RunSolver<FMC_MP_BOTH_SIDES,SqliteVisitor>(TORRESANI_MP_BOTH_SIDES_INPUT,graphMatchingHassanDatasets,options,"house","AMP-B");

   RunSolver<FMC_MP_LEFT,SqliteVisitor>(TORRESANI_MP_LEFT_INPUT,graphMatchingHotelDatasets,options,"hotel","AMP-O");
   RunSolver<FMC_MP_LEFT,SqliteVisitor>(TORRESANI_MP_LEFT_INPUT,graphMatchingHouseDatasets,options,"house","AMP-O");
   //RunSolver<FMC_MP_LEFT,SqliteVisitor>(TORRESANI_MP_LEFT_INPUT,graphMatchingHassanDatasets,options,"house","AMP-O");

   RunSolver<FMC_MP_RIGHT,SqliteVisitor>(TORRESANI_MP_RIGHT_INPUT,graphMatchingHotelDatasets,options,"hotel","AMP-I");
   RunSolver<FMC_MP_RIGHT,SqliteVisitor>(TORRESANI_MP_RIGHT_INPUT,graphMatchingHouseDatasets,options,"house","AMP-I");
   //RunSolver<FMC_MP_RIGHT,SqliteVisitor>(TORRESANI_MP_RIGHT_INPUT,graphMatchingHassanDatasets,options,"house","AMP-I");

   RunSolver<FMC_MCF_BOTH_SIDES,SqliteVisitor>(TORRESANI_MCF_BOTH_SIDES_INPUT,graphMatchingHotelDatasets,options,"hotel","AMCF-B");
   RunSolver<FMC_MCF_BOTH_SIDES,SqliteVisitor>(TORRESANI_MCF_BOTH_SIDES_INPUT,graphMatchingHouseDatasets,options,"house","AMCF-B");
   //RunSolver<FMC_MCF_BOTH_SIDES,SqliteVisitor>(TORRESANI_MCF_BOTH_SIDES_INPUT,graphMatchingHassanDatasets,options,"house","AMCF-B");

   RunSolver<FMC_MCF_LEFT,SqliteVisitor>(TORRESANI_MCF_LEFT_INPUT,graphMatchingHotelDatasets,options,"hotel","AMCF-O");
   RunSolver<FMC_MCF_LEFT,SqliteVisitor>(TORRESANI_MCF_LEFT_INPUT,graphMatchingHouseDatasets,options,"house","AMCF-O");
   //RunSolver<FMC_MCF_LEFT,SqliteVisitor>(TORRESANI_MCF_LEFT_INPUT,graphMatchingHassanDatasets,options,"house","AMCF-O");

   RunSolver<FMC_MCF_RIGHT,SqliteVisitor>(TORRESANI_MCF_RIGHT_INPUT,graphMatchingHotelDatasets,options,"hotel","AMCF-I");
   RunSolver<FMC_MCF_RIGHT,SqliteVisitor>(TORRESANI_MCF_RIGHT_INPUT,graphMatchingHouseDatasets,options,"house","AMCF-I");
   //RunSolver<FMC_MCF_RIGHT,SqliteVisitor>(TORRESANI_MCF_RIGHT_INPUT,graphMatchingHassanDatasets,options,"house","AMCF-I");
   
   RunSolver<FMC_GM_LEFT,SqliteVisitor>(TORRESANI_GM_LEFT_INPUT,graphMatchingHotelDatasets,options,"hotel","GM-O");
   RunSolver<FMC_GM_LEFT,SqliteVisitor>(TORRESANI_GM_LEFT_INPUT,graphMatchingHouseDatasets,options,"house","GM-O");
   //RunSolver<FMC_GM_LEFT,SqliteVisitor>(TORRESANI_GM_LEFT_INPUT,graphMatchingHassanDatasets,options,"house","GM-O");
   
   RunSolver<FMC_GM_RIGHT,SqliteVisitor>(TORRESANI_GM_RIGHT_INPUT,graphMatchingHotelDatasets,options,"hotel","GM-I");
   RunSolver<FMC_GM_RIGHT,SqliteVisitor>(TORRESANI_GM_RIGHT_INPUT,graphMatchingHouseDatasets,options,"house","GM-I");
   //RunSolver<FMC_GM_RIGHT,SqliteVisitor>(TORRESANI_GM_RIGHT_INPUT,graphMatchingHassanDatasets,options,"house","GM-I");

   //RunSolver<FMC_MP_LEFT,SqliteVisitor>(UAI_MP_INPUT,graphMatchingWormsDatasets,options,"worms","AMP-O");
   //RunSolver<FMC_MCF_LEFT,SqliteVisitor>(UAI_MCF_INPUT,graphMatchingWormsDatasets,options,"worms","AMCF-o");
   //RunSolver<FMC_GM_LEFT,SqliteVisitor>(UAI_GM_INPUT,graphMatchingWormsDatasets,options,"worms","GM-O");

   //FinishTikzFiles(graphMatchingDatasets);
}
