/*
 * this meta-solver collects several solvers together with runtime options and a list of problem categories consisting of problem instances.
 * All problem categories are solved by all solvers applicable to them and iteration information is stored in a database
 *
 * additional ideas: distribute computations over several computers and fetch problems over ssh-server and write to one database, also acting as synchronization
 */

#ifndef LP_MP_EVALUATE_HXX
#define LP_MP_EVALUATE_HXX

#include "LP_MP.h"
#include "help_functions.hxx"
#include "parse_rules.h"
#include "solver.hxx"


namespace LP_MP {

template<typename FMC, class VISITOR, typename SOLVER, typename INPUT_FUNCTION>
void RunSolver(INPUT_FUNCTION f, const std::vector<std::string>& datasets, std::vector<std::string> options, const std::string& datasetName, const std::string& algorithmName)
{
   options.push_back("--datasetName");
   options.push_back(datasetName);
   options.push_back("--algorithmName");
   options.push_back(algorithmName);
   options.push_back("--algorithmFMC");
   options.push_back(FMC::name);

   for(auto& dataset : datasets) {

      auto solverOptions = options;
      solverOptions.insert(solverOptions.begin(), std::string("LP_MP evaluation tool"));
      solverOptions.push_back("-i");
      solverOptions.push_back(dataset);
   
      // do zrobienia: fetch version number automatically from CMakeLists.txt
      //TCLAP::CmdLine cmd(std::string("Command line options for evaluation run"), ' ', "0.0.1");
      VisitorSolver<SOLVER,VISITOR> s(solverOptions);
      std::cout << "run solver " << FMC::name << " on problem " << dataset << "\n";
      s.ReadProblem(f);
      s.Solve();
   }
}



} // end namespace LP_MP

#endif // LP_MP_EVALUATE_HXX
