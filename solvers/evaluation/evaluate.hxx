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

template<typename FMC, typename SOLVER, typename INPUT_FUNCTION>
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
      solverOptions.push_back("-o");
      std::string output_file = ExtractFilename(dataset);
      output_file.append("_solution_");
      output_file.append(algorithmName);
      output_file.append(".txt");
      solverOptions.push_back(output_file);
   
      //SOLVER s(solverOptions);

      // convert std::vector<std::string> to char**
      std::vector<char*> solverOptionsRaw;
      std::transform(solverOptions.begin(), solverOptions.end(), std::back_inserter(solverOptionsRaw), 
            [](const std::string& s) {
            char *pc = new char[s.size()+1];
            std::strcpy(pc, s.c_str());
            return pc;
            });

      SOLVER s(solverOptionsRaw.size(), &solverOptionsRaw[0]);
      std::cout << "run solver " << FMC::name << " on problem " << dataset << "\n";
      s.ReadProblem(f);
      s.Solve();

      for (INDEX i=0; i<solverOptionsRaw.size(); i++) {
         delete [] solverOptionsRaw[i];
      }
   }
}



} // end namespace LP_MP

#endif // LP_MP_EVALUATE_HXX
