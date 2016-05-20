/*
 * this meta-solver collects several solvers together with runtime options and a list of problem categories consisting of problem instances.
 * All problem categories are solved by all solvers applicable to them and convergence plots are generated.
 * Also a table summarising final primal-dual gaps and runtimes is generated.
 *
 * additional ideas: distribute computations over several computers and fetch problems over ssh-server and write back convergence plots.
 */

#ifndef LP_MP_EVALUATE_HXX
#define LP_MP_EVALUATE_HXX

#include "LP_MP.h"
#include "problem_decomposition.hxx"
#include "help_functions.hxx"
#include "parse_rules.h"


namespace LP_MP {

/*
// this is not the nicest thing with its mix of variables and typenames below. However input function needs to be stored via value
// evaluate expects a meta::list of Solver structs parametrizing which solvers to take
template<typename FMC_SOLVER, typename INPUT, typename VISITOR = TikzVisitor>
struct Solver {
   INPUT& InputFunction;
   using FMC = FMC_SOLVER;
   using VisitorType = VISITOR;
};

template<typename FMC_SOLVER, typename INPUT, typename VISITOR = TikzVisitor>
Solver<FMC_SOLVER,INPUT,VISITOR> CreateSolver(FMC_SOLVER, INPUT f)
{
   return Solver<FMC_SOLVER,INPUT,VISITOR>{f};
}
*/

//template<typename... SOLVER_LIST>
//void RunSolver(std::tuple<SOLVER_LIST...>, std::vector<std::string>, const std::string& ) {}
//template<typename SOLVER, typename... SOLVER_LIST>
//void RunSolver(std::tuple<SOLVER, SOLVER_LIST...> solvers, const std::string& dataset, const std::vector<std::string>& options)


/*
std::string GetRuntimeFilename(const std::string& inputFile)
{
   std::string tikzFile = ExtractFilename(inputFile);
   tikzFile.append("_runtime.tikz");
   return tikzFile;
}
std::string GetIterationFilename(const std::string& inputFile)
{
   // we want to write to current directory, hence only retain the current file name, excluding ending
   std::string tikzFile = ExtractFilename(inputFile);
   tikzFile.append("_iteration.tikz");
   return tikzFile;
}

void InitStatisticsFile(const std::string& statisticsFileName)
{
   // clear file
   std::cout << "Empty file " << statisticsFileName << "\n";
   std::ofstream f;
   f.open(statisticsFileName, std::ofstream::out | std::ofstream::trunc);
   f.close();
}
void InitTikzFiles(const std::vector<std::string>& datasets)
{
   const std::string preamble = 
      std::string("\\begin{figure}\\begin{tikzpicture}\n\\begin{axis}[legend style={at={(0.5,-0.1)}, anchor=north}]\n\n");

   for(auto& dataset : datasets) {
      std::string runtimePlotFile = GetRuntimeFilename(dataset);
      std::string iterationPlotFile = GetIterationFilename(dataset);
      std::cout << "remove files " << runtimePlotFile << " and " << iterationPlotFile << "\n";
      std::remove(runtimePlotFile.c_str());
      std::remove(iterationPlotFile.c_str());

      std::ofstream runtimePlotStream({runtimePlotFile});
      std::ofstream iterationPlotStream({iterationPlotFile});

      runtimePlotStream << preamble;
      iterationPlotStream << preamble;

      runtimePlotStream.close();
      iterationPlotStream.close();
   }
}

void FinishTikzFiles(const std::vector<std::string>& datasets)
{

   for(auto& dataset : datasets) {
      const std::string ending = std::string("\n\n\\end{axis}\n\\end{tikzpicture}\n\\caption{" + LatexEscape(ExtractFilename(dataset)) + "}\n\\end{figure}");
      std::string runtimePlotFile = GetRuntimeFilename(dataset);
      std::string iterationPlotFile = GetIterationFilename(dataset);
      std::cout << "finish file " << runtimePlotFile << "\n";

      std::ofstream runtimePlotStream(runtimePlotFile, std::ios_base::app);
      std::ofstream iterationPlotStream(iterationPlotFile, std::ios_base::app);
      
      runtimePlotStream << ending;
      iterationPlotStream << ending;

      runtimePlotStream.close();
      iterationPlotStream.close();
   }

   // now collect all results in one large tex file
   // do zrobienia: get name for experiments and use here
   std::string runtimePlotFile("runtime_results.tex");
   std::ofstream runtimePlotStream(runtimePlotFile);

   std::string iterationPlotFile("iteration_results.tex");
   std::ofstream iterationPlotStream(iterationPlotFile);

   const std::string preamble = 
"\\documentclass[a4paper,11pt,twoside]{scrbook}\n\\usepackage{tikz}\n\\usepackage{pgfplots}\n\\begin{document}\n\n";
   runtimePlotStream << preamble;
   iterationPlotStream << preamble;

   for(auto& dataset : datasets) {
      std::string runtimePlotFile = GetRuntimeFilename(dataset);
      std::string iterationPlotFile = GetIterationFilename(dataset);
      runtimePlotStream << "\\input{" << runtimePlotFile << "}\n\\newpage\n";
      iterationPlotStream << "\\input{" << iterationPlotFile << "}\n\\newpage\n";
   }

   const std::string ending = "\\end{document}";
   runtimePlotStream << ending;
   iterationPlotStream << ending;
}
*/

//template<typename FMC, typename VISITOR , typename INPUT_FUNCTION>
template<typename FMC, template<typename> class VISITOR , typename INPUT_FUNCTION>
void RunSolver(INPUT_FUNCTION f, const std::vector<std::string>& datasets, std::vector<std::string> options, const std::string& datasetName, const std::string& algorithmName)
{
   options.push_back("--datasetName");
   options.push_back(datasetName);
   options.push_back("--algorithmName");
   options.push_back(algorithmName);
   //if(lineStyleCounter >= lineStyles.size()-1) {
   //   throw std::runtime_error("not enough line styles defined for tikz");
   //}
   // nicer syntax: push option and argument simultaneously
   for(auto& dataset : datasets) {

      auto solverOptions = options;
      solverOptions.insert(solverOptions.begin(), std::string("LP_MP evaluation tool"));
      solverOptions.push_back("-i");
      solverOptions.push_back(dataset);
      // better output in current directory
      //solverOptions.push_back("--runtimePlotFile");
      //solverOptions.push_back(GetRuntimeFilename(dataset));
      //solverOptions.push_back("--iterationPlotFile");
      //solverOptions.push_back(GetIterationFilename(dataset));
      //solverOptions.push_back("--lowerBoundLegend");
      //solverOptions.push_back(FMC::name);
      //solverOptions.push_back("--lowerBoundColor");
      //solverOptions.push_back(lineStyles[lineStyleCounter]);
      //solverOptions.push_back("--upperBoundColor");
      //solverOptions.push_back(lineStyles[lineStyleCounter+1]);
   
      // do zrobienia: fetch version number automatically from CMakeLists.txt
      TCLAP::CmdLine cmd(std::string("Command line options for evaluation run"), ' ', "0.0.1");
      ProblemDecomposition<FMC> pd(cmd);
      VISITOR<ProblemDecomposition<FMC>> visitor(cmd,pd);
      cmd.parse(solverOptions);
      std::cout << "run solver " << FMC::name << " on problem " << dataset << "\n";
      pd.ParseProblem(f);
      pd.Run(visitor);
   }
   //lineStyleCounter += 2;
}


/*
// SOLVER_LIST is a meta::list of Solver structs
template<typename SOLVER_LIST>
void Run(SOLVER_LIST SolverList, const std::vector<std::string>& datasets, std::vector<std::string> options)
{
   // first entry of options must be the program name
   options.insert(options.begin(), std::string("LP_MP evaluation tool"));
   // run each solver on each dataset
   for(auto dataset : datasets) {
      //meta::for_each(SOLVER_LIST{},Run); // would be nice, write own
      RunSolver(SolverList, dataset, options);
   }
}
*/

} // end namespace LP_MP

#endif // LP_MP_EVALUATE_HXX
