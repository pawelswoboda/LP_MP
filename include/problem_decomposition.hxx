#ifndef LP_MP_PROBLEM_DECOMPOSITION_HXX
#define LP_MP_PROBLEM_DECOMPOSITION_HXX

#include "LP_MP.h"
#include "meta/meta.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <regex>

// construct problem from factor-message network with problem decomposition by invoking the correct problem decomposition routines of the respective classes

namespace LP_MP {

// do zrobienia: better name for ProblemDecomposition?
template<typename FMC>
class ProblemDecomposition {
   using problem_decomposition = typename FMC::problem_decomposition;
public:
   ProblemDecomposition() : lp_(nullptr) {}
   ~ProblemDecomposition() { if(lp_ != nullptr) delete lp_; }

   ProblemDecomposition& ReadFile(const std::string filename) 
   {
      std::fstream fs{filename};
      if(!fs.is_open()) { 
         std::cout << "Cannot open " << filename << "\n";
         throw std::runtime_error("Could not open file");
      }

      // read each line and pass it to the correct problem constructor
      INDEX current_problem_no = 0;
      for(std::string line; std::getline(fs,line);) {
         if(line.find("problem ") != std::string::npos) {
            line.erase(0,8); // remove string "problem "
            current_problem_no = std::stol(line);
            if(current_problem_no > std::tuple_size<decltype(problem_constructor_)>()) throw std::runtime_error("current problem number must be smaller than " + std::to_string(std::tuple_size<decltype(problem_constructor_)>()));
         }
         /*
         std::regex problem_no_regex {"problem\\s*(\\d+)"};
         std::cout << "built regex\n";
         std::smatch problem_no_match;
         if(std::regex_search(line, problem_no_match, problem_no_regex)) {// get problem number
            current_problem_no = std::stol(problem_no_match[0]);
            std::cout << "read problem no = " << current_problem_no << "\n";
            assert( current_problem_no < std::tuple_size<decltype(problem_constructor_)>() );
         }
         */
         else {
            ReadLine(problem_decomposition{},current_problem_no, line);
         }
      }

      fs.close();
      return *this;
   }

   template<class... PROBLEM_CONSTRUCTOR_REST>
   void ReadLine(meta::list<PROBLEM_CONSTRUCTOR_REST...> t, const INDEX problem_no, const std::string& line) 
   {}

   template<class PROBLEM_CONSTRUCTOR, class... PROBLEM_CONSTRUCTOR_REST>
   void ReadLine(meta::list<PROBLEM_CONSTRUCTOR, PROBLEM_CONSTRUCTOR_REST...> t, const INDEX problem_no, const std::string& line) 
   {
      if(problem_no == 0) {
         PROBLEM_CONSTRUCTOR& pc = GetProblemConstructor<PROBLEM_CONSTRUCTOR>();
         pc.ReadLine(line);
         return;
      } else {
         return ReadLine(meta::list<PROBLEM_CONSTRUCTOR_REST...>{}, problem_no - 1, line);
      }
   }

   template<class PROBLEM_CONSTRUCTOR>
   PROBLEM_CONSTRUCTOR& GetProblemConstructor() 
   {
      constexpr INDEX n = meta::find_index<problem_decomposition, PROBLEM_CONSTRUCTOR>::value;
      static_assert(n < std::tuple_size<decltype(problem_constructor_)>(),"");
      return std::get<n>(problem_constructor_);
   }

   template<class... PROBLEM_CONSTRUCTOR_REST>
   void Construct(meta::list<PROBLEM_CONSTRUCTOR_REST...> t) 
   {}

   template<class PROBLEM_CONSTRUCTOR, class... PROBLEM_CONSTRUCTOR_REST>
   void Construct(meta::list<PROBLEM_CONSTRUCTOR, PROBLEM_CONSTRUCTOR_REST...> t)
   {
      GetProblemConstructor<PROBLEM_CONSTRUCTOR>().Construct(*this);
      return Construct(meta::list<PROBLEM_CONSTRUCTOR_REST...>{});
   }

   // call the Construct from each problem decomposition
   ProblemDecomposition& Construct()
   {
      if(lp_ == nullptr) delete lp_;
      lp_ = new LP();
      Construct(problem_decomposition{});
      return *this;
   }

   LP* GetLP() const { return lp_; }
   
private:
   LP* lp_;
   tuple_from_list<problem_decomposition> problem_constructor_;
};

// Macro for generating main function reading in file and then solving problem
#define LP_MP_CONSTRUCT_SOLVER(FMC) \
int main(int argc, char * argv[]) \
{ \
   ProblemDecomposition<typename FMC> pd; \
\
   if(argc != 2) { \
      std::cout << "1 argument must be supplied" << std::endl; \
      exit(1); \
   } \
   std::string filename(argv[1]); \
   std::cout.precision(9); \
\
   pd.ReadFile(filename); \
   pd.Construct(); \
   pd.GetLP()->Solve(750); \
}

} // end namespace LP_MP


#endif // LP_MP_PROBLEM_DECOMPOSITION_HXX
