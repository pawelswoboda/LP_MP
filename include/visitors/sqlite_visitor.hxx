#ifndef LP_MP_SQLITE_VISITOR_HXX
#define LP_MP_SQLITE_VISITOR_HXX

#include <sqlite3.h>
#include <iostream>
#include <fstream>

#include "standard_visitor.hxx"
//#include "solvers/evaluation/evaluate.hxx" // do zrobienia: for extract file name. Put this in special file

namespace LP_MP {

struct IterationStatistics {
   INDEX iteration_;
   INDEX timeElapsed_; // in milliseconds
   REAL lowerBound_;
   REAL upperBound_;
};

// table layout of sqlite database is
// table ProblemSolverTable (UniqueId id, TEXT dataset, TEXT instance, TEXT algorithmName, TEXT algorithmFMC) // the three fields after id must be unique . algorithmFMC is the name given in the FMC.
// table Iterations (INTEGER pst_id, INTEGER iteration, TIME runtime, REAL lowerBound, REAL upperBound) // id refers to id of ProblemSolverTable, pair (id,iteration) must be unique.

// this visitor connects to given sqlite database and writes or updates the runtime and iteration data of the algorithm.
// do zrobienia: transaction support to avoid concurrent writing to database?
template<class PROBLEM_DECOMPOSITION>
class SqliteVisitor : public StandardVisitor<PROBLEM_DECOMPOSITION> {

public:
   SqliteVisitor(TCLAP::CmdLine& cmd, PROBLEM_DECOMPOSITION& pd) 
      :
         StandardVisitor<PROBLEM_DECOMPOSITION>(cmd,pd),
         pd_(pd),
         databaseFileArg_("","databaseFile","sqlite database into which to protocolate runtime/iteration vs. primal/dual energy",true,"","file name",cmd),
         datasetNameArg_("","datasetName","name of dataset the input file belongs to",true,"","string",cmd),
         algorithmNameArg_("","algorithmName","name of algorithm",true,"","string",cmd),
         overwriteDbRecordArg_("","overwriteDbRecord","overwrite previous record, if not, and record is present, abort optimization",cmd,false)
   {}

   LPVisitorReturnType begin(const LP* lp) // called, after problem is constructed. 
   {
      const auto ret = StandardVisitor<PROBLEM_DECOMPOSITION>::begin(lp);
      try {
         databaseFile_ = databaseFileArg_.getValue();
         datasetName_ = datasetNameArg_.getValue();
         algorithmName_ = algorithmNameArg_.getValue();
         overwriteDbRecord_ = overwriteDbRecordArg_.getValue();
      } catch (TCLAP::ArgException &e) {
         std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
         exit(1);
      }
      if(sqlite3_open(databaseFile.c_str(), &database_) != SQLITE_OK) {
         std::cout << "Could not open database file " << databaseFile_ << ", error: " << sqlite3_errmsg(database_) << "\n";
         return LPVisitorReturnType::Error;
      }
      // if record is already present, and overwrite is set to false, return; otherwise insert empty record
      const std::string inputFile = ExtractFileName(pd_.GetInputFileName());
      const std::string algorithmFMC_ = PROBLEM_DECOMPOSITION::FactorMessageConnection::name;
      std::string checkRecordPresentStatement = "SELECT (id) from ProblemSolverTable WHERE dataset = " << datasetName_ << " AND instance = " << pd_.GetInputFileName() << " AND algorithmName = " << algorithmName << ";\n";
      // lambda records id of entry
      char* zErrMsg = nullptr;
      int rc;
      // do zrobienia: begin transaction here
      rc = sqlite3_exec(database_, checkRecordPresentStatement, [&dbRecordId_](void* int argc, char **argv, char **azColName) ->int { return 0; },0, &zErrMsg);
      if(rc || (!overwriteDbRecordArg_ && dbRecordId_ >= 0)) {
         sqlite3_close(database_);
         return LPVisitorReturnType::Error; // immediately abort;
      }
      if(dbRecordId_ < 0) { // insert new record
         const std::string insertIdStatement = "INSERT INTO ProblemSolverTable (dataset, instance, algorithmName, algorithmFMC) VALUES (" << datasetName_ << "," << inputFile << "," << algorithmName_ << "," << algorithmFMC_ << ");";
         rc = sqlite3_exec(database_, insertIdStatement, [&dbRecordId_](void* int argc, char **argv, char **azColName) ->int { return 0; },0, &zErrMsg);
         assert(rc == 0);
      }
      // do zrobienia: end transaction here

      return ret;
   }

   template<LPVisitorReturnType LP_STATE>
   LPVisitorReturnType visit(const LP* lp)
   {
      auto ret_state = this->template StandardVisitor<PROBLEM_DECOMPOSITION>::template visit<LP_STATE>(lp);
      switch(LP_STATE) {
         case LPVisitorReturnType::ReparametrizeAndComputePrimal:
            {
            const REAL lowerBound = this->GetLowerBound();
            const REAL upperBound = lp->BestPrimalBound();
            const INDEX timeElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - StandardVisitor<PROBLEM_DECOMPOSITION>::GetBeginTime()).count();
            const INDEX curIter = StandardVisitor<PROBLEM_DECOMPOSITION>::GetIter();
            iterationStatistics_.push_back({curIter,timeElapsed,lowerBound,upperBound});
            break;
            }
         case LPVisitorReturnType::Break:
            {
            WriteBounds(iterationStatistics_);
            sqlite3_close(database_);
            break;
            }
      }
      return ret_state;
   }
   void WriteBounds(const std::string& filename, const std::string datasetName, const std::string algorithmName, const std::vector<IterationStatistics>& iterStats)
   {
      // connection to database is open and we have a corresponding id in ProblemSolverTable
      // Write into Iterations all iterations
      // first remove all iterations previously in this table
      std::string rmIterStmt = "DELETE FROM Iterations WHERE pst_id = " << dbRecordId_ << "\n";
      int rc = sqlite3_exec(database_, rmIterStmt, nullptr, 0, nullptr);
      assert(rc == 0);
      for(const auto& it : iterStats) {
         std::string stmt = "INSERT INTO Iterations (pst_id,iteration,runtime,lowerBound,upperBound) VALUES (" << dbRecordId_ << "," << it.iteration_ << "," << it.timeElapsed_ << "," << it.lowerBound_ << "," << it.upperBound_ << ");";
         rc = sqlite3_exec(database_, stmt, [](...) -> int { return 0; },0, nullptr);
         assert(rc);
      }
   }

private:
   TCLAP::ValueArg<std::string> databaseFileArg_;
   TCLAP::ValueArg<std::string> datasetNameArg_;
   TCLAP::ValueArg<std::string> algorithmNameArg_; // custom name given for algorithm, to differentiate between same algorithm with differing options
   TCLAP::SwitchArg overwriteDbRecordArg_;

   std::string databaseFile_;
   std::string datasetName_;
   std::string algorithmName_;
   bool overwriteDbRecord_;
   SIGNED_INDEX dbRecordId_;

   std::vector<IterationStatistics> iterationStatistics_;
   PROBLEM_DECOMPOSITION& pd_;

   sqlite* database_;
};

} // end namespace LP_MP

#endif // LP_MP_SQLITE_VISITOR_HXX

