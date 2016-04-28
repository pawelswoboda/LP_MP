#ifndef LP_MP_SQLITE_VISITOR_HXX
#define LP_MP_SQLITE_VISITOR_HXX

#include <sqlite3.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "standard_visitor.hxx"
#include "help_functions.hxx"
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
         overwriteDbRecordArg_("","overwriteDbRecord","if true: overwrite previous record. if false: if record is present, abort optimization",cmd,false)
   {}

   static int CountCallback(void* c_ptr, int argc, char** argv, char**)
   {
      assert(argc == 1);
      *static_cast<int*>(c_ptr) = atoi(argv[0]);
      return 0;
   }

   bool TableExists(const std::string& tableName)
   {
      int rc, ret;
      const std::string existsTable = "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='" + tableName + "';";
      rc = sqlite3_exec(database_, existsTable.c_str(), &CountCallback, &ret, NULL);
      assert(rc == 0);
      assert(ret == 0 || ret == 1);
      return ret == 1;
   }

   void CreateTable(const std::string tableName, const std::vector<std::string>& columns)
   {
      std::string createTable = "CREATE TABLE " + tableName + "(";
      for(INDEX i=0; i<columns.size(); ++i) {
         createTable += columns[i];
         if(i < columns.size()-1) {
            createTable += ", ";
         }
      }
      createTable += ");";
      std::cout << createTable << "\n";
      int rc = sqlite3_exec(database_, createTable.c_str(), nullptr, nullptr, nullptr);
      assert(rc == 0);
   }
   void ConditionallyCreateTable(const std::string tableName, const std::vector<std::string>& columns)
   {
      if(!TableExists(tableName)) { CreateTable(tableName, columns); }
   }

   void BuildDb() 
   {
      ConditionallyCreateTable("Solvers", 
            {"id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL",
            "algorithmName TEXT NOT NULL",
            "algorithmFMC TEXT NOT NULL",
            "UNIQUE(algorithmName, algorithmFMC)"});
      ConditionallyCreateTable("Datasets", 
            {"id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL",
            "name TEXT NOT NULL",
            "UNIQUE(name)"});
      ConditionallyCreateTable("Instances", 
            {"id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL",
            "dataset_id INT",
            "name TEXT NOT NULL",
            "FOREIGN KEY(dataset_id) REFERENCES Datasets(id)",
            "UNIQUE(dataset_id, name)"});
      ConditionallyCreateTable("Iterations", 
            {"instance_id INTEGER NOT NULL",
            "solver_id INTEGER NOT NULL",
            "iteration INTEGER NOT NULL",
            "runtime INT",
            "lower_bound DOUBLE",
            "upper_bound DOUBLE",
            "FOREIGN KEY(instance_id) REFERENCES Instances(id)",
            "FOREIGN KEY(solver_id) REFERENCES Solvers(id)",
            "UNIQUE(instance_id, solver_id, iteration)"});
   }

   // assume that conditionSQL returns an id if record is present, otherwise insert record and retrieve id again
   struct IdRet {bool recordPresent; int id;};
   static int RecordIdCallback(void* i, int argc, char** argv, char**) 
   {
      assert(argc <= 1);
      if(argc == 0) {
         *static_cast<IdRet*>(i) = IdRet({false,0});
      } else {
         *static_cast<IdRet*>(i) = IdRet({true,atoi(argv[0])});
      }
      return 0;
   }
      
   int ConditionallyInsertById(const std::string& conditionSQL, const std::string& insertSQL)
   {
      int rc;
      IdRet i;
      rc = sqlite3_exec(database_, conditionSQL.c_str(), &RecordIdCallback, &i, nullptr);
      assert(rc == 0);
      if(i.recordPresent) {
         return i.id;
      } else {
         rc = sqlite3_exec(database_, insertSQL.c_str(), nullptr, nullptr, nullptr);
         assert(rc == 0);
         rc = sqlite3_exec(database_, conditionSQL.c_str(), &RecordIdCallback, &i, nullptr);
         assert(rc == 0);
         assert(i.recordPresent == true);
         return i.id;
      }
   }

   // return id of solver
   int GetSolverId(const std::string& algorithmName, const std::string& algorithmFMC)
   {
      const std::string getSolverId = "SELECT (id) FROM Solvers WHERE algorithmName = '" + algorithmName + "' and algorithmFMC = '" + algorithmFMC + "';";
      const std::string insertSolver = "INSERT INTO Solvers (algorithmName, algorithmFMC) VALUES ('" + algorithmName + "', '" + algorithmFMC + "');";
      return ConditionallyInsertById(getSolverId, insertSolver);
   }
   int GetDatasetId(const std::string& dataset)
   {
      const std::string getDatasetId = "SELECT (id) FROM Datasets WHERE name = '" + dataset + "';";
      const std::string insertDataset = "INSERT INTO Datasets (name) VALUES ('" + dataset + "');";
      return ConditionallyInsertById(getDatasetId, insertDataset);
   }

   int GetInstanceId(const std::string& instance, const int dataset_id)
   {
      const std::string getInstanceId = "SELECT (id) FROM Instances WHERE name = '" + instance + "' AND dataset_id = " + std::to_string(dataset_id) + ";";
      const std::string insertInstance = "INSERT INTO Instances (name, dataset_id) VALUES ('" + instance + "', " + std::to_string(dataset_id) + ");";
      std::cout << getInstanceId << "\n";
      std::cout << insertInstance << "\n";
      return ConditionallyInsertById(getInstanceId, insertInstance);
   }

   static int GetIdCallback(void* id_ptr, const int argc, char** argv, char**) 
   { 
      assert(argc == 1);
      *static_cast<int*>(id_ptr) = atoi(argv[0]);
      return 0; 
   }


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

      if(sqlite3_open(databaseFile_.c_str(), &database_) != SQLITE_OK) {
         std::cout << "Could not open database file " << databaseFile_ << ", error: " << sqlite3_errmsg(database_) << "\n";
         return LPVisitorReturnType::Error;
      }

      BuildDb();

      const std::string inputFile = ExtractFilename(pd_.GetInputFileName());
      const std::string algorithmFMC_(PROBLEM_DECOMPOSITION::FactorMessageConnection::name);

      solver_id_ = GetSolverId(algorithmName_, algorithmFMC_);
      dataset_id_ = GetDatasetId(datasetName_);
      instance_id_ = GetInstanceId(inputFile, dataset_id_);

      return ret;
   }

   template<LPVisitorReturnType LP_STATE>
   LPVisitorReturnType visit(LP* lp)
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
               std::cout << "Write iteration statistics\n";
               WriteBounds(iterationStatistics_);
               sqlite3_close(database_);
               break;
            }
      }
      return ret_state;
   }

   void WriteBounds(const std::vector<IterationStatistics>& iterStats)
   {
      // connection to database is open and we have a corresponding id in ProblemSolverTable
      // Write into Iterations all iterations
      // first remove all iterations previously in this table
      std::cout << "overwrite previous iterations\n";
      const std::string rmIterStmt = "DELETE FROM Iterations WHERE solver_id = " + std::to_string(solver_id_) + " AND instance_id = " + std::to_string(instance_id_) + "\n";
      int rc = sqlite3_exec(database_, rmIterStmt.c_str(), nullptr, nullptr, nullptr);
      assert(rc == 0);
      for(const auto& it : iterStats) {
         const std::string stmt = "INSERT INTO Iterations (solver_id, instance_id, iteration, runtime, lower_bound, upper_bound) VALUES ('" + std::to_string(solver_id_) + "', '" + std::to_string(instance_id_) + "', '" + std::to_string(it.iteration_) + "', '" + std::to_string(it.timeElapsed_) + "', '" + std::to_string(it.lowerBound_) + "', '" + std::to_string(it.upperBound_) + "');";
         std::cout << stmt << "\n";
         rc = sqlite3_exec(database_, stmt.c_str(), nullptr, 0, nullptr);
         assert(rc == 0);
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

   std::vector<IterationStatistics> iterationStatistics_;
   PROBLEM_DECOMPOSITION& pd_;

   sqlite3* database_;
   int solver_id_;
   int dataset_id_;
   int instance_id_;
};

} // end namespace LP_MP

#endif // LP_MP_SQLITE_VISITOR_HXX

