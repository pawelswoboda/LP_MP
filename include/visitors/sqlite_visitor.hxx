#ifndef LP_MP_SQLITE_VISITOR_HXX
#define LP_MP_SQLITE_VISITOR_HXX

#include <sqlite3.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "standard_visitor.hxx"
#include "help_functions.hxx"

namespace LP_MP {

struct IterationStatistics {
   INDEX iteration_;
   INDEX timeElapsed_; // in milliseconds
   REAL lowerBound_;
   REAL upperBound_;
};

// this visitor connects to given sqlite database and writes or updates the runtime and iteration data of the algorithm.
// do zrobienia: transaction support to avoid concurrent writing to database?
// do zrobienia: error handling
template<class BASE_VISITOR = StandardVisitor>


class SqliteVisitor : public BASE_VISITOR {

const std::string sql = 
R"(
CREATE TABLE IF NOT EXISTS SOLVERS (
id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
algorithmName TEXT NOT NULL,
algorithmFMC TEXT NOT NULL,
UNIQUE(algorithmName, algorithmFMC) 
);

CREATE TABLE IF NOT EXISTS Datasets (
id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
name TEXT NOT NULL,
UNIQUE(name) 
);

CREATE TABLE IF NOT EXISTS Instances (
id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
dataset_id INT,
name TEXT NOT NULL,
ground_truth TEXT,
FOREIGN KEY(dataset_id) REFERENCES Datasets(id),
UNIQUE(dataset_id, name) 
);

CREATE INDEX IF NOT EXISTS InstancesIndex ON Instances(id,dataset_id);

CREATE TABLE IF NOT EXISTS Iterations (
instance_id INTEGER NOT NULL,
solver_id INTEGER NOT NULL,
iteration INTEGER NOT NULL,
runtime INT,
lowerBound DOUBLE PRECISION,
upperBound DOUBLE PRECISION,
FOREIGN KEY(instance_id) REFERENCES Instances(id)
FOREIGN KEY(solver_id) REFERENCES Solvers(id)
UNIQUE(instance_id, solver_id, iteration) 
);

CREATE TABLE IF NOT EXISTS Solutions (
instance_id INTEGER NOT NULL,
solver_id INTEGER NOT NULL,
solution TEXT NOT NULL,
FOREIGN_KEY(instance_id) REFERENCES Instances(id)
FOREIGN_KEY(solver_id) REFERENCES Solvers(id)
UNIQUE(instance_id, solver_id)
);

CREATE INDEX IF NOT EXISTS IterationsIndex ON Iterations(instance_id,solver_id,lowerBound,upperBound,runtime);

CREATE VIEW IF NOT EXISTS LowerBoundView AS
SELECT MAX(lowerBound),solver_id,instance_id FROM Iterations GROUP BY solver_id,instance_id;

CREATE VIEW IF NOT EXISTS MaxLowerBoundView AS
SELECT MAX(lowerBound),instance_id FROM Iterations GROUP BY instance_id;

CREATE VIEW IF NOT EXISTS UpperBoundView AS
SELECT MIN(upperBound),solver_id,instance_id FROM Iterations GROUP BY solver_id,instance_id;

CREATE VIEW IF NOT EXISTS MinUpperBoundView AS
SELECT MIN(upperBound),instance_id FROM Iterations GROUP BY instance_id;

CREATE VIEW IF NOT EXISTS RuntimeView AS 
SELECT MAX(runtime),solver_id,instance_id FROM Iterations GROUP BY solver_id,instance_id;

CREATE VIEW IF NOT EXISTS AggregateIterationsHelper AS
SELECT LowerBoundView.lowerBound AS lowerBound, UpperBoundView.upperBound AS upperBound, RuntimeView.runtime AS runtime, LowerBoundView.solver_id AS solver_id, LowerBoundView.instance_id AS instance_id
FROM LowerBoundView
INNER JOIN UpperBoundView ON (LowerBoundView.solver_id = UpperBoundView.solver_id AND LowerBoundView.instance_id = UpperBoundView.instance_id)
INNER JOIN RuntimeView ON (LowerBoundView.solver_id = RuntimeView.solver_id AND LowerBoundView.instance_id = RuntimeView.instance_id);

CREATE VIEW IF NOT EXISTS MinMaxBoundInstancesView AS
SELECT MAX(lowerBound),MIN(upperBound),instance_id FROM Iterations GROUP BY instance_id;

CREATE VIEW IF NOT EXISTS AggregateIterations AS
SELECT MAX(lowerBound) AS lowerBound, MIN(upperBound) AS upperBound, MAX(runtime) as runtime,instance_id, Datasets.id AS dataset_id, Datasets.name AS datasetName, Solvers.id AS solver_id, Solvers.algorithmName AS algorithmName, Solvers.algorithmFMC AS algorithmFmc FROM Iterations INNER JOIN Instances ON (Instances.id = instance_id) INNER JOIN Datasets ON (Instances.dataset_id = Datasets.id) INNER JOIN Solvers ON (Solvers.id = solver_id) GROUP BY instance_id,solver_id;

CREATE VIEW IF NOT EXISTS AggregateInstances AS
SELECT AVG(ai.lowerBound) AS lowerBound, AVG(ai.upperBound) AS upperBound, AVG(ai.runtime) AS runtime, Instances.dataset_id AS dataset_id, ai.datasetName AS datasetName, ai.solver_id AS solver_id FROM AggregateIterations AS ai INNER JOIN Instances ON (ai.instance_id = Instances.id)  GROUP BY Instances.dataset_id, ai.solver_id;

CREATE VIEW IF NOT EXISTS MinMaxBoundDatasetsView AS
SELECT MAX(ai.lowerBound) as lowerBound, MIN(ai.upperBound) as upperBound, MIN(ai.runtime) AS runtime, ai.dataset_id AS dataset_id, ai.datasetName AS datasetName FROM AggregateInstances AS ai GROUP BY ai.dataset_id;
)";

   using BaseVisitor = BASE_VISITOR;

public:
   SqliteVisitor(TCLAP::CmdLine& cmd) 
      :
         BaseVisitor(cmd),
         databaseFileArg_("","databaseFile","sqlite database into which to protocolate runtime/iteration vs. primal/dual energy",true,"","file name",cmd),
         datasetNameArg_("","datasetName","name of dataset the input file belongs to",true,"","string",cmd),
         algorithmNameArg_("","algorithmName","name of algorithm",true,"","string",cmd),
         algorithmFMCArg_("","algorithmFMC","FMC of algorithm", true, "", "string", cmd),
         overwriteDbRecordArg_("","overwriteDbRecord","if true: overwrite previous record. if false: if record is present, abort optimization",cmd,false),
         database_(nullptr)
   {
      // get inputFile argument from cmd
      auto argList = cmd.getArgList();
      inputFileArg_ = nullptr;
      for(auto arg : argList) {
         if(arg->getName() == "inputFile") {
            inputFileArg_ = arg;
         }
      }
      if(inputFileArg_ == nullptr) {
         throw std::runtime_error("input file argument must have been specified before");
      }
   }

   ~SqliteVisitor()
   {
      if(database_) {
         sqlite3_close(database_);
      }
   }

   static int CountCallback(void* c_ptr, int argc, char** argv, char**)
   {
      assert(argc == 1);
      *static_cast<int*>(c_ptr) = atoi(argv[0]);
      return 0;
   }

   void BuildDb() 
   {
      int rc = sqlite3_exec(database_, sql.c_str(), nullptr, nullptr, nullptr);
      if(rc != SQLITE_OK) {
         throw std::runtime_error(std::string("Could not create schema: ") + sqlite3_errmsg(database_));
      }
   }

   // assume that conditionSQL returns an id if record is present, otherwise insert record and retrieve id again
   struct IdRet {bool recordPresent; int id;};
   static int RecordIdCallback(void* i, int argc, char** argv, char**) 
   {
      assert(argc >= 1);
      if (argc == 1) {
         *static_cast<IdRet*>(i) = IdRet({true,atoi(argv[0])});
      } else {
         throw std::runtime_error("Only one record should be present");
      }
      return 0;
   }
      
   int ConditionallyInsertById(const std::string& conditionSQL, const std::string& insertSQL)
   {
      int rc;
      IdRet i = {false,0}; // note: callback will not be called if condition is not met
      rc = sqlite3_exec(database_, conditionSQL.c_str(), &RecordIdCallback, &i, nullptr);
      if(rc != SQLITE_OK) {
         throw std::runtime_error("Could not execute " + conditionSQL);
      }
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
      return ConditionallyInsertById(getInstanceId, insertInstance);
   }

   bool CheckIterationsPresent(const int solver_id, const int instance_id)
   {
      const std::string checkIterationsPresent = "SELECT COUNT(*) FROM Iterations WHERE solver_id = '" + std::to_string(solver_id) + "' AND instance_id = '" + std::to_string(instance_id) + "';";
      int i = 0;
      int rc = sqlite3_exec(database_, checkIterationsPresent.c_str(), &CountCallback, &i, nullptr);
      return i>0;
   }

   LpControl begin(LP& lp) // called, after problem is constructed. 
   {
      auto ret = BaseVisitor::begin(lp);
      try {
         databaseFile_ = databaseFileArg_.getValue();
         datasetName_ = datasetNameArg_.getValue();
         algorithmName_ = algorithmNameArg_.getValue();
         algorithmFMC_ = algorithmFMCArg_.getValue();
         overwriteDbRecord_ = overwriteDbRecordArg_.getValue();
      } catch (TCLAP::ArgException &e) {
         std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
         exit(1);
      }

      if(sqlite3_open(databaseFile_.c_str(), &database_) != SQLITE_OK) {
         std::cout << "Could not open database file " << databaseFile_ << ", error: " << sqlite3_errmsg(database_) << "\n";
         ret.error = true;
         return ret;
      }

      BuildDb();

      const std::string inputFile = ExtractFilename( dynamic_cast<TCLAP::ValueArg<std::string>*>(inputFileArg_)->getValue() ); // this is very bad design!
      std::cout << "input file = " << inputFile << "\n";

      solver_id_ = GetSolverId(algorithmName_, algorithmFMC_);
      dataset_id_ = GetDatasetId(datasetName_);
      instance_id_ = GetInstanceId(inputFile, dataset_id_);

      // do zrobienia: functions not in parallel! -> in parallel more than one optimization can take place. make lock in database
      if(!overwriteDbRecord_ && CheckIterationsPresent(solver_id_, instance_id_)) { 
         std::cout << "Not performing optimization, as instance was already optimized with same algorithm\n";
         ret.error = true;
      }

      return ret;
   }

   LpControl visit(LpControl c, const REAL lowerBound, const REAL upperBound)
   {
      auto ret_state = this->BaseVisitor::visit(c, lowerBound, upperBound);
      
      const INDEX timeElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - BaseVisitor::GetBeginTime()).count();
      const INDEX curIter = BaseVisitor::GetIter();
      iterationStatistics_.push_back({curIter,timeElapsed,lowerBound,upperBound});

      return ret_state;
   }

   void end(const REAL lowerBound, const REAL upperBound)
   {
      const INDEX timeElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - BaseVisitor::GetBeginTime()).count();
      const INDEX curIter = BaseVisitor::GetIter();
      iterationStatistics_.push_back({curIter+1,timeElapsed,lowerBound,upperBound}); // additional fake iteration, e.g. for post-processing, collecting primal rounding by external rounding routines etc.

      std::cout << "write bounds to database\n";
      WriteBounds(iterationStatistics_);
   }

   void WriteBounds(const std::vector<IterationStatistics>& iterStats)
   {
      // connection to database is open and we have a corresponding id in ProblemSolverTable
      // Write into Iterations all iterations
      // first remove all iterations previously in this table
      int rc;
      rc = sqlite3_exec(database_,"BEGIN TRANSACTION;", nullptr, nullptr, nullptr);
      assert(rc == SQLITE_OK);
      const std::string rmIterStmt = "DELETE FROM Iterations WHERE solver_id = " + std::to_string(solver_id_) + " AND instance_id = " + std::to_string(instance_id_) + "\n";
      rc = sqlite3_exec(database_, rmIterStmt.c_str(), nullptr, nullptr, nullptr);
      assert(rc == 0);
      for(const auto& it : iterStats) {
         const std::string stmt = "INSERT INTO Iterations (solver_id, instance_id, iteration, runtime, lowerBound, upperBound) VALUES ('" + std::to_string(solver_id_) + "', '" + std::to_string(instance_id_) + "', '" + std::to_string(it.iteration_) + "', '" + std::to_string(it.timeElapsed_) + "', '" + std::to_string(it.lowerBound_) + "', '" + std::to_string(it.upperBound_) + "');";
         rc = sqlite3_exec(database_, stmt.c_str(), nullptr, nullptr, nullptr);
         if(rc != SQLITE_OK) {
            sqlite3_exec(database_,"ROLLBACK TRANSACTION;", nullptr, nullptr, nullptr);
            throw std::runtime_error("Could not insert iteration information: " + stmt);
         }
         assert(rc == 0);
      }
      if(sqlite3_exec(database_,"END TRANSACTION;", nullptr, nullptr, nullptr) != SQLITE_OK) {
         throw std::runtime_error(std::string("Could not commit transaction: ") + sqlite3_errmsg(database_) );
      }
      assert(rc == SQLITE_OK);
   }

private:
   TCLAP::ValueArg<std::string> databaseFileArg_;
   TCLAP::ValueArg<std::string> datasetNameArg_;
   TCLAP::ValueArg<std::string> algorithmNameArg_; // custom name given for algorithm, to differentiate between same algorithm with differing options
   TCLAP::ValueArg<std::string> algorithmFMCArg_; 
   TCLAP::SwitchArg overwriteDbRecordArg_;

   std::string databaseFile_;
   std::string datasetName_;
   std::string algorithmName_;
   std::string algorithmFMC_;
   bool overwriteDbRecord_;
   TCLAP::Arg* inputFileArg_;

   std::vector<IterationStatistics> iterationStatistics_;

   sqlite3* database_;
   int solver_id_;
   int dataset_id_;
   int instance_id_;
};

} // end namespace LP_MP

#endif // LP_MP_SQLITE_VISITOR_HXX

