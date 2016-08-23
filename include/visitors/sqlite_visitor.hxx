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

   bool TableExists(const std::string& tableName)
   {
      int rc, ret = 0; // callback is not called when select does not return anything
      const std::string existsTable = "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name='" + tableName + "';";
      rc = sqlite3_exec(database_, existsTable.c_str(), &CountCallback, &ret, NULL);
      if(rc != SQLITE_OK || !(rc == 0 || rc == 1)) {
         throw std::runtime_error("Could not determine whether table " + tableName + " exists in database");
      }  
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
      if(rc != SQLITE_OK) {
         throw std::runtime_error("Could not create table " + tableName);
      }
      assert(rc == 0);
   }
   void ConditionallyCreateTable(const std::string tableName, const std::vector<std::string>& columns)
   {
      if(!TableExists(tableName)) { CreateTable(tableName, columns); }
   }

   bool ViewExists(const std::string& viewName) 
   {
      int rc, ret = 0;
      const std::string existsView = "SELECT COUNT(*) FROM sqlite_master WHERE type='view' AND name='" + viewName + "';";
      rc = sqlite3_exec(database_, existsView.c_str(), &CountCallback, &ret, NULL);
      if(rc != SQLITE_OK || !(rc == 0 || rc == 1)) {
         throw std::runtime_error("Could not determine whether view " + viewName + " exists in database");
      }  
      assert(rc == 0);
      assert(ret == 0 || ret == 1);
      return ret == 1;
   }
   void CreateView(const std::string& viewName, const std::string command)
   {
      const std::string sql = "CREATE VIEW " + viewName + " AS " + command;
      int rc = sqlite3_exec(database_, sql.c_str(), nullptr, nullptr, nullptr);
      if(rc != SQLITE_OK) {
         throw std::runtime_error("Could not create view " + viewName + ": " + sql);
      }
      assert(rc == 0);
   }
   void ConditionallyCreateView(const std::string& viewName, const std::string command)
   {
      if(!ViewExists(viewName)) { CreateView(viewName, command); }
   }

   bool IndexExists(const std::string& indexName) 
   {
      int rc, ret = 0;
      const std::string existsIndex = "SELECT COUNT(*) FROM sqlite_master WHERE type='index' AND name='" + indexName + "';";
      rc = sqlite3_exec(database_, existsIndex.c_str(), &CountCallback, &ret, NULL);
      if(rc != SQLITE_OK || !(rc == 0 || rc == 1)) {
         throw std::runtime_error("Could not determine whether index " + indexName + " exists in database");
      }  
      assert(rc == 0);
      assert(ret == 0 || ret == 1);
      return ret == 1;
   }
   void CreateIndex(const std::string& indexName, const std::string command)
   {
      const std::string sql = "CREATE INDEX " + indexName + " ON " + command;
      int rc = sqlite3_exec(database_, sql.c_str(), nullptr, nullptr, nullptr);
      if(rc != SQLITE_OK) {
         throw std::runtime_error("Could not create index " + indexName + ": " + sql);
      }
   }
   void ConditionallyCreateIndex(const std::string& indexName, const std::string command)
   {
      if(!IndexExists(indexName)) { CreateIndex(indexName, command); }
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
      ConditionallyCreateIndex("InstancesIndex","Instances(id,dataset_id)");
      ConditionallyCreateTable("Iterations", 
            {"instance_id INTEGER NOT NULL",
            "solver_id INTEGER NOT NULL",
            "iteration INTEGER NOT NULL",
            "runtime INT",
            "lowerBound DOUBLE PRECISION", // or REAL
            "upperBound DOUBLE PRECISION",
            "FOREIGN KEY(instance_id) REFERENCES Instances(id)",
            "FOREIGN KEY(solver_id) REFERENCES Solvers(id)",
            "UNIQUE(instance_id, solver_id, iteration)"});
      // do zrobienia: or should one make a covering index?
      ConditionallyCreateIndex("IterationsIndex1","Iterations(instance_id,solver_id,lowerBound,upperBound,runtime)");
      ConditionallyCreateIndex("InstancesIndex1","Instances(id,dataset_id)");
      //ConditionallyCreateIndex("IterationsIndex2","Iterations(instance_id,solver_id,lowerBound)");
      //ConditionallyCreateIndex("IterationsIndex3","Iterations(instance_id,solver_id,upperBound)");
      //ConditionallyCreateIndex("IterationsIndex4","Iterations(instance_id,solver_id,runtime)");

      ConditionallyCreateView("LowerBoundView", 
"SELECT lowerBound,iteration,solver_id,instance_id FROM Iterations AS a WHERE NOT EXISTS (SELECT 1 FROM Iterations as b WHERE a.solver_id = b.solver_id AND a.instance_id = b.instance_id AND a.lowerBound < b.lowerBound) GROUP BY solver_id, instance_id;");
      // maximum lower bound over all solvers
      ConditionallyCreateView("MaxLowerBoundView",
"SELECT MAX(lowerBound) as lowerBound, instance_id FROM Iterations AS a WHERE NOT EXISTS (SELECT 1 FROM Iterations as b WHERE a.instance_id = b.instance_id AND a.lowerBound < b.lowerBound) GROUP BY instance_id;");
      ConditionallyCreateView("UpperBoundView", 
"SELECT upperBound,iteration,solver_id,instance_id FROM Iterations AS a WHERE NOT EXISTS (SELECT 1 FROM Iterations AS b WHERE a.solver_id = b.solver_id AND a.instance_id = b.instance_id AND a.upperBound > b.upperBound) GROUP BY solver_id, instance_id");
      ConditionallyCreateView("MinUpperBoundView",
"SELECT MIN(upperBound) as upperBound, instance_id FROM Iterations AS a WHERE NOT EXISTS (SELECT 1 FROM Iterations as b WHERE a.instance_id = b.instance_id AND a.upperBound > b.upperBound) GROUP BY instance_id;");
      ConditionallyCreateView("RuntimeView", 
"SELECT runtime,iteration,solver_id,instance_id FROM Iterations AS a WHERE NOT EXISTS (SELECT 1 FROM Iterations AS b WHERE a.solver_id = b.solver_id AND a.instance_id = b.instance_id AND a.runtime < b.runtime) GROUP BY solver_id, instance_id;");
      ConditionallyCreateView("AggregateIterationsHelper",
" SELECT LowerBoundView.lowerBound AS lowerBound, UpperBoundView.upperBound AS upperBound, RuntimeView.runtime AS runtime, LowerBoundView.solver_id AS solver_id, LowerBoundView.instance_id AS instance_id"
" FROM LowerBoundView "
" INNER JOIN UpperBoundView ON (LowerBoundView.solver_id = UpperBoundView.solver_id AND LowerBoundView.instance_id = UpperBoundView.instance_id)"
" INNER JOIN RuntimeView ON (LowerBoundView.solver_id = RuntimeView.solver_id AND LowerBoundView.instance_id = RuntimeView.instance_id)"
";"); 
      ConditionallyCreateView("MinMaxBoundInstancesView",
" SELECT MaxLowerBoundView.lowerBound, MinUpperBoundView.upperBound, RuntimeView.runtime, RuntimeView.iteration, MaxLowerBoundView.instance_id"
" FROM MaxLowerBoundView "
" INNER JOIN MinUpperBoundView ON (MaxLowerBoundView.instance_id = MinUpperBoundView.instance_id)"
" INNER JOIN RuntimeView ON (MaxLowerBoundView.instance_id = RuntimeView.instance_id)"
";");
      ConditionallyCreateView("AggregateIterations", 
" SELECT AggregateIterationsHelper.lowerBound, AggregateIterationsHelper.upperBound, AggregateIterationsHelper.runtime, AggregateIterationsHelper.instance_id, AggregateIterationsHelper.solver_id, Instances.name AS instanceName, Datasets.id AS dataset_id, Datasets.name AS datasetName, Solvers.algorithmName, Solvers.algorithmFmc"
" FROM AggregateIterationsHelper"
" INNER JOIN Instances ON (Instances.id = AggregateIterationsHelper.instance_id)"
" INNER JOIN Datasets ON (Instances.dataset_id = Datasets.id)"
" INNER JOIN Solvers ON (Solvers.id = AggregateIterationsHelper.solver_id)"
";");
      ConditionallyCreateView("AggregateInstances",
" SELECT AVG(ai.lowerBound) AS lowerBound, AVG(ai.upperBound) AS upperBound, AVG(ai.runtime) AS runtime, ai.dataset_id AS dataset_id, ai.datasetName AS datasetName, ai.solver_id AS solver_id FROM AggregateIterations AS ai "
" WHERE (SELECT 1 FROM Instances WHERE Instances.dataset_id = ai.dataset_id AND Instances.id = ai.instance_id) "
//" INNER JOIN Instances ON (Instances.dataset_id = ai.dataset_id AND Instances.id = ai.instance_id)"
" GROUP BY ai.dataset_id, ai.solver_id"
";");
      ConditionallyCreateView("MinMaxBoundDatasetsView",
"SELECT MAX(ai.lowerBound) AS lowerBound, MIN(ai.upperBound) AS upperBound, MIN(ai.runtime) AS runtime, ai.dataset_id AS dataset_id, ai.datasetName AS datasetName FROM AggregateInstances AS ai "
" WHERE (SELECT 1 FROM AggregateInstances WHERE AggregateInstances.dataset_id = ai.dataset_id)"
" GROUP BY ai.dataset_id"
";");

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

      if(!overwriteDbRecord_ && CheckIterationsPresent(solver_id_, instance_id_)) { 
         std::cout << "Not performing optimization, as instance was already optimized with same algorithm\n";
         ret.end = true;
      }

      return ret;
   }

   LpControl visit(LpControl c, const REAL lowerBound, const REAL upperBound)
   {
      auto ret_state = this->BaseVisitor::visit(c, lowerBound, upperBound);
      
      if(!(c.error || c.end)) {
         const INDEX timeElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - BaseVisitor::GetBeginTime()).count();
         const INDEX curIter = BaseVisitor::GetIter();
         iterationStatistics_.push_back({curIter,timeElapsed,lowerBound,upperBound});
      }
      if(c.end) {
         std::cout << "write bounds to database\n";
         WriteBounds(iterationStatistics_);
      }

      return ret_state;
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

