#ifndef LP_MP_POSTGRESQL_VISITOR_HXX
#define LP_MP_POSTGRESQL_VISITOR_HXX

#include <pqxx/pqxx>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "standard_visitor.hxx"
#include "help_functions.hxx"

namespace LP_MP {
  struct Credentials {
    std::string user = "";
    std::string pwd = "";
  
    Credentials& operator=(const std::string &str)
    {
      size_t pos = str.find(":");
      if( pos == std::string::npos ){
        throw TCLAP::ArgParseException(str + ", please type username:password");
      }
      if( pos == str.size()-1 || pos == 0 ){
        throw TCLAP::ArgParseException(str + ", user or password is missing");
      }
      user = str.substr(0,pos);
      pwd = str.substr(pos+1);
    
      return *this;
    }
  };
}

namespace TCLAP {
  template<>
  struct ArgTraits<LP_MP::Credentials> {
    typedef StringLike ValueCategory;
  };
}

namespace LP_MP { 

  struct IterationStatistics {
    INDEX iteration_;
    INDEX timeElapsed_; // in milliseconds
    REAL lowerBound_;
    REAL upperBound_;
  };
  
  // this visitor connects to given sqlite database and writes or updates the runtime and iteration data of the algorithm.

  template<class BASE_VISITOR = StandardVisitor>
  class PostgresqlVisitor : public BASE_VISITOR {

    using BaseVisitor = BASE_VISITOR;

  public:
    PostgresqlVisitor(TCLAP::CmdLine& cmd) 
      :
      BaseVisitor(cmd),
      databaseIpArg_("","databaseIp","address to postgresql database into which to protocolate runtime/iteration vs. primal/dual energy",true,"","ip or dns",cmd),
      databaseNameArg_("","databaseName","name of postgres database into which to protocolate runtime/iteration vs. primal/dual energy",true,"","db name",cmd),
      databaseAuthArg_("","databaseAuth","credentials (username:password) to write in postgresql database into which to protocolate runtime/iteration vs. primal/dual energy",true,Credentials(),"string:string",cmd),
      datasetNameArg_("","datasetName","name of dataset the input file belongs to",true,"","string",cmd),
      algorithmNameArg_("","algorithmName","name of algorithm",true,"","string",cmd),
      algorithmFMCArg_("","algorithmFMC","FMC of algorithm", true, "", "string", cmd),
      overwriteDbRecordArg_("","overwriteDbRecord","if true: overwrite previous record. if false: if record is present, abort optimization",cmd,false)
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

      ~PostgresqlVisitor()
    {
      if(database_) { 
        try{
          database_->disconnect();
        }catch(const std::exception &e){
          std::cerr << e.what() << std::endl;
        } 
      }
    }

    static int CountCallback(void* c_ptr, int argc, char** argv, char**)
    {
      assert(argc == 1);
      *static_cast<int*>(c_ptr) = atoi(argv[0]);
      return 0;
    }

    void ExecuteStatement(std::vector<std::string> stmts){
      pqxx::work transaction(*database_);
      for( auto stmt : stmts ){
        transaction.exec(stmt.c_str());
      }
      transaction.commit();
    }

    void ExecuteStatement(std::string stmt){
      pqxx::work transaction(*database_);
      transaction.exec(stmt.c_str());
      transaction.commit();
    }
    
    void CreateTable(const std::string tableName, const std::vector<std::string>& columns)
    {
      std::string createTable = "CREATE TABLE IF NOT EXISTS " + tableName + "(";
      for(INDEX i=0; i<columns.size(); ++i) {
        createTable += columns[i];
        if(i < columns.size()-1) {
          createTable += ", ";
        }
      }
      createTable += ");";
      //std::cout << createTable << "\n";
      ExecuteStatement(createTable.c_str());
    }
    void ConditionallyCreateTable(const std::string tableName, const std::vector<std::string>& columns)
    {
      CreateTable(tableName, columns);
    }

    bool ViewExists(const std::string& viewName) 
    {
      return false;
    }

    void CreateView(const std::string& viewName, const std::string command)
    {
      const std::string sql = "CREATE OR REPLACE VIEW " + viewName + " AS " + command;
      ExecuteStatement(sql.c_str());
    }
    void ConditionallyCreateView(const std::string& viewName, const std::string command)
    {
      CreateView(viewName, command);
    }

    bool IndexExists(const std::string& indexName) 
    {
      return false;
    }
    void CreateIndex(const std::string& indexName, const std::string command)
    {
      const std::string sql = "CREATE INDEX IF NOT EXISTS " + indexName + " ON " + command;
      ExecuteStatement(sql.c_str());
    }
    void ConditionallyCreateIndex(const std::string& indexName, const std::string command)
    {
      CreateIndex(indexName, command);
    }
    void BuildDb() 
    {
      ConditionallyCreateTable("Solvers", 
                               {"id SERIAL PRIMARY KEY NOT NULL",
                                   "algorithmName TEXT NOT NULL",
                                   "algorithmFMC TEXT NOT NULL",
                                   "UNIQUE(algorithmName, algorithmFMC)"});
      ConditionallyCreateTable("Datasets", 
                               {"id SERIAL PRIMARY KEY NOT NULL",
                                   "name TEXT NOT NULL",
                                   "UNIQUE(name)"});
      ConditionallyCreateTable("Instances", 
                               {"id SERIAL PRIMARY KEY NOT NULL",
                                   "dataset_id INT",
                                   "name TEXT NOT NULL",
                                   "FOREIGN KEY(dataset_id) REFERENCES Datasets(id)",
                                   "UNIQUE(dataset_id, name)"});

      ConditionallyCreateTable("Instances_Status", 
                               {   "solver_id INT NOT NULL",
                                   "instance_id INT NOT NULL",
                                   "status TEXT NOT NULL",
                                   "FOREIGN KEY(solver_id) REFERENCES Solvers(id)",
                                   "FOREIGN KEY(instance_id) REFERENCES Instances(id)",
                                   "UNIQUE(solver_id,instance_id)"});


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
 
      ConditionallyCreateView("LowerBoundView", 
                              "SELECT MAX(lowerBound) as lowerbound,solver_id,instance_id FROM Iterations GROUP BY solver_id,instance_id;"
                              );
      // maximum lower bound over all solvers
      ConditionallyCreateView("MaxLowerBoundView",
                              "SELECT MAX(lowerBound) as lowerbound,instance_id FROM Iterations GROUP BY instance_id;"
                              );

      ConditionallyCreateView("UpperBoundView", 
                              "SELECT MIN(upperBound) as upperbound,solver_id,instance_id FROM Iterations GROUP BY solver_id,instance_id;"
                              );
      ConditionallyCreateView("MinUpperBoundView",
                              "SELECT MIN(upperBound) as upperbound,instance_id FROM Iterations GROUP BY instance_id;"
                              );
      ConditionallyCreateView("RuntimeView", 
                              "SELECT MAX(runtime) as runtime,solver_id,instance_id FROM Iterations GROUP BY solver_id,instance_id;"
                              );
      ConditionallyCreateView("AggregateIterationsHelper",
                              " SELECT LowerBoundView.lowerBound AS lowerBound, UpperBoundView.upperBound AS upperBound," \
                              " RuntimeView.runtime AS runtime, LowerBoundView.solver_id AS solver_id, LowerBoundView.instance_id AS instance_id"
                              " FROM LowerBoundView"
                              " INNER JOIN UpperBoundView ON (LowerBoundView.solver_id = UpperBoundView.solver_id AND LowerBoundView.instance_id = UpperBoundView.instance_id)"
                              " INNER JOIN RuntimeView ON (LowerBoundView.solver_id = RuntimeView.solver_id AND LowerBoundView.instance_id = RuntimeView.instance_id)"
                              ";"); 
      ConditionallyCreateView("MinMaxBoundInstancesView",
                              "SELECT MAX(lowerBound),MIN(upperBound),instance_id FROM Iterations GROUP BY instance_id;"
                              );
      ConditionallyCreateView("AggregateIterations", 
                              "SELECT MAX(lowerBound) AS lowerBound, MIN(upperBound) AS upperBound, MAX(runtime) as runtime,instance_id, Datasets.id AS dataset_id, Datasets.name AS datasetName, Solvers.id AS solver_id, Solvers.algorithmName AS algorithmName, Solvers.algorithmFMC AS algorithmFmc " \
                              "FROM Iterations INNER JOIN Instances ON (Instances.id = instance_id) INNER JOIN Datasets ON (Instances.dataset_id = Datasets.id) INNER JOIN Solvers ON (Solvers.id = solver_id) GROUP BY instance_id,solver_id,datasets.id,solvers.id;");
      ConditionallyCreateView("AggregateInstances",
                              "SELECT AVG(ai.lowerBound) AS lowerBound, AVG(ai.upperBound) AS upperBound, AVG(ai.runtime) AS runtime, Instances.dataset_id AS dataset_id, ai.datasetName AS datasetName, ai.solver_id AS solver_id FROM AggregateIterations AS ai INNER JOIN Instances ON (ai.instance_id = Instances.id)  GROUP BY Instances.dataset_id, ai.solver_id,ai.datasetname;"); 
      ConditionallyCreateView("MinMaxBoundDatasetsView",
                              "SELECT MAX(ai.lowerBound) as lowerBound, MIN(ai.upperBound) as upperBound, MIN(ai.runtime) AS runtime, ai.dataset_id AS dataset_id, ai.datasetName AS datasetName FROM AggregateInstances AS ai GROUP BY ai.dataset_id,ai.datasetname;");
    }

    int ConditionallyInsertById(const std::string& select,const std::string& insert){
      pqxx::work transaction(*database_);
      
      pqxx::result res_insert = transaction.exec(select.c_str());
      if( res_insert.size() == 0 ){
        transaction.exec(insert.c_str()); 
      }
      
      pqxx::result res = transaction.exec(select.c_str());
      if( res.size() == 0 ){
        throw std::runtime_error("Insert command '" + insert + " did not work!");
      }

      auto rows = res.begin();
      int id = rows[0].as<int>();
      transaction.commit();
      return id;
    }
    
    // return id of solver
    int GetSolverId(const std::string& algorithmName, const std::string& algorithmFMC)
    {
      const std::string getSolverId = "SELECT (id) FROM Solvers WHERE algorithmName = '" + algorithmName + "' and algorithmFMC = '" + algorithmFMC + "';";
      const std::string insertSolver = "INSERT INTO Solvers (algorithmName, algorithmFMC) VALUES ('" + algorithmName + "', '" + algorithmFMC + "')" \
        " ON CONFLICT (algorithmName, algorithmFMC) DO NOTHING;";
      return ConditionallyInsertById(getSolverId,insertSolver);
    }
    int GetDatasetId(const std::string& dataset)
    {
      const std::string getDatasetId = "SELECT (id) FROM Datasets WHERE name = '" + dataset + "';";
      const std::string insertDataset = "INSERT INTO Datasets (name) VALUES ('" + dataset + "')" \
        " ON CONFLICT (name) DO NOTHING;";
      return ConditionallyInsertById(getDatasetId, insertDataset);
    }
    int GetInstanceId(const std::string& instance, const int dataset_id)
    {
      const std::string getInstanceId = "SELECT (id) FROM Instances WHERE name = '" + instance + "' AND dataset_id = " + std::to_string(dataset_id) + ";";
      const std::string insertInstance = "INSERT INTO Instances (name, dataset_id) VALUES ('" + instance + "', " + std::to_string(dataset_id) + ")" \
        " ON CONFLICT (name,dataset_id) DO NOTHING;";
      return ConditionallyInsertById(getInstanceId, insertInstance);
    }

    bool CheckIterationsPresent(const int solver_id, const int instance_id)
    {
      pqxx::nontransaction transaction(*database_);
      const std::string checkIterationsPresent = "SELECT COUNT(*) FROM Iterations WHERE solver_id = '" + \
        std::to_string(solver_id) + "' AND instance_id = '" + std::to_string(instance_id) + "';";
      pqxx::result res = transaction.exec(checkIterationsPresent.c_str());
      if( res.size() != 1 ){ throw std::runtime_error("CheckIterations received nothing!");  }
      auto rows = res.begin();
      return rows[0].as<int>() > 0;
    }
  
    LpControl begin(LP& lp) // called, after problem is constructed. 
    {
      auto ret = BaseVisitor::begin(lp);
      try {
        databaseIp_ = databaseIpArg_.getValue();
        databaseName_ = databaseNameArg_.getValue();
        datasetName_ = datasetNameArg_.getValue();
        databaseAuth_ = databaseAuthArg_.getValue();
        algorithmName_ = algorithmNameArg_.getValue();
        algorithmFMC_ = algorithmFMCArg_.getValue();
        overwriteDbRecord_ = overwriteDbRecordArg_.getValue();
      } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
        exit(1);
      }

      std::string con = "dbname=" + databaseName_ + " user=" + databaseAuth_.user + " " \
        "password=" + databaseAuth_.pwd + " hostaddr=" + databaseIp_ + " port=5432 sslmode=require";
      std::unique_ptr<pqxx::connection> database (new pqxx::connection(con.c_str()));
      database_ = std::move(database);
      if( !database_->is_open() ) {
        std::cout << "Could not open database " << databaseName_ << " on ";
        std::cout << databaseIp_ << " with " << databaseAuth_.user << ":" << databaseAuth_.pwd << std::endl;
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
        ret.error = true;
      } else {
        const std::string InsertStatus = "INSERT INTO Instances_Status (solver_id,instance_id,status) VALUES " \
          "(" + std::to_string(solver_id_) + "," + std::to_string(instance_id_) + ",'processing') " \
          "ON CONFLICT (solver_id,instance_id) DO NOTHING RETURNING *;";

        pqxx::work transaction(*database_);
        pqxx::result res = transaction.exec(InsertStatus.c_str());
        transaction.commit();        

        if(res.size() == 0){
          std::cout << "Not performing optimization, as instance is processing or processed\n";
          ret.error = true;
        }
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
      iterationStatistics_.push_back({curIter+1,timeElapsed,lowerBound,upperBound}); // additional fake iteration, e.g. for post-processing

      std::cout << "write bounds to database\n";
      WriteBounds(iterationStatistics_);
    }

    void WriteBounds(const std::vector<IterationStatistics>& iterStats)
    {
      // connection to database is open and we have a corresponding id in ProblemSolverTable
      // Write into Iterations all iterations
      // first remove all iterations previously in this table
    
      pqxx::work transaction(*database_);
      const std::string rmIterStmt = "DELETE FROM Iterations WHERE solver_id = " + std::to_string(solver_id_) + \
        " AND instance_id = " + std::to_string(instance_id_) + "\n"; 

      transaction.exec(rmIterStmt.c_str());
      
      for(const auto& it : iterStats) {
        const std::string stmt = "INSERT INTO Iterations (solver_id, instance_id, iteration, runtime, lowerBound, upperBound) VALUES ('" + \
          std::to_string(solver_id_) + "', '" + std::to_string(instance_id_) + "', '" + std::to_string(it.iteration_) + "', '" + \
          std::to_string(it.timeElapsed_) + "', '" + std::to_string(it.lowerBound_) + "', '" + std::to_string(it.upperBound_) + "');";

        transaction.exec(stmt.c_str());
      }
 
      const std::string InsertStatus= "UPDATE instances_status "\
        "SET status = 'processed' " \
        "WHERE solver_id = " + std::to_string(solver_id_) + " " + \
        "AND instance_id = " + std::to_string(instance_id_) + " " + \
        "RETURNING * ;";
      pqxx::result res = transaction.exec(InsertStatus.c_str());
      if( res.size() != 1 ){
        throw std::runtime_error("solver_id " + std::to_string(solver_id_) + ", " + "instance_id" + std::to_string(instance_id_) + ": " \
                                 "It occurs a problem by updating the instance status after optimizing -> " + std::to_string(res.size()));
      }      
      transaction.commit();
    }

  private:
    TCLAP::ValueArg<std::string> databaseIpArg_;
    TCLAP::ValueArg<std::string> databaseNameArg_;
    TCLAP::ValueArg<Credentials> databaseAuthArg_;
    TCLAP::SwitchArg overwriteDbRecordArg_;
    
    TCLAP::ValueArg<std::string> datasetNameArg_;
    TCLAP::ValueArg<std::string> algorithmNameArg_; // custom name given for algorithm, to differentiate between same algorithm with differing options
    TCLAP::ValueArg<std::string> algorithmFMCArg_; 
    
    std::string databaseIp_;
    std::string databaseName_;
    Credentials databaseAuth_;
    bool overwriteDbRecord_;
    
    std::string datasetName_;
    std::string algorithmName_;
    std::string algorithmFMC_;
    
    TCLAP::Arg* inputFileArg_;

    std::vector<IterationStatistics> iterationStatistics_;

    std::unique_ptr<pqxx::connection> database_;
       
    int solver_id_;
    int dataset_id_;
    int instance_id_;
  };

} // end namespace LP_MP

#endif // LP_MP_POSTGRESQL_VISITOR_HXX

