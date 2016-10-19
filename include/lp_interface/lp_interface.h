#ifndef LP_MP_LP_INTERFACE_HXX
#define LP_MP_LP_INTERFACE_HXX

#include "config.hxx"
#include "primal_solution_storage.hxx"

#ifdef USE_GUROBI
#include <gurobi_c++.h>
#elif USE_CPLEX
#include <ilcplex/ilocplex.h>
#elif USE_LOCALSOLVER
#include "localsolver.h" 
#else
#include "lp_auxiliary.hxx"
#endif

namespace LP_MP {

#ifdef USE_GUROBI
  typedef GRBVar LpVariable;
  typedef GRBLinExpr LinExpr;
#elif USE_CPLEX
  typedef IloNumVar LpVariable;
  typedef IloNumExpr LinExpr;
#elif USE_LOCALSOLVER
  typedef localsolver::LSExpression LpVariable;
  typedef localsolver::LSExpression LinExpr;
#endif

    
  class LpInterfaceAdapter {
  public:
 
    LpInterfaceAdapter(){ } 
    
    virtual ~LpInterfaceAdapter(){ }
    
    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    LpInterfaceAdapter( FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd, 
                        MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd,
                        REAL epsilon = std::numeric_limits<REAL>::infinity(),
                        bool integer = true)
      : epsilon_(epsilon),primal_(PrimalSolutionStorage(factorBegin,factorEnd)),integer_(integer),
        LpInterfaceAdapter(factorBegin,factorEnd,messageBegin,messageEnd,primal_,epsilon_,integer_)
    { } 
    
    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    LpInterfaceAdapter( FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd, 
                        MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd,
                        PrimalSolutionStorage primal,
                        REAL epsilon = std::numeric_limits<REAL>::infinity(),
                        bool integer = true)
      : epsilon_(epsilon),primal_(primal),integer_(integer)
    { }
    
    virtual LinExpr CreateLinExpr() = 0;
    
    virtual void CreateMainVariable(REAL lb, REAL ub,bool integer = false) = 0;
    virtual void CreateMainVariables(INDEX n,REAL lb, REAL ub,bool integer = false) = 0;

    virtual void CreateAuxVariable(REAL lb, REAL ub,bool integer = false) = 0;
    virtual void CreateAuxVariables(INDEX n,REAL lb, REAL ub,bool integer = false) = 0;
    
    virtual REAL GetVariableValue(const INDEX i) const = 0;
    virtual REAL GetObjectiveValue() const = 0;
    virtual REAL GetBestBound() const = 0;
    
    virtual void SetTimeLimit(REAL t) = 0;
    virtual void SetNumberOfThreads(INDEX t) = 0;
    virtual void SetDisplayLevel(INDEX t) = 0;
    
    virtual void addLinearEquality(LinExpr lhs,LinExpr rhs) = 0;
    virtual void addLinearInequality(LinExpr lhs,LinExpr rhs) = 0;
    
    virtual INDEX solve() = 0;
    
    virtual void AddObjective(INDEX i,REAL value){ /* Have to be reimplemented by the subclass */
      primal_[Offset_ + i] = unknownState;
    }

    bool IsObjective(INDEX i) const { // do not use it in messages!
      assert(primal_.size() > 0);
      assert(Offset_ + i < primal_.size());
      if(primal_[Offset_ + i] == false){
        return false;
      } else {
        return true;
      }
    }

    bool IsLeftObjective(INDEX i){
      assert(primal_.size() > 0);
      assert(OffsetLeft_ + i < primal_.size());
      if(primal_[OffsetLeft_ + i] == false){
        return false;
      } else {
        return true;
      }
    }

    bool IsRightObjective(INDEX i){
      assert(primal_.size() > 0);
      assert(OffsetRight_ + i < primal_.size());
      if(primal_[OffsetRight_ + i] == false){
        return false;
      } else {
        return true;
      }
    }

    unsigned char IsPrimal(INDEX i){ // do not use it in messages!
      return primal_[Offset_ + i];
    }

    unsigned char IsLeftPrimal(INDEX i){
      return primal_[OffsetLeft_ + i];
    }

    unsigned char IsRightPrimal(INDEX i){
      return primal_[OffsetRight_ + i];
    }
    
    INDEX GetFactorSize() const { return size_; }
    INDEX GetLeftFactorSize() const { return leftSize_; }
    INDEX GetRightFactorSize() const { return rightSize_; }
    
    LpVariable GetVariable(const INDEX i) const 
    { assert(i < size_);  assert(Offset_ + i < noVars_); return MainVars_[Offset_ + i]; }
    
    LpVariable GetLeftVariable(const INDEX i) const 
    { assert(i < leftSize_); assert(OffsetLeft_ + i < noVars_); return MainVars_[OffsetLeft_ + i]; }
    
    LpVariable GetRightVariable(const INDEX i) const 
    { assert(i < rightSize_); assert(OffsetRight_ + i < noVars_); return MainVars_[OffsetRight_ + i]; }
    
    LpVariable GetAuxVariable(const INDEX i) const { assert(i < noAuxVars_); return MainAuxVars_[OffsetAux_ + i]; }

    REAL GetEpsilon() const { return epsilon_; }
    
    std::vector<REAL>& GetRepam(){ return pot_; };
    
  protected:
    std::vector<LpVariable> MainVars_;
    std::vector<LpVariable> MainAuxVars_;
    std::vector<REAL> pot_;
    PrimalSolutionStorage primal_;
    bool integer_;
  
    INDEX Offset_,OffsetAux_,OffsetLeft_,OffsetRight_;
    INDEX size_,sizeAux_,leftSize_,rightSize_;
    INDEX noVars_,noAuxVars_;    
    REAL epsilon_;
    
    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    void CreateVariables(FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd, MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd){
      
      /* Create Variables */
      MainVars_ = std::vector<LpVariable>();
      for(auto factorIt = factorBegin; factorIt != factorEnd; ++factorIt) {
        Offset_ = factorIt->GetPrimalOffset();
        for(INDEX i=0;i<factorIt->size();i++){
          if(primal_[Offset_ + i] == true){
            CreateMainVariable(1.0,1.0,integer_);
          }
          else if(primal_[Offset_ + i] == false){
            CreateMainVariable(0.0,0.0,integer_);
          } 
          else {
            CreateMainVariable(0.0,1.0,integer_);
            primal_[Offset_ + i] = false;
          }
        }
        
        factorIt->SetAuxOffset(MainAuxVars_.size()); // <-- factoradapter provides this
        factorIt->CreateAuxVariables(this); // <-- factor provides this
      }
      noVars_ = MainVars_.size();
      noAuxVars_ = MainAuxVars_.size();
    }
    
    template<typename FACTOR_ITERATOR, typename MESSAGE_ITERATOR>
    void Build(FACTOR_ITERATOR factorBegin, FACTOR_ITERATOR factorEnd, MESSAGE_ITERATOR messageBegin, MESSAGE_ITERATOR messageEnd){
      
      /* Process Factors */
      //printf("Process Factors\n");
      for(auto factorIt = factorBegin; factorIt != factorEnd; ++factorIt) {
        Offset_ = factorIt->GetPrimalOffset();          // <-- factoradapter provides this
        size_ = factorIt->size();                       // <-- factoradapter provides this
        OffsetAux_ = factorIt->GetAuxOffset();          // <-- factoradapter provides this
        //sizeAux_ = factorIt->GetNumberOfAuxVariables(); // <-- factor provides this
        pot_ = factorIt->GetReparametrizedPotential();  // <-- factoradapter provides this
        //printf("Factor\n");
        factorIt->CreateConstraints(this);              // <-- factor provides this
      }
            
      /* Process Messages */
      for(auto messageIt = messageBegin; messageIt != messageEnd; ++messageIt) {
        leftSize_ = messageIt->GetLeftFactor()->size();               // <-- messageadapter provides this
        rightSize_ = messageIt->GetRightFactor()->size();             // <-- messageadapter provides this
        OffsetLeft_ = messageIt->GetLeftFactor()->GetPrimalOffset();  // <-- messageadapter provides this
        OffsetRight_ = messageIt->GetRightFactor()->GetPrimalOffset();// <-- messageadapter provides this
        messageIt->CreateConstraints(this);                           // <-- message provides this
      }
    }
  };
  
}

#endif // LP_MP_LP_INTERFACE_HXX
