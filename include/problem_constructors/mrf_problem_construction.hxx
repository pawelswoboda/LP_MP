#ifndef LP_MP_MRF_PROBLEM_CONSTRUCTION_HXX
#define LP_MP_MRF_PROBLEM_CONSTRUCTION_HXX

#include "solver.hxx"
#include "cycle_inequalities.hxx"
#include "parse_rules.h"
#include "pegtl/parse.hh"
#include "tree_decomposition.hxx"
#include "arboricity.h"

#include <string>

namespace LP_MP {

// expects simplex factor as unary and pairwise factors and marg message such that unary factor is on the left side and pairwise factor is on the right side
// possibly use static inheritance instead of virtual functions
template<class FACTOR_MESSAGE_CONNECTION, INDEX UNARY_FACTOR_NO, INDEX PAIRWISE_FACTOR_NO, INDEX LEFT_MESSAGE_NO, INDEX RIGHT_MESSAGE_NO>
class MRFProblemConstructor {
public:
   using FMC = FACTOR_MESSAGE_CONNECTION;

   using UnaryFactorContainer = meta::at_c<typename FMC::FactorList, UNARY_FACTOR_NO>;
   using UnaryFactorType = typename UnaryFactorContainer::FactorType;
   using PairwiseFactorContainer = meta::at_c<typename FMC::FactorList, PAIRWISE_FACTOR_NO>;
   using PairwiseFactorType = typename PairwiseFactorContainer::FactorType;
   using LeftMessageContainer = typename meta::at_c<typename FMC::MessageList, LEFT_MESSAGE_NO>::MessageContainerType;
   using LeftMessageType = typename LeftMessageContainer::MessageType;
   using RightMessageContainer = typename meta::at_c<typename FMC::MessageList, RIGHT_MESSAGE_NO>::MessageContainerType;
   using RightMessageType = typename RightMessageContainer::MessageType;


   template<typename SOLVER>
   MRFProblemConstructor(SOLVER& solver) : lp_(&solver.GetLP()) {}

   virtual void ConstructUnaryFactor(UnaryFactorType& u, const std::vector<REAL>& cost) = 0;
   virtual void ConstructPairwiseFactor(PairwiseFactorType& p, const INDEX leftDim, const INDEX rightDim) = 0;
   virtual RightMessageType ConstructRightUnaryPairwiseMessage(UnaryFactorContainer* const right, PairwiseFactorContainer* const p) = 0;
   virtual LeftMessageType ConstructLeftUnaryPairwiseMessage(UnaryFactorContainer* const right, PairwiseFactorContainer* const p) = 0;

   UnaryFactorContainer* AddUnaryFactor(const std::vector<REAL>& cost)
   {
      return AddUnaryFactor(unaryFactor_.size(), cost);
   }
   UnaryFactorContainer* AddUnaryFactor(const INDEX node_number, const std::vector<REAL>& cost)
   {
      //UnaryFactorContainer* u = new UnaryFactorContainer(UnaryFactorType(cost), cost);
      auto* u = new UnaryFactorContainer( cost.size() );
      ConstructUnaryFactor( *(u->GetFactor()), cost );
      if(node_number >= unaryFactor_.size()) {
         unaryFactor_.resize(node_number+1,nullptr);
      } else {
         if(unaryFactor_[node_number] != nullptr) { throw std::runtime_error("unary factor " + std::to_string(node_number) + " already present"); }
      }
      unaryFactor_[node_number] = u;
      if(node_number > 0 && unaryFactor_[node_number-1]) { // fails for non-contiguous access
         lp_->AddFactorRelation(unaryFactor_[node_number-1], unaryFactor_[node_number]);
      }
      
      lp_->AddFactor(u);;

      return u;
   }
   // unary factor was created elsewhere, let mrf know it
   void RegisterUnaryFactor(const INDEX node_number, UnaryFactorContainer* u)
   {
      if(node_number >= unaryFactor_.size()) {
         unaryFactor_.resize(node_number+1,nullptr);
      } else {
         if(unaryFactor_[node_number] != nullptr) { throw std::runtime_error("unary factor " + std::to_string(node_number) + " already present"); }
      }
      unaryFactor_[node_number] = u;
   }
   template<typename COST>
   PairwiseFactorContainer* AddPairwiseFactor(INDEX var1, INDEX var2, const COST& cost)
   { 
      //if(var1 > var2) std::swap(var1,var2);
      assert(var1<var2);
      assert(!HasPairwiseFactor(var1,var2));
      //assert(cost.size() == GetNumberOfLabels(var1) * GetNumberOfLabels(var2));
      //assert(pairwiseMap_.find(std::make_tuple(var1,var2)) == pairwiseMap_.end());
      //PairwiseFactorContainer* p = new PairwiseFactorContainer(PairwiseFactor(cost), cost);
      auto* p = new PairwiseFactorContainer(GetNumberOfLabels(var1), GetNumberOfLabels(var2), cost);
      lp_->AddFactor(p);
      ConstructPairwiseFactor(*(p->GetFactor()), var1, var2);
      pairwiseFactor_.push_back(p);
      pairwiseIndices_.push_back(std::array<INDEX,2>({var1,var2}));
      const INDEX factorId = pairwiseFactor_.size()-1;
      pairwiseMap_.insert(std::make_pair(std::make_tuple(var1,var2), factorId));
      LinkUnaryPairwiseFactor(unaryFactor_[var1], p, unaryFactor_[var2]);

      lp_->AddFactorRelation(unaryFactor_[var1], p);
      lp_->AddFactorRelation(p, unaryFactor_[var2]);

      return p;
   }

   void LinkUnaryPairwiseFactor(UnaryFactorContainer* const left, PairwiseFactorContainer* const p, UnaryFactorContainer* right)
   {
      auto* l = new LeftMessageContainer(ConstructLeftUnaryPairwiseMessage(left, p), left, p);
      lp_->AddMessage(l);
      auto* r = new RightMessageContainer(ConstructRightUnaryPairwiseMessage(right, p), right, p);
      lp_->AddMessage(r);
   }


   UnaryFactorContainer* GetUnaryFactor(const INDEX i) const { assert(i<unaryFactor_.size()); return unaryFactor_[i]; }
   PairwiseFactorContainer* GetPairwiseFactor(const INDEX i) const { assert(i<pairwiseFactor_.size()); return pairwiseFactor_[i]; }
   PairwiseFactorContainer* GetPairwiseFactor(const INDEX i, const INDEX j) const { 
      assert(i<j);    
      assert(j<unaryFactor_.size());
      const INDEX factor_id = GetPairwiseFactorId(i,j);
      return pairwiseFactor_[factor_id]; 
   }

   LeftMessageContainer* get_left_message(const INDEX i, const INDEX j) const
   {
      assert(i < j);
      auto* f = GetPairwiseFactor(i,j);
      auto msgs = f->template get_messages<LeftMessageContainer>();
      assert(msgs.size() == 1);
      return msgs[0];
   }

   RightMessageContainer* get_right_message(const INDEX i, const INDEX j) const
   {
      assert(i < j);
      auto* f = GetPairwiseFactor(i,j);
      auto msgs = f->template get_messages<RightMessageContainer>();
      assert(msgs.size() == 1);
      return msgs[0];
   }

   INDEX GetNumberOfVariables() const 
   { 
      assert(!(unaryFactor_.size() == 0 && pairwiseFactor_.size() > 0)); 
      return unaryFactor_.size(); 
   } // note: this is not a good idea, if unaryFactors are populated elsewhere: take maximum in pairwise factor indices then.
   INDEX GetPairwiseFactorId(const INDEX var1, const INDEX var2) const 
   {
      assert(var1<var2);
      assert(pairwiseMap_.find(std::make_tuple(var1,var2)) != pairwiseMap_.end());
      return pairwiseMap_.find(std::make_tuple(var1,var2))->second; 
   }
   bool HasPairwiseFactor(const INDEX var1, const INDEX var2) const
   {
      assert(var1<var2);
      if(pairwiseMap_.find(std::make_tuple(var1,var2)) != pairwiseMap_.end()) { 
         return true;
      } else {
         return false;
      }
   }
   INDEX GetNumberOfPairwiseFactors() const { return pairwiseFactor_.size(); }
   std::array<INDEX,2> GetPairwiseVariables(const INDEX factorNo) const { return pairwiseIndices_[factorNo]; }
   INDEX GetNumberOfLabels(const INDEX i) const { return unaryFactor_[i]->GetFactor()->size(); }
   REAL GetPairwiseValue(const INDEX factorId, const INDEX i1, const INDEX i2) const
   {
      assert(i1 < GetNumberOfLabels( GetPairwiseVariables(factorId)[0] ));
      assert(i2 < GetNumberOfLabels( GetPairwiseVariables(factorId)[1] ));
      return (*pairwiseFactor_[factorId]->GetFactor())(i1,i2);
   }


   template<typename SOLVER>
   void Construct(SOLVER& pd) 
   {
      std::cout << "Construct MRF problem with " << unaryFactor_.size() << " unary factors and " << pairwiseFactor_.size() << " pairwise factors\n";

      // add order relations. These are important for the anisotropic weights computation to work.
      if(unaryFactor_.size() > 1) {
         for(auto it=unaryFactor_.begin(); it+1 != unaryFactor_.end(); ++it) {
            assert(*it != nullptr);
            lp_->AddFactorRelation(*it,*(it+1));
         }
      }

      assert(pairwiseIndices_.size() == pairwiseFactor_.size());
   }

   template<typename STREAM>
   void WritePrimal(STREAM& s) const 
   {
      if(unaryFactor_.size() > 0) {
         for(INDEX i=0; i<unaryFactor_.size()-1; ++i) {
            s << unaryFactor_[i]->GetFactor()->primal() << ", ";
         }
         auto* f = unaryFactor_.back();
         s << f->GetFactor()->primal() << "\n";
      }
   }

  // build tree of unary and pairwise factors
  LP_tree add_tree(std::vector<PairwiseFactorContainer*> p)
  {
     LP_tree t;
     assert(p.size() > 0);

     // extract messages joining unaries and pairwise
     std::set<UnaryFactorContainer*> u_set; // u_set is not strictly needed, but helps in checking correctness of algorithm
     std::set<PairwiseFactorContainer*> p_set;
     std::vector<std::tuple<Chirality, MessageTypeAdapter*>> ordered_msgs;
     UnaryFactorContainer* root = nullptr;
     for(auto* f : p) {
        p_set.insert(f); 

        auto left_msgs = f->template get_messages<LeftMessageContainer>();
        assert(left_msgs.size() == 1);
        UnaryFactorContainer* left = (*left_msgs.begin())->GetLeftFactor();
        u_set.insert(left);

        auto right_msgs = f->template get_messages<RightMessageContainer>();
        assert(right_msgs.size() == 1);
        UnaryFactorContainer* right = (*right_msgs.begin())->GetLeftFactor(); 
        u_set.insert(right);

        assert(left != right);

        root = left;
     }

     // assume some arbitrary root. build tree recursively from root
     assert(u_set.size() == p_set.size() + 1);

     std::deque<UnaryFactorContainer*> u_stack;
     u_stack.push_back(root);

     while(!u_stack.empty()) {
        auto* f = u_stack.front();
        u_stack.pop_front();
        u_set.erase(f);

        {
           auto msgs = f->template get_messages<LeftMessageContainer>();
           for(auto it=msgs.begin(); it!= msgs.end(); ++it) {
              auto* p_cand = (*it)->GetRightFactor();
              if(p_set.find(p_cand) != p_set.end()) {
                 p_set.erase(p_cand);
                 //t.AddMessage((*it), Chirality::left); // or right?
                 ordered_msgs.push_back(std::make_tuple(Chirality::left, *it));
                 // search for the other unary connected to p_cand
                 auto msgs_other = p_cand->template get_messages<RightMessageContainer>();
                 assert(msgs_other.size() == 1);
                 auto* u_other = (*msgs_other.begin())->GetLeftFactor();
                 assert(u_set.find(u_other) != u_set.end());
                 //t.AddMessage(*(msgs_other.begin()), Chirality::right);
                 ordered_msgs.push_back(std::make_tuple(Chirality::right, *(msgs_other.begin())));
                 u_stack.push_back(u_other);
              }
           }
        }

        {
           auto msgs = f->template get_messages<RightMessageContainer>();
           for(auto it=msgs.begin(); it!= msgs.end(); ++it) {
              auto* p_cand = (*it)->GetRightFactor();
              if(p_set.find(p_cand) != p_set.end()) {
                 p_set.erase(p_cand);
                 //t.AddMessage((*it), Chirality::left); // or right?
                 ordered_msgs.push_back(std::make_tuple(Chirality::left, *it));
                 // search for the other unary connected to p_cand
                 auto msgs_other = p_cand->template get_messages<LeftMessageContainer>();
                 assert(msgs_other.size() == 1);
                 auto* u_other = (*msgs_other.begin())->GetLeftFactor();
                 assert(u_set.find(u_other) != u_set.end());
                 //t.AddMessage(*(msgs_other.begin()), Chirality::right);
                 ordered_msgs.push_back(std::make_tuple(Chirality::right, *(msgs_other.begin())));
                 u_stack.push_back(u_other);
              }
           }
        }
     }

     assert(p_set.empty());

     std::reverse(ordered_msgs.begin(), ordered_msgs.end());
     for(auto e : ordered_msgs) {
        t.AddMessage(std::get<1>(e), std::get<0>(e));
     }

     t.init();

     t.compute_subgradient();
     assert(std::abs(t.lower_bound() - t.primal_cost()) <= eps);

     return std::move(t);
  }

  // compute forest cover of MRF and add each resulting tree

  auto compute_forest_cover()
  {
     return compute_forest_cover(pairwiseIndices_);
  }
  std::vector<LP_tree> compute_forest_cover(const std::vector<std::array<INDEX,2>>& pairwiseIndices)
  {
     UndirectedGraph g = UndirectedGraph(unaryFactor_.size(), pairwiseFactor_.size());
     for(auto e : pairwiseIndices) {
        const auto i = e[0];
        const auto j = e[1];
        g.AddEdge(i,j,1);
     }
     
     const INDEX forest_num = g.Solve();
     std::cout << "decomposed mrf into " << forest_num << " trees\n";
     std::vector<LP_tree> trees;

     std::vector<int> parents(unaryFactor_.size());
     for(INDEX k=0; k<forest_num; ++k) {
        g.GetForestParents(k, parents.data()); // possible GetForestEdges will return pairwise ids, hence GetPairwiseFactor(id) then can be called, which is faster

        UnionFind uf(unaryFactor_.size());
        for(INDEX i=0; i<parents.size(); ++i) {
           if(parents[i] != -1) {
              uf.merge(i, parents[i]);
           }
        }
        auto contiguous_ids = uf.get_contiguous_ids();
        
        const INDEX no_trees = uf.count();
        std::vector<std::vector<PairwiseFactorContainer*>> pairwise(no_trees);

        for(INDEX i=0; i<parents.size(); ++i) {
           if(parents[i] != -1) {
              const INDEX tree_id = contiguous_ids[uf.find(i)];
              const INDEX j = parents[i];
              assert(tree_id == contiguous_ids[uf.find(j)]);
              PairwiseFactorContainer* p = GetPairwiseFactor(std::min(i,j), std::max(i,j));
              pairwise[tree_id].push_back(p);
           } 
        }
        for(INDEX t=0; t<no_trees; ++t) {
           if(pairwise[t].size() > 0) {
              trees.push_back(add_tree(pairwise[t])); 
           }
        }
     }

     auto check_pairwise_factors_present = [&trees]() -> INDEX {
        INDEX pairwise_factors_in_trees = 0;
        for(auto& tree : trees) {
           for(auto* f : tree.factors_) {
              if(dynamic_cast<PairwiseFactorContainer*>(f)) {
                 pairwise_factors_in_trees++;
              }
           }
        }
        return pairwise_factors_in_trees;
     };
     assert(check_pairwise_factors_present() == pairwiseIndices.size());

     return std::move(trees);
  }

protected:
   std::vector<UnaryFactorContainer*> unaryFactor_;
   std::vector<PairwiseFactorContainer*> pairwiseFactor_;
   
   std::vector<std::array<INDEX,2>> pairwiseIndices_;

   std::map<std::tuple<INDEX,INDEX>, INDEX> pairwiseMap_; // given two sorted indices, return factorId belonging to that index.

   //INDEX unaryFactorIndexBegin_, unaryFactorIndexEnd_; // do zrobienia: not needed anymore

   LP* lp_;
};

// overloads virtual functions above for standard SimplexFactor and SimplexMarginalizationMessage
template<class FACTOR_MESSAGE_CONNECTION, INDEX UNARY_FACTOR_NO, INDEX PAIRWISE_FACTOR_NO, INDEX LEFT_MESSAGE_NO, INDEX RIGHT_MESSAGE_NO>
class StandardMrfConstructor : public MRFProblemConstructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR_NO, PAIRWISE_FACTOR_NO, LEFT_MESSAGE_NO, RIGHT_MESSAGE_NO>
{
public:
   using MRFProblemConstructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR_NO, PAIRWISE_FACTOR_NO, LEFT_MESSAGE_NO, RIGHT_MESSAGE_NO>::MRFProblemConstructor;
   // this is not nice: how to import all aliases directly here? Is it possible at all, while base class is parametrized by templates?
   using BaseConstructor = MRFProblemConstructor<FACTOR_MESSAGE_CONNECTION, UNARY_FACTOR_NO, PAIRWISE_FACTOR_NO, LEFT_MESSAGE_NO, RIGHT_MESSAGE_NO>;
   using UnaryFactorType = typename BaseConstructor::UnaryFactorType;
   using UnaryFactorContainer = typename BaseConstructor::UnaryFactorContainer;
   using PairwiseFactorType = typename BaseConstructor::PairwiseFactorType;
   using PairwiseFactorContainer = typename BaseConstructor::PairwiseFactorContainer;
   using RightMessageType = typename BaseConstructor::RightMessageType;
   using LeftMessageType = typename BaseConstructor::LeftMessageType;

   void ConstructUnaryFactor(UnaryFactorType& u, const std::vector<REAL>& cost) 
   { 
      for(INDEX i=0; i<cost.size(); ++i) {
         u[i] = cost[i];
      }
   }

   virtual void ConstructPairwiseFactor(PairwiseFactorType& p, const INDEX i1, const INDEX i2) {} // nothing needs to be done

   RightMessageType ConstructRightUnaryPairwiseMessage(UnaryFactorContainer* const right, PairwiseFactorContainer* const p)
   { 
      //using RightUnaryLoopType = typename RightMessageType::LeftLoopType;
      //using RightPairwiseLoopType = typename RightMessageType::RightLoopType;

      const INDEX rightDim = right->GetFactor()->size();
      const INDEX leftDim = p->GetFactor()->dim1();

      //RightUnaryLoopType rightUnaryLoop(rightDim);
      //std::array<INDEX,2> pairwiseDim = {{leftDim, rightDim}};
      //RightPairwiseLoopType rightPairwiseLoop( pairwiseDim );
      return RightMessageType(leftDim, rightDim);
   }

   LeftMessageType ConstructLeftUnaryPairwiseMessage(UnaryFactorContainer* const left, PairwiseFactorContainer* const p)
   {
      //using LeftUnaryLoopType = typename LeftMessageType::LeftLoopType;
      //using LeftPairwiseLoopType = typename LeftMessageType::RightLoopType;

      const INDEX leftDim = left->GetFactor()->size();
      const INDEX rightDim = p->GetFactor()->dim2();

      //LeftUnaryLoopType leftUnaryLoop(leftDim);
      //std::array<INDEX,2> pairwiseDim = {{leftDim, rightDim}};
      //LeftPairwiseLoopType leftPairwiseLoop( pairwiseDim );
      //return LeftMessageType(leftUnaryLoop, leftPairwiseLoop);
      return LeftMessageType(leftDim, rightDim);
   }
};

// assumes underlying label space comes from an assignment problem on a bipartite graph, where certain nodes corresponds to unaries in the graphical model. Should inherit from StandardMrfConstructor
template<class MRF_PROBLEM_CONSTRUCTOR>
class AssignmentGmConstructor : public MRF_PROBLEM_CONSTRUCTOR
{
public:
   using PairwiseFactorType = typename MRF_PROBLEM_CONSTRUCTOR::PairwiseFactorType;

   template<typename SOLVER>
   AssignmentGmConstructor(SOLVER& pd) : MRF_PROBLEM_CONSTRUCTOR(pd) {}

   void SetGraph(const std::vector<std::vector<INDEX>> graph) { graph_ = graph; }

   virtual void ConstructPairwiseFactor(PairwiseFactorType& p, const INDEX i1, const INDEX i2) 
   { 
      assert(i1 < graph_.size()+1 && i2 < graph_.size()+1);
      assert(i1 < i2);

      const INDEX dim1 = this->GetNumberOfLabels(i1);
      const INDEX dim2 = this->GetNumberOfLabels(i2);
      assert(dim1 == p.dim1() && dim2 == p.dim2());

      if(i1 < graph_.size() && i2 < graph_.size()) {
         // put infinities on diagonal
         assert(this->GetNumberOfLabels(i1) == graph_[i1].size() + 1);
         assert(this->GetNumberOfLabels(i2) == graph_[i2].size() + 1);
         for(INDEX x1=0; x1<this->GetNumberOfLabels(i1)-1; ++x1) { // last label is non-assignment
            for(INDEX x2=0; x2<this->GetNumberOfLabels(i2)-1; ++x2) { // last label is non-assignment
               if(graph_[i1][x1] == graph_[i2][x2]) {
                  p.cost(x1,x2) = std::numeric_limits<REAL>::infinity();
               }
            }
         }
      }
   }

   bool CheckPrimalConsistency() const // this function should not be needed, as gm model and mp model for graph matching have equality messages and mcf and hungarian have mcf factor which takes care of primal feasibility. Only enable in debug build
   {
      std::cout << "check assignment\n";
      INDEX no_labels = 0;
      for(auto&v : graph_) {
         for(auto l : v) {
            no_labels = std::max(no_labels, l);
         }
      }
      ++no_labels;

      std::vector<bool> labels_taken(no_labels,false);
      for(INDEX i=0; i<this->unaryFactor_.size(); ++i) {
         const INDEX state = this->unaryFactor_[i]->GetFactor()->primal(); 
         if(state < graph_[i].size()) {
            const INDEX label = graph_[i][state];
            if(labels_taken[label]) { 
               std::cout << "var " << i << ", state " << state << ", label " << label << " conflict\n";
               return false; 
            }
            labels_taken[ label ] = true;
         }
      }

      return true;
   }

   const std::vector<INDEX>& assignment_nodes(const INDEX i) const
   {
      assert(i<graph_.size());
      return graph_[i];
   }

private:
   std::vector<std::vector<INDEX>> graph_;

};

// derives from a given mrf problem constructor and adds tightening capabilities on top of it, as implemented in cycle_inequalities and proposed by David Sontag
template<class MRF_PROBLEM_CONSTRUCTOR,
   INDEX TERNARY_FACTOR_NO, INDEX PAIRWISE_TRIPLET_MESSAGE12_NO, INDEX PAIRWISE_TRIPLET_MESSAGE13_NO, INDEX PAIRWISE_TRIPLET_MESSAGE23_NO> // the last indices indicate triplet factor and possible messages
class TighteningMRFProblemConstructor : public MRF_PROBLEM_CONSTRUCTOR
{
protected:
   using MRFPC = MRF_PROBLEM_CONSTRUCTOR;
   using FMC = typename MRFPC::FMC;
   using MrfConstructorType = TighteningMRFProblemConstructor<MRFPC, TERNARY_FACTOR_NO, PAIRWISE_TRIPLET_MESSAGE12_NO, PAIRWISE_TRIPLET_MESSAGE13_NO, PAIRWISE_TRIPLET_MESSAGE23_NO>;

   using TripletFactorContainer = meta::at_c<typename FMC::FactorList, TERNARY_FACTOR_NO>; 
   using TripletFactor = typename TripletFactorContainer::FactorType;
   using PairwiseTripletMessage12Container = typename meta::at_c<typename FMC::MessageList, PAIRWISE_TRIPLET_MESSAGE12_NO>::MessageContainerType;
   using PairwiseTripletMessage13Container = typename meta::at_c<typename FMC::MessageList, PAIRWISE_TRIPLET_MESSAGE13_NO>::MessageContainerType;
   using PairwiseTripletMessage23Container = typename meta::at_c<typename FMC::MessageList, PAIRWISE_TRIPLET_MESSAGE23_NO>::MessageContainerType;

public:
   template<typename SOLVER>
   TighteningMRFProblemConstructor(SOLVER& pd)
      : MRF_PROBLEM_CONSTRUCTOR(pd)
   {}

   TripletFactorContainer* AddTripletFactor(const INDEX var1, const INDEX var2, const INDEX var3, const std::vector<REAL>& cost)
   {
      assert(var1<var2 && var2<var3);
      assert(var3<this->GetNumberOfVariables());
      assert(tripletMap_.find(std::array<INDEX,3>({var1,var2,var3})) == tripletMap_.end());
      
      assert(this->pairwiseMap_.find(std::make_tuple(var1,var2)) != this->pairwiseMap_.end());
      assert(this->pairwiseMap_.find(std::make_tuple(var1,var3)) != this->pairwiseMap_.end());
      assert(this->pairwiseMap_.find(std::make_tuple(var2,var3)) != this->pairwiseMap_.end());

      const INDEX factor12Id = this->pairwiseMap_.find(std::make_tuple(var1,var2))->second;
      const INDEX factor13Id = this->pairwiseMap_.find(std::make_tuple(var1,var3))->second;
      const INDEX factor23Id = this->pairwiseMap_.find(std::make_tuple(var2,var3))->second;

      TripletFactorContainer* t = new TripletFactorContainer(this->GetNumberOfLabels(var1), this->GetNumberOfLabels(var2), this->GetNumberOfLabels(var3));
      this->lp_->AddFactor(t);
      tripletFactor_.push_back(t);
      tripletIndices_.push_back(std::array<INDEX,3>({var1,var2,var3}));
      const INDEX factorId = tripletFactor_.size()-1;
      tripletMap_.insert(std::make_pair(std::array<INDEX,3>({var1,var2,var3}), factorId));

      LinkPairwiseTripletFactor<PairwiseTripletMessage12Container>(factor12Id,factorId);
      LinkPairwiseTripletFactor<PairwiseTripletMessage13Container>(factor13Id,factorId);
      LinkPairwiseTripletFactor<PairwiseTripletMessage23Container>(factor23Id,factorId); 

      //this->lp_->ForwardPassFactorRelation(this->GetPairwiseFactor(factor12Id),t);
      //this->lp_->ForwardPassFactorRelation(this->GetPairwiseFactor(factor13Id),t);
      //this->lp_->ForwardPassFactorRelation(t,this->GetPairwiseFactor(factor23Id));

      //this->lp_->BackwardPassFactorRelation(this->GetPairwiseFactor(factor23Id),t);
      //this->lp_->BackwardPassFactorRelation(this->GetPairwiseFactor(factor13Id),t);
      //this->lp_->BackwardPassFactorRelation(t,this->GetPairwiseFactor(factor12Id));
      return t;
   }
   template<typename PAIRWISE_TRIPLET_MESSAGE_CONTAINER>
   void LinkPairwiseTripletFactor(const INDEX pairwiseFactorId, const INDEX tripletFactorId)
   {
      using PairwiseTripletMessageType = typename PAIRWISE_TRIPLET_MESSAGE_CONTAINER::MessageType;

      typename MRFPC::PairwiseFactorContainer* const p = this->pairwiseFactor_[pairwiseFactorId];
      const INDEX pairwiseVar1 = this->pairwiseIndices_[pairwiseFactorId][0];
      const INDEX pairwiseVar2 = this->pairwiseIndices_[pairwiseFactorId][1];
      assert(pairwiseVar1 < pairwiseVar2);
      const INDEX pairwiseDim1 = this->GetNumberOfLabels(pairwiseVar1);
      const INDEX pairwiseDim2 = this->GetNumberOfLabels(pairwiseVar2);

      TripletFactorContainer* const t = tripletFactor_[tripletFactorId];
      const INDEX tripletVar1 = tripletIndices_[tripletFactorId][0];
      const INDEX tripletVar2 = tripletIndices_[tripletFactorId][1];
      const INDEX tripletVar3 = tripletIndices_[tripletFactorId][2];
      assert(tripletVar1 < tripletVar2 && tripletVar2 < tripletVar3);
      const INDEX tripletDim1 = this->GetNumberOfLabels(tripletVar1);
      const INDEX tripletDim2 = this->GetNumberOfLabels(tripletVar2);
      const INDEX tripletDim3 = this->GetNumberOfLabels(tripletVar3);
         
      assert(pairwiseDim1*pairwiseDim2 == p->GetFactor()->size());

      using MessageType = typename PAIRWISE_TRIPLET_MESSAGE_CONTAINER::MessageType;
      MessageType m = MessageType(tripletDim1, tripletDim2, tripletDim3);
      PAIRWISE_TRIPLET_MESSAGE_CONTAINER* mc = new PAIRWISE_TRIPLET_MESSAGE_CONTAINER(m, p, t);
      this->lp_->AddMessage(mc);
   }
   INDEX GetNumberOfTripletFactors() const { return tripletFactor_.size(); }

   std::array<INDEX,3> GetTripletIndices(const INDEX factor_id)
   {
      assert(factor_id < GetNumberOfTripletFactors());
      return tripletIndices_[factor_id];
   }

   void AddEmptyPairwiseFactor(const INDEX var1, const INDEX var2)
   {
      assert(this->pairwiseMap_.find(std::make_tuple(var1,var2)) == this->pairwiseMap_.end()); 
      this->AddPairwiseFactor(var1, var2, matrix<REAL>(this->GetNumberOfLabels(var1), this->GetNumberOfLabels(var2), 0));
   }

   // do zrobienia: use references for pi
   bool AddTighteningTriplet(const INDEX var1, const INDEX var2, const INDEX var3)//, const std::vector<INDEX> pi1, const std::vector<INDEX> pi2, const std::vector<INDEX> pi3)
   {
      assert(var1 < var2 && var2 < var3 && var3 < this->GetNumberOfVariables());
      if(tripletMap_.find(std::array<INDEX,3>({var1,var2,var3})) == tripletMap_.end()) {
         // first check whether necessary pairwise factors are present. If not, add them.
         if(this->pairwiseMap_.find(std::make_tuple(var1,var2)) == this->pairwiseMap_.end()) {
            AddEmptyPairwiseFactor(var1,var2);
         }
         if(this->pairwiseMap_.find(std::make_tuple(var1,var3)) == this->pairwiseMap_.end()) {
            AddEmptyPairwiseFactor(var1,var3);
         }
         if(this->pairwiseMap_.find(std::make_tuple(var2,var3)) == this->pairwiseMap_.end()) {
            AddEmptyPairwiseFactor(var2,var3);
         }

         const INDEX tripletSize = this->GetNumberOfLabels(var1) *  this->GetNumberOfLabels(var2) * this->GetNumberOfLabels(var3);
         AddTripletFactor(var1,var2,var3, std::vector<REAL>(tripletSize,0.0));
         return true;
      } else {
         return false;
      }
   }

   INDEX add_triplets(const std::vector<triplet_candidate>& tc, const INDEX max_triplets_to_add = std::numeric_limits<INDEX>::max())
   {
      INDEX no_triplets_added = 0;
      for(auto t : tc) {
         if(AddTighteningTriplet(t.i, t.j, t.k)) {
            ++no_triplets_added;
            if(no_triplets_added >= max_triplets_to_add) {
               break;
            }
         }
      }
      return no_triplets_added; 
   }

   INDEX Tighten(const INDEX noTripletsToAdd)
   {
      assert(noTripletsToAdd > 0);
      if(debug()) {
         std::cout << "Tighten mrf with cycle inequalities, no triplets to add = " << noTripletsToAdd << "\n";
      }

      //auto fp = [this](const INDEX v1, const INDEX v2, const INDEX v3) { return this->AddTighteningTriplet(v1,v2,v3); }; // do zrobienia: do not give this via template, as Cycle already has gm_ object.

      triplet_search<typename std::remove_reference<decltype(*this)>::type> triplets(*this, eps);
      if(debug()) { std::cout << "search for triplets\n"; }
      auto triplet_candidates = triplets.search();
      if(debug()) { std::cout << "done\n"; }
      INDEX no_triplets_added = add_triplets(triplet_candidates, noTripletsToAdd);
      if(diagnostics()) { std::cout << "added " << no_triplets_added << " by triplet search\n"; }

      if(no_triplets_added < 0.2*noTripletsToAdd) {
         if(debug()) {
            std::cout << "----------------- do cycle search with k-projection graph---------------\n";
         }
         k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(*this)>::type, false> cycle_search(*this, eps);
         triplet_candidates = cycle_search.search();
         if(debug()) { std::cout << "... done\n"; }
         const INDEX no_triplet_k_projection_graph = add_triplets(triplet_candidates, noTripletsToAdd-no_triplets_added);
         if(diagnostics()) {std::cout << "added " << no_triplet_k_projection_graph << " by cycle search in k-projection graph\n"; }
         no_triplets_added += no_triplet_k_projection_graph;

         // search in expanded projection graph
         if(no_triplets_added < 0.2*noTripletsToAdd) {
            if(debug()) {
               std::cout << "----------------- do cycle search with expanded projection graph---------------\n";
            }
            k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(*this)>::type, true> cycle_search(*this, eps);
            triplet_candidates = cycle_search.search();
            if(debug()) { std::cout << "... done\n"; }
            const INDEX no_triplet_k_projection_graph = add_triplets(triplet_candidates, noTripletsToAdd-no_triplets_added);
            if(diagnostics()) { std::cout << "added " << no_triplet_k_projection_graph << " by cycle search in expanded projection graph\n";
            }
            no_triplets_added += no_triplet_k_projection_graph;
         }
      }

      if(debug()) {
         std::cout << "no triplets found, now search for triplets with guaranteed increase also smaller than " << eps << "\n";
      }
       
      assert(eps >= std::numeric_limits<REAL>::epsilon());

      if(no_triplets_added == 0) {
         triplet_search<typename std::remove_reference<decltype(*this)>::type> triplets(*this, std::numeric_limits<REAL>::epsilon());
         if(debug()) {
            std::cout << "search for triplets (any will do)\n";
         }
         auto triplet_candidates = triplets.search();
         if(debug()) {
            std::cout << "done\n";
         }
         no_triplets_added = add_triplets(triplet_candidates, noTripletsToAdd);
         if(diagnostics()) {
            std::cout << "added " << no_triplets_added << " by triplet search with no increase\n";
         }
      }

      if(no_triplets_added == 0) {
         if(debug()) {
            std::cout << "----------------- do cycle search with k-projection graph, add all---------------\n";
         }
         k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(*this)>::type, false> cycle_search(*this, std::numeric_limits<REAL>::epsilon());
         triplet_candidates = cycle_search.search();
         if(debug()) { std::cout << "... done\n"; }
         const INDEX no_triplet_k_projection_graph = add_triplets(triplet_candidates, noTripletsToAdd-no_triplets_added);
         if(diagnostics()) {
            std::cout << "added " << no_triplet_k_projection_graph << " by cycle search in k-projection grapa with no increase\n";
         }
         no_triplets_added += no_triplet_k_projection_graph;
      }

      if(no_triplets_added == 0) {
         if(debug()) {
            std::cout << "----------------- do cycle search with expanded projection graph, add all---------------\n";
         }
         k_ary_cycle_inequalities_search<typename std::remove_reference<decltype(*this)>::type, true> cycle_search(*this, std::numeric_limits<REAL>::epsilon());
         triplet_candidates = cycle_search.search();
         if(debug()) { std::cout << "... done\n"; }
         const INDEX no_triplet_k_projection_graph = add_triplets(triplet_candidates, noTripletsToAdd-no_triplets_added);
         if(diagnostics()) {
            std::cout << "added " << no_triplet_k_projection_graph << " by cycle search in expanded projection graph with no increase\n";
         }
         no_triplets_added += no_triplet_k_projection_graph;
      }

      return no_triplets_added;
   }


protected:
   std::vector<TripletFactorContainer*> tripletFactor_;
   std::vector<std::array<INDEX,3>> tripletIndices_;
   std::map<std::array<INDEX,3>, INDEX> tripletMap_; // given two sorted indices, return factorId belonging to that index.
};


///////////////////////////////////////////////////////////////////
//
// input classes. Given an input (filename, string) construct factors 
//
///////////////////////////////////////////////////////////////////

// file format described in http://www.cs.huji.ac.il/project/PASCAL/fileFormat.php
namespace UaiMrfInput {

   struct MrfInput {
      INDEX number_of_variables_;
      INDEX number_of_cliques_;
      std::vector<INDEX> cardinality_;
      std::vector<std::vector<INDEX>> clique_scopes_;
      std::vector<std::vector<REAL>> function_tables_;
   };

   // import basic parsers
   using Parsing::opt_whitespace;
   using Parsing::mand_whitespace;
   using Parsing::opt_invisible;
   using Parsing::mand_invisible;
   using Parsing::positive_integer;
   using Parsing::real_number;


   struct init_line : pegtl::seq< opt_whitespace, pegtl::string<'M','A','R','K','O','V'>, opt_whitespace > {};
   struct number_of_variables : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace > {};
   // vector of integers denoting how many labels each variable has
   struct cardinality : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace > {};
   struct number_of_cliques : pegtl::seq< opt_whitespace, positive_integer, opt_whitespace> {};
   // first is the number of variables in the clique, then the actual variables.
   // the clique_scopes should match number_of_clique_lines, each line consisting of a sequence of integers
   struct new_clique_scope : pegtl::seq< positive_integer > {};
   struct clique_scope : pegtl::seq< positive_integer > {};
   struct clique_scope_line : pegtl::seq< opt_whitespace, new_clique_scope, pegtl::plus< opt_whitespace, clique_scope >, opt_whitespace, pegtl::eol > {};
   struct clique_scopes_end
   {
      //template< pegtl::apply_mode A, pegtl::rewind_mode M, template< typename ... > class Action, template< typename ... > class Control, typename Input >
      template< pegtl::apply_mode A, template< typename ... > class Action, template< typename ... > class Control, typename Input >
         static bool match( Input &, MrfInput& input )
         {
            return input.number_of_cliques_ == input.clique_scopes_.size();
         }
   };
   struct clique_scopes : pegtl::until< clique_scopes_end, clique_scope_line > {};
   // a function table is begun by number of entries and then a list of real numbers. Here we record all the values in the real stack
   // do zrobienia: treat whitespace
   struct new_function_table : pegtl::seq< positive_integer > {};
   struct function_table_entry : pegtl::seq< real_number > {};
   struct function_tables_end
   {
      //template< pegtl::apply_mode A, pegtl::rewind_mode M, template< typename ... > class Action, template< typename ... > class Control, typename Input >
      template< pegtl::apply_mode A, template< typename ... > class Action, template< typename ... > class Control, typename Input >
         static bool match( Input &, MrfInput& input )
         {
            return input.number_of_cliques_ == input.function_tables_.size();
         }
   };
   struct function_table_end
   {
      //template< pegtl::apply_mode A, pegtl::rewind_mode M, template< typename ... > class Action, template< typename ... > class Control, typename Input >
      template< pegtl::apply_mode A, template< typename ... > class Action, template< typename ... > class Control, typename Input >
         static bool match( Input &, MrfInput& input )
         {
            auto& table = input.function_tables_.back();
            if(table.back() + 1 == table.size()) {
               table.resize(table.size()-1); // remove last entry which holds the end size of the table
               return true;
            } else {
               return false;
            }
         }
   };
   struct function_table : pegtl::seq< new_function_table, opt_invisible, pegtl::until< function_table_end, opt_invisible, function_table_entry >, opt_invisible > {};
   struct function_tables : pegtl::seq< opt_invisible, pegtl::until< function_tables_end, function_table >, opt_invisible > {};//,pegtl::seq<pegtl::star<pegtl::sor<mand_whitespace, pegtl::eol>>, real_number, pegtl::plus<pegtl::star<pegtl::sor<mand_whitespace, pegtl::eol>>, real_number>> {};
   //template<> struct control< function_table > : pegtl::change_state_and_action< function_table, ..., object_action > {};


   struct grammar : pegtl::seq<
                    init_line, pegtl::eol,
                    number_of_variables, pegtl::eol,
                    pegtl::plus< cardinality >, pegtl::eol,
                    number_of_cliques, pegtl::eol,
                    clique_scopes,
                    opt_invisible,
                    function_tables
                   > {};


   template< typename Rule >
      struct action
      : pegtl::nothing< Rule > {};

   template<> struct action< number_of_variables > {
      template<typename Input>
      static void apply(const Input& in, MrfInput& input) 
      {
         input.number_of_variables_ = std::stoul(in.string());
      }
   };

   template<> struct action< number_of_cliques > {
      template<typename Input>
      static void apply(const Input & in, MrfInput& input)
      {
         input.number_of_cliques_ = std::stoul(in.string()); 
      }
   };

   template<> struct action< cardinality > {
      template<typename Input>
      static void apply(const Input & in, MrfInput& input)
      {
         input.cardinality_.push_back(std::stoul(in.string()));
      }
   };

   template<> struct action< new_clique_scope > {
      template<typename Input>
      static void apply(const Input &, MrfInput& input)
      {
         input.clique_scopes_.push_back(std::vector<INDEX>(0));
      }
   };
   template<> struct action< clique_scope > {
      template<typename Input>
      static void apply(const Input & in, MrfInput& input)
      {
         input.clique_scopes_.back().push_back(std::stoul(in.string()));
         assert(input.clique_scopes_.back().back() < input.number_of_variables_);
      }
   };
   template<> struct action< new_function_table > {
      template<typename Input>
      static void apply(const Input & in, MrfInput& input)
      {
         const INDEX no_entries = std::stoul(in.string());
         std::vector<REAL> entries;
         entries.reserve(no_entries+1);
         entries.push_back(no_entries);
         input.function_tables_.push_back(std::move(entries));
      }
   };
   template<> struct action< function_table_entry > {
      template<typename Input>
      static void apply(const Input & in, MrfInput& input)
      {
         auto& table = input.function_tables_.back();
         table.push_back(std::stod(in.string()));
         std::swap(table.back(), *(table.rbegin()+1)); // exchange last element, which always denotes the final number of entries in the function table
      }
   };

   template<typename MRF_CONSTRUCTOR>
      void build_mrf(MRF_CONSTRUCTOR& mrf, const MrfInput& input)
      {
         assert(input.number_of_cliques_ == input.clique_scopes_.size());
         assert(input.number_of_cliques_ == input.function_tables_.size());

         // only unary and pairwise potentials supported right now
         for(INDEX i=0; i<input.number_of_cliques_; ++i) {
            assert(input.clique_scopes_[i].size() < 3);
         }

         // first input the unaries, as pairwise potentials need them to be able to link to them
         // add unary factors with cost zero for each variables. There are models where unaries are not explicitly added.
         for(INDEX i=0; i<input.number_of_variables_; ++i) {
            const INDEX noLabels = input.cardinality_[i];
            mrf.AddUnaryFactor(i,std::vector<REAL>(noLabels,0.0));
         }

         for(INDEX i=0; i<input.number_of_cliques_; ++i) {
            if(input.clique_scopes_[i].size() == 1) {
               const INDEX var = input.clique_scopes_[i][0];
               auto* f = mrf.GetUnaryFactor(var);
               assert(input.function_tables_[i].size() == input.cardinality_[var]);
               for(INDEX x=0; x<input.function_tables_[i].size(); ++x) {
                  assert( (*f->GetFactor())[x] == 0.0);
                  (*f->GetFactor())[x] = input.function_tables_[i][x];
               }
            }
         }
         // now the pairwise potentials. 
         for(INDEX i=0; i<input.number_of_cliques_; ++i) {
            if(input.clique_scopes_[i].size() == 2) {
               const INDEX var1 = input.clique_scopes_[i][0];
               const INDEX var2 = input.clique_scopes_[i][1];
               const INDEX dim1 = mrf.GetNumberOfLabels(var1);
               const INDEX dim2 = mrf.GetNumberOfLabels(var2);
               assert(var1<var2 && var2 < input.number_of_variables_);
               assert(input.function_tables_[i].size() == input.cardinality_[var1]*input.cardinality_[var2]);
               matrix<REAL> pairwise_cost(dim1,dim2);
               for(INDEX l1=0; l1<dim1; ++l1) {
                  for(INDEX l2=0; l2<dim2; ++l2) {
                     pairwise_cost(l1,l2) = input.function_tables_[i][l2*dim1 + l1];
                  }
               }
               mrf.AddPairwiseFactor(var1,var2,pairwise_cost); // or do we have to transpose the values?
            }
         }
      }


   template<typename SOLVER, INDEX PROBLEM_CONSTRUCTOR_NO>
   bool ParseString(const std::string& instance, SOLVER& s)
   {
      std::cout << "parsing string\n";
      MrfInput input;
      //bool read_suc = pegtl::parse_string<grammar, action>(instance,"",input);
      bool read_suc = pegtl::parse<grammar, action>(instance,"",input);
      if(read_suc) {
         auto& mrf_constructor = s.template GetProblemConstructor<PROBLEM_CONSTRUCTOR_NO>();
         build_mrf(mrf_constructor, input);
      }
      return read_suc;
   }

   template<typename SOLVER>
   bool ParseProblem(const std::string& filename, SOLVER& s)
   {
      std::cout << "parsing " << filename << "\n";
      pegtl::file_parser problem(filename);
      MrfInput input;
      bool read_suc = problem.parse< grammar, action >(input);
      if(read_suc) {
         auto& mrf_constructor = s.template GetProblemConstructor<0>();
         build_mrf(mrf_constructor, input);
      }
      return read_suc;
   }
}

// for graphical models in opengm's hdf5 format and with explicit function tables, function-id-16000
namespace HDF5Input {
   template<typename SOLVER>
   bool ParseProblem(const std::string filename, SOLVER& s)
   {
      auto& mrf = s.template GetProblemConstructor<0>();
      return ParseGM(filename, mrf);
   }
}

} // end namespace LP_MP

#endif // LP_MP_MRF_PROBLEM_CONSTRUCTION_HXX

