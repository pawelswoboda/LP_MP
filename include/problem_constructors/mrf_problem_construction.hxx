#ifndef LP_MP_MRF_PROBLEM_CONSTRUCTION_HXX
#define LP_MP_MRF_PROBLEM_CONSTRUCTION_HXX

#include "solver.hxx"
#include "cycle_inequalities.hxx"

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


   MRFProblemConstructor(Solver<FMC>& solver) : lp_(&solver.GetLP()) {}
   // do zrobienia: this object is not movable
   /*
   MRFProblemConstructor(const MRFProblemConstructor& other)
   {
      assert(false);
   }
   */
   virtual UnaryFactorType ConstructUnaryFactor(const std::vector<REAL>& cost) = 0;
   virtual PairwiseFactorType ConstructPairwiseFactor(const std::vector<REAL>& cost, const INDEX leftDim, const INDEX rightDim) = 0;
   virtual RightMessageType ConstructRightUnaryPairwiseMessage(UnaryFactorContainer* const right, PairwiseFactorContainer* const p) = 0;
   virtual LeftMessageType ConstructLeftUnaryPairwiseMessage(UnaryFactorContainer* const right, PairwiseFactorContainer* const p) = 0;

   /*
   INDEX AddUnaryFactor(const std::vector<REAL>& cost)
   {
      UnaryFactorContainer* u = new UnaryFactorContainer(UnaryFactorType(cost), cost);
      unaryFactor_.push_back(u);
      lp_->AddFactor(u);;
      return unaryFactor_.size()-1;
   }
   */
   UnaryFactorContainer* AddUnaryFactor(const INDEX node_number, const std::vector<REAL>& cost)
   {
      //UnaryFactorContainer* u = new UnaryFactorContainer(UnaryFactorType(cost), cost);
      auto* u = new UnaryFactorContainer( ConstructUnaryFactor(cost), cost);
      if(node_number >= unaryFactor_.size()) {
         unaryFactor_.resize(node_number+1,nullptr);
      } else {
         if(unaryFactor_[node_number] != nullptr) { throw std::runtime_error("unary factor " + std::to_string(node_number) + " already present"); }
      }
      unaryFactor_[node_number] = u;
      
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
   PairwiseFactorContainer* AddPairwiseFactor(INDEX var1, INDEX var2, const std::vector<REAL>& cost)
   { 
      //if(var1 > var2) std::swap(var1,var2);
      assert(var1<var2);
      assert(!HasPairwiseFactor(var1,var2));
      assert(cost.size() == GetNumberOfLabels(var1) * GetNumberOfLabels(var2));
      //assert(pairwiseMap_.find(std::make_tuple(var1,var2)) == pairwiseMap_.end());
      //PairwiseFactorContainer* p = new PairwiseFactorContainer(PairwiseFactor(cost), cost);
      auto* p = new PairwiseFactorContainer(ConstructPairwiseFactor(cost, GetNumberOfLabels(var1), GetNumberOfLabels(var2)), cost);
      pairwiseFactor_.push_back(p);
      pairwiseIndices_.push_back(std::make_tuple(var1,var2));
      const INDEX factorId = pairwiseFactor_.size()-1;
      pairwiseMap_.insert(std::make_pair(std::make_tuple(var1,var2), factorId));
      LinkUnaryPairwiseFactor(unaryFactor_[var1], p, unaryFactor_[var2]);

      lp_->AddFactor(p);
      lp_->AddFactorRelation(unaryFactor_[var1], p);
      lp_->AddFactorRelation(p, unaryFactor_[var2]);

      return p;
   }
   /*
   void LinkLeftUnaryPairwiseFactor(UnaryFactorContainer* const left, PairwiseFactorContainer* const p, LeftMessageType msg)
   {
      //assert(false); // left->size need not be msg size. Use instead message size
      auto* m = new LeftMessageContainer(msg, left, p, left->size());
      leftMessage_.push_back(m);
      lp_->AddMessage(m);
   }
   void LinkRightUnaryPairwiseFactor(UnaryFactorContainer* const right, PairwiseFactorContainer* const p, RightMessageType msg)
   {
      //assert(false); // left->size need not be msg size. Use instead message size
      auto* m = new RightMessageContainer(msg, right, p, right->size());
      rightMessage_.push_back(m);
      lp_->AddMessage(m);
   }
   */
   void LinkUnaryPairwiseFactor(UnaryFactorContainer* const left, PairwiseFactorContainer* const p, UnaryFactorContainer* right)
   {
      auto* l = new LeftMessageContainer(ConstructLeftUnaryPairwiseMessage(left, p), left, p, left->size());
      leftMessage_.push_back(l);
      lp_->AddMessage(l);
      auto* r = new RightMessageContainer(ConstructRightUnaryPairwiseMessage(right, p), right, p, right->size());
      rightMessage_.push_back(r);
      lp_->AddMessage(r);
      /*
      using LeftUnaryLoopType = typename LeftMessageType::LeftLoopType;
      using LeftPairwiseLoopType = typename LeftMessageType::RightLoopType;
      using RightUnaryLoopType = typename RightMessageType::LeftLoopType;
      using RightPairwiseLoopType = typename RightMessageType::RightLoopType;

      const INDEX leftDim = left->size();
      const INDEX rightDim = right->size();
      if(leftDim*rightDim != p->size()) throw std::runtime_error("dimensions for pairwise potential do not match");

      LeftUnaryLoopType leftUnaryLoop(leftDim);
      RightUnaryLoopType rightUnaryLoop(rightDim);
      std::array<INDEX,2> pairwiseDim = {{leftDim, rightDim}};
      LeftPairwiseLoopType leftPairwiseLoop( pairwiseDim );
      RightPairwiseLoopType rightPairwiseLoop( pairwiseDim );

      LinkLeftUnaryPairwiseFactor(left, p, LeftMessageType(leftUnaryLoop, leftPairwiseLoop));
      LinkRightUnaryPairwiseFactor(right, p, RightMessageType(rightUnaryLoop, rightPairwiseLoop));
      */
   }


   UnaryFactorContainer* GetUnaryFactor(const INDEX i) const { assert(i<unaryFactor_.size()); return unaryFactor_[i]; }
   PairwiseFactorContainer* GetPairwiseFactor(const INDEX i) const { assert(i<pairwiseFactor_.size()); return pairwiseFactor_[i]; }

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
   std::tuple<INDEX,INDEX> GetPairwiseVariables(const INDEX factorNo) const { return pairwiseIndices_[factorNo]; }
   INDEX GetNumberOfLabels(const INDEX i) const { return unaryFactor_[i]->size(); }
   REAL GetPairwiseValue(const INDEX factorId, const INDEX i1, const INDEX i2) const
   {
      assert(i1 < GetNumberOfLabels( std::get<0>(GetPairwiseVariables(factorId)) ));
      assert(i2 < GetNumberOfLabels( std::get<1>(GetPairwiseVariables(factorId)) ));
      const INDEX var1 = std::get<0>(GetPairwiseVariables(factorId));
      const INDEX label = i1 + i2*GetNumberOfLabels(var1);
      //const INDEX var2 = std::get<1>(GetPairwiseVariables(factorId));
      //const INDEX label = i2 + i1*GetNumberOfLabels(var2);
      return pairwiseFactor_[factorId]->operator[](label);
   }


   void Construct(Solver<FMC>& pd) 
   {
      std::cout << "Construct MRF problem with " << unaryFactor_.size() << " unary factors and " << pairwiseFactor_.size() << " pairwise factors\n";

      // add order relations. These are important for the anisotropic weights computation to work.
      unaryFactorIndexBegin_ = lp_->GetNumberOfFactors();
      if(unaryFactor_.size() > 1) {
         for(auto it=unaryFactor_.begin(); it+1 != unaryFactor_.end(); ++it) {
            assert(*it != nullptr);
            lp_->AddFactorRelation(*it,*(it+1));
         }
      }

      assert(pairwiseIndices_.size() == pairwiseFactor_.size());
      /*
      auto pairwiseFactorIt = pairwiseFactor_.begin();
      auto pairwiseIndicesIt = pairwiseIndices_.begin();
      for(; pairwiseFactorIt != pairwiseFactor_.end(); ++pairwiseFactorIt, ++pairwiseIndicesIt) {
         const INDEX i = std::get<0>(*pairwiseIndicesIt);
         const INDEX j = std::get<1>(*pairwiseIndicesIt);
         lp->AddFactorRelation(unaryFactor_[i], *pairwiseFactorIt);
         lp->AddFactorRelation(*pairwiseFactorIt, unaryFactor_[j]);
      }
      */
      /*
      for(LeftMessageContainer* messageIt : leftMessage_) {
         lp->AddMessage(messageIt);
      }
      for(RightMessageContainer* messageIt : rightMessage_) {
         lp->AddMessage(messageIt);
      }
      */
   }

protected:
   std::vector<UnaryFactorContainer*> unaryFactor_;
   std::vector<PairwiseFactorContainer*> pairwiseFactor_;
   
   std::vector<std::tuple<INDEX,INDEX>> pairwiseIndices_;

   std::vector<LeftMessageContainer*> leftMessage_;
   std::vector<RightMessageContainer*> rightMessage_;

   std::map<std::tuple<INDEX,INDEX>, INDEX> pairwiseMap_; // given two sorted indices, return factorId belonging to that index.

   INDEX unaryFactorIndexBegin_, unaryFactorIndexEnd_; 

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

   UnaryFactorType ConstructUnaryFactor(const std::vector<REAL>& cost) 
   { return UnaryFactorType(cost); }

   PairwiseFactorType ConstructPairwiseFactor(const std::vector<REAL>& cost, const INDEX leftDim, const INDEX rightDim) 
   { return PairwiseFactorType(cost); }

   RightMessageType ConstructRightUnaryPairwiseMessage(UnaryFactorContainer* const right, PairwiseFactorContainer* const p)
   { 
      using RightUnaryLoopType = typename RightMessageType::LeftLoopType;
      using RightPairwiseLoopType = typename RightMessageType::RightLoopType;

      const INDEX rightDim = right->size();
      const INDEX leftDim = p->size() / rightDim;

      RightUnaryLoopType rightUnaryLoop(rightDim);
      std::array<INDEX,2> pairwiseDim = {{leftDim, rightDim}};
      RightPairwiseLoopType rightPairwiseLoop( pairwiseDim );
      return RightMessageType(rightUnaryLoop, rightPairwiseLoop);
   }

   LeftMessageType ConstructLeftUnaryPairwiseMessage(UnaryFactorContainer* const left, PairwiseFactorContainer* const p)
   {
      using LeftUnaryLoopType = typename LeftMessageType::LeftLoopType;
      using LeftPairwiseLoopType = typename LeftMessageType::RightLoopType;

      const INDEX leftDim = left->size();
      const INDEX rightDim = p->size() / leftDim;

      LeftUnaryLoopType leftUnaryLoop(leftDim);
      std::array<INDEX,2> pairwiseDim = {{leftDim, rightDim}};
      LeftPairwiseLoopType leftPairwiseLoop( pairwiseDim );
      return LeftMessageType(leftUnaryLoop, leftPairwiseLoop);
   }
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
   TighteningMRFProblemConstructor(Solver<FMC>& pd)
      : MRF_PROBLEM_CONSTRUCTOR(pd)
   {}

   TripletFactorContainer* AddTripletFactor(const INDEX var1, const INDEX var2, const INDEX var3, const std::vector<REAL>& cost)
   {
      assert(var1<var2 && var2<var3);
      assert(var3<this->GetNumberOfVariables());
      assert(tripletMap_.find(std::make_tuple(var1,var2,var3)) == tripletMap_.end());
      
      assert(this->pairwiseMap_.find(std::make_tuple(var1,var2)) != this->pairwiseMap_.end());
      assert(this->pairwiseMap_.find(std::make_tuple(var1,var3)) != this->pairwiseMap_.end());
      assert(this->pairwiseMap_.find(std::make_tuple(var2,var3)) != this->pairwiseMap_.end());

      const INDEX factor12Id = this->pairwiseMap_.find(std::make_tuple(var1,var2))->second;
      const INDEX factor13Id = this->pairwiseMap_.find(std::make_tuple(var1,var3))->second;
      const INDEX factor23Id = this->pairwiseMap_.find(std::make_tuple(var2,var3))->second;

      TripletFactorContainer* t = new TripletFactorContainer(TripletFactor(cost), cost);
      tripletFactor_.push_back(t);
      tripletIndices_.push_back(std::make_tuple(var1,var2,var3));
      const INDEX factorId = tripletFactor_.size()-1;
      tripletMap_.insert(std::make_pair(std::make_tuple(var1,var2,var3), factorId));

      LinkPairwiseTripletFactor<PairwiseTripletMessage12Container>(factor12Id,factorId);
      LinkPairwiseTripletFactor<PairwiseTripletMessage13Container>(factor13Id,factorId);
      LinkPairwiseTripletFactor<PairwiseTripletMessage23Container>(factor23Id,factorId);

      this->lp_->AddFactor(t);

      this->lp_->AddFactorRelation(this->GetPairwiseFactor(factor12Id),t);
      this->lp_->AddFactorRelation(this->GetPairwiseFactor(factor13Id),t);
      this->lp_->AddFactorRelation(t,this->GetPairwiseFactor(factor23Id));
      return t;
   }
   template<typename PAIRWISE_TRIPLET_MESSAGE_CONTAINER>
   void LinkPairwiseTripletFactor(const INDEX pairwiseFactorId, const INDEX tripletFactorId)
   {
      using PairwiseTripletMessageType = typename PAIRWISE_TRIPLET_MESSAGE_CONTAINER::MessageType;

      using PairwiseLoopType = typename PairwiseTripletMessageType::LeftLoopType;
      using TripletLoopType = typename PairwiseTripletMessageType::RightLoopType;

      typename MRFPC::PairwiseFactorContainer* const p = this->pairwiseFactor_[pairwiseFactorId];
      const INDEX pairwiseVar1 = std::get<0>(this->pairwiseIndices_[pairwiseFactorId]);
      const INDEX pairwiseVar2 = std::get<1>(this->pairwiseIndices_[pairwiseFactorId]);
      const INDEX pairwiseDim1 = this->GetNumberOfLabels(pairwiseVar1);
      const INDEX pairwiseDim2 = this->GetNumberOfLabels(pairwiseVar2);

      TripletFactorContainer* const t = tripletFactor_[tripletFactorId];
      const INDEX tripletVar1 = std::get<0>(tripletIndices_[tripletFactorId]);
      const INDEX tripletVar2 = std::get<1>(tripletIndices_[tripletFactorId]);
      const INDEX tripletVar3 = std::get<2>(tripletIndices_[tripletFactorId]);
      const INDEX tripletDim1 = this->GetNumberOfLabels(tripletVar1);
      const INDEX tripletDim2 = this->GetNumberOfLabels(tripletVar2);
      const INDEX tripletDim3 = this->GetNumberOfLabels(tripletVar3);
         
      assert(pairwiseDim1*pairwiseDim2 == p->size());
      assert(tripletDim1*tripletDim2*tripletDim3 == t->size());

      PairwiseLoopType pairwiseLoop( pairwiseDim1*pairwiseDim2 );

      std::array<INDEX,3> tripletDim = {{tripletDim1, tripletDim2, tripletDim3}};
      TripletLoopType tripletLoop( tripletDim );

      using MessageType = typename PAIRWISE_TRIPLET_MESSAGE_CONTAINER::MessageType;
      MessageType m = MessageType(pairwiseLoop, tripletLoop);
      PAIRWISE_TRIPLET_MESSAGE_CONTAINER* mc = new PAIRWISE_TRIPLET_MESSAGE_CONTAINER(m, p, t, p->size());
      tripletMessage_.push_back( mc );
      this->lp_->AddMessage(mc);
   }
   INDEX GetNumberOfTripletFactors() const { return tripletFactor_.size(); }

   void AddEmptyPairwiseFactor(const INDEX var1, const INDEX var2)
   {
      assert(this->pairwiseMap_.find(std::make_tuple(var1,var2)) == this->pairwiseMap_.end()); 
      const INDEX dim = this->GetNumberOfLabels(var1) * this->GetNumberOfLabels(var2);
      this->AddPairwiseFactor(var1,var2,std::vector<REAL>(dim,0));
   }

   // do zrobienia: use references for pi
   bool AddTighteningTriplet(const SIGNED_INDEX var1, const SIGNED_INDEX var2, const SIGNED_INDEX var3)//, const std::vector<SIGNED_INDEX> pi1, const std::vector<SIGNED_INDEX> pi2, const std::vector<SIGNED_INDEX> pi3)
   {
      assert(var1 < var2 && var2 < var3 && var3 < this->GetNumberOfVariables());
      if(tripletMap_.find(std::make_tuple(var1,var2,var3)) == tripletMap_.end()) {
         std::cout << "Add tightening triplet, do zrobienia: add infinity on diagonals for matching problem\n";
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
         AddTripletFactor(var1,var2,var3, std::vector<REAL>(tripletSize));
         return true;
      } else {
         return false;
      }
   }

   INDEX Tighten(const REAL epsilon, const INDEX noTripletsToAdd)
   {
      assert(noTripletsToAdd > 0);
      std::cout << "Tighten mrf with cycle inequalities, epsilon = " << epsilon << ", no triplets to add = " << noTripletsToAdd << "\n";
      // do zrobienia: templatize addTriplet function to avoid this declaration
      //std::function<void(const SIGNED_INDEX,const SIGNED_INDEX, const SIGNED_INDEX, const std::vector<SIGNED_INDEX>, const std::vector<SIGNED_INDEX>, const std::vector<SIGNED_INDEX>)> addTriplet = &MrfConstructorType::AddTighteningTriplet;

      std::map<std::vector<int>, bool > tripletSet;
      Cycle<decltype(*this)> cycle(*this);

      std::cout << this->GetNumberOfVariables() << "\n";
      std::cout << this->GetNumberOfPairwiseFactors() << "\n";
      std::cout << this << "\n";

      auto fp = [this](const INDEX v1, const INDEX v2, const INDEX v3) { return this->AddTighteningTriplet(v1,v2,v3); }; // do zrobienia: do not give this via template, as Cycle already has gm_ object.
      INDEX noTripletsAdded = cycle.TightenTriplet(fp, noTripletsToAdd,epsilon,tripletSet);
      std::cout << "Added " << noTripletsAdded << " < " << noTripletsToAdd << " triplets by triplet searching\n";
      //std::cout << "return already\n";
      //return noTripletsAdded;
      if(noTripletsAdded < noTripletsToAdd) {
         std::vector<int> projection_imap;
         std::vector<std::vector<int> > partition_imap;
         std::vector<std::list<int> > cycle_set;
         const INDEX noTripletsAddedByCycle = cycle.TightenCycle(fp, noTripletsToAdd - noTripletsAdded, epsilon, projection_imap, partition_imap, cycle_set, 1);
         std::cout << "Added " << noTripletsAddedByCycle << " triplets by cycle searching (k projection graph)\n";
         noTripletsAdded += noTripletsAddedByCycle;

         if(noTripletsAdded < noTripletsToAdd) {
            const INDEX noTripletsAddedByCycle = cycle.TightenCycle(fp, noTripletsToAdd - noTripletsAdded, epsilon, projection_imap, partition_imap, cycle_set, 2);
            std::cout << "Added " << noTripletsAddedByCycle << " triplets by cycle searching (full projection graph)\n";
            noTripletsAdded += noTripletsAddedByCycle;
         }

         std::cout << this->GetNumberOfVariables() << "\n";
         std::cout << this->GetNumberOfPairwiseFactors() << "\n";
         std::cout << this->GetNumberOfTripletFactors() << "\n";
         std::cout << this << "\n";
      }

      return noTripletsAdded;
   }


protected:
   std::vector<TripletFactorContainer*> tripletFactor_;
   std::vector<std::tuple<INDEX,INDEX,INDEX>> tripletIndices_;
   std::vector<MessageTypeAdapter*> tripletMessage_;
   std::map<std::tuple<INDEX,INDEX,INDEX>, INDEX> tripletMap_; // given two sorted indices, return factorId belonging to that index.
};

} // end namespace LP_MP

#endif // LP_MP_MRF_PROBLEM_CONSTRUCTION_HXX

