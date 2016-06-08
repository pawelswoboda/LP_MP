#ifndef LP_MP_xxx
#define LP_MP_xxx

namespace LP_MP {

  template<class FMC, INDEX MRF_PROBLEM_CONSTRUCTOR_NO, INDEX DISCRETE_TOMOGRAPHY_COUNTING_FACTOR_NO,
	   INDEX DISCRETE_TOMOGRAPHY_MESSAGE_LEFT, INDEX DISCRETE_TOMOGRAPHY_MESSAGE_RIGHT, INDEX DISCRETE_TOMOGRAPHY_COUNTING_PAIRWISE_MESSAGE_NO>
  class DiscreteTomographyTreeConstructor {
  public:
    using MrfConstructorType = typename meta::at_c<FMC::ProblemDecompositionList,MRF_PROBLEM_CONSTRUCTOR_NO>;
    using DiscreteTomographyCountingFactorContainer = typename meta::at_c<FMC::FactorList,DISCRETE_TOMOGRAPHY_COUNTING_FACTOR_NO>;
    using DiscreteTomographyCountingMessageLeft = typename meta::at_c<typename FMC::MessageList, DISCRETE_TOMOGRAPHY_MESSAGE_LEFT>::MessageContainerType;
    using DiscreteTomographyCountingMessageRight = typename meta::at_c<typename FMC::MessageList, DISCRETE_TOMOGRAPHY_MESSAGE_RIGHT>::MessageContainerType;
    using DiscreteTomographyCountingPairwiseMessageContainer = typename meta::at_c<typename FMC::MessageList, DISCRETE_TOMOGRAPHY_COUNTGIN_PAIRWISE_MESSAGE_NO>::MessageContainerType;

    struct Tree {
       std::vector<INDEX> projectionVar_;
    };

    void AddProjection(const std::vector<INDEX>& projectionVar, const std::vector<REAL>& summationCost)
    {
       auto& mrfConstructor = pd_.GetProblemConstructor<MRF_PROBLEM_CONSTRUCTOR_NO>();

       for(INDEX i=0;i<projectionVar.size()-1;++i) {
          const INDEX i1 = std::min(projectionVar[i],projectionVar[i+1]);
          const INDEX i2 = std::max(projectionVar[i],projectionVar[i+1]);

          if(!mrfConstructor.HasPairwiseFactor(i1,i2)) {
             mrfConstructor.AddPairwiseFactor(i1,i2,std::vector<REAL>(noLabels_*noLabels,0.0));
          }
       }

       //auto *f = mrfConstructor.GetPairwiseFactor(i1,i2);

       struct treeIdx{
          treeIdx(INDEX i,INDEX j,INDEX no,INDEX ido) : a(i),b(j),n(no),id(ido){}
          const INDEX a,b,n,id;	
       };

       auto stack = std::vector<treeIdx>();
       auto queue = std::vector<treeIdx>();
       auto factors = std::vector<DiscreteCountingFactorContainer*>();

       for( INDEX i=0;i<projectionVar.size();++i){
          queue.push_back(treeIdx(projectionVar[i],projectionVar[i],1));
       };

       std::vector<treeIdx> *PointerToStack;
       std::vector<treeIdx> *PointerToQueue;
       PointerToStack = &stack;
       PointerToQueue = &queue;

       while( PointerToQueue->size() > 0 ){

          if( PointerToQueue->size() == 1 ){
             PointerToStack->push_back(PointerToQueue->back());
             PointerToQueue->pop_back();
          }
          treeIdx idxL = (*PointerToQueue)[0];
          treeIdx idxR = (*PointerToQueue)[1];
          PointerToStack->erase( PointerToStack->begin(),
                PointerToStack->begin() + 2 );
          PointerToStack->push_back(treeIdx(idxL.b,
                   idxR.a,
                   idxL.n+idxR.n,
                   factors.size()
                   )
                );
          std::vector<REAL> repam;

          auto *f = new DiscreteTomographyCountingFactorContainer(DiscreteTomographyFactorCounting(noLabels_,idxL.n,idxR.n),repam);
          pd_.GetLP()->AddFactor(f);

          auto *reg = mrfConstructor.GetPairwiseFactor(idxL.b,idxR.a);
          auto *m_pairwise = new DiscreteTomographyCountingPairwiseMessageContainer(DiscreteTomographyMessageCountingPairwise(noLabels_),
                reg,factors.back(),pow(noLabels_,2));
          pd_.GetMessage()->AddMessage(m_pairwise);






          std::swap(PointerToStack,PointerToQueue);	
       }

       //auto *f = new DiscreteTomographyCountingFactorContainer(DiscreteTomographyCountingFactor(...),std::vector<REAL>(,0.0)); // <-- Woher kennt er (i1.i2)
       //pd_.GetLP()->AddFactor(f);


       auto *m_left = new DiscreteTomographyCountingMessageLeftContainer(DiscreteTomographyMessage<Direction::Left>(), leftFactorContPtr, rightFactorContPtr,...);
       pd_.GetMessage()->AddMessage(m_left);
       auto *m_right = new DiscreteTomographyCountingMessageLeftContainer(DiscreteTomographyMessage<Direction::Right>(),...);
       pd_.GetMessage()->AddMessage(m_right);
       auto *m_pairwise = new DiscreteTomographyCountingMessageLeftContainer(DiscreteTomographyMessage,...);
       pd_.GetMessage()->AddMessage(m_pairwise);


    }

  private:
    std::vector<Tree> treeIndices_;
    const INDEX noLabels_;
    ... pd_;
  }

}

#endif // LP_MP_xxx
