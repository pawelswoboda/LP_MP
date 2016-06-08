#ifndef LP_MP_DT_TREE_CONSTRUCTOR
#define LP_MP_DT_TREE_CONSTRUCTOR

namespace LP_MP {

  template<class FMC,
	   INDEX MRF_PROBLEM_CONSTRUCTOR_NO,
	   INDEX DISCRETE_TOMOGRAPHY_COUNTING_FACTOR_NO,
	   INDEX DISCRETE_TOMOGRAPHY_MESSAGE_LEFT,
	   INDEX DISCRETE_TOMOGRAPHY_MESSAGE_RIGHT,
	   INDEX DISCRETE_TOMOGRAPHY_COUNTING_PAIRWISE_MESSAGE_NO>
  class DiscreteTomographyTreeConstructor {
  public:
    using MrfConstructorType =
      typename meta::at_c<typename FMC::ProblemDecompositionList,MRF_PROBLEM_CONSTRUCTOR_NO>;
    using DiscreteTomographyCountingFactorContainer =
      typename meta::at_c<typename FMC::FactorList,DISCRETE_TOMOGRAPHY_COUNTING_FACTOR_NO>;
    using DiscreteTomographyCountingMessageLeft =
      typename meta::at_c<typename FMC::MessageList, DISCRETE_TOMOGRAPHY_MESSAGE_LEFT>::MessageContainerType;
    using DiscreteTomographyCountingMessageRight =
      typename meta::at_c<typename FMC::MessageList, DISCRETE_TOMOGRAPHY_MESSAGE_RIGHT>::MessageContainerType;
    using DiscreteTomographyCountingPairwiseMessageContainer =
      typename meta::at_c<typename FMC::MessageList, DISCRETE_TOMOGRAPHY_COUNTING_PAIRWISE_MESSAGE_NO>::MessageContainerType;

    DiscreteTomographyTreeConstructor(ProblemDecomposition<FMC>& pd) : pd_(pd) {}
    
    void SetNumberOfLabels(const INDEX noLabels) { noLabels_ = noLabels; }

    void AddProjection(const std::vector<INDEX>& projectionVar, const std::vector<REAL>& summationCost)
    {      
      auto& mrfConstructor = pd_.template GetProblemConstructor<MRF_PROBLEM_CONSTRUCTOR_NO>();

      for(INDEX i=0;i<projectionVar.size()-1;++i) {
	const INDEX i1 = std::min(projectionVar[i],projectionVar[i+1]);
	const INDEX i2 = std::max(projectionVar[i],projectionVar[i+1]);

	if(!mrfConstructor.HasPairwiseFactor(i1,i2)) {
	  mrfConstructor.AddPairwiseFactor(i1,i2,std::vector<REAL>(pow(noLabels_,2),0.0));
	}
      }
      
      struct treeIdx{
	treeIdx(INDEX i,INDEX j,INDEX no,INDEX ido) : a(i),b(j),n(no),id(ido){}
	const INDEX a,b,n,id;	
      };

      auto stack = std::deque<treeIdx>();
      auto queue = std::deque<treeIdx>();
      auto factors = std::vector<DiscreteTomographyCountingFactorContainer*>();
    
      for( INDEX i=0;i<projectionVar.size();++i){
	queue.push_back(treeIdx(projectionVar[i],projectionVar[i],1,0));
      };
    
      std::deque<treeIdx> *PointerToStack;
      std::deque<treeIdx> *PointerToQueue;
      PointerToStack = &stack;
      PointerToQueue = &queue;

      while( PointerToQueue->size() > 0 ){

	if( PointerToQueue->size() == 1 ){
	  PointerToStack->push_back(PointerToQueue->back());
	  PointerToQueue->pop_back();
	}
	else {
	  treeIdx idxL = (*PointerToQueue)[0];
	  treeIdx idxR = (*PointerToQueue)[1];
	  PointerToStack->pop_front();
	  PointerToStack->pop_front();
	
	  PointerToStack->push_back(treeIdx(idxL.a,
					    idxR.b,
					    idxL.n+idxR.n,
					    factors.size()
					    )
				    );

	
	  // factors
	  DiscreteTomographyFactorCounting t(noLabels_,idxL.n,idxR.n,summationCost.size());
	  INDEX repamSize =
	    t.getSize(DiscreteTomographyFactorCounting::NODE::up) +
	    t.getSize(DiscreteTomographyFactorCounting::NODE::left) +
	    t.getSize(DiscreteTomographyFactorCounting::NODE::right) +
	    t.getSize(DiscreteTomographyFactorCounting::NODE::reg);
	  std::vector<REAL> repam(repamSize,0.0);
	   
	  if( idxL.n + idxR.n != projectionVar.size() ){
	    assert(PointerToQueue->size() > 0);
	  }
	  else{
	    assert(PointerToQueue->size() == 0);
	    
	    for(INDEX i=0;i<summationCost.size();i++){
	      for(INDEX j=0;j<pow(noLabels_,2);j++){
		repam[j+i*pow(noLabels_,2)]=summationCost[i];
	      }
	    }
	  }

	  if( idxL.n == 1 ){
	    for(INDEX i=0;i<t.getSize(DiscreteTomographyFactorCounting::NODE::left);i++){
	      if( (i % (1 + noLabels_ + pow(noLabels_,2))) != 0 ){
		repam[t.getSize(DiscreteTomographyFactorCounting::NODE::up) + i] = std::numeric_limits<REAL>::max();
	      }
	    }
	  }
	  if( idxR.n == 1 ){
	    for(INDEX i=0;i<t.getSize(DiscreteTomographyFactorCounting::NODE::right);i++){
	      if( (i % (1 + noLabels_ + pow(noLabels_,2))) != 0 ){
		repam[t.getSize(DiscreteTomographyFactorCounting::NODE::up) +
		  t.getSize(DiscreteTomographyFactorCounting::NODE::left) + i] = std::numeric_limits<REAL>::max();
	      }
	    }
	  }
		  
	  auto *f = new DiscreteTomographyCountingFactorContainer(t,repam);
	  pd_.GetLP()->AddFactor(f);
	  factors.push_back(f);
	  	  
	  auto *reg = mrfConstructor.GetPairwiseFactor(idxL.b,idxR.a);

	  // messages
	  auto *m_pairwise = new DiscreteTomographyCountingPairwiseMessageContainer(DiscreteTomographyMessageCountingPairwise(noLabels_),
										    reg,factors.back(),pow(noLabels_,2));
	  pd_.GetLP()->AddMessage(m_pairwise);

	  if(idxL.n != 1){
	    auto *m_left = new DiscreteTomographyCountingMessageLeft(DiscreteTomographyMessageCounting<DIRECTION::left>(noLabels_,idxL.n),
								     factors[idxL.id],factors.back(),
								     pow(noLabels_,2)*std::min(summationCost.size(),idxL.n*(noLabels_-1)+1));
	    pd_.GetLP()->AddMessage(m_left);
	  }
	  if(idxR.n != 1){
	    auto *m_right = new DiscreteTomographyCountingMessageRight(DiscreteTomographyMessageCounting<DIRECTION::right>(noLabels_,idxR.n),
								       factors[idxR.id],factors.back(),
								       pow(noLabels_,2)*std::min(summationCost.size(),(idxR.n*(noLabels_-1)+1)));
	    pd_.GetLP()->AddMessage(m_right);
	  }
	}
	std::swap(PointerToStack,PointerToQueue);	
      }
	
      
    }

  private:
    //std::vector<Tree> treeIndices_;
    INDEX noLabels_ = 0;
    ProblemDecomposition<FMC>& pd_;
  };

}

#endif // LP_MP_xxx
