#ifndef LP_MP_DT_TREE_CONSTRUCTOR_HXX
#define LP_MP_DT_TREE_CONSTRUCTOR_HXX

namespace LP_MP {

  template<typename FMC,
           INDEX MRF_PROBLEM_CONSTRUCTOR_NO,
           typename COUNTING_FACTOR,
           typename COUNTING_MESSAGE_LEFT,
           typename COUNTING_MESSAGE_RIGHT,
           typename PAIRWISE_COUNTING_CENTER,
           typename PAIRWISE_COUNTING_LEFT,
           typename PAIRWISE_COUNTING_RIGHT,
           typename LEFT_PAIRWISE_COUNTING_LEFT,
           typename RIGHT_PAIRWISE_COUNTING_LEFT,
           typename LEFT_PAIRWISE_COUNTING_RIGHT,
           typename RIGHT_PAIRWISE_COUNTING_RIGHT
          > 
  class DiscreteTomographyTreeConstructor {
  public:
    using MrfConstructorType =
      typename meta::at_c<typename FMC::ProblemDecompositionList,MRF_PROBLEM_CONSTRUCTOR_NO>;
    using PairwiseFactorType = typename MrfConstructorType::PairwiseFactorContainer;
    //using DiscreteTomographyCountingFactorContainer =
    //  typename meta::at_c<typename FMC::FactorList,DISCRETE_TOMOGRAPHY_COUNTING_FACTOR_NO>;
    //using DiscreteTomographyCountingMessageLeft =
    //  typename meta::at_c<typename FMC::MessageList, DISCRETE_TOMOGRAPHY_MESSAGE_LEFT>::MessageContainerType;
    //using DiscreteTomographyCountingMessageRight =
    //  typename meta::at_c<typename FMC::MessageList, DISCRETE_TOMOGRAPHY_MESSAGE_RIGHT>::MessageContainerType;
    //using DiscreteTomographyCountingPairwiseMessageContainer =
    //  typename meta::at_c<typename FMC::MessageList, DISCRETE_TOMOGRAPHY_COUNTING_PAIRWISE_MESSAGE_NO>::MessageContainerType;

    DiscreteTomographyTreeConstructor(Solver<FMC>& pd) : pd_(pd) {}
    
    void SetNumberOfLabels(const INDEX noLabels) { noLabels_ = noLabels; }

    void connect_pairwise_counting_center(PairwiseFactorType* p, COUNTING_FACTOR* c) {
       auto* m = new PAIRWISE_COUNTING_CENTER(false, p, c);
       pd_.GetLP().AddMessage(m);
       pd_.GetLP().AddFactorRelation(p,c);
    }
    void connect_pairwise_counting_left(PairwiseFactorType* p, COUNTING_FACTOR* c) {
       auto* m = new PAIRWISE_COUNTING_LEFT(false, p, c);
       pd_.GetLP().AddMessage(m);
       pd_.GetLP().AddFactorRelation(p,c);
    }
    void connect_pairwise_counting_right(PairwiseFactorType* p, COUNTING_FACTOR* c) {
       auto* m = new PAIRWISE_COUNTING_RIGHT(false, p, c);
       pd_.GetLP().AddMessage(m);
       pd_.GetLP().AddFactorRelation(p,c);
    }
    void connect_pairwise_counting_left(PairwiseFactorType* p1, PairwiseFactorType* p2, COUNTING_FACTOR* c) {
       auto* m1 = new LEFT_PAIRWISE_COUNTING_LEFT(false,p1,c);
       pd_.GetLP().AddMessage(m1);
       auto* m2 = new RIGHT_PAIRWISE_COUNTING_LEFT(false,p2,c);
       pd_.GetLP().AddMessage(m2);
       pd_.GetLP().AddFactorRelation(p2,c);
    }
    void connect_pairwise_counting_right(PairwiseFactorType* p1, PairwiseFactorType* p2, COUNTING_FACTOR* c) {
       auto* m1 = new LEFT_PAIRWISE_COUNTING_RIGHT(false,p1,c);
       pd_.GetLP().AddMessage(m1);
       auto* m2 = new RIGHT_PAIRWISE_COUNTING_RIGHT(false,p2,c);
       pd_.GetLP().AddMessage(m2);
       pd_.GetLP().AddFactorRelation(p2,c);
    }

    void connect_counting_factors_left(COUNTING_FACTOR* c_left, COUNTING_FACTOR* c_top) {
       pd_.GetLP().AddFactorRelation(c_left,c_top);
       auto* m = new COUNTING_MESSAGE_LEFT(c_left, c_top);
       pd_.GetLP().AddMessage(m);
    }
    void connect_counting_factors_right(COUNTING_FACTOR* c_right, COUNTING_FACTOR* c_top) {
       pd_.GetLP().AddFactorRelation(c_right,c_top);
       auto* m = new COUNTING_MESSAGE_RIGHT(c_right, c_top);
       pd_.GetLP().AddMessage(m);
    }

    void AddProjection(const std::vector<INDEX>& projectionVar, const std::vector<REAL>& summationCost)
    {      
       const INDEX max_sum = std::max(noLabels_,INDEX(summationCost.size()));
       auto& mrf = pd_.template GetProblemConstructor<MRF_PROBLEM_CONSTRUCTOR_NO>();
       assert(projectionVar.size() > 3); // otherwise just use a ternary factor to enforce consistency
       assert(std::is_sorted(projectionVar.begin(), projectionVar.end())); // support unsorted projectionVar (transpose in messages) later

      auto& mrfConstructor = pd_.template GetProblemConstructor<MRF_PROBLEM_CONSTRUCTOR_NO>();

      for(INDEX i=0;i<projectionVar.size()-1;++i) {
        const INDEX i1 = std::min(projectionVar[i],projectionVar[i+1]);
        const INDEX i2 = std::max(projectionVar[i],projectionVar[i+1]);

        if(!mrfConstructor.HasPairwiseFactor(i1,i2)) {
          mrfConstructor.AddPairwiseFactor(i1,i2,std::vector<REAL>(pow(noLabels_,2),0.0));
        }
      }
      
      struct counting_factor_rec {
         COUNTING_FACTOR* f;
         INDEX rightmost_var;
      };
      std::deque<counting_factor_rec> queue;

      // partition projectionVar into subsets of size at least 4 and at most 6
      INDEX begin_rest;
      // do zrobienia: make test for one-dimensional tomography problem on 4,5,6,7,8,9,10,11 and 12 variables
      switch(projectionVar.size() % 6) {
         case 0: begin_rest = 0; break;// nothing to do
         case 1: {// use two counting factors. Right counting subfactor of left counting factor is connected to top subcounting factor of right counting factor 
                 auto* f1 = new COUNTING_FACTOR(noLabels_, 1, 1, std::min(2*noLabels_,max_sum));
                 pd_.GetLP().AddFactor(f1);
                 connect_pairwise_counting_left(mrf.GetPairwiseFactor(projectionVar[0],projectionVar[1]),f1);
                 connect_pairwise_counting_center(mrf.GetPairwiseFactor(projectionVar[1],projectionVar[2]),f1);
                 connect_pairwise_counting_right(mrf.GetPairwiseFactor(projectionVar[2],projectionVar[3]),f1);
                 pd_.GetLP().AddFactorRelation(f1, mrf.GetUnaryFactor(projectionVar[3]));

                 auto* f2 = new COUNTING_FACTOR(noLabels_, f1->GetFactor()->up_sum_size(), 1, std::min(f1->GetFactor()->up_sum_size() + 2*noLabels_, max_sum));
                 pd_.GetLP().AddFactor(f2);
                 connect_counting_factors_left(f1,f2);
                 connect_pairwise_counting_center(mrf.GetPairwiseFactor(projectionVar[3],projectionVar[4]),f2);
                 connect_pairwise_counting_right(mrf.GetPairwiseFactor(projectionVar[4],projectionVar[5]),f2);
                 pd_.GetLP().AddFactorRelation(f2, mrf.GetUnaryFactor(projectionVar[5]));

                 begin_rest = 7;
                 queue.push_back({f2,6});
                 
                 break;
                 }
         case 2: {// for the remaining cases we use two separate counting factors.
                 auto* f1 = new COUNTING_FACTOR(noLabels_, 1, 1, std::min(2*noLabels_, max_sum));
                 pd_.GetLP().AddFactor(f1);
                 connect_pairwise_counting_left(mrf.GetPairwiseFactor(projectionVar[0],projectionVar[1]), f1);
                 connect_pairwise_counting_center(mrf.GetPairwiseFactor(projectionVar[1],projectionVar[2]),f1);
                 connect_pairwise_counting_right(mrf.GetPairwiseFactor(projectionVar[2],projectionVar[3]),f1);
                 pd_.GetLP().AddFactorRelation(f1, mrf.GetUnaryFactor(projectionVar[3]));
                 queue.push_back({f1,3});

                 auto* f2 = new COUNTING_FACTOR(noLabels_, 1, 1, std::min(2*noLabels_, max_sum));
                 pd_.GetLP().AddFactor(f2);
                 connect_pairwise_counting_left(mrf.GetPairwiseFactor(projectionVar[4],projectionVar[5]), f2);
                 connect_pairwise_counting_center(mrf.GetPairwiseFactor(projectionVar[5],projectionVar[6]),f2);
                 connect_pairwise_counting_right(mrf.GetPairwiseFactor(projectionVar[6],projectionVar[7]),f2);
                 pd_.GetLP().AddFactorRelation(f2, mrf.GetUnaryFactor(projectionVar[7]));
                 queue.push_back({f2,7});

                 begin_rest = 8;
                 break;
                 }
         case 3: {
                 auto* f1 = new COUNTING_FACTOR(noLabels_, noLabels_, 1, std::min(3*noLabels_, max_sum));
                 pd_.GetLP().AddFactor(f1);
                 connect_pairwise_counting_left(mrf.GetPairwiseFactor(projectionVar[0],projectionVar[1]), mrf.GetPairwiseFactor(projectionVar[1], projectionVar[2]), f1);
                 connect_pairwise_counting_center(mrf.GetPairwiseFactor(projectionVar[2],projectionVar[3]),f1);
                 connect_pairwise_counting_right(mrf.GetPairwiseFactor(projectionVar[3],projectionVar[4]),f1);
                 pd_.GetLP().AddFactorRelation(f1, mrf.GetUnaryFactor(projectionVar[4]));
                 queue.push_back({f1,4});

                 auto* f2 = new COUNTING_FACTOR(noLabels_, 1, 1, std::min(2*noLabels_, max_sum));
                 pd_.GetLP().AddFactor(f2);
                 connect_pairwise_counting_left(mrf.GetPairwiseFactor(projectionVar[5],projectionVar[6]), f2);
                 connect_pairwise_counting_center(mrf.GetPairwiseFactor(projectionVar[6],projectionVar[7]),f2);
                 connect_pairwise_counting_right(mrf.GetPairwiseFactor(projectionVar[7],projectionVar[8]),f2);
                 pd_.GetLP().AddFactorRelation(f2, mrf.GetUnaryFactor(projectionVar[8]));
                 queue.push_back({f2,8});

                 begin_rest = 9;
                 break;
                 }
         case 4: {
                 auto* f = new COUNTING_FACTOR(noLabels_, 1, 1, std::min(2*noLabels_, max_sum));
                 pd_.GetLP().AddFactor(f);
                 connect_pairwise_counting_left(mrf.GetPairwiseFactor(projectionVar[0],projectionVar[1]), f);
                 connect_pairwise_counting_center(mrf.GetPairwiseFactor(projectionVar[1],projectionVar[2]),f);
                 connect_pairwise_counting_right(mrf.GetPairwiseFactor(projectionVar[2],projectionVar[3]), f);
                 pd_.GetLP().AddFactorRelation(f, mrf.GetUnaryFactor(projectionVar[3]));
                 queue.push_back({f,3});

                 begin_rest = 4;
                 break;
                 }
         case 5: {
                 auto* f = new COUNTING_FACTOR(noLabels_, noLabels_, 1, std::min(3*noLabels_, max_sum));
                 pd_.GetLP().AddFactor(f);
                 connect_pairwise_counting_left(mrf.GetPairwiseFactor(projectionVar[0],projectionVar[1]), mrf.GetPairwiseFactor(projectionVar[1], projectionVar[2]), f);
                 connect_pairwise_counting_center(mrf.GetPairwiseFactor(projectionVar[2],projectionVar[3]),f);
                 connect_pairwise_counting_right(mrf.GetPairwiseFactor(projectionVar[3],projectionVar[4]), f);
                 pd_.GetLP().AddFactorRelation(f, mrf.GetUnaryFactor(projectionVar[4]));
                 queue.push_back({f,4});

                 begin_rest = 5;
                 break; 
                 }
      }
      assert((projectionVar.size() - begin_rest) % 6 == 0);

      for(INDEX i=begin_rest;i<projectionVar.size(); i+=6){ // a counting factor on base level can take care of six unaries
         auto* f = new COUNTING_FACTOR(noLabels_, noLabels_, noLabels_, std::min(4*noLabels_, max_sum));
         pd_.GetLP().AddFactor(f);
         connect_pairwise_counting_left(mrf.GetPairwiseFactor(projectionVar[i],projectionVar[i+1]), mrf.GetPairwiseFactor(projectionVar[i+1], projectionVar[i+2]), f);
         connect_pairwise_counting_center(mrf.GetPairwiseFactor(projectionVar[i+2],projectionVar[i+3]),f);
         connect_pairwise_counting_right(mrf.GetPairwiseFactor(projectionVar[i+3],projectionVar[i+4]), mrf.GetPairwiseFactor(projectionVar[i+4], projectionVar[i+5]), f);
         pd_.GetLP().AddFactorRelation(f, mrf.GetUnaryFactor(projectionVar[i+5]));
         queue.push_back({f,i+5});
      }

      // recursively join factors on lower levels
      while(queue.size() > 1) {
         const INDEX q_size = queue.size();
         counting_factor_rec first_counting = queue.front();
         if(q_size%2 == 1) { 
            queue.pop_front(); 
            assert(queue.size() == q_size-1);
         }
         assert(queue.size()%2 == 0);
         for(INDEX i=q_size%2; i<q_size; i+=2) {
            auto left = queue.front();
            queue.pop_front();
            auto right = queue.front();
            queue.pop_front();
            assert(left.rightmost_var < right.rightmost_var);

            const INDEX left_sum_size = left.f->GetFactor()->up_sum_size();
            const INDEX right_sum_size = right.f->GetFactor()->up_sum_size();
            auto* f = new COUNTING_FACTOR(noLabels_,  left_sum_size, right_sum_size, std::min(left_sum_size + 2*noLabels_ + right_sum_size, max_sum));
            pd_.GetLP().AddFactor(f);
            connect_counting_factors_left(left.f, f);
            assert(left.rightmost_var+1 < projectionVar.size());
            connect_pairwise_counting_center( mrf.GetPairwiseFactor(projectionVar[left.rightmost_var], projectionVar[left.rightmost_var+1]), f);
            connect_counting_factors_right(right.f,f);

            queue.push_back({f,right.rightmost_var}); 
         }
         if(q_size%2 == 1) { 
            queue.push_front(first_counting);
         } 
      }

      assert(queue.size() == 1);
      // set projection constraints in top factor
      auto* f = queue.front().f;
      f->GetFactor()->summation_cost(summationCost); 
      assert(f->GetFactor()->LowerBound() < std::numeric_limits<REAL>::infinity());
    }

  private:
    //std::vector<Tree> treeIndices_;
    INDEX noLabels_ = 0;
    Solver<FMC>& pd_;
  };

// construct one-dimensional discrete tomography problem with sequential factors
template<typename FMC,
     INDEX MRF_PROBLEM_CONSTRUCTOR_NO,
     typename SUM_FACTOR,
     typename SUM_PAIRWISE_FACTOR,
     typename SUM_PAIRWISE_MESSAGE_LEFT,
     typename SUM_PAIRWISE_MESSAGE_RIGHT,
     typename SUM_PAIRWISE_PAIRWISE_MESSAGE>
class dt_sequential_constructor {
public:
   using MrfConstructorType =
      typename meta::at_c<typename FMC::ProblemDecompositionList,MRF_PROBLEM_CONSTRUCTOR_NO>;
   using PairwiseFactorType = typename MrfConstructorType::PairwiseFactorContainer;

   dt_sequential_constructor(Solver<FMC>& pd) : pd_(pd) {}

   void SetNumberOfLabels(const INDEX noLabels) { noLabels_ = noLabels; }

   void AddProjection(const std::vector<INDEX>& projectionVar, const std::vector<REAL>& summationCost)
   { 
      assert(summationCost.size() > 0);
      const INDEX max_sum = std::max(noLabels_,INDEX(summationCost.size()));
      auto& mrf = pd_.template GetProblemConstructor<MRF_PROBLEM_CONSTRUCTOR_NO>();
      assert(std::is_sorted(projectionVar.begin(), projectionVar.end())); // support unsorted projectionVar (transpose in messages) later

      auto& mrfConstructor = pd_.template GetProblemConstructor<MRF_PROBLEM_CONSTRUCTOR_NO>();

      for(INDEX i=0;i<projectionVar.size()-1;++i) {
         const INDEX i1 = std::min(projectionVar[i],projectionVar[i+1]);
         const INDEX i2 = std::max(projectionVar[i],projectionVar[i+1]);

         if(!mrfConstructor.HasPairwiseFactor(i1,i2)) {
            mrfConstructor.AddPairwiseFactor(i1,i2,std::vector<REAL>(pow(noLabels_,2),0.0));
         }
      }

      auto* f_prev = new SUM_FACTOR(noLabels_, 1);
      pd_.GetLP().AddFactor(f_prev);
      for(INDEX i=1; i<projectionVar.size(); ++i) {
         const INDEX sum_size = std::min(i*noLabels_, max_sum);
         auto* f = new SUM_FACTOR(noLabels_, sum_size);
         pd_.GetLP().AddFactor(f);
         auto* f_p = new SUM_PAIRWISE_FACTOR(noLabels_, f_prev->GetFactor()->sum_size(), sum_size);
         pd_.GetLP().AddFactor(f_p);
         auto* m_l = new SUM_PAIRWISE_MESSAGE_LEFT(false,f_prev,f_p);
         pd_.GetLP().AddMessage(m_l);
         auto* m_r = new SUM_PAIRWISE_MESSAGE_RIGHT(false,f,f_p);
         pd_.GetLP().AddMessage(m_r);
         auto* m_c = new SUM_PAIRWISE_PAIRWISE_MESSAGE(false,mrf.GetPairwiseFactor(projectionVar[i-1], projectionVar[i]),f_p);
         pd_.GetLP().AddMessage(m_c);
         f_prev = f;

         pd_.GetLP().AddFactorRelation(f_prev,f_p);
         pd_.GetLP().AddFactorRelation(f_p,f);
         pd_.GetLP().AddFactorRelation(f_p,mrf.GetUnaryFactor(projectionVar[i]));
         pd_.GetLP().AddFactorRelation(mrf.GetUnaryFactor(projectionVar[i]), f_p);
      }

      f_prev->GetFactor()->summation_cost(summationCost);
   }
private:
   INDEX noLabels_ = 0;
   Solver<FMC>& pd_;


};

}

#endif // LP_MP_xxx
