#ifndef LP_MP_TOMOGRAPHY_H
#define LP_MP_TOMOGRAPHY_H

/*
  This file defines the FMC (factors,messages,...) for the tomography solver
*/

namespace LP_MP{

  typedef UnaryLoop<> UnaryLoopType;
  typedef PairwiseLoop<0> LeftLoopType;
  typedef PairwiseLoop<1> RightLoopType;

  typedef MultiplexMargMessage<UnaryLoopType,LeftLoopType,true,false,false,true> LeftMargMessage;
  typedef MultiplexMargMessage<UnaryLoopType,RightLoopType,true,false,false,true> RightMargMessage;

  typedef MultiplexFactor<std::vector<REAL>, const_ones_array, const_one> Simplex;


  struct FMC_DT {

    typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_MP_PARAM, 0, true, true > UnaryFactor;
    typedef FactorContainer<Simplex, ExplicitRepamStorage, FMC_MP_PARAM, 1, false, false > PairwiseFactor;
   
    typedef MessageContainer<LeftMargMessage, StandardMessageStorage, FMC_MP_PARAM, 1 > UnaryPairwiseMessageLeft;
    typedef MessageContainer<RightMargMessage, StandardMessageStorage, FMC_MP_PARAM, 2 > UnaryPairwiseMessageRight;
    using FactorList = meta::list< UnaryFactor, PairwiseFactor >;
    using MessageList = meta::list<
      MessageListItem< UnaryPairwiseMessageLeft,  0, 1, std::vector, FixedSizeMessageContainer<1>::type >,
      MessageListItem< UnaryPairwiseMessageRight, 0, 1, std::vector, FixedSizeMessageContainer<1>::type >
      >;

    using mrf = MRFProblemConstructor<FMC_MP_PARAM,0,1,0,1>;
    using ProblemDecompositionList = meta::list<mrf>;
	  
  };  

}
#endif // LP_MP_TOMOGRAPHY_H
