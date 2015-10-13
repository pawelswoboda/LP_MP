#ifndef LP_MP_PROBLEM_CONSTRUCTION_HELPER
#define LP_MP_PROBLEM_CONSTRUCTION_HELPER

#include "messages/equality_message.hxx"

namespace LP_MP {

template<class MULTIPLEX_CONTAINER>
MULTIPLEX_CONTAINER* AddMultiplexFactor(LP* lp, const std::vector<REAL>& cost, const std::vector<INDEX>& varCapacity, const INDEX sum)
{
   assert(cost.size() == varCapacity.size());
   typedef typename MULTIPLEX_CONTAINER::FactorType MultiplexType;
   MULTIPLEX_CONTAINER * f = new MULTIPLEX_CONTAINER(MultiplexType(cost, varCapacity, sum), cost);
   lp->AddFactor(f);
   return f;
}

template<class MULTIPLEX_CONTAINER>
MULTIPLEX_CONTAINER* AddSimplexFactor(LP* lp, const std::vector<REAL>& cost)
{
   std::vector<INDEX> varCapacity(cost.size(), 1);
   INDEX sum = 1;
   return AddMultiplexFactor<MULTIPLEX_CONTAINER>(lp, cost, varCapacity, sum);
}

template<class MULTIPLEX_CONTAINER>
MULTIPLEX_CONTAINER* AddPairwiseMultiplexFactor(LP* lp, const std::vector<REAL>& cost, 
      const std::vector<INDEX>& capacity1, const std::vector<INDEX>& capacity2, 
      const INDEX sum1, const INDEX sum2)
{
   std::vector<INDEX> pairwise_capacity(cost.size(),1);
   assert(cost.size() == capacity1.size()*capacity2.size());

   for(INDEX i1_index=0; i1_index<capacity1.size(); i1_index++)
      for(INDEX i2_index=0; i2_index<capacity2.size(); i2_index++)
         pairwise_capacity[i1_index + i2_index*capacity1.size()] = capacity1[i1_index]*capacity2[i2_index];

   const INDEX pairwise_sum = sum1*sum2; // number of possible assignments
   
   // for non SIMD
   auto f = AddMultiplexFactor<MULTIPLEX_CONTAINER>(lp, cost, pairwise_capacity, pairwise_sum); // for non SIMD
   // for SIMD
   //typedef typename MULTIPLEX_CONTAINER::FactorType MultiplexType;
   //MULTIPLEX_CONTAINER* f = new MULTIPLEX_CONTAINER(MultiplexType(capacity1.size(), capacity2.size(), cost, pairwise_capacity, pairwise_sum), cost);
   //lp->AddFactor(f);
   return f;
}

/*
template<class MESSAGE_CONTAINER, class LEFT_FACTOR_TYPE, class RIGHT_FACTOR_TYPE, INDEX COMMON_VAR_LEFT, INDEX COMMON_VAR_RIGHT>
void LinkPairwiseMultiplexFactorsImpl(
      LP* lp, 
      LEFT_FACTOR_TYPE* left, RIGHT_FACTOR_TYPE* right, 
      const INDEX lX, const INDEX lY, 
      const INDEX rX, const INDEX rY)
{
   static_assert(0 <= COMMON_VAR_LEFT && COMMON_VAR_LEFT <= 1 && 0 <= COMMON_VAR_RIGHT && COMMON_VAR_RIGHT <= 1, "template parameters must be 0/1");
   //assert(typeid(left) == typeid(MultiplexFactor*));
   //assert(typeid(right) == typeid(MultiplexFactor*));

   typedef typename MESSAGE_CONTAINER::MessageType::LeftLoopType LeftLoopType;
   typedef typename MESSAGE_CONTAINER::MessageType::RightLoopType RightLoopType;
   //typedef PairwiseLoop<COMMON_VAR_LEFT> LeftLoopType;
   //typedef PairwiseLoop<COMMON_VAR_RIGHT> RightLoopType;
   const std::array<INDEX,2> l = {{lX,lY}};
   LeftLoopType ll( l );
   std::array<INDEX,2> r = {{rX,rY}};
   RightLoopType lr( r );

   INDEX commonVarDim;
   if(COMMON_VAR_LEFT == 0) {
      commonVarDim = lX;
   } else {
      commonVarDim = lY;
   } 
   if(COMMON_VAR_RIGHT == 0) {
      assert( rX == commonVarDim);
   } else {
      assert( rY == commonVarDim);
   } 

   // do zrobienia: prawidlowo (left|right)->GetSum()
   MESSAGE_CONTAINER* m = new MESSAGE_CONTAINER(MultiplexMargMessage<LeftLoopType,RightLoopType>(ll,lr,right->GetFactor()->GetSum(), left->GetFactor()->GetSum()), left, right, commonVarDim);

   lp->AddMessage(m);
}

// MESSAGE_CONTAINER must contain MultiplexMargMessage
template<class MESSAGE_CONTAINER, class LEFT_FACTOR_TYPE, class RIGHT_FACTOR_TYPE>
void LinkPairwiseMultiplexFactors(
      LP* lp, 
      LEFT_FACTOR_TYPE* left, RIGHT_FACTOR_TYPE* right, 
      const INDEX lX, const INDEX lY, const INDEX commonVarLeft, 
      const INDEX rX, const INDEX rY, const INDEX commonVarRight)
{
   if(commonVarLeft == 1 && commonVarRight == 1) 
      LinkPairwiseMultiplexFactorsImpl<MESSAGE_CONTAINER,LEFT_FACTOR_TYPE,RIGHT_FACTOR_TYPE,1,1>(lp, left, right, lX, lY, rX, rY);
   else if(commonVarLeft == 1 && commonVarRight == 0) 
      LinkPairwiseMultiplexFactorsImpl<MESSAGE_CONTAINER,LEFT_FACTOR_TYPE,RIGHT_FACTOR_TYPE,1,0>(lp, left, right, lX, lY, rX, rY);
   else if(commonVarLeft == 0 && commonVarRight == 1) 
      LinkPairwiseMultiplexFactorsImpl<MESSAGE_CONTAINER,LEFT_FACTOR_TYPE,RIGHT_FACTOR_TYPE,0,1>(lp, left, right, lX, lY, rX, rY);
   else if(commonVarLeft == 0 && commonVarRight == 0) 
      LinkPairwiseMultiplexFactorsImpl<MESSAGE_CONTAINER,LEFT_FACTOR_TYPE,RIGHT_FACTOR_TYPE,0,0>(lp, left, right, lX, lY, rX, rY);
   else throw std::runtime_error("common variable indicators must be 0/1");
}
*/

template<class MESSAGE_CONTAINER, class UNARY_FACTOR_TYPE, class PAIRWISE_FACTOR_TYPE>
void LinkUnaryPairwiseMultiplexFactors(
      LP* lp,
      UNARY_FACTOR_TYPE* unary, PAIRWISE_FACTOR_TYPE* pairwise,
      const INDEX left_dim, const INDEX right_dim, const INDEX common_var)
{
   assert(common_var == 0 || common_var == 1);
   typedef typename MESSAGE_CONTAINER::MessageType::LeftLoopType UnaryLoopType;
   typedef typename MESSAGE_CONTAINER::MessageType::RightLoopType PairwiseLoopType;

   const INDEX common_var_dim = common_var == 0 ? left_dim : right_dim;

   UnaryLoopType ll(common_var_dim);
   std::array<INDEX,2> dim = {{left_dim, right_dim}};
   PairwiseLoopType lr( dim );

   typedef typename MESSAGE_CONTAINER::MessageType MultiplexMargMessageType;

   MESSAGE_CONTAINER* m = new MESSAGE_CONTAINER(MultiplexMargMessageType(ll, lr, unary->GetFactor()->GetSum(), pairwise->GetFactor()->GetSum()), unary, pairwise, common_var_dim);
   lp->AddMessage(m);
}

// MESSAGE_CONTAINER must contain EqualityMessage
template<class MESSAGE_CONTAINER, class LEFT_FACTOR_TYPE, class RIGHT_FACTOR_TYPE>
void LinkSingleVariableEqualityMessage(LP* lp, LEFT_FACTOR_TYPE* left_fac, const INDEX left_var, RIGHT_FACTOR_TYPE* right_fac, const INDEX right_var)
{
   MESSAGE_CONTAINER* m = new MESSAGE_CONTAINER(EqualityMessage(left_var,right_var),left_fac,right_fac,1);
   lp->AddMessage(m);
}


} // end namespace LP_MP

#endif // LP_MP_PROBLEM_CONSTRUCTION_HELPER
