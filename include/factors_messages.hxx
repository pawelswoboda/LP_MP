#ifndef LP_MP_FACTORS_MESSAGES_HXX
#define LP_MP_FACTORS_MESSAGES_HXX

#include <tuple>
#include <array>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <functional>
#include <utility>
#include <memory>
#include <limits>
#include <exception>
#include <typeinfo>
#include <type_traits>
#include <variant>
#include <assert.h>
#include <cxxabi.h>

#include <mutex>

#include "template_utilities.hxx"
#include "function_existence.hxx"
#include "meta/meta.hpp"
#include "static_if.hxx" // remove
#include "MemoryPool.h"

#include "memory_allocator.hxx"

#include "LP_MP.h"

// do zrobienia: remove these
//#include "factors/reparametrization_storage.hxx"  // also delete file
//#include "messages/message_storage.hxx"
#include <fstream>
#include <sstream>
#include <vector>
#include "vector"

#include "DD_ILP.hxx"

// this file provides message and factor containers. The factors and messages are plugged into the container and then every method call is dispatched correctly with static polymorphism and template tricks.

// do zrobienia: Introduce MessageConstraint and FactorConstraint for templates
// cleanup name inconsistencies: MessageType, MessageDispatcher etc

namespace LP_MP {

// we must check existence of functions in message classes. The necessary test code is concentrated here. 
namespace FunctionExistence {

// Macros to construct help functions for checking existence of member functions of classes
LP_MP_FUNCTION_EXISTENCE_CLASS(HasReceiveMessageFromRight,ReceiveMessageFromRight)
LP_MP_FUNCTION_EXISTENCE_CLASS(HasReceiveMessageFromLeft, ReceiveMessageFromLeft)
   
LP_MP_FUNCTION_EXISTENCE_CLASS(HasReceiveRestrictedMessageFromRight,ReceiveRestrictedMessageFromRight)
LP_MP_FUNCTION_EXISTENCE_CLASS(HasReceiveRestrictedMessageFromLeft, ReceiveRestrictedMessageFromLeft)

LP_MP_FUNCTION_EXISTENCE_CLASS(HasSendMessagesToRight,SendMessagesToRight)
LP_MP_FUNCTION_EXISTENCE_CLASS(HasSendMessagesToLeft, SendMessagesToLeft)

LP_MP_FUNCTION_EXISTENCE_CLASS(has_send_message_to_right_improvement, send_message_to_right_improvement)
LP_MP_FUNCTION_EXISTENCE_CLASS(has_send_message_to_left_improvement, send_message_to_left_improvement)
LP_MP_FUNCTION_EXISTENCE_CLASS(has_send_messages_to_right_improvement, send_messages_to_right_improvement)
LP_MP_FUNCTION_EXISTENCE_CLASS(has_send_messages_to_left_improvement, send_messages_to_left_improvement)

LP_MP_FUNCTION_EXISTENCE_CLASS(HasRepamRight, RepamRight)
LP_MP_FUNCTION_EXISTENCE_CLASS(HasRepamLeft, RepamLeft)

LP_MP_FUNCTION_EXISTENCE_CLASS(HasComputeLeftFromRightPrimal, ComputeLeftFromRightPrimal)
LP_MP_FUNCTION_EXISTENCE_CLASS(HasComputeRightFromLeftPrimal, ComputeRightFromLeftPrimal)

LP_MP_FUNCTION_EXISTENCE_CLASS(HasCheckPrimalConsistency, CheckPrimalConsistency)

LP_MP_FUNCTION_EXISTENCE_CLASS(HasPropagatePrimal, PropagatePrimal)
LP_MP_FUNCTION_EXISTENCE_CLASS(HasMaximizePotential, MaximizePotential)
LP_MP_FUNCTION_EXISTENCE_CLASS(HasMaximizePotentialAndComputePrimal, MaximizePotentialAndComputePrimal)

LP_MP_FUNCTION_EXISTENCE_CLASS(has_apply, apply)

LP_MP_FUNCTION_EXISTENCE_CLASS(has_create_constraints, create_constraints)

LP_MP_ASSIGNMENT_FUNCTION_EXISTENCE_CLASS(IsAssignable, operator[])
}

template<typename FMC>
struct empty_message_fmc
{
using FactorList = typename FMC::FactorList;
using MessageList = meta::list<>;
};

// function getters for statically dispatching ReceiveMessage and SendMessage to left and right side correctly, used in FactorContainer
template<typename MSG_CONTAINER>
struct LeftMessageFuncGetter
{
   using ConnectedFactorType = typename MSG_CONTAINER::RightFactorContainer;

   //constexpr static decltype(&MSG_CONTAINER::GetLeftMessage) GetMessageFunc() { return &MSG_CONTAINER::GetLeftMessage; }

   constexpr static decltype(&MSG_CONTAINER::ReceiveMessageFromRightContainer) GetReceiveFunc() { return &MSG_CONTAINER::ReceiveMessageFromRightContainer; }

   constexpr static decltype(&MSG_CONTAINER::ReceiveRestrictedMessageFromRightContainer) GetReceiveRestrictedFunc() { return &MSG_CONTAINER::ReceiveRestrictedMessageFromRightContainer; }

   constexpr static decltype(&MSG_CONTAINER::SendMessageToRightContainer) GetSendFunc() { return &MSG_CONTAINER::SendMessageToRightContainer; }

   constexpr static decltype(&MSG_CONTAINER::test_send_message_to_right) get_test_send_message_func() { return &MSG_CONTAINER::test_send_message_to_right; }

   template<typename LEFT_FACTOR, typename MSG_ITERATOR>
   constexpr static decltype(&MSG_CONTAINER::template SendMessagesToRightContainer<LEFT_FACTOR, MSG_ITERATOR>) GetSendMessagesFunc() 
   { return &MSG_CONTAINER::template SendMessagesToRightContainer<LEFT_FACTOR, MSG_ITERATOR>; }

   template<typename MSG_ITERATOR>
   constexpr static decltype(&MSG_CONTAINER::template test_send_messages_to_right<MSG_ITERATOR>) get_test_send_messages_func() 
   { return &MSG_CONTAINER::template test_send_messages_to_right<MSG_ITERATOR>; }

   constexpr static decltype(&MSG_CONTAINER::send_message_to_right_improvement) get_send_message_improvement_func()
   { return &MSG_CONTAINER::send_message_to_right_improvement; }
   template<typename MSG_ITERATOR>
   constexpr static decltype(&MSG_CONTAINER::template send_messages_to_right_improvement<MSG_ITERATOR>) get_send_messages_improvement_func()
   { return &MSG_CONTAINER::template send_messages_to_right_improvement<MSG_ITERATOR>; }

#ifdef LP_MP_PARALLEL
   constexpr static decltype(&MSG_CONTAINER::ReceiveMessageFromRightContainerSynchronized) GetReceiveSynchronizedFunc() { return &MSG_CONTAINER::ReceiveMessageFromRightContainerSynchronized; }
   constexpr static decltype(&MSG_CONTAINER::SendMessageToRightContainerSynchronized) GetSendSynchronizedFunc() { return &MSG_CONTAINER::SendMessageToRightContainerSynchronized; }
   template<typename LEFT_FACTOR, typename MSG_ARRAY, typename ITERATOR>
   constexpr static decltype(&MSG_CONTAINER::template SendMessagesToRightContainerSynchronized<LEFT_FACTOR, MSG_ARRAY, ITERATOR>) GetSendMessagesSynchronizedFunc() 
   { return &MSG_CONTAINER::template SendMessagesToRightContainerSynchronized<LEFT_FACTOR, MSG_ARRAY, ITERATOR>; }
#endif

   constexpr static bool sends_message_to_adjacent_factor() { return MSG_CONTAINER::sends_message_to_right_constexpr(); }
   constexpr static bool receives_message_from_adjacent_factor() { return MSG_CONTAINER::receives_message_from_right_constexpr(); }
   constexpr static bool adjacent_factor_sends_message() { return MSG_CONTAINER::sends_message_to_left_constexpr(); }
   constexpr static bool adjacent_factor_receives_message() { return MSG_CONTAINER::receives_message_from_left_constexpr(); }

   constexpr static bool
   CanCallReceiveRestrictedMessage()
   { return MSG_CONTAINER::CanCallReceiveRestrictedMessageFromRightContainer(); }

   constexpr static bool 
   CanCallSendMessages()
   { return MSG_CONTAINER::CanCallSendMessagesToRightContainer(); }

   // do zrobienia: rename CanPropagatePrimalThroughMessage
   constexpr static bool CanComputePrimalThroughMessage()
   { return MSG_CONTAINER::CanComputeRightFromLeftPrimal(); }
   constexpr static decltype(&MSG_CONTAINER::ComputeRightFromLeftPrimal) GetComputePrimalThroughMessageFunc()
   { return &MSG_CONTAINER::ComputeRightFromLeftPrimal; }

   constexpr static Chirality get_chirality() { return Chirality::left; }

   static auto* get_adjacent_factor(const MSG_CONTAINER& msg) { return msg.GetRightFactor(); }
};

template<typename MSG_CONTAINER>
struct RightMessageFuncGetter
{
   using ConnectedFactorType = typename MSG_CONTAINER::LeftFactorContainer;

   //constexpr static decltype(&MSG_CONTAINER::GetRightMessage) GetMessageFunc() { return &MSG_CONTAINER::GetRightMessage; }

   constexpr static decltype(&MSG_CONTAINER::ReceiveMessageFromLeftContainer) GetReceiveFunc() { return &MSG_CONTAINER::ReceiveMessageFromLeftContainer; }

   constexpr static decltype(&MSG_CONTAINER::ReceiveRestrictedMessageFromLeftContainer) GetReceiveRestrictedFunc() { return &MSG_CONTAINER::ReceiveRestrictedMessageFromLeftContainer; }

   constexpr static decltype(&MSG_CONTAINER::SendMessageToLeftContainer) GetSendFunc() { return &MSG_CONTAINER::SendMessageToLeftContainer; }

   constexpr static decltype(&MSG_CONTAINER::test_send_message_to_left) get_test_send_message_func() { return &MSG_CONTAINER::test_send_message_to_left; }

   template<typename RIGHT_FACTOR, typename MSG_ITERATOR>
   constexpr static decltype(&MSG_CONTAINER::template SendMessagesToLeftContainer<RIGHT_FACTOR, MSG_ITERATOR>) GetSendMessagesFunc() 
   { return &MSG_CONTAINER::template SendMessagesToLeftContainer<RIGHT_FACTOR, MSG_ITERATOR>; }

   template<typename MSG_ITERATOR>
   constexpr static decltype(&MSG_CONTAINER::template test_send_messages_to_left<MSG_ITERATOR>) get_test_send_messages_func() 
   { return &MSG_CONTAINER::template test_send_messages_to_left<MSG_ITERATOR>; }

   constexpr static decltype(&MSG_CONTAINER::send_message_to_left_improvement) get_send_message_improvement_func()
   { return &MSG_CONTAINER::send_message_to_left_improvement; }

   template<typename MSG_ITERATOR>
   constexpr static decltype(&MSG_CONTAINER::template send_messages_to_left_improvement<MSG_ITERATOR>) get_send_messages_improvement_func()
   { return &MSG_CONTAINER::template send_messages_to_left_improvement<MSG_ITERATOR>; }

#ifdef LP_MP_PARALLEL
   constexpr static decltype(&MSG_CONTAINER::ReceiveMessageFromLeftContainerSynchronized) GetReceiveSynchronizedFunc() { return &MSG_CONTAINER::ReceiveMessageFromLeftContainerSynchronized; }
   constexpr static decltype(&MSG_CONTAINER::SendMessageToLeftContainerSynchronized) GetSendSynchronizedFunc() { return &MSG_CONTAINER::SendMessageToLeftContainerSynchronized; }
   template<typename RIGHT_FACTOR, typename MSG_ARRAY, typename ITERATOR>
   constexpr static decltype(&MSG_CONTAINER::template SendMessagesToLeftContainerSynchronized<RIGHT_FACTOR, MSG_ARRAY, ITERATOR>) GetSendMessagesSynchronizedFunc() 
   { return &MSG_CONTAINER::template SendMessagesToLeftContainerSynchronized<RIGHT_FACTOR, MSG_ARRAY, ITERATOR>; }
#endif

   constexpr static bool sends_message_to_adjacent_factor() { return MSG_CONTAINER::sends_message_to_left_constexpr(); }
   constexpr static bool receives_message_from_adjacent_factor() { return MSG_CONTAINER::receives_message_from_left_constexpr(); }
   constexpr static bool adjacent_factor_sends_message() { return MSG_CONTAINER::sends_message_to_right_constexpr(); }
   constexpr static bool adjacent_factor_receives_message() { return MSG_CONTAINER::receives_message_from_right_constexpr(); }

   constexpr static bool CanCallReceiveRestrictedMessage() 
   { return MSG_CONTAINER::CanCallReceiveRestrictedMessageFromLeftContainer(); }

   constexpr static bool
   CanCallSendMessages()
   { return MSG_CONTAINER::CanCallSendMessagesToLeftContainer(); }

   constexpr static bool CanComputePrimalThroughMessage()
   { return MSG_CONTAINER::CanComputeLeftFromRightPrimal(); }
   constexpr static decltype(&MSG_CONTAINER::ComputeLeftFromRightPrimal) GetComputePrimalThroughMessageFunc()
   { return &MSG_CONTAINER::ComputeLeftFromRightPrimal; }

   constexpr static Chirality get_chirality() { return Chirality::right; }

   static auto* get_adjacent_factor(const MSG_CONTAINER& msg) { return msg.GetLeftFactor(); }
};

template<class MSG_CONTAINER, template<typename> class FuncGetter>
struct MessageDispatcher
{
   using ConnectedFactorType = typename FuncGetter<MSG_CONTAINER>::ConnectedFactorType; // this is the type of factor container to which the message is connected

   static void ReceiveMessage(MSG_CONTAINER& t)
   {
      auto staticMemberFunc = FuncGetter<MSG_CONTAINER>::GetReceiveFunc();
      (t.*staticMemberFunc)();
   }
#ifdef LP_MP_PARALLEL
   static void ReceiveMessageSynchronized(MSG_CONTAINER& t)
   {
      auto staticMemberFunc = FuncGetter<MSG_CONTAINER>::GetReceiveSynchronizedFunc();
      (t.*staticMemberFunc)();
   }
#endif

   constexpr static bool CanCallReceiveRestrictedMessage() { return FuncGetter<MSG_CONTAINER>::CanCallReceiveRestrictedMessage(); }
   static void ReceiveRestrictedMessage(MSG_CONTAINER& t)
   {
      auto staticMemberFunc = FuncGetter<MSG_CONTAINER>::GetReceiveRestrictedFunc();
      (t.*staticMemberFunc)();
   }

   constexpr static bool sends_message_to_adjacent_factor() { return FuncGetter<MSG_CONTAINER>::sends_message_to_adjacent_factor(); }
   constexpr static bool receives_message_from_adjacent_factor() { return FuncGetter<MSG_CONTAINER>::receives_message_from_adjacent_factor(); }
   constexpr static bool adjacent_factor_sends_message() { return FuncGetter<MSG_CONTAINER>::adjacent_factor_sends_message(); }
   constexpr static bool adjacent_factor_receives_message() { return FuncGetter<MSG_CONTAINER>::adjacent_factor_receives_message(); }

   template<typename FACTOR_TYPE>
   static void SendMessage(FACTOR_TYPE* f, MSG_CONTAINER& t, const REAL omega)
   {
      auto staticMemberFunc = FuncGetter<MSG_CONTAINER>::GetSendFunc();
      (t.*staticMemberFunc)(f, omega);
   }

   static void test_send_message(MSG_CONTAINER& t)
   {
      auto fun = FuncGetter<MSG_CONTAINER>::get_test_send_message_func();
      (t.*fun)(); 
   }

   static REAL send_message_improvement(MSG_CONTAINER& t)
   {
       auto static_member_func = FuncGetter<MSG_CONTAINER>::get_send_message_improvement_func();
       return (t.*static_member_func)();
   }


#ifdef LP_MP_PARALLEL
   template<typename FACTOR_TYPE>
   static void SendMessageSynchronized(FACTOR_TYPE* f, MSG_CONTAINER& t, const REAL omega)
   {
      auto staticMemberFunc = FuncGetter<MSG_CONTAINER>::GetSendSynchronizedFunc();
      (t.*staticMemberFunc)(f, omega);
   }
#endif

   // batch message sending
   constexpr static bool CanCallSendMessages() { return FuncGetter<MSG_CONTAINER>::CanCallSendMessages(); }

   template<typename FACTOR, typename MSG_ITERATOR>
   static void SendMessages(const FACTOR& f, MSG_ITERATOR msgs_begin, MSG_ITERATOR msgs_end, const REAL omega)
   {
      auto staticMemberFunc = FuncGetter<MSG_CONTAINER>::template GetSendMessagesFunc<FACTOR, MSG_ITERATOR>();
      (*staticMemberFunc)(f, msgs_begin, msgs_end, omega);
   }

   template<typename MSG_ITERATOR>
   static void test_send_messages(MSG_ITERATOR msg_begin, MSG_ITERATOR msg_end)
   {
      auto fun = FuncGetter<MSG_CONTAINER>::template get_test_send_messages_func<MSG_ITERATOR>();
      (*fun)(msg_begin, msg_end); 
   }


   template<typename MSG_ITERATOR>
   static REAL send_messages_improvement(MSG_ITERATOR msg_begin, MSG_ITERATOR msg_end)
   {
       auto static_member_func = FuncGetter<MSG_CONTAINER>::template get_send_messages_improvement_func<MSG_ITERATOR>();
       return (*static_member_func)(msg_begin, msg_end);
   }

#ifdef LP_MP_PARALLEL
   template<typename FACTOR, typename MSG_ARRAY, typename ITERATOR>
   static void SendMessagesSynchronized(const FACTOR& f, const MSG_ARRAY& msgs, ITERATOR omegaBegin)
   {
      auto staticMemberFunc = FuncGetter<MSG_CONTAINER>::template GetSendMessagesSynchronizedFunc<FACTOR, MSG_ARRAY, ITERATOR>();
      (*staticMemberFunc)(f, msgs, omegaBegin);
   }
#endif

   constexpr static bool CanComputePrimalThroughMessage() // do zrobienia: return false, if the factor from which this is called computes its own primal already
   {
      return FuncGetter<MSG_CONTAINER>::CanComputePrimalThroughMessage();
   }

   static void ComputePrimalThroughMessage(MSG_CONTAINER& t) 
   {
      auto staticMemberFunc = FuncGetter<MSG_CONTAINER>::GetComputePrimalThroughMessageFunc();
      return (t.*staticMemberFunc)();
   }
   constexpr static Chirality get_chirality() { return FuncGetter<MSG_CONTAINER>::get_chirality(); }

   static auto* get_adjacent_factor(const MSG_CONTAINER& msg) { return FuncGetter<MSG_CONTAINER>::get_adjacent_factor(msg); }
};

template<typename MSG_CONTAINER>
struct empty_next_left_message_container {
    MSG_CONTAINER* next_left_msg() const { assert(false); return nullptr; }
    void set_next_left_msg(MSG_CONTAINER* m) { assert(false); }
};

template<typename MSG_CONTAINER>
struct next_left_message_container {
    void set_next_left_msg(MSG_CONTAINER* m) { next = m; }
    MSG_CONTAINER* next_left_msg() const { return next; }
    MSG_CONTAINER* next = nullptr;
};
   
template<typename MSG_CONTAINER>
struct empty_next_right_message_container {
    MSG_CONTAINER* next_right_msg() const { assert(false); return nullptr; }
    void set_next_right_msg(MSG_CONTAINER*) { assert(false); }
};

template<typename MSG_CONTAINER>
struct next_right_message_container {
   void set_next_right_msg(MSG_CONTAINER* m) { next = m; }
   MSG_CONTAINER* next_right_msg() const { return next; }
   MSG_CONTAINER* next = nullptr;
};

// holds a pointer to next msg if messages are not linked through message storage held by factors
template<typename MESSAGE_CONTAINER_TYPE, SIGNED_INDEX NO_LEFT_FACTORS, SIGNED_INDEX NO_RIGHT_FACTORS, Chirality CHIRALITY>
struct next_msg_container_selector {
    class empty{};

   using type = 
   std::conditional_t<(CHIRALITY == Chirality::left && NO_LEFT_FACTORS < 0 && NO_RIGHT_FACTORS < 0), empty_next_left_message_container<MESSAGE_CONTAINER_TYPE>,
   std::conditional_t<(CHIRALITY == Chirality::left && NO_LEFT_FACTORS == 0 && NO_RIGHT_FACTORS < 0), next_left_message_container<MESSAGE_CONTAINER_TYPE>,
   std::conditional_t<(CHIRALITY == Chirality::left && NO_LEFT_FACTORS > 0 && NO_RIGHT_FACTORS < 0), empty_next_left_message_container<MESSAGE_CONTAINER_TYPE>,

   std::conditional_t<(CHIRALITY == Chirality::left && NO_LEFT_FACTORS < 0 && NO_RIGHT_FACTORS == 0), empty_next_left_message_container<MESSAGE_CONTAINER_TYPE>,
   std::conditional_t<(CHIRALITY == Chirality::left && NO_LEFT_FACTORS == 0 && NO_RIGHT_FACTORS == 0), empty_next_left_message_container<MESSAGE_CONTAINER_TYPE>,
   std::conditional_t<(CHIRALITY == Chirality::left && NO_LEFT_FACTORS > 0 && NO_RIGHT_FACTORS == 0), empty_next_left_message_container<MESSAGE_CONTAINER_TYPE>,

   std::conditional_t<(CHIRALITY == Chirality::left && NO_LEFT_FACTORS < 0 && NO_RIGHT_FACTORS > 0), empty_next_left_message_container<MESSAGE_CONTAINER_TYPE>,
   std::conditional_t<(CHIRALITY == Chirality::left && NO_LEFT_FACTORS == 0 && NO_RIGHT_FACTORS > 0), next_left_message_container<MESSAGE_CONTAINER_TYPE>,
   std::conditional_t<(CHIRALITY == Chirality::left && NO_LEFT_FACTORS > 0 && NO_RIGHT_FACTORS > 0), empty_next_left_message_container<MESSAGE_CONTAINER_TYPE>,
   

   std::conditional_t<(CHIRALITY == Chirality::right && NO_RIGHT_FACTORS > 0 && NO_LEFT_FACTORS < 0), empty_next_right_message_container<MESSAGE_CONTAINER_TYPE>,
   std::conditional_t<(CHIRALITY == Chirality::right && NO_RIGHT_FACTORS == 0 && NO_LEFT_FACTORS < 0), next_right_message_container<MESSAGE_CONTAINER_TYPE>,
   std::conditional_t<(CHIRALITY == Chirality::right && NO_RIGHT_FACTORS > 0 && NO_LEFT_FACTORS < 0), empty_next_right_message_container<MESSAGE_CONTAINER_TYPE>,

   std::conditional_t<(CHIRALITY == Chirality::right && NO_RIGHT_FACTORS > 0 && NO_LEFT_FACTORS == 0), empty_next_right_message_container<MESSAGE_CONTAINER_TYPE>,
   std::conditional_t<(CHIRALITY == Chirality::right && NO_RIGHT_FACTORS == 0 && NO_LEFT_FACTORS == 0), empty_next_right_message_container<MESSAGE_CONTAINER_TYPE>,
   std::conditional_t<(CHIRALITY == Chirality::right && NO_RIGHT_FACTORS > 0 && NO_LEFT_FACTORS == 0), empty_next_right_message_container<MESSAGE_CONTAINER_TYPE>,

   std::conditional_t<(CHIRALITY == Chirality::right && NO_RIGHT_FACTORS > 0 && NO_LEFT_FACTORS > 0), empty_next_right_message_container<MESSAGE_CONTAINER_TYPE>,
   std::conditional_t<(CHIRALITY == Chirality::right && NO_RIGHT_FACTORS == 0 && NO_LEFT_FACTORS > 0), next_right_message_container<MESSAGE_CONTAINER_TYPE>,
   std::conditional_t<(CHIRALITY == Chirality::right && NO_RIGHT_FACTORS > 0 && NO_LEFT_FACTORS > 0), empty_next_right_message_container<MESSAGE_CONTAINER_TYPE>,  
   empty
   >>> >>> >>>  >>> >>> >>>;

   ~next_msg_container_selector()
   {
       static_assert(!std::is_same_v<type, empty>, "");
   }
       
};


// Class holding message and left and right factor
// do zrobienia: possibly replace {LEFT|RIGHT}_FACTOR_NO by their type
template<typename MESSAGE_TYPE, 
         INDEX LEFT_FACTOR_NO, INDEX RIGHT_FACTOR_NO,
         message_passing_schedule MPS,
         SIGNED_INDEX NO_OF_LEFT_FACTORS, SIGNED_INDEX NO_OF_RIGHT_FACTORS,
         typename FACTOR_MESSAGE_TRAIT, 
         INDEX MESSAGE_NO,
         INDEX ESTIMATED_NO_OF_LEFT_FACTORS = 4, INDEX ESTIMATED_NO_OF_RIGHT_FACTORS = 4
         >
class MessageContainer : 
    public next_msg_container_selector< 
    MessageContainer<MESSAGE_TYPE,LEFT_FACTOR_NO,RIGHT_FACTOR_NO,MPS,NO_OF_LEFT_FACTORS,NO_OF_RIGHT_FACTORS,FACTOR_MESSAGE_TRAIT,MESSAGE_NO,ESTIMATED_NO_OF_LEFT_FACTORS,ESTIMATED_NO_OF_RIGHT_FACTORS>,
    NO_OF_LEFT_FACTORS, NO_OF_RIGHT_FACTORS, Chirality::left>::type,

    public next_msg_container_selector< 
    MessageContainer<MESSAGE_TYPE,LEFT_FACTOR_NO,RIGHT_FACTOR_NO,MPS,NO_OF_LEFT_FACTORS,NO_OF_RIGHT_FACTORS,FACTOR_MESSAGE_TRAIT,MESSAGE_NO,ESTIMATED_NO_OF_LEFT_FACTORS,ESTIMATED_NO_OF_RIGHT_FACTORS>, 
    NO_OF_LEFT_FACTORS, NO_OF_RIGHT_FACTORS, Chirality::right>::type
{
public:
   using FMC = FACTOR_MESSAGE_TRAIT;

   static constexpr INDEX leftFactorNumber = LEFT_FACTOR_NO;
   static constexpr INDEX rightFactorNumber = RIGHT_FACTOR_NO;

   static constexpr INDEX no_left_factors() { return NO_OF_LEFT_FACTORS; }
   static constexpr INDEX no_right_factors() { return NO_OF_RIGHT_FACTORS; }

   static constexpr INDEX estimated_no_left_factors() { return ESTIMATED_NO_OF_LEFT_FACTORS; }
   static constexpr INDEX estimated_no_right_factors() { return ESTIMATED_NO_OF_RIGHT_FACTORS; }

   using MessageContainerType = MessageContainer<MESSAGE_TYPE, LEFT_FACTOR_NO, RIGHT_FACTOR_NO, MPS, NO_OF_LEFT_FACTORS, NO_OF_RIGHT_FACTORS, FACTOR_MESSAGE_TRAIT, MESSAGE_NO, ESTIMATED_NO_OF_LEFT_FACTORS,ESTIMATED_NO_OF_RIGHT_FACTORS>;
   using MessageType = MESSAGE_TYPE;
   using next_left_message_container_type = next_msg_container_selector<MessageContainerType, NO_OF_LEFT_FACTORS, NO_OF_RIGHT_FACTORS, Chirality::left>;
   using next_right_message_container_type = next_msg_container_selector<MessageContainerType, NO_OF_LEFT_FACTORS, NO_OF_RIGHT_FACTORS, Chirality::left>;

   // FactorContainer
   using LeftFactorContainer = meta::at_c<typename FACTOR_MESSAGE_TRAIT::FactorList, leftFactorNumber>;
   using RightFactorContainer = meta::at_c<typename FACTOR_MESSAGE_TRAIT::FactorList, rightFactorNumber>;
   // Factor
   using LeftFactorType = typename LeftFactorContainer::FactorType;
   using RightFactorType = typename RightFactorContainer::FactorType;

   template<typename ...ARGS>
   MessageContainer(LeftFactorContainer* const l, RightFactorContainer* const r, ARGS... args) 
   : msg_op_(args...),
   leftFactor_(l),
   rightFactor_(r)
   {
       assert(l != nullptr && r != nullptr);
      //leftFactor_->template AddMessage<MessageDispatcher<MessageContainerType, LeftMessageFuncGetter>, MessageContainerType>(this);
      //rightFactor_->template AddMessage<MessageDispatcher<MessageContainerType, RightMessageFuncGetter>, MessageContainerType>(this);
   }

   /* seems not to work, as arguments are matched greedily???
   template<typename ...ARGS>
   MessageContainer(ARGS... args, LeftFactorContainer* const l, RightFactorContainer* const r) 
   : msg_op_(args...),
   leftFactor_(l),
   rightFactor_(r)
   {
      leftFactor_->template AddMessage<MessageDispatcher<MessageContainerType, LeftMessageFuncGetter>, MessageContainerType>(this);
      rightFactor_->template AddMessage<MessageDispatcher<MessageContainerType, RightMessageFuncGetter>, MessageContainerType>(this);
   }
   */

   MessageContainer(LeftFactorContainer* const l, RightFactorContainer* const r, MESSAGE_TYPE msg_op) 
       : msg_op_(msg_op),
       leftFactor_(l),
       rightFactor_(r)
   {
      //leftFactor_->template AddMessage<MessageDispatcher<MessageContainerType, LeftMessageFuncGetter>, MessageContainerType>(this);
      //rightFactor_->template AddMessage<MessageDispatcher<MessageContainerType, RightMessageFuncGetter>, MessageContainerType>(this);
      //int status;
      //std::cout << "msg holding type = " << abi::__cxa_demangle(typeid(*this).name(),0,0,&status) << "\n";
      //std::cout << FunctionExistence::IsAssignable<RightFactorContainer,REAL,INDEX>() << "\n";
      //std::cout << "msg holding type = " << abi::__cxa_demangle(typeid(msg_op_).name(),0,0,&status) << "\n";
      //std::cout << "left factor number = " << leftFactorNumber << "\n";
      //std::cout << "right factor number = " << rightFactorNumber << "\n";
      //std::cout << "left factor type = " << abi::__cxa_demangle(typeid(LeftFactorContainer).name(),0,0,&status) << "\n";
      //std::cout << "right factor type = " << abi::__cxa_demangle(typeid(RightFactorContainer).name(),0,0,&status) << "\n";
      // register messages in factors
   }


   ~MessageContainer() {
      static_assert(meta::unique<typename FACTOR_MESSAGE_TRAIT::MessageList>::size() == FACTOR_MESSAGE_TRAIT::MessageList::size(), 
            "Message list must have unique elements");
      //static_assert(MESSAGE_NO >= 0 && MESSAGE_NO < FACTOR_MESSAGE_TRAIT::MessageList::size(), "message number must be smaller than length of message list");
      static_assert(leftFactorNumber < FACTOR_MESSAGE_TRAIT::FactorList::size(), "left factor number out of bound");
      static_assert(rightFactorNumber < FACTOR_MESSAGE_TRAIT::FactorList::size(), "right factor number out of bound");
      // do zrobienia: put message constraint here, i.e. which methods MESSAGE_TYPE must minimally implement
   } 
   
   // overloaded new so that factor containers are allocated by global block allocator consecutively
   //void* operator new(std::size_t size)
   //{
   //   assert(size == sizeof(MessageContainerType));
   //   //INDEX s = size/sizeof(REAL);
   //   //if(size % sizeof(REAL) != 0) { s++; }
   //   //return (void*) global_real_block_allocator.allocate(s,1);
   //   return Allocator::get().allocate(1);
   //}

   //void operator delete(void* mem)
   //{
   //   Allocator::get().deallocate((MessageContainerType*) mem);
   //   //assert(false);
   //   //global_real_block_allocator.deallocate(mem,sizeof(FactorContainerType));
   //}

   using free_message_container_type = MessageContainer<MESSAGE_TYPE, LEFT_FACTOR_NO, RIGHT_FACTOR_NO, MPS, 0,0, empty_message_fmc<FMC>, MESSAGE_NO>;
   // return message container not embedded in message storage
   free_message_container_type free_message() const
   {
       free_message_container_type m(leftFactor_, rightFactor_, msg_op_);
       return m; 
   }

   // this view of the message container is given to left and right factor respectively when receiving or sending messages
   class test_zero_message_val {
   public:
      test_zero_message_val()
      {}

      test_zero_message_val& operator-=(const REAL x)
      {
          assert( std::abs(x) <= eps);
          return *this;
      }
      test_zero_message_val& operator+=(const REAL x)
      {
         assert(false);
         return *this;
      }
   };

   class test_zero_message : public MessageContainerType {
   public:
      test_zero_message_val operator[](const INDEX i) 
      {
         return test_zero_message_val();
      }

      template<typename ARRAY>
      test_zero_message& operator-=(const ARRAY& diff) {
        for(std::size_t i=0; i<diff.size(); ++i) {
            assert( std::abs(diff[i]) <= eps );
        }
        return *this;
      }
   };

   void send_message_to_left(const REAL omega = 1.0) 
   {
      send_message_to_left(rightFactor_->GetFactor(), omega);
   }
   void send_message_to_left(RightFactorType* r, const REAL omega)
   {
#ifndef NDEBUG // check whether double application of sending messages will result in zero update (for omega=1)
     test_send_message_to_left();
#endif

     msg_op_.send_message_to_left(*r, *static_cast<MessageContainerView<Chirality::right>*>(this), omega); 
   }

   void test_send_message_to_left()
   {
       RightFactorContainer right_factor_copy(*rightFactor_);
       // rewire rightFactor_, so that a message to left will reparametrize the copy
       RightFactorContainer* r = rightFactor_;
       rightFactor_ = &right_factor_copy;
       msg_op_.send_message_to_left(*right_factor_copy.GetFactor(), *static_cast<OneSideMessageContainerView<Chirality::right>*>(this), 1.0);
       // send again the message. It should be zero now.
       msg_op_.send_message_to_left(*right_factor_copy.GetFactor(), *static_cast<test_zero_message*>(this), 1.0); 
       rightFactor_ = r;
   }

#ifdef LP_MP_PARALLEL
   void send_message_to_left_synchronized(const REAL omega = 1.0) 
   {
      send_message_to_left_synchronized(rightFactor_->GetFactor(), omega);
   }
   void send_message_to_left_synchronized(RightFactorType* r, const REAL omega)
   {
     auto& mtx = GetLeftFactor()->mutex_;
     std::unique_lock<std::recursive_mutex> lck(mtx,std::defer_lock);
     if(lck.try_lock()) {
       msg_op_.send_message_to_left(*r, *static_cast<MessageContainerView<Chirality::right>*>(this), omega); 
     } else {
#ifndef NDEBUG
       if(debug()) {
         std::cout << "locking failed for " << GetLeftFactor() << "\n";
       }
#endif
     }
   }
#endif

   constexpr static bool can_compute_send_message_to_left_improvement()
   {
      return FunctionExistence::has_send_message_to_left_improvement<MessageType, REAL, LeftFactorType, RightFactorType>(); 
   }

   REAL send_message_to_left_improvement()
   {
       if constexpr(can_compute_send_message_to_left_improvement()) {
           return msg_op_.send_message_to_left_improvement(*leftFactor_->GetFactor(), *rightFactor_->GetFactor());
       } else {
           assert(false);
           return 0.0;
       }
   }

   void send_message_to_right(const REAL omega = 1.0) 
   {
      send_message_to_right(leftFactor_->GetFactor(), omega);
   }
   void send_message_to_right(LeftFactorType* l, const REAL omega)
   {
#ifndef NDEBUG // check whether double application of sending messages will result in zero update (for omega=1)
     test_send_message_to_right();
#endif

     msg_op_.send_message_to_right(*l, *static_cast<MessageContainerView<Chirality::left>*>(this), omega); 
   }

   void test_send_message_to_right()
   {
       LeftFactorContainer left_factor_copy(*leftFactor_);
       // rewire rightFactor_, so that a message to right will reparametrize the copy
       LeftFactorContainer* r = leftFactor_;
       leftFactor_ = &left_factor_copy;
       msg_op_.send_message_to_right(*left_factor_copy.GetFactor(), *static_cast<OneSideMessageContainerView<Chirality::left>*>(this), 1.0);
       // send again the message. It should be zero now.
       msg_op_.send_message_to_right(*left_factor_copy.GetFactor(), *static_cast<test_zero_message*>(this), 1.0); 
       leftFactor_ = r;
   }

#ifdef LP_MP_PARALLEL
   void send_message_to_right_synchronized(const REAL omega = 1.0) 
   {s
      send_message_to_right_synchronized(leftFactor_->GetFactor(), omega);
   }
   void send_message_to_right_synchronized(LeftFactorType* l, const REAL omega)
   {
     auto& mtx = GetRightFactor()->mutex_;
     std::unique_lock<std::recursive_mutex> lck(mtx,std::defer_lock);
     if(lck.try_lock()) {
       msg_op_.send_message_to_right(*l, *static_cast<MessageContainerView<Chirality::left>*>(this), omega); 
     } else {
#ifndef NDEBUG
       if(debug()) {
         std::cout << "locking failed for " << GetRightFactor() << "\n";
       }
#endif
     }
   }
#endif

   constexpr static bool can_compute_send_message_to_right_improvement()
   {
      return FunctionExistence::has_send_message_to_right_improvement<MessageType, REAL, LeftFactorType, RightFactorType>(); 
   }

   REAL send_message_to_right_improvement()
   {
       if constexpr(can_compute_send_message_to_right_improvement()) {
           return msg_op_.send_message_to_right_improvement(*leftFactor_->GetFactor(), *rightFactor_->GetFactor());
       } else {
           assert(false);
           return 0.0;
       }
   }

   constexpr static bool
   CanCallReceiveMessageFromRightContainer()
   { 
      return MPS == message_passing_schedule::left || MPS==message_passing_schedule::full;
      // obsolete
      return FunctionExistence::HasReceiveMessageFromRight<MessageType, void, 
      RightFactorType, MessageContainerType>(); 
   }
   void ReceiveMessageFromRightContainer()
   {
#ifndef NDEBUG
      const REAL before_left_lb = leftFactor_->LowerBound();
      const REAL before_right_lb = rightFactor_->LowerBound();
#endif

      send_message_to_left();

#ifndef NDEBUG
      const REAL after_left_lb = leftFactor_->LowerBound();
      const REAL after_right_lb = rightFactor_->LowerBound();
      assert(before_left_lb + before_right_lb <= after_left_lb + after_right_lb + eps); 
#endif
   }

#ifdef LP_MP_PARALLEL
   void ReceiveMessageFromRightContainerSynchronized()
   {
      send_message_to_left_synchronized();
   }
#endif

   constexpr static bool
   CanCallReceiveRestrictedMessageFromRightContainer()
   { 
      return FunctionExistence::HasReceiveRestrictedMessageFromRight<MessageType, void, 
      RightFactorType, MessageContainerType>(); // do zrobienia: signature is slighly different: MessageContainerType is not actually used
   }
   void ReceiveRestrictedMessageFromRightContainer()
   {
      rightFactor_->conditionally_init_primal(leftFactor_->primal_access_);
      msg_op_.ReceiveRestrictedMessageFromRight(*(rightFactor_->GetFactor()), *static_cast<OneSideMessageContainerView<Chirality::left>*>(this));
   }

   constexpr static bool 
   CanCallReceiveMessageFromLeftContainer()
   { 
      return MPS == message_passing_schedule::right || MPS == message_passing_schedule::full;
   }
   void ReceiveMessageFromLeftContainer()
   { 
#ifndef NDEBUG
      const REAL before_left_lb = leftFactor_->LowerBound();
      const REAL before_right_lb = rightFactor_->LowerBound();
#endif

      send_message_to_right();

#ifndef NDEBUG
      const REAL after_left_lb = leftFactor_->LowerBound();
      const REAL after_right_lb = rightFactor_->LowerBound();
      assert(before_left_lb + before_right_lb <= after_left_lb + after_right_lb + eps); 
#endif
   }

#ifdef LP_MP_PARALLEL
   void ReceiveMessageFromLeftContainerSynchronized()
   { 
      send_message_to_right_synchronized();
   }
#endif

   constexpr static bool
   CanCallReceiveRestrictedMessageFromLeftContainer()
   { 
      return FunctionExistence::HasReceiveRestrictedMessageFromLeft<MessageType, void, LeftFactorType, MessageContainerType>(); 
   }
   void ReceiveRestrictedMessageFromLeftContainer()
   {
      leftFactor_->conditionally_init_primal(rightFactor_->primal_access_);
      msg_op_.ReceiveRestrictedMessageFromLeft(*(leftFactor_->GetFactor()), *static_cast<OneSideMessageContainerView<Chirality::right>*>(this));
   }


   constexpr static bool 
   CanCallSendMessageToRightContainer()
   { 
      return MPS == message_passing_schedule::left || MPS == message_passing_schedule::full || MPS == message_passing_schedule::only_send;
   }

   void SendMessageToRightContainer(LeftFactorType* l, const REAL omega)
   {
      send_message_to_right(l, omega);
   }

#ifdef LP_MP_PARALLEL
   void SendMessageToRightContainerSynchronized(LeftFactorType* l, const REAL omega)
   {
      send_message_to_right_synchronized(l, omega);
   }
#endif

   constexpr static bool
   CanCallSendMessageToLeftContainer()
   { 
      return MPS == message_passing_schedule::right || MPS == message_passing_schedule::full || MPS == message_passing_schedule::only_send;
   }

   void SendMessageToLeftContainer(RightFactorType* r, const REAL omega)
   {
      send_message_to_left(r, omega);
   }

#ifdef LP_MP_PARALLEL
   void SendMessageToLeftContainerSynchronized(RightFactorType* r, const REAL omega)
   {
      send_message_to_left_synchronized(r, omega);
   }
#endif

   constexpr static bool CanCallSendMessagesToLeftContainer()
   {
      // possibly the below is to complicated. meta::find will be easier
      constexpr INDEX msg_array_number = RightFactorContainer::template FindMessageDispatcherTypeIndex<MessageDispatcher<MessageContainerType,RightMessageFuncGetter>>();
      using msg_container_type = meta::at_c<typename RightFactorContainer::msg_container_type_list, msg_array_number>;
      using MSG_ARRAY_ITERATOR = decltype(std::declval<msg_container_type>().begin());
      return FunctionExistence::HasSendMessagesToLeft<MessageType, void, RightFactorType, MSG_ARRAY_ITERATOR, MSG_ARRAY_ITERATOR, REAL>();
   }

   template<Chirality C> class MessageContainerView; // forward declaration. Put MessageIteratorView after definition of MessageContainerView

   template<Chirality CHIRALITY, typename MESSAGE_ITERATOR>
   class MessageIteratorView {
       public:
           MessageIteratorView(MESSAGE_ITERATOR it) : it_(it) {} 
           const MessageContainerView<CHIRALITY>& operator*() const {
               return *(static_cast<MessageContainerView<CHIRALITY>*>( &*it_ )); 
           }
           MessageContainerView<CHIRALITY>& operator*() {
               return *(static_cast<MessageContainerView<CHIRALITY>*>( &*it_ )); 
           }
           MessageIteratorView<CHIRALITY,MESSAGE_ITERATOR>& operator++() {
               ++it_;
               return *this;
           }
           bool operator==(const MessageIteratorView<CHIRALITY,MESSAGE_ITERATOR>& o) const {
               return it_ == o.it_; 
           }
           bool operator!=(const MessageIteratorView<CHIRALITY,MESSAGE_ITERATOR>& o) const {
               return it_ != o.it_; 
           }
       private:
           MESSAGE_ITERATOR it_;
   };

   template<Chirality CHIRALITY> class OneSideMessageContainerView; // forward declaration

   template<Chirality CHIRALITY, typename MESSAGE_ITERATOR>
   class one_sided_message_iterator_view {
       public:
           one_sided_message_iterator_view(MESSAGE_ITERATOR it) : it_(it) {} 
           const OneSideMessageContainerView<CHIRALITY>& operator*() const {
               return *(static_cast<OneSideMessageContainerView<CHIRALITY>*>( &*it_ )); 
           }
           OneSideMessageContainerView<CHIRALITY>& operator*() {
               return *(static_cast<OneSideMessageContainerView<CHIRALITY>*>( &*it_ )); 
           }
           one_sided_message_iterator_view<CHIRALITY,MESSAGE_ITERATOR>& operator++() {
               ++it_;
               return *this;
           }
           bool operator==(const one_sided_message_iterator_view<CHIRALITY,MESSAGE_ITERATOR>& o) const {
               return it_ == o.it_; 
           }
           bool operator!=(const one_sided_message_iterator_view<CHIRALITY,MESSAGE_ITERATOR>& o) const {
               return it_ != o.it_; 
           }
       private:
           MESSAGE_ITERATOR it_;
   };

   // for testing whether zero messages are computed
   template<typename MESSAGE_ITERATOR>
   class test_zero_message_iterator_view {
       public:
           test_zero_message_iterator_view(MESSAGE_ITERATOR it) : it_(it) {} 

           const test_zero_message& operator*() const {
               return *static_cast<test_zero_message*>(&*it_);
           }
           test_zero_message& operator*() {
               return *static_cast<test_zero_message*>(&*it_);
           }
           test_zero_message_iterator_view<MESSAGE_ITERATOR>& operator++() {
               ++it_;
               return *this;
           }
           bool operator==(const test_zero_message_iterator_view<MESSAGE_ITERATOR>& o) const {
               return it_ == o.it_; 
           }
           bool operator!=(const test_zero_message_iterator_view<MESSAGE_ITERATOR>& o) const {
               return it_ != o.it_; 
           }
       private:
           MESSAGE_ITERATOR it_;
   };

#ifdef LP_MP_PARALLEL
   template<Chirality CHIRALITY, typename MESSAGE_ITERATOR, typename LOCK_ITERATOR>
   struct MessageIteratorViewSynchronized {
     MessageIteratorViewSynchronized(MESSAGE_ITERATOR it, LOCK_ITERATOR lock_it) : it_(it), lock_it_(lock_it) {}
     MessageContainerView<CHIRALITY>& operator*() const {
       return *(static_cast<MessageContainerView<CHIRALITY>*>( *it_ )); 
     }
     MessageIteratorViewSynchronized<CHIRALITY,MESSAGE_ITERATOR,LOCK_ITERATOR>& operator++() {
       ++it_;
       ++lock_it_;
       while(*lock_it_ == false) { // this will always terminate: the lock_rec has one more entry than there are msgs and last entry is always true
         ++it_;
         ++lock_it_; 
       }
       return *this;
     }
     bool operator==(const MessageIteratorViewSynchronized<CHIRALITY,MESSAGE_ITERATOR,LOCK_ITERATOR>& o) const {
       return it_ == o.it_; 
     }
     bool operator!=(const MessageIteratorViewSynchronized<CHIRALITY,MESSAGE_ITERATOR,LOCK_ITERATOR>& o) const {
       return it_ != o.it_; 
     }
     private:
     MESSAGE_ITERATOR it_;
     LOCK_ITERATOR lock_it_;
   };
#endif

   template<typename IT, typename LOCK_ITERATOR>
   struct omega_iterator_with_lock {
     omega_iterator_with_lock(IT it, LOCK_ITERATOR lock_it) : it_(it), lock_it_(lock_it) {}
     omega_iterator_with_lock(const omega_iterator_with_lock& o) : it_(o.it_), lock_it_(o.lock_it_) {}
     omega_iterator_with_lock& operator++() {
       ++it_;
       ++lock_it_;
       while(*lock_it_ == false) { // this will always terminate: the lock_it_ has one more entry than there are elements pointed to by it_ and last entry is always true
         ++it_;
         ++lock_it_;
       }
       return *this;
     }
     omega_iterator_with_lock operator+(const INDEX s) {
        omega_iterator_with_lock o(*this);
        for(INDEX i=0; i<s; ++i) {
           ++o;
        }
        return o;
     }
     auto operator*() const { return *it_; }
     bool operator==(const omega_iterator_with_lock<IT,LOCK_ITERATOR>& o) const { return it_ == o.it_; }
     bool operator!=(const omega_iterator_with_lock<IT,LOCK_ITERATOR>& o) const { return it_ != o.it_; }

     private:
     IT it_;
     LOCK_ITERATOR lock_it_;
   };

   // rename to send_messages_to_left_container
   template<typename RIGHT_FACTOR, typename MSG_ITERATOR>
   static void SendMessagesToLeftContainer(const RIGHT_FACTOR& rightFactor, MSG_ITERATOR msgs_begin, MSG_ITERATOR msgs_end, const REAL omega) 
   {
#ifndef NDEBUG
       test_send_messages_to_left(msgs_begin, msgs_end);
#endif
      using MessageIteratorType = MessageIteratorView<Chirality::right, MSG_ITERATOR>;
      return MessageType::SendMessagesToLeft(rightFactor, MessageIteratorType(msgs_begin), MessageIteratorType(msgs_end), omega);
   }

   template<typename MSG_ITERATOR>
   static void test_send_messages_to_left(MSG_ITERATOR msg_begin, MSG_ITERATOR msg_end)
   {
       auto* right_factor = (*msg_begin).GetRightFactor();
       RightFactorContainer right_factor_copy(*right_factor);
       // rewire right factor_, so that SendMessagesToLeft right will reparametrize the copy
       for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
           (*msg_it).SetRightFactor(&right_factor_copy);
       }

       using message_iterator_type = one_sided_message_iterator_view<Chirality::right, MSG_ITERATOR>;
       MessageType::SendMessagesToLeft(*right_factor_copy.GetFactor(), message_iterator_type(msg_begin), message_iterator_type(msg_end), 1.0);
       // send again the message. It should be zero now.
       using zero_message_iterator_type = test_zero_message_iterator_view<MSG_ITERATOR>;
       MessageType::SendMessagesToLeft(*right_factor_copy.GetFactor(), zero_message_iterator_type(msg_begin), zero_message_iterator_type(msg_end), 1.0);

       for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
           (*msg_it).SetRightFactor(right_factor);
       }
   }

#ifdef LP_MP_PARALLEL
   template<typename RIGHT_FACTOR, typename MSG_ARRAY, typename ITERATOR, typename LOCK_VECTOR>
   static void SendMessagesToLeftContainerSynchronized_impl(const RIGHT_FACTOR& rightFactor, const MSG_ARRAY& msgs, ITERATOR omegaBegin, LOCK_VECTOR& lock_rec) 
   {
      // first lock as many adjacent factors as possible.
      auto lock_it = lock_rec.begin();
      for(auto it=msgs.begin(); it!=msgs.end(); ++it, ++lock_it) {
        auto& mtx = (*it)->GetLeftFactor()->mutex_;
        if(mtx.try_lock()) { // mark that factor was locked by this process
          *lock_it = true;
        } else {
          *lock_it = false; 
#ifndef NDEBUG
          if(debug()) {
            std::cout << "locking failed for " << (*it)->GetLeftFactor() << "\n";
          }
#endif
        }
      }
      assert(lock_it+1 == lock_rec.end());
      //std::fill(lock_rec.begin(), lock_rec.end(), true);

      using MessageIteratorType = MessageIteratorViewSynchronized<Chirality::right, decltype(msgs.begin()), decltype(lock_it)>;
      omega_iterator_with_lock<decltype(omegaBegin), decltype(lock_rec.cbegin())> omega_it(omegaBegin, lock_rec.cbegin()) ;
      MessageType::SendMessagesToLeft(rightFactor, MessageIteratorType(msgs.begin(), lock_rec.begin()), MessageIteratorType(msgs.end(), lock_rec.end()-1), omega_it);

      // unlock those factors which were locked above
      lock_it = lock_rec.begin();
      for(auto it=msgs.begin(); it!=msgs.end(); ++it, ++lock_it) {
        if(*lock_it) {
          (*it)->GetLeftFactor()->mutex_.unlock();
        }
      }
      assert(lock_it+1 == lock_rec.end());
   }
#endif
   template<typename RIGHT_FACTOR, typename MSG_ARRAY, typename ITERATOR>
   static void SendMessagesToLeftContainerSynchronized(const RIGHT_FACTOR& rightFactor, const MSG_ARRAY& msgs, ITERATOR omegaBegin) 
   {
     // record which factors were locked here
     if(msgs.size() <= 63) {
       std::array<bool,64> lock_rec;
       lock_rec[msgs.size()] = true;
       SendMessagesToLeftContainerSynchronized_impl(rightFactor, msgs, omegaBegin, lock_rec);
     } else {
       std::vector<bool> lock_rec(msgs.size()+1); // replace with own vector
       lock_rec[msgs.size()] = true;
       SendMessagesToLeftContainerSynchronized_impl(rightFactor, msgs, omegaBegin, lock_rec);
     }
   }

   template<typename MSG_ITERATOR>
   constexpr static bool can_compute_send_messages_to_left_improvement()
   {
      return FunctionExistence::has_send_messages_to_left_improvement<MessageType, REAL, RightFactorType, MSG_ITERATOR>(); 
   }

   template<typename MSG_ITERATOR>
   static REAL send_messages_to_left_improvement(MSG_ITERATOR msg_begin, MSG_ITERATOR msg_end) 
   {
       if constexpr(can_compute_send_messages_to_left_improvement<MSG_ITERATOR>()) {
           auto* r = (*msg_begin).GetRightFactor()->GetFactor();
           auto& msg_op = (*msg_begin).GetMessageOp();
           return msg_op.send_messages_to_left_improvement(*r, msg_begin, msg_end);
       } else {
           assert(false); // not implemented yet
       }
   }

   constexpr static bool CanCallSendMessagesToRightContainer()
   {
      // possibly the below is to complicated. meta::find will be easier
      constexpr INDEX msg_array_number = LeftFactorContainer::template FindMessageDispatcherTypeIndex<MessageDispatcher<MessageContainerType,LeftMessageFuncGetter>>();
      using msg_container_type = meta::at_c<typename LeftFactorContainer::msg_container_type_list, msg_array_number>;
      using MSG_ARRAY_ITERATOR = decltype(std::declval<msg_container_type>().begin());
      return FunctionExistence::HasSendMessagesToRight<MessageType, void, LeftFactorType, MSG_ARRAY_ITERATOR, MSG_ARRAY_ITERATOR, REAL>();
   }

   // rename send_messages_to_right_container
   template<typename LEFT_FACTOR, typename MSG_ITERATOR>
   static void SendMessagesToRightContainer(const LEFT_FACTOR& leftFactor, MSG_ITERATOR msgs_begin, MSG_ITERATOR msgs_end, const REAL omega) 
   {
#ifndef NDEBUG
       test_send_messages_to_right(msgs_begin, msgs_end);
#endif
      using MessageIteratorType = MessageIteratorView<Chirality::left, MSG_ITERATOR>;
      return MessageType::SendMessagesToRight(leftFactor, MessageIteratorType(msgs_begin), MessageIteratorType(msgs_end), omega);
   }

   template<typename MSG_ITERATOR>
   static void test_send_messages_to_right(MSG_ITERATOR msg_begin, MSG_ITERATOR msg_end)
   {
       auto* left_factor = (*msg_begin).GetLeftFactor();
       LeftFactorContainer left_factor_copy(*left_factor);
       // rewire left factor_, so that SendMessagesToRight right will reparametrize the copy
       for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
           (*msg_it).SetLeftFactor(&left_factor_copy);
       }

       using message_iterator_type = one_sided_message_iterator_view<Chirality::left, MSG_ITERATOR>;
       MessageType::SendMessagesToRight(*left_factor_copy.GetFactor(), message_iterator_type(msg_begin), message_iterator_type(msg_end), 1.0);
       // send again the message. It should be zero now.
       using zero_message_iterator_type = test_zero_message_iterator_view<MSG_ITERATOR>;
       MessageType::SendMessagesToRight(*left_factor_copy.GetFactor(), zero_message_iterator_type(msg_begin), zero_message_iterator_type(msg_end), 1.0);

       for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
           (*msg_it).SetLeftFactor(left_factor);
       }
   }

#ifdef LP_MP_PARALLEL
   template<typename LEFT_FACTOR, typename MSG_ARRAY, typename ITERATOR, typename LOCK_VECTOR>
   static void SendMessagesToRightContainerSynchronized_impl(const LEFT_FACTOR& leftFactor, const MSG_ARRAY& msgs, ITERATOR omegaBegin, LOCK_VECTOR& lock_rec) 
   {
      // first lock as many adjacent factors as possible.
      auto lock_it = lock_rec.begin();
      for(auto it=msgs.begin(); it!=msgs.end(); ++it, ++lock_it) {
        if((*it)->GetRightFactor()->mutex_.try_lock()) { // mark that factor was locked by this process
          *lock_it = true;
        } else {
          *lock_it = false; 
#ifndef NDEBUG
          if(debug()) {
            std::cout << "locking failed for " << (*it)->GetRightFactor() << "\n";
          }
#endif
        }
      }
      assert(lock_it+1 == lock_rec.end());
      //std::fill(lock_rec.begin(), lock_rec.end(), true);

      using MessageIteratorType = MessageIteratorViewSynchronized<Chirality::left, decltype(msgs.begin()), decltype(lock_it)>;
      omega_iterator_with_lock<decltype(omegaBegin), decltype(lock_rec.cbegin())> omega_it(omegaBegin, lock_rec.cbegin()) ;
      MessageType::SendMessagesToRight(leftFactor, MessageIteratorType(msgs.begin(), lock_rec.begin()), MessageIteratorType(msgs.end(), lock_rec.end()-1), omega_it);

      // unlock those factors which were locked above
      lock_it = lock_rec.begin();
      for(auto it=msgs.begin(); it!=msgs.end(); ++it, ++lock_it) {
        if(*lock_it) {
          (*it)->GetRightFactor()->mutex_.unlock();
        }
      }
      assert(lock_it+1 == lock_rec.end());
   }
   template<typename LEFT_FACTOR, typename MSG_ARRAY, typename ITERATOR>
   static void SendMessagesToRightContainerSynchronized(const LEFT_FACTOR& leftFactor, const MSG_ARRAY& msgs, ITERATOR omegaBegin) 
   {
     // record which factors were locked here
     if(msgs.size() <= 63) {
       std::array<bool,64> lock_rec;
       lock_rec[msgs.size()] = true;
       SendMessagesToRightContainerSynchronized_impl(leftFactor, msgs, omegaBegin, lock_rec);
     } else {
       std::vector<bool> lock_rec(msgs.size()+1); // replace with own vector
       lock_rec[msgs.size()] = true;
       SendMessagesToRightContainerSynchronized_impl(leftFactor, msgs, omegaBegin, lock_rec);
     }
   }
#endif

   template<typename MSG_ITERATOR>
   constexpr static bool can_compute_send_messages_to_right_improvement()
   {
      return FunctionExistence::has_send_messages_to_right_improvement<MessageType, REAL, LeftFactorType, MSG_ITERATOR>(); 
   }

   template<typename MSG_ITERATOR>
   static REAL send_messages_to_right_improvement(MSG_ITERATOR msg_begin, MSG_ITERATOR msg_end) 
   {
       if constexpr(can_compute_send_messages_to_right_improvement<MSG_ITERATOR>()) {
           assert(msg_begin != msg_end);
           auto* l = (*msg_begin).GetLeftFactor()->GetFactor();
           auto& msg_op = (*msg_begin).GetMessageOp();
           return msg_op.send_messages_to_right_improvement(*l, msg_begin, msg_end);
       } else {
           assert(false); // not implemented yet
       }
   }

   constexpr static bool
   CanComputeRightFromLeftPrimal()
   {
      return CanComputeRightFromLeftPrimalWithoutReturn() || CanComputeRightFromLeftPrimalWithReturn();
   }
   constexpr static bool
   CanComputeLeftFromRightPrimal()
   {
      return CanComputeLeftFromRightPrimalWithoutReturn() || CanComputeLeftFromRightPrimalWithReturn();
   } 

   constexpr static bool
   CanComputeRightFromLeftPrimalWithoutReturn()
   {
      return FunctionExistence::HasComputeRightFromLeftPrimal<MessageType, void, LeftFactorType, RightFactorType>();
   }
   constexpr static bool
   CanComputeLeftFromRightPrimalWithoutReturn()
   {
      return FunctionExistence::HasComputeLeftFromRightPrimal<MessageType, void, LeftFactorType, RightFactorType>();
   }

   constexpr static bool
   CanComputeRightFromLeftPrimalWithReturn()
   {
      return FunctionExistence::HasComputeRightFromLeftPrimal<MessageType, bool, LeftFactorType, RightFactorType>();
   }
   constexpr static bool
   CanComputeLeftFromRightPrimalWithReturn()
   {
      return FunctionExistence::HasComputeLeftFromRightPrimal<MessageType, bool, LeftFactorType, RightFactorType>();
   }

   void ComputeRightFromLeftPrimal() 
   {
      rightFactor_->conditionally_init_primal(leftFactor_->primal_access_);

      if constexpr (CanComputeRightFromLeftPrimalWithoutReturn()) {
        msg_op_.ComputeRightFromLeftPrimal(*leftFactor_->GetFactor(), *rightFactor_->GetFactor());
        rightFactor_->PropagatePrimal();
        rightFactor_->propagate_primal_through_messages(); 
      } else if constexpr (MessageContainerType::CanComputeRightFromLeftPrimalWithReturn()) {
        const bool changed = msg_op_.ComputeRightFromLeftPrimal(*leftFactor_->GetFactor(), *rightFactor_->GetFactor());
        if(changed) {
          rightFactor_->PropagatePrimal();
          rightFactor_->propagate_primal_through_messages();
        }
      }

      /*
      static_if<CanComputeRightFromLeftPrimalWithoutReturn()>([&](auto f) {
        f(msg_op_).ComputeRightFromLeftPrimal(*leftFactor_->GetFactor(), *rightFactor_->GetFactor());
        rightFactor_->PropagatePrimal();
        rightFactor_->propagate_primal_through_messages();
      }).else_([&](auto) {
         static_if<MessageContainerType::CanComputeRightFromLeftPrimalWithReturn()>([&](auto f) {
               const bool changed = f(msg_op_).ComputeRightFromLeftPrimal(*leftFactor_->GetFactor(), *rightFactor_->GetFactor());
               if(changed) {
                  rightFactor_->PropagatePrimal();
                  rightFactor_->propagate_primal_through_messages();
               }
         });
      });
      */
   }

   void ComputeLeftFromRightPrimal()
   {
      leftFactor_->conditionally_init_primal(rightFactor_->primal_access_);
      if constexpr(CanComputeLeftFromRightPrimalWithoutReturn()) {
        msg_op_.ComputeLeftFromRightPrimal(*leftFactor_->GetFactor(), *rightFactor_->GetFactor());
        leftFactor_->PropagatePrimal();
        leftFactor_->propagate_primal_through_messages();
      } else if constexpr(CanComputeLeftFromRightPrimalWithReturn()) {
          const bool changed = msg_op_.ComputeLeftFromRightPrimal(*leftFactor_->GetFactor(), *rightFactor_->GetFactor());
          if(changed) {
              leftFactor_->PropagatePrimal();
              leftFactor_->propagate_primal_through_messages();
          }
      }
   }

   constexpr static bool
   CanCheckPrimalConsistency()
   {
      return FunctionExistence::HasCheckPrimalConsistency<MessageType,bool,
          typename LeftFactorContainer::FactorType,
          typename RightFactorContainer::FactorType>();
   }

   bool CheckPrimalConsistency() const
   { 
      bool ret;
      if constexpr(CanCheckPrimalConsistency()) {
          return msg_op_.CheckPrimalConsistency(*leftFactor_->GetFactor(), *rightFactor_->GetFactor());
      } else {
          return true;
      }
   }

   // do zrobienia: not needed anymore
   // do zrobienia: this does not capture write back functions not returning REAL&
   constexpr static bool IsAssignableLeft() {
      return FunctionExistence::IsAssignable<LeftFactorType, REAL, INDEX>();
   }
   constexpr static bool IsAssignableRight() {
      return FunctionExistence::IsAssignable<RightFactorType, REAL, INDEX>();
   }

   template<typename ARRAY, bool IsAssignable = IsAssignableLeft()>
   constexpr static bool CanBatchRepamLeft()
   {
      return FunctionExistence::HasRepamLeft<MessageType,void,LeftFactorType,ARRAY>();
   }
   template<typename ARRAY, bool IsAssignable = IsAssignableLeft()>
   //typename std::enable_if<CanBatchRepamLeft<ARRAY>() == true && IsAssignable == true>::type
   void
   RepamLeft(const ARRAY& m)
   { 
      //assert(false); // no -+ distinguishing
      if constexpr(CanBatchRepamLeft<ARRAY>()) {
            msg_op_.RepamLeft(*(leftFactor_->GetFactor()), m);
      } else {
         const auto s = m.size();
         for(INDEX i=0; i<s; ++i) {
            msg_op_.RepamLeft(*(leftFactor_->GetFactor()), m[i], i);
         }
      }
   }
   /*
   template<typename ARRAY, bool IsAssignable = IsAssignableLeft()>
   typename std::enable_if<CanBatchRepamLeft<ARRAY>() == false && IsAssignable == true>::type
   RepamLeft(const ARRAY& m)
   { 
      //assert(false); // no -+ distinguishing
      for(INDEX i=0; i<m.size(); ++i) {
         msg_op_.RepamLeft(*(leftFactor_->GetFactor()), m[i], i);
      }
   }
   template<typename ARRAY, bool IsAssignable = IsAssignableLeft()>
   typename std::enable_if<IsAssignable == false>::type
   RepamLeft(const ARRAY& m)
   {
   assert(false);
   }
   */

   template<bool IsAssignable = IsAssignableLeft()>
   //typename std::enable_if<IsAssignable == true>::type
   void
   RepamLeft(const REAL diff, const INDEX dim) {
      msg_op_.RepamLeft(*(leftFactor_->GetFactor()), diff, dim); // note: in right, we reparametrize by +diff, here by -diff
   }
   /*
   template<bool IsAssignable = IsAssignableLeft()>
   typename std::enable_if<IsAssignable == false>::type
   RepamLeft(const REAL diff, const INDEX dim)
   {
   assert(false);
   }
   */

   template<typename ARRAY>
   constexpr static bool CanBatchRepamRight()
   {
      // do zrobienia: replace Container by actual factor
      return FunctionExistence::HasRepamRight<MessageType,void,RightFactorType,ARRAY>();
   }
   template<typename ARRAY, bool IsAssignable = IsAssignableRight()>
   //typename std::enable_if<CanBatchRepamRight<ARRAY>() == true && IsAssignable == true>::type
   void
   RepamRight(const ARRAY& m)
   { 
      //assert(false); // no -+ distinguishing
      if constexpr(CanBatchRepamRight<ARRAY>()) {
            msg_op_.RepamRight(*(rightFactor_->GetFactor()), m);
      } else {
          const auto s = m.size();
         for(INDEX i=0; i<s; ++i) {
            msg_op_.RepamRight(*(rightFactor_->GetFactor()), m[i], i);
         }
      }
   }
   /*
   template<typename ARRAY, bool IsAssignable = IsAssignableRight()>
   typename std::enable_if<CanBatchRepamRight<ARRAY>() == false && IsAssignable == true>::type
   RepamRight(const ARRAY& m)
   {
      //assert(false); // no -+ distinguishing
      for(INDEX i=0; i<m.size(); ++i) {
         msg_op_.RepamRight(*(rightFactor_->GetFactor()), m[i], i);
      }
   }
   template<typename ARRAY, bool IsAssignable = IsAssignableRight()>
   typename std::enable_if<IsAssignable == false>::type
   RepamRight(const ARRAY& m)
   {
   assert(false);
   }
   */

   template<bool IsAssignable = IsAssignableRight()>
   //typename std::enable_if<IsAssignable == true>::type
   void
   RepamRight(const REAL diff, const INDEX dim) {
      msg_op_.RepamRight(*(rightFactor_->GetFactor()), diff, dim);
   }
   /*
   template<bool IsAssignable = IsAssignableRight()>
   typename std::enable_if<IsAssignable == false>::type
   RepamRight(const REAL diff, const INDEX dim)
   {
   assert(false);
   }
   */

   // do zrobienia: better name?
   //REAL GetLeftMessage(const INDEX i) const { return msg_op_.GetLeftMessage(i,*this); }
   //REAL GetRightMessage(const INDEX i) const { return msg_op_.GetRightMessage(i,*this);  }

   FactorTypeAdapter* GetLeftFactorTypeAdapter() const { return leftFactor_; }
   FactorTypeAdapter* GetRightFactorTypeAdapter() const { return rightFactor_; }
   // do zrobienia: Rename Get{Left|Right}FactorContainer
   auto* GetLeftFactor() const { return leftFactor_; }
   auto* GetRightFactor() const { return rightFactor_; }

   void SetLeftFactor(FactorTypeAdapter* l) 
   {
      assert(dynamic_cast<LeftFactorContainer*>(l));
      leftFactor_ = static_cast<LeftFactorContainer*>(l); 
   }
   void SetRightFactor(FactorTypeAdapter* r)
   {
      assert(dynamic_cast<RightFactorContainer*>(r));
      rightFactor_ = static_cast<RightFactorContainer*>(r); 
   }

   // class for storing a callback upon new assignment of message: update left and right factors
   // convention is as follows: original message is for right factor. Inverted message is for left one
   template<Chirality CHIRALITY>
   class MsgVal {
   public:
      MsgVal(MessageContainerType* msg, const INDEX dim) : 
         msg_(msg), 
         dim_(dim)
      {}
      // do zrobienia: do not support this operation! Goal is to not hold messages anymore, except for implicitly held reparametrizations.
      /*
      MsgVal& operator=(const REAL x) __attribute__ ((always_inline))
      {
         assert(false);
         const REAL diff = x - msg_->operator[](dim_);
         // set new message
         static_cast<typename MessageContainerType::MessageStorageType*>(msg_)->operator[](dim_) = x;
         // propagate difference to left and right factor
         msg_->RepamLeft( diff, dim_);
         msg_->RepamRight( diff, dim_);
         return *this;
      }
      */
      MsgVal& operator-=(const REAL x) __attribute__ ((always_inline))
      {
         if(CHIRALITY == Chirality::right) { // message is computed by right factor
            msg_->RepamLeft( +x, dim_);
            msg_->RepamRight(-x, dim_);
         } else if (CHIRALITY == Chirality::left) { // message is computed by left factor
            msg_->RepamLeft(  -x, dim_);
            msg_->RepamRight( +x, dim_);
            //msg_->RepamLeft( +x, dim_);
            //msg_->RepamRight( +x, dim_);
         } else {
            assert(false);
         }
         return *this;
      }
      MsgVal& operator+=(const REAL x) __attribute__ ((always_inline))
      {
         assert(false);
         if(CHIRALITY == Chirality::right) { // message is computed by right factor
            msg_->RepamLeft( x, dim_);
            msg_->RepamRight( x, dim_);
         } else if(CHIRALITY == Chirality::left) { // message is computed by left factor
            msg_->RepamLeft( x, dim_);
            msg_->RepamRight( x, dim_);
            //msg_->RepamLeft( -x, dim_);
            //msg_->RepamRight( -x, dim_);
         } else {
            assert(false);
         }
         return *this;
      }
      // do zrobienia: this value should never be used. Remove function
      //operator REAL() const __attribute__ ((always_inline)) { return static_cast<typename MessageContainerType::MessageStorageType*>(msg_)->operator[](dim_); }
   private:
      MessageContainerType* const msg_;
      const INDEX dim_;
   };

   // this view of the message container is given to left and right factor respectively when receiving or sending messages
   template<Chirality CHIRALITY>
   class MessageContainerView : public MessageContainerType {
   public:
      //using MessageContainerType;
      MsgVal<CHIRALITY> operator[](const INDEX i) 
      {
         return MsgVal<CHIRALITY>(this,i);
      }

      template<typename ARRAY>
      MessageContainerType& operator-=(const ARRAY& diff) {
        // note: order of below operations is important: When the message is e.g. just the potential, we must reparametrize the other side first!
        if(CHIRALITY == Chirality::right) {
          RepamLeft(diff);
          RepamRight(-diff);
        } else if(CHIRALITY == Chirality::left) {
          RepamRight(diff);
          RepamLeft(-diff);
        } else {
          assert(false);
        }
        return *this;
      }

   };

   // for primal computation: record message change only in one side and into a special array
   template<Chirality CHIRALITY>
   class OneSideMsgVal
   {
   public:
      OneSideMsgVal(MessageContainerType* msg, const INDEX dim) : 
         msg_(msg), 
         dim_(dim)
      {}

      OneSideMsgVal& operator-=(const REAL x) __attribute__ ((always_inline))
      {
         if(CHIRALITY == Chirality::right) { // message is received by right factor
            msg_->RepamRight(-x, dim_);
         } else if (CHIRALITY == Chirality::left) { // message is received by left factor
            msg_->RepamLeft(-x, dim_);
         } else {
            assert(false);
         }
         return *this;
      }

      OneSideMsgVal& operator+=(const REAL x) __attribute__ ((always_inline))
      {
         assert(false);
         return *this;
      }

   private:
      MessageContainerType* const msg_;
      const INDEX dim_;
   };

   // this view is given to receive restricted message operations. 
   // Reparametrization is recorded only on one side
   template<Chirality CHIRALITY>
   class OneSideMessageContainerView : public MessageContainerType {
   public:
      //using MessageContainerType;
      OneSideMsgVal<CHIRALITY> operator[](const INDEX i) 
      {
         return OneSideMsgVal<CHIRALITY>(this,i);
      }

      template<typename ARRAY>
      MessageContainerType& operator-=(const ARRAY& diff) {
        if(CHIRALITY == Chirality::right) {
          RepamRight(-diff);
        } else if(CHIRALITY == Chirality::left) {
          RepamLeft(-diff);
        } else {
          assert(false);
        }
        return *this;
      }

   };


   // there must be four different implementations of msg updating with SIMD: 
   // (i) If parallel reparametrization is not supported by either left and right factor
   // If (ii) left or (iii) right factor supports reparametrization but not both
   // If (iv) both left and right factor support reparametrization


   template<typename ARRAY>
   MessageContainerType& operator-=(const ARRAY& diff) {
      assert(false); // update to left right -+
      RepamLeft(-diff);
      RepamRight(-diff);
      return *this;
   }

   template<typename ARRAY>
   MessageContainerType& operator+=(const ARRAY& diff) {
      assert(false); // update to left right -+
      RepamLeft(diff);
      RepamRight(diff);
      return *this;
   }

   // possibly not the best choice: Sometimes msg_op_ needs access to this class
   const MessageType& GetMessageOp() const { return msg_op_; }

   MessageType& GetMessageOp() { return msg_op_; }

   // for weight computations these functions are necessary
   static constexpr bool sends_message_to_left_constexpr()
   {
      return MPS == message_passing_schedule::right || MPS == message_passing_schedule::full || MPS == message_passing_schedule::only_send;
   }
   static constexpr bool sends_message_to_right_constexpr() 
   {
      return MPS == message_passing_schedule::left || MPS == message_passing_schedule::full || MPS == message_passing_schedule::only_send;
   }
   static constexpr bool receives_message_from_left_constexpr() 
   {
     return MPS == message_passing_schedule::right || MPS == message_passing_schedule::full;
   }
   static constexpr bool receives_message_from_right_constexpr() 
   {
     return MPS == message_passing_schedule::left || MPS == message_passing_schedule::full;
   }

   bool SendsMessageToLeft() const { return sends_message_to_left_constexpr(); }
   bool SendsMessageToRight() const { return sends_message_to_right_constexpr(); }
   bool ReceivesMessageFromLeft() const { return receives_message_from_left_constexpr(); }
   bool ReceivesMessageFromRight() const { return receives_message_from_right_constexpr(); }

   /*
   constexpr static bool
   can_construct_constraints()
   {
      return has_construct_sat_clauses<MessageType,void,sat_solver&, LeftFactorType, RightFactorType, sat_var, sat_var>();
   }

   void construct_sat_clauses(sat_solver& s, sat_var left_var, sat_var right_var) 
   {
      static_if<can_construct_sat_clauses()>([&](auto f) {
            f(msg_op_).construct_sat_clauses(s, *leftFactor_->GetFactor(), *rightFactor_->GetFactor(), left_var, right_var);
            });
      if(!can_construct_sat_clauses()) {
         assert(false);
      }
   }
   */ 
   
   // for traversing a tree
   virtual void send_message_up(Chirality c) 
   {
      if(c == Chirality::right) { // right factor is top one
         leftFactor_->GetFactor()->init_primal();
         //if constexpr(LeftFactorContainer::CanMaximizePotentalAndComputePrimal()) {
         //    leftFactor_->MaximizePotentialAndComputePrimal();
         //}
        static_if<LeftFactorContainer::CanMaximizePotentialAndComputePrimal()>([&](auto f) {
               f(leftFactor_)->MaximizePotentialAndComputePrimal(); 
         });
         this->send_message_to_right();
         //static_if<CanCallReceiveMessageFromLeftContainer()>([&](auto f) {
         //      f(this)->ReceiveMessageFromLeftContainer();
         //}).else_([&](auto) {
         //      static_if<MessageContainerType::CanCallSendMessageToRightContainer()>([&](auto f) {
         //               f(this)->SendMessageToRightContainer(leftFactor_->GetFactor(),1.0);
         //      }).else_([](auto) {
         //         assert(false); // possibly try to call SendMessagesToRightContainer with exactly one message
         //      });
         //});
      } else {
         rightFactor_->GetFactor()->init_primal();
         //if constexpr(RightFactorContainer::CanMaximizePotentalAndComputePrimal()) {
         //    rightFactor_->MaximizePotentialAndComputePrimal();
         //}
         static_if<RightFactorContainer::CanMaximizePotentialAndComputePrimal()>([&](auto f) {
               f(rightFactor_)->MaximizePotentialAndComputePrimal(); 
         });
         this->send_message_to_left();
         //static_if<CanCallReceiveMessageFromRightContainer()>([&](auto f) {
         //      f(this)->ReceiveMessageFromRightContainer();
         //}).else_([&](auto) {
         //      static_if<MessageContainerType::CanCallSendMessageToLeftContainer()>([&](auto f) {
         //               f(this)->SendMessageToLeftContainer(rightFactor_->GetFactor(),1.0);
         //      }).else_([](auto) {
         //         assert(false); // possibly try to call SendMessagesToRightContainer with exactly one message
         //      });
         //});
      }
   }

   
   virtual void track_solution_down(Chirality c) 
   {
      // we can assume that upper factor has already (partially) computed primal.
      // we check whether we can receive restricted messages from upper and compute primal in lower. If yes, we receive restricted message, compute primal in lower factor and propagate it back to upper factor.
      if(c == Chirality::right) {
         //leftFactor_->init_primal(); // initialization is already done in upward pass
         static_if<MessageContainerType::CanComputeLeftFromRightPrimal()>([&](auto f) {
               f(msg_op_).ComputeLeftFromRightPrimal(*leftFactor_->GetFactor(), *rightFactor_->GetFactor());
         });
         static_if<LeftFactorContainer::CanMaximizePotentialAndComputePrimal()>([&](auto f2) {
               f2(leftFactor_)->MaximizePotentialAndComputePrimal(); 
         });
      } else {
         assert(c == Chirality::left);
         //rightFactor_->init_primal();
         static_if<MessageContainerType::CanComputeRightFromLeftPrimal()>([&](auto f) {
               f(msg_op_).ComputeRightFromLeftPrimal(*leftFactor_->GetFactor(), *rightFactor_->GetFactor());
         });
         static_if<RightFactorContainer::CanMaximizePotentialAndComputePrimal()>([&](auto f2) {
               f2(rightFactor_)->MaximizePotentialAndComputePrimal(); 
         });
      }
      return;


      // if this is not possible, we propagate primal labeling of upper to lower
      if(c == Chirality::right) { // right factor is upper
         static_if<LeftFactorContainer::CanMaximizePotentialAndComputePrimal() && CanCallReceiveRestrictedMessageFromRightContainer()>([&](auto f) {
                  // receive restricted messages 
                  serialization_archive ar(leftFactor_->GetFactor(), leftFactor_->GetFactor()+1, [](auto& f, auto& ar) { f.serialize_dual(ar); });
                  save_archive s_ar(ar);
                  leftFactor_->GetFactor()->serialize_dual( s_ar );

                  f(this)->ReceiveRestrictedMessageFromRightContainer();

                  // compute primal in lower
                  leftFactor_->MaximizePotentialAndComputePrimal();

                  // restore dual reparametrization to before restricted messages were sent.
                  load_archive l_ar(ar);
                  leftFactor_->GetFactor()->serialize_dual( l_ar );

                  // propagate back to upper
                  f(this)->ComputeRightFromLeftPrimal(); 

         }).else_([&](auto) {
            static_if<MessageContainerType::CanComputeLeftFromRightPrimal()>([&](auto f) {
                  f(this)->ComputeLeftFromRightPrimal();
                  static_if<LeftFactorContainer::CanMaximizePotentialAndComputePrimal()>([&](auto f2) {
                        f2(leftFactor_)->MaximizePotentialAndComputePrimal(); 
                  });
            }).else_([&](auto) {
               assert(false);
            });
         });

      } else if(c == Chirality::left) { // left factor is upper
         static_if<RightFactorContainer::CanMaximizePotentialAndComputePrimal() && CanCallReceiveRestrictedMessageFromLeftContainer()>([&](auto f) {
                  std::stringstream ss;
                  // receive restricted messages 
                  serialization_archive ar(rightFactor_->GetFactor(), rightFactor_->GetFactor()+1, [](auto& f, auto& ar) { f.serialize_dual(ar); });
                  save_archive s_ar(ar);
                  rightFactor_->GetFactor()->serialize_dual( s_ar );

                  f(this)->ReceiveRestrictedMessageFromLeftContainer();

                  // compute primal in lower
                  rightFactor_->MaximizePotentialAndComputePrimal();

                  // restore dual reparametrization to before restricted messages were sent.
                  load_archive l_ar(ar);
                  rightFactor_->GetFactor()->serialize_dual( l_ar );

                  // propagate back to upper
                  f(this)->ComputeLeftFromRightPrimal(); 

         }).else_([&](auto) {
            static_if<MessageContainerType::CanComputeRightFromLeftPrimal()>([&](auto f) {
                  f(this)->ComputeRightFromLeftPrimal();
                  static_if<RightFactorContainer::CanMaximizePotentialAndComputePrimal()>([&](auto f2) {
                        f2(rightFactor_)->MaximizePotentialAndComputePrimal(); 
                  });
            }).else_([&](auto) {
               assert(false);
            });
         });

      } else {
         assert(false);
      }
   }

   // construct constraints
   template<typename SOLVER>
   void construct_constraints_impl(SOLVER& s, const typename DD_ILP::variable_counters& left_variable_counters, const typename DD_ILP::variable_counters& right_variable_counters)
   {
      auto current_variable_counters = s.get_variable_counters();

      auto left_vars = leftFactor_->GetFactor()->export_variables();
      s.set_variable_counters(left_variable_counters);
      auto left_external_vars = std::apply([this,&s](auto... x){ return std::make_tuple(this->leftFactor_->load_external_variables(s, x)...); }, left_vars);

      auto right_vars = rightFactor_->GetFactor()->export_variables();
      s.set_variable_counters(right_variable_counters);
      auto right_external_vars = std::apply([this,&s](auto... x){ return std::make_tuple(this->rightFactor_->load_external_variables(s, x)...); }, right_vars);

      auto t = std::tuple_cat(std::tie(*leftFactor_->GetFactor()), left_external_vars, std::tie(*rightFactor_->GetFactor()), right_external_vars);
      auto construct_constraints_fun = [this,&s](auto... x) { this->msg_op_.construct_constraints(s, x...); };
      std::apply(construct_constraints_fun, t);

      s.set_variable_counters(current_variable_counters);
   }

   virtual void construct_constraints(
       DD_ILP::external_solver_interface<DD_ILP::sat_solver>& s, 
       const typename DD_ILP::variable_counters& left_variable_counters,
       const typename DD_ILP::variable_counters& right_variable_counters 
       ) 
   { construct_constraints_impl(s, left_variable_counters, right_variable_counters); }
   virtual void construct_constraints(
       DD_ILP::external_solver_interface<DD_ILP::problem_export>& s, 
       const typename DD_ILP::variable_counters& left_variable_counters,
       const typename DD_ILP::variable_counters& right_variable_counters 
       ) 
   { construct_constraints_impl(s, left_variable_counters, right_variable_counters); }
#ifdef DD_ILP_WITH_GUROBI
   virtual void construct_constraints(
       DD_ILP::external_solver_interface<DD_ILP::gurobi_interface>& s, 
       const typename DD_ILP::variable_counters& left_variable_counters,
       const typename DD_ILP::variable_counters& right_variable_counters 
       ) 
   { construct_constraints_impl(s, left_variable_counters, right_variable_counters); }
#endif


protected:
   MessageType msg_op_; // possibly inherit privately from MessageType to apply empty base optimization when applicable
   LeftFactorContainer* leftFactor_;
   RightFactorContainer* rightFactor_;
};


// message container storage options for holding messages in factor containers.
// explicitly hold N message containers
template<typename MESSAGE_CONTAINER_TYPE, std::size_t N>
class message_container_storage_array {
public:
    using message_type = typename MESSAGE_CONTAINER_TYPE::MessageType;
    using iterator = MESSAGE_CONTAINER_TYPE*;

    static constexpr std::size_t message_storage_byte_size() 
    {
        // encountered bug: We can possibly overwrite the message_container_storage_array with factor 5*sizeof(void*), but static_assert in destructor does not raise this!
        return N*(sizeof(message_type) + 6*sizeof(void*)); 

        /*
        constexpr auto no_left_factors = MESSAGE_CONTAINER_TYPE::no_left_factors();
        constexpr auto no_right_factors = MESSAGE_CONTAINER_TYPE::no_right_factors();
        if constexpr(no_left_factors == 0 && no_right_factors != 0) {
            return N*(sizeof(message_type) + 4*sizeof(void*));
        } else if constexpr(no_left_factors != 0 && no_right_factors == 0) {
            return N*(sizeof(message_type) + 4*sizeof(void*));
        } else {
            return N*(sizeof(message_type) + 4*sizeof(void*)); // for left and right factor + vtbl
        }
        */
        //return sizeof(MESSAGE_CONTAINER_TYPE) * N; 
    }
    //static constexpr std::size_t message_storage_byte_size_ = sizeof(MESSAGE_CONTAINER_TYPE)*N;

    message_container_storage_array()
    {
        std::fill(storage_.begin(), storage_.end(), 0);
        assert(occupied() == 0);
    }

    ~message_container_storage_array()
    {
        static_assert(message_storage_byte_size() >= N*sizeof(MESSAGE_CONTAINER_TYPE));
        static_assert(N > 0);
        //std::string("size for message storage = ") + std::to_string(message_storage_byte_size()/N) + " must be equal to size of message container = " + std::to_string(sizeof(MESSAGE_CONTAINER_TYPE));
    }

    static constexpr std::size_t capacity() { return N; }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename ...ARGS>
    MESSAGE_CONTAINER_TYPE* push_back(LEFT_FACTOR* l, RIGHT_FACTOR* r, ARGS... args) {
        const std::size_t i = occupied();
        assert(i<N);
        auto* ptr = &storage_[i*sizeof(MESSAGE_CONTAINER_TYPE)];
        new(ptr) MESSAGE_CONTAINER_TYPE(l, r, args...); // placement new
        return reinterpret_cast<MESSAGE_CONTAINER_TYPE*>(ptr);
    }

    void set_next_message(MESSAGE_CONTAINER_TYPE* m) {}

    const MESSAGE_CONTAINER_TYPE& operator[](const std::size_t i) const {
        assert(i<capacity());
        char* ptr = &storage_[i*sizeof(MESSAGE_CONTAINER_TYPE)];
        return *reinterpret_cast<MESSAGE_CONTAINER_TYPE*>(ptr);
    }

    MESSAGE_CONTAINER_TYPE& operator[](const std::size_t i) {
        assert(i<capacity());
        auto* ptr = &storage_[i*sizeof(MESSAGE_CONTAINER_TYPE)];
        return *reinterpret_cast<MESSAGE_CONTAINER_TYPE*>(ptr);
    }


    std::size_t occupied() const 
    {
        // check storage until encountering a slot with all zeros
        auto cap = capacity();
        for(std::size_t i=0; i<cap; ++i) {
            bool occupied = false;
            for(std::size_t j=0; j<sizeof(MESSAGE_CONTAINER_TYPE); ++j) {
                if(storage_[i*sizeof(MESSAGE_CONTAINER_TYPE) + j] != 0) {
                    occupied = true;
                }
            }
            if(!occupied) { return i; }
        }
        return capacity(); 
    }

protected:
    std::array<unsigned char, message_storage_byte_size()> storage_;
};

// hold list of chunks of N messages each
template<typename MESSAGE_CONTAINER_TYPE, std::size_t N>
class variable_message_container_storage {

    class variable_message_container_storage_chunk : public message_container_storage_array<MESSAGE_CONTAINER_TYPE,N> {
        friend class variable_message_container_storage<MESSAGE_CONTAINER_TYPE,N>;
        public:
        using storage_type = variable_message_container_storage_chunk;
        using message_type = typename MESSAGE_CONTAINER_TYPE::MessageType;

        //add new operator that uses block allocator here

        variable_message_container_storage_chunk()
            : message_container_storage_array<MESSAGE_CONTAINER_TYPE,N>(),
            next_(nullptr)
        {}

        ~variable_message_container_storage_chunk()
        {
            static_assert(N > 0);
        }

        // overloaded new so that factor containers are allocated by global block allocator consecutively
        void* operator new(std::size_t size)
        {
            assert(size == sizeof(storage_type));
            return Allocator::get().allocate(1);
        }

        void operator delete(void* mem)
        {
            Allocator::get().deallocate((storage_type*) mem);
        }


        private:
        variable_message_container_storage_chunk* next_;

        struct Allocator {
            using type = MemoryPool<storage_type,4096*(sizeof(storage_type)+sizeof(void*))>; 
            static type& get() {
                static type allocator;
                return allocator;
            }
        };
    };

public:
    using chunk_type = variable_message_container_storage_chunk;

    variable_message_container_storage() :
        end_(&first_chunk_,0),
        size_(0)
    {} 

    chunk_type* last_chunk() const
    {
        auto* cur_chunk = &first_chunk_;

        while(cur_chunk->next_ != nullptr) {
            cur_chunk = cur_chunk->next_;
        } 
        return cur_chunk; 
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename ...ARGS>
    MESSAGE_CONTAINER_TYPE* push_back(LEFT_FACTOR* l, RIGHT_FACTOR* r, ARGS... args) 
    {
        auto* cur_chunk = last_chunk();
#ifndef NDEBUG
        const auto no_messages = size();
#endif
        if(cur_chunk->occupied() == N) {
            auto* new_chunk = new chunk_type();
            cur_chunk->next_ = new_chunk;
            cur_chunk = new_chunk;
        }

        auto* m = cur_chunk->push_back(l, r, args...);

        end_ = iterator(cur_chunk, cur_chunk->occupied());
        ++size_;
#ifndef NDEBUG
        assert(no_messages + 1 == size());
#endif 
        return m;
    }

    void set_next_message(MESSAGE_CONTAINER_TYPE* m) {}

    std::size_t size() const
    {
        return size_;
        // can be optimized using end_
        std::size_t c = 0;
        auto* cur_chunk = &first_chunk_;
        while(cur_chunk->next_ != nullptr) {
            cur_chunk = cur_chunk->next_;
            c += N;
        }
        c += cur_chunk->occupied(); 
        return c;
    }

   class iterator {
      public:
         iterator(chunk_type* c) : chunk_(c), i_(0) {}
         iterator(chunk_type* c, std::size_t i) : chunk_(c), i_(i) {
             if(i_ == N) {
                 chunk_ = chunk_->next_;
                 i_ = 0;
             }
         }
         iterator operator++() {
             assert(chunk_ != nullptr);
             ++i_;
             if(i_ == N) {
                 chunk_ = chunk_->next_;
                 i_ = 0;
             }
             return *this;
         }
         MESSAGE_CONTAINER_TYPE& operator*() { return (*chunk_)[i_]; } 
         const MESSAGE_CONTAINER_TYPE& operator*() const { return (*chunk_)[i_]; } 
         bool operator==(const iterator& o) const { return chunk_ == o.chunk_ && i_ == o.i_; }
         bool operator!=(const iterator& o) const { return !(*this == o); }
         chunk_type* chunk() const { return chunk_; }
         std::size_t chunk_size() const { return i_; }
      private:
         chunk_type* chunk_;
         std::size_t i_;
   };

   const iterator begin() const {
      return iterator(&first_chunk_);
   }
   const iterator end() const {
       assert(end_ == iterator(last_chunk(), last_chunk()->occupied()));
       return end_;
      //auto* c = last_chunk();
      //return iterator(c, c->occupied());
   }

   iterator begin() {
      return iterator(&first_chunk_);
   }
   iterator end() {
       assert(end_ == iterator(last_chunk(), last_chunk()->occupied()));
       return end_;
      //auto* c = last_chunk();
      //return iterator(c, c->occupied());
   }

    bool empty() const { return begin() == end(); }

private:
    mutable chunk_type first_chunk_;
    iterator end_;
    std::size_t size_;
};


// hold up to N messages
template<typename MESSAGE_CONTAINER_TYPE, std::size_t N>
class up_to_message_container_storage : public message_container_storage_array<MESSAGE_CONTAINER_TYPE,N> {
public:
    up_to_message_container_storage()
        : message_container_storage_array<MESSAGE_CONTAINER_TYPE,N>(),
        end_(this->begin())
    {
        std::cout << &end_ << "\n";
        assert(end_ != nullptr);
    }

    ~up_to_message_container_storage()
    {
        static_assert(N > 0);
    } 

    std::size_t size() const 
    { 
        const std::size_t n = std::distance(begin(), end()); 
        assert(n == this->occupied());
        return n;
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename ...ARGS>
    MESSAGE_CONTAINER_TYPE* push_back(LEFT_FACTOR* l, RIGHT_FACTOR* r, ARGS... args) 
    {
        assert(size() < this->capacity());
#ifndef NDEBUG
        const auto n = size();
#endif
        new(end_) MESSAGE_CONTAINER_TYPE(l, r, args...); // placement new
        ++end_;
        assert(n+1 == size());
        return end_;
    }

    const MESSAGE_CONTAINER_TYPE* begin() const { return reinterpret_cast<const MESSAGE_CONTAINER_TYPE*>(&this->storage_[0]); } 
    const MESSAGE_CONTAINER_TYPE* end() const { assert(end_ != nullptr); assert(this->occupied() == std::distance(begin(), reinterpret_cast<const MESSAGE_CONTAINER_TYPE*>(end_))); return end_; }

    MESSAGE_CONTAINER_TYPE* begin() { return reinterpret_cast<MESSAGE_CONTAINER_TYPE*>(&this->storage_[0]); } 
    MESSAGE_CONTAINER_TYPE* end() { assert(end_ != nullptr); assert(this->occupied() == std::distance(begin(), end_)); return end_; }

    bool empty() const { assert(end_ != nullptr); assert(this->occupied() == std::distance(begin(), end())); return begin() == end(); }
private:
    MESSAGE_CONTAINER_TYPE* end_;
};

// hold exactly N messages
template<typename MESSAGE_CONTAINER_TYPE, std::size_t N>
class fixed_message_container_storage : public message_container_storage_array<MESSAGE_CONTAINER_TYPE,N> {
public:
    fixed_message_container_storage()
        : message_container_storage_array<MESSAGE_CONTAINER_TYPE,N>()
    {}

    ~fixed_message_container_storage()
    {
        static_assert(N > 0);
    } 

    std::size_t size() const { assert(this->capacity() == this->occupied()); return this->capacity(); }

    const auto* begin() const { assert(this->capacity() == this->occupied()); return reinterpret_cast<const MESSAGE_CONTAINER_TYPE*>(&this->storage_[0]); } 
    const auto* end() const { assert(this->capacity() == this->occupied()); return reinterpret_cast<const MESSAGE_CONTAINER_TYPE*>(&this->storage_[0]) + this->capacity(); }

    auto* begin() { assert(this->capacity() == this->occupied()); return reinterpret_cast<MESSAGE_CONTAINER_TYPE*>(&this->storage_[0]); } 
    auto* end() { assert(this->capacity() == this->occupied()); return reinterpret_cast<MESSAGE_CONTAINER_TYPE*>(&this->storage_[0]) + this->capacity(); }

    static constexpr bool empty() { return false; }
};

template<typename MESSAGE_CONTAINER_TYPE, Chirality CHIRALITY>
class pointer_to_message_container_storage {
public:

    std::size_t size() const {
        std::size_t s=0;
        auto* p = ptr;
        while(p != nullptr) {
            ++s;
            if(CHIRALITY == Chirality::left) {
                p = p->next_left_msg();
            } else {
                assert(CHIRALITY == Chirality::right);
                p = p->next_right_msg();
            }
        }
        return s;
    }

    template<typename LEFT_FACTOR, typename RIGHT_FACTOR, typename ...ARGS>
    MESSAGE_CONTAINER_TYPE* push_back(LEFT_FACTOR* l, RIGHT_FACTOR* r, ARGS... args)
    {
        return nullptr;
    }

    void set_next_message(MESSAGE_CONTAINER_TYPE* m)
    {
        assert(m != nullptr);
        if(CHIRALITY == Chirality::left) {
            m->set_next_left_msg(ptr);
            ptr = m;
        } else {
            assert(CHIRALITY == Chirality::right);
            m->set_next_right_msg(ptr);
            ptr = m; 
        }
    }


    class iterator {
      public:
         iterator(MESSAGE_CONTAINER_TYPE* m) : m_(m) {}
         iterator operator++() {
             // next_msg is only called when one side has fixed number of messages and the other a variable number.
             if(MESSAGE_CONTAINER_TYPE::no_left_factors() == 0 && MESSAGE_CONTAINER_TYPE::no_right_factors() != 0) {
                 m_ = m_->next_left_msg();
             } else if(MESSAGE_CONTAINER_TYPE::no_right_factors() == 0 && MESSAGE_CONTAINER_TYPE::no_left_factors() != 0) {
                 m_ = m_->next_right_msg();
             } else {
                 assert(false);
             }
             return *this;
         }

         MESSAGE_CONTAINER_TYPE& operator*() const { return *m_; } 
         bool operator==(const iterator& o) const { return m_ == o.m_; }
         bool operator!=(const iterator& o) const { return m_ != o.m_; }
      private:
         MESSAGE_CONTAINER_TYPE* m_;
   };

    auto begin() const { return iterator(ptr); }
    auto end() const { return iterator(nullptr); }

    bool empty() const 
    { 
        if(ptr == nullptr) { assert(size() == 0); }
        if(ptr != nullptr) { assert(size() > 0); }
        return ptr == nullptr; 
    }
private:
    MESSAGE_CONTAINER_TYPE* ptr;
};

// N=0 means variable number of messages, > 0 means compile time fixed number of messages and <0 means at most compile time number of messages
// see config.hxx for shortcuts
template<typename MESSAGE_CONTAINER_TYPE, SIGNED_INDEX NO_LEFT_FACTORS, SIGNED_INDEX NO_RIGHT_FACTORS, INDEX ESTIMATED_NO_OF_LEFT_FACTORS, INDEX ESTIMATED_NO_OF_RIGHT_FACTORS, Chirality CHIRALITY>
struct message_container_selector {
    // the following cases may arise:
   // exactly one side has fixed number of messages, hold it in fixed_message_storage in factor of corresponding side.  The other side holds pointer to first message and messages holds next pointer.
   // both sides have variable number of messages: Hold two copies in of message container in variable_message_container
   // both sides have fixed number of messages: not implemented yet.

   // number of entries in variable_message_container_storage chunks
    static constexpr std::size_t variable_message_container_storage_size = 4;

    class empty {};

   using type = 
   std::conditional_t<(CHIRALITY == Chirality::left && NO_LEFT_FACTORS < 0 && NO_RIGHT_FACTORS < 0), up_to_message_container_storage<MESSAGE_CONTAINER_TYPE, -NO_LEFT_FACTORS>,
   std::conditional_t<(CHIRALITY == Chirality::left && NO_LEFT_FACTORS == 0 && NO_RIGHT_FACTORS < 0), pointer_to_message_container_storage<MESSAGE_CONTAINER_TYPE,CHIRALITY>,
   std::conditional_t<(CHIRALITY == Chirality::left && NO_LEFT_FACTORS > 0 && NO_RIGHT_FACTORS < 0), fixed_message_container_storage<MESSAGE_CONTAINER_TYPE, NO_LEFT_FACTORS>,

   std::conditional_t<(CHIRALITY == Chirality::left && NO_LEFT_FACTORS < 0 && NO_RIGHT_FACTORS == 0), up_to_message_container_storage<MESSAGE_CONTAINER_TYPE, -NO_LEFT_FACTORS>,
   std::conditional_t<(CHIRALITY == Chirality::left && NO_LEFT_FACTORS == 0 && NO_RIGHT_FACTORS == 0), variable_message_container_storage<MESSAGE_CONTAINER_TYPE, ESTIMATED_NO_OF_LEFT_FACTORS>,
   std::conditional_t<(CHIRALITY == Chirality::left && NO_LEFT_FACTORS > 0 && NO_RIGHT_FACTORS == 0), fixed_message_container_storage<MESSAGE_CONTAINER_TYPE, NO_LEFT_FACTORS>,

   std::conditional_t<(CHIRALITY == Chirality::left && NO_LEFT_FACTORS < 0 && NO_RIGHT_FACTORS > 0), up_to_message_container_storage<MESSAGE_CONTAINER_TYPE, -NO_LEFT_FACTORS>,
   std::conditional_t<(CHIRALITY == Chirality::left && NO_LEFT_FACTORS == 0 && NO_RIGHT_FACTORS > 0), pointer_to_message_container_storage<MESSAGE_CONTAINER_TYPE,CHIRALITY>,
   std::conditional_t<(CHIRALITY == Chirality::left && NO_LEFT_FACTORS > 0 && NO_RIGHT_FACTORS > 0), fixed_message_container_storage<MESSAGE_CONTAINER_TYPE, NO_LEFT_FACTORS>,
   

   std::conditional_t<(CHIRALITY == Chirality::right && NO_RIGHT_FACTORS < 0 && NO_LEFT_FACTORS < 0), up_to_message_container_storage<MESSAGE_CONTAINER_TYPE, -NO_RIGHT_FACTORS>,
   std::conditional_t<(CHIRALITY == Chirality::right && NO_RIGHT_FACTORS == 0 && NO_LEFT_FACTORS < 0), pointer_to_message_container_storage<MESSAGE_CONTAINER_TYPE,CHIRALITY>,
   std::conditional_t<(CHIRALITY == Chirality::right && NO_RIGHT_FACTORS > 0 && NO_LEFT_FACTORS < 0), fixed_message_container_storage<MESSAGE_CONTAINER_TYPE, NO_RIGHT_FACTORS>,

   std::conditional_t<(CHIRALITY == Chirality::right && NO_RIGHT_FACTORS < 0 && NO_LEFT_FACTORS == 0), up_to_message_container_storage<MESSAGE_CONTAINER_TYPE, -NO_RIGHT_FACTORS>,
   std::conditional_t<(CHIRALITY == Chirality::right && NO_RIGHT_FACTORS == 0 && NO_LEFT_FACTORS == 0), variable_message_container_storage<MESSAGE_CONTAINER_TYPE, ESTIMATED_NO_OF_RIGHT_FACTORS>,
   std::conditional_t<(CHIRALITY == Chirality::right && NO_RIGHT_FACTORS > 0 && NO_LEFT_FACTORS == 0), fixed_message_container_storage<MESSAGE_CONTAINER_TYPE, NO_RIGHT_FACTORS>,

   std::conditional_t<(CHIRALITY == Chirality::right && NO_RIGHT_FACTORS < 0 && NO_LEFT_FACTORS > 0), up_to_message_container_storage<MESSAGE_CONTAINER_TYPE, -NO_RIGHT_FACTORS>,
   std::conditional_t<(CHIRALITY == Chirality::right && NO_RIGHT_FACTORS == 0 && NO_LEFT_FACTORS > 0), pointer_to_message_container_storage<MESSAGE_CONTAINER_TYPE,CHIRALITY>,
   std::conditional_t<(CHIRALITY == Chirality::right && NO_RIGHT_FACTORS > 0 && NO_LEFT_FACTORS > 0), fixed_message_container_storage<MESSAGE_CONTAINER_TYPE, NO_RIGHT_FACTORS>,
   empty
   >>> >>> >>>  >>> >>> >>>;

   ~message_container_selector()
   {
       static_assert(!std::is_same_v<type, empty>, "");
       static_assert(ESTIMATED_NO_OF_LEFT_FACTORS > 0);
       static_assert(ESTIMATED_NO_OF_RIGHT_FACTORS > 0);
   }
};



// container class for factors. Here we hold the factor, all connected messages, reparametrization storage and perform reparametrization and coordination for sending and receiving messages.
// derives from REPAM_STORAGE_TYPE to mixin a class for storing the reparametrized potential
// implements the interface from FactorTypeAdapter for access from LP_MP
// if COMPUTE_PRIMAL_SOLUTION is true, MaximizePotential is expected to return either an integer of type INDEX or a std::vector<INDEX>
// if WRITE_PRIMAL_SOLUTION is false, WritePrimal will not output anything
// do zrobienia: introduce enum classes for COMPUTE_PRIMAL_SOLUTION and WRITE_PRIMAL_SOLUTION
template<typename FACTOR_TYPE, 
         class FACTOR_MESSAGE_TRAIT,
         INDEX FACTOR_NO,
         bool COMPUTE_PRIMAL_SOLUTION = false> 
class FactorContainer : public FactorTypeAdapter
{
public:
   using FactorContainerType = FactorContainer<FACTOR_TYPE, FACTOR_MESSAGE_TRAIT, FACTOR_NO, COMPUTE_PRIMAL_SOLUTION>;
   using FactorType = FACTOR_TYPE;
   using FMC = FACTOR_MESSAGE_TRAIT;

   // do zrobienia: templatize cosntructor to allow for more general initialization of reparametrization storage and factor
   template<typename ...ARGS>
   FactorContainer(ARGS... args) : factor_(args...) {}

   FactorContainer(const FactorType&& factor) : factor_(std::move(factor)) {
      //INDEX status;
      //std::cout << "msg_ type= "  << abi::__cxa_demangle(typeid(msg_).name(),0,0,&status) << "\n";
      //std::cout << "dispatcher list = "  << abi::__cxa_demangle(typeid(MESSAGE_DISPATCHER_TYPELIST).name(),0,0,&status) << "\n";
      //std::cout << "msg_ type= "  << abi::__cxa_demangle(typeid(msg_).name(),0,0,&status) << "\n";
      //std::cout << "left message list = " << abi::__cxa_demangle(typeid(left_message_list_).name(),0,0,&status) << "\n";
      //std::cout << "left message list = " << abi::__cxa_demangle(typeid(left_message_list_1).name(),0,0,&status) << "\n";
   
   }
   FactorContainer(const FactorType& factor) : factor_(factor) 
   {}
   ~FactorContainer() { 
      static_assert(meta::unique<MESSAGE_DISPATCHER_TYPELIST>::size() == MESSAGE_DISPATCHER_TYPELIST::size(), 
            "Message dispatcher typelist must have unique elements");
      static_assert(FACTOR_NO >= 0 && FACTOR_NO < FACTOR_MESSAGE_TRAIT::FactorList::size(), "factor number must be smaller than length of factor list");
   }

   // overloaded new so that factor containers are allocated by global block allocator consecutively
   void* operator new(std::size_t size)
   {
      assert(size == sizeof(FactorContainerType));
      //INDEX s = size/sizeof(REAL);
      //if(size % sizeof(REAL) != 0) { s++; }
      //return (void*) global_real_block_allocator.allocate(s,1);
      return Allocator::get().allocate(1);
   }
   void operator delete(void* mem)
   {
      Allocator::get().deallocate((FactorContainerType*) mem);
      //assert(false);
      //global_real_block_allocator.deallocate((double*)mem,sizeof(FactorContainerType)/sizeof(REAL)+1);
   }

   using empty_message_storage_factor_container = FactorContainer<FACTOR_TYPE, empty_message_fmc<FMC>, FACTOR_NO, COMPUTE_PRIMAL_SOLUTION>;
   empty_message_storage_factor_container* no_message_factor_clone()
   {
       auto* c = new empty_message_storage_factor_container(factor_); 
       return c;
   }

   virtual FactorTypeAdapter* clone() const final
   {
      auto* c = new FactorContainer(factor_);
      return c;
   }

   template<typename MESSAGE_CONTAINER_TYPE, Chirality CHIRALITY, typename ADJACENT_FACTOR, typename... ARGS>
   auto add_message(ADJACENT_FACTOR* a_f, ARGS... args)
   {
       if constexpr(CHIRALITY == Chirality::left) {
           constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<MessageDispatcher<MESSAGE_CONTAINER_TYPE,LeftMessageFuncGetter>>();
           return std::get<n>(msg_).push_back(this, a_f, args...); 
       } else {
           static_assert(CHIRALITY == Chirality::right);
           constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<MessageDispatcher<MESSAGE_CONTAINER_TYPE,RightMessageFuncGetter>>();
           return std::get<n>(msg_).push_back(a_f, this, args...); 
       }
   }

   template<typename MESSAGE_CONTAINER_TYPE>
   void set_left_msg(MESSAGE_CONTAINER_TYPE* m)
   {
       using message_dispatcher = MessageDispatcher<MESSAGE_CONTAINER_TYPE, LeftMessageFuncGetter>;
       constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<message_dispatcher>();
       std::get<n>(msg_).set_next_message(m); 
   }

   template<typename MESSAGE_CONTAINER_TYPE>
   void set_right_msg(MESSAGE_CONTAINER_TYPE* m)
   {
       using message_dispatcher = MessageDispatcher<MESSAGE_CONTAINER_TYPE, RightMessageFuncGetter>;
       constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<message_dispatcher>();
       std::get<n>(msg_).set_next_message(m); 
   }

   //template<typename MESSAGE_DISPATCHER_TYPE, typename MESSAGE_TYPE> 
   //void AddMessage(MESSAGE_TYPE* m) { 
   //   constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
   //   static_assert( n < meta::size<MESSAGE_DISPATCHER_TYPELIST>(), "message dispatcher not supported");
   //   static_assert( n < std::tuple_size<decltype(msg_)>(), "message dispatcher not supported");
   //   //INDEX status;
   //   //std::cout << "msg dispatcher list =\n" << abi::__cxa_demangle(typeid(MESSAGE_DISPATCHER_TYPELIST).name(),0,0,&status) << "\n";
   //   //std::cout << "dispatcher  type =\n" << abi::__cxa_demangle(typeid(MESSAGE_DISPATCHER_TYPE).name(),0,0,&status) << "\n";
   //   //std::cout << " number = " << n << "\n" ;
   //   //std::cout << "message type = " << abi::__cxa_demangle(typeid(MESSAGE_TYPE).name(),0,0,&status) << "\n";
   //
   //      std::get<n>(msg_).push_back(m);
   //}

   void update_factor_uniform(const REAL leave_weight) final
   {
       receive_messages();
       MaximizePotential();
       send_messages(leave_weight);
   }
   void UpdateFactor(const weight_slice omega, const receive_slice receive_mask) final
   {
      ReceiveMessages(receive_mask);
      MaximizePotential();
      SendMessages(omega);
   }

   void update_factor_adaptive(const weight_slice omega, const receive_slice receive_mask) final
   {
      ReceiveMessages(receive_mask);
      MaximizePotential();
      send_messages_with_adaptive_weights(omega); 
   }

   void update_factor_residual(const weight_slice omega, const receive_slice receive_mask) final
   {
      assert(*std::min_element(omega.begin(), omega.end()) >= 0.0);
      assert(*std::max_element(omega.begin(), omega.end()) <= 1.0+eps);
      assert(std::distance(omega.begin(), omega.end()) == no_send_messages());
      assert(receive_mask.size() == no_receive_messages());
      ReceiveMessages(receive_mask);
      MaximizePotential();
      send_messages_residual(omega); // other message passing type shall be called "shared"
   }

#ifdef LP_MP_PARALLEL
   void UpdateFactorSynchronized(const weight_slice& omega) final
   {
      assert(*std::min_element(omega.begin(), omega.end()) >= 0.0);
      assert(std::accumulate(omega.begin(), omega.end(), 0.0) <= 1.0 + eps);
      assert(std::distance(omega.begin(), omega.end()) == no_send_messages());
      std::lock_guard<std::recursive_mutex> lock(mutex_); // only here do we wait for the mutex. In all other places try_lock is allowed only
      ReceiveMessagesSynchronized(omega);
      MaximizePotential();
      SendMessagesSynchronized(omega);
   }

   void UpdateFactorPrimalSynchronized(const weight_slice& omega, const INDEX iteration) final
   {
     //std::cout << "not implemented\n";
     //assert(false);
   }
#endif

   // do zrobienia: possibly also check if method present
   constexpr static bool
   CanComputePrimal()
   {
      return COMPUTE_PRIMAL_SOLUTION;
   }

   constexpr static bool
   CanMaximizePotentialAndComputePrimal()
   {
      return FunctionExistence::HasMaximizePotentialAndComputePrimal<FactorType,void>();
   }

   constexpr static bool
   CanPropagatePrimal()
   {
      return FunctionExistence::HasPropagatePrimal<FactorType,void>();
   }

   void PropagatePrimal() 
   {
      if constexpr(CanPropagatePrimal()) {
          factor_.PropagatePrimal();
      }
   }

   constexpr static bool
   CanMaximizePotential()
   {
      return FunctionExistence::HasMaximizePotential<FactorType,void>();
   }


   void UpdateFactorPrimal(const weight_slice& omega, const receive_slice& receive_mask, INDEX primal_access) final
   {
#ifdef LP_MP_PARALLEL
     std::lock_guard<std::recursive_mutex> lock(mutex_); // only here do we wait for the mutex. In all other places try_lock is allowed only
#endif
      assert(primal_access > 0); // otherwise primal is not initialized in first iteration
      conditionally_init_primal(primal_access);
      if(CanComputePrimal()) { // do zrobienia: for now
         primal_access_ = primal_access;
         if(false && CanReceiveRestrictedMessage() && ReceivesRestrictedMessage()) {

            serialization_archive ar(GetFactor(), GetFactor()+1, [](auto& f, auto& ar) { f.serialize_dual(ar); });
            save_archive s_ar(ar);
            factor_.serialize_dual( s_ar );

            // now we change the dual information
            // first we compute restricted incoming messages, on which to compute the primal
            ReceiveRestrictedMessages();

            // now we compute primal w.r.t. the changed dual information!
            MaximizePotentialAndComputePrimal();

            // restore dual reparametrization to before restricted messages were sent.
            load_archive l_ar(ar);
            factor_.serialize_dual( l_ar );

            ReceiveMessages(receive_mask);
            MaximizePotential();
            SendMessages(omega);
         } else {
            ReceiveMessages(receive_mask);
            MaximizePotentialAndComputePrimal();
            SendMessages(omega);
         }
         // now propagate primal to adjacent factors
         propagate_primal_through_messages();
      } else {
         ReceiveMessages(receive_mask);
         MaximizePotential();
         SendMessages(omega);
      }  
   }

   void MaximizePotential()
   {
       if constexpr(CanMaximizePotential()) {
           factor_.MaximizePotential();
       }
   }

   virtual void MaximizePotentialAndComputePrimal() final
   {
       if constexpr(CanMaximizePotentialAndComputePrimal()) {
           factor_.MaximizePotentialAndComputePrimal();
       } else {
           assert(false);
       }
   }

   // do zrobienia: rename PropagatePrimalThroughMessages
   virtual void propagate_primal_through_messages() final
   {
      meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [this](auto l) {
            if constexpr(l.CanComputePrimalThroughMessage()) {
                  constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
                  auto msg_begin = std::get<n>(msg_).begin();
                  auto msg_end = std::get<n>(msg_).end();
                  for(auto it = msg_begin; it != msg_end; ++it) {
                     l.ComputePrimalThroughMessage(*it);
                  }
            }
      });
   }

   virtual bool check_primal_consistency() final
   {
       // to do: possibly only go through messages where current factor is on left side. Otherwise, messages are checked twice
       bool consistent = true;
       meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [this, &consistent](auto l) {
              constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
              if(consistent) {
                auto msg_begin = std::get<n>(msg_).begin();
                auto msg_end = std::get<n>(msg_).end();
                for(auto it = msg_begin; it != msg_end; ++it) {
                const bool message_consistent = (*it).CheckPrimalConsistency();
                if(!message_consistent) {
                    consistent = false;
                    break;
                }
               }
              }
      });
      return consistent;
   }

   void receive_messages()
   {
      meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [this](auto l) {
            constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
            if constexpr(l.receives_message_from_adjacent_factor()) {
                auto msg_begin = std::get<n>(msg_).begin();
                auto msg_end = std::get<n>(msg_).end();
                for(auto it = msg_begin; it != msg_end; ++it) {
                    l.ReceiveMessage(*it);
                }
            }
      });
   }


   template<typename WEIGHT_VEC>
   void ReceiveMessages(const WEIGHT_VEC& receive_mask) 
   {
      assert(std::distance(receive_mask.begin(), receive_mask.end()) == no_receive_messages()); 
      assert(receive_mask.size() == 0 || *std::max_element(receive_mask.begin(), receive_mask.end()) <= 1);
      assert(receive_mask.size() == 0 || *std::min_element(receive_mask.begin(), receive_mask.end()) >= 0);

      assert(receive_mask.size() == no_receive_messages());
      auto receive_it = receive_mask.begin();

      meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [this,&receive_it](auto l) {
            constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
            if constexpr(l.receives_message_from_adjacent_factor()) {
                auto msg_begin = std::get<n>(msg_).begin();
                auto msg_end = std::get<n>(msg_).end();
                  for(auto it = msg_begin; it != msg_end; ++it, ++receive_it) {
                     if(*receive_it) {
#ifndef NDEBUG
                        const REAL before_lb = LowerBound() + l.get_adjacent_factor(*it)->LowerBound();
#endif
                        l.ReceiveMessage(*it);
#ifndef NDEBUG
                        const REAL after_lb = LowerBound() + l.get_adjacent_factor(*it)->LowerBound();
                        assert(before_lb <= after_lb + eps);
#endif
                     }
                  }
            }
            
      });

   }

#ifdef LP_MP_PARALLEL
   template<typename WEIGHT_VEC>
   void ReceiveMessagesSynchronized(const WEIGHT_VEC& omega) 
   {
       assert(false); // introduce receive mask
      // note: currently all messages are received, even if not needed. Change this again.
      auto omegaIt = omega.begin();
      meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [this,&omegaIt](auto l) {
            constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
            if constexpr(l.receives_message_from_adjacent_factor() {
                  for(auto it = std::get<n>(msg_).begin(); it != std::get<n>(msg_).end(); ++it, ++omegaIt) {
                     //if(*omegaIt == 0.0) { // makes large difference for cosegmentation_bins, why?
                     l.ReceiveMessageSynchronized(*it);
                     //}
                  }
            }
            
            //std::advance(omegaIt, std::get<n>(msg_).size());
      });
   }
#endif

   // we write message change not into original reparametrization, but into temporary one named pot
   void ReceiveRestrictedMessages() 
   {
      meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [this](auto l) {
            constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
            if constexpr(l.CanCallReceiveRestrictedMessage()) {
                auto msg_begin = std::get<n>(msg_).begin();
                auto msg_end = std::get<n>(msg_).end();
                  for(auto it = msg_begin; it != msg_end; ++it) {
                     l.ReceiveRestrictedMessage(*it); 
                  }
            }
      });
   }

   struct can_receive_restricted_message {
      template<typename MESSAGE_DISPATCHER_TYPE>
         using invoke = typename std::is_same<std::integral_constant<bool,MESSAGE_DISPATCHER_TYPE::CanCallReceiveRestrictedMessage()>, std::integral_constant<bool,true> >::type;
   };
   constexpr static bool CanReceiveRestrictedMessage() 
   {
      return meta::any_of<MESSAGE_DISPATCHER_TYPELIST, can_receive_restricted_message>{};
   }

   // methods used by MessageIterator
   INDEX no_messages() const final
   {
      INDEX no_msgs = 0;
      meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [this,&no_msgs](auto l) {
            constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
            no_msgs += std::get<n>(msg_).size();
            } );
      return no_msgs;
   }

   // counts number of messages for which messages are sent
   INDEX no_send_messages() const final
   {
      INDEX no_calls = 0;
      meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [this,&no_calls](auto l) {
            constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
            if(l.sends_message_to_adjacent_factor()) {
               no_calls += std::get<n>(msg_).size();
            }
            } );
      return no_calls;
   }

   // as above, but if batch messages sending is enabled, such messages are counted only once.
   INDEX no_send_messages_calls() const 
   {
      INDEX no_calls = 0;
      meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [this,&no_calls](auto l) {
            constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
            if(l.CanCallSendMessages() && l.sends_message_to_adjacent_factor()) {
               if(!std::get<n>(msg_).empty()) {
                  ++no_calls;
               }
            } else if(l.sends_message_to_adjacent_factor()) {
               no_calls += std::get<n>(msg_).size();
            }
            } );
      return no_calls;
   }

   INDEX no_receive_messages() const
   {
       INDEX no_messages = 0;
       meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [this,&no_messages](auto l) {
            constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
            if(l.receives_message_from_adjacent_factor()) {
               no_messages += std::get<n>(msg_).size();
            }
       } );
       return no_messages;
   }

   // given some sending mask, skip over those iterators that have zero weight
   template<typename MSG_ITERATOR, typename WEIGHT_ITERATOR>
   struct conditional_message_iterator {
   public:
       conditional_message_iterator(MSG_ITERATOR msg_it, WEIGHT_ITERATOR weight_it)
           : msg_it_(msg_it),
           weight_it_(weight_it)
       {}

       const auto& operator*() const
       {
           while(*weight_it_ == 0) {
               ++msg_it_;
               ++weight_it_;
           }
           return *msg_it_;
       }
       auto& operator*() 
       {
           while(*weight_it_ == 0) {
               ++msg_it_;
               ++weight_it_;
           }
           return *msg_it_;
       }
       conditional_message_iterator& operator++()
       {
           ++msg_it_;
           ++weight_it_;
           return *this; 
       }
       bool operator==(const conditional_message_iterator o) const { return msg_it_ == o.msg_it_; }
       bool operator!=(const conditional_message_iterator o) const { return !(*this == o); } 
   private:
       mutable MSG_ITERATOR msg_it_;
       mutable WEIGHT_ITERATOR weight_it_; 
   };

   template<typename MSG_ITERATOR, typename WEIGHT_ITERATOR>
   conditional_message_iterator<MSG_ITERATOR,WEIGHT_ITERATOR> conditional_message_iterator_begin(MSG_ITERATOR msg_begin, WEIGHT_ITERATOR weight_it)
   {
       return conditional_message_iterator<MSG_ITERATOR,WEIGHT_ITERATOR>(msg_begin, weight_it); 
   }

   template<typename MSG_ITERATOR, typename WEIGHT_ITERATOR>
   conditional_message_iterator<MSG_ITERATOR,WEIGHT_ITERATOR> conditional_message_iterator_end(MSG_ITERATOR msg_begin, MSG_ITERATOR msg_end, WEIGHT_ITERATOR weight_it)
   {
       // find the last message
       auto last_active_msg = msg_begin;
       for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it, ++weight_it) {
           if(*weight_it > 0.0) {
               last_active_msg = msg_it;
           } 
       }
       ++last_active_msg;
       return conditional_message_iterator<MSG_ITERATOR,WEIGHT_ITERATOR>(last_active_msg, weight_it); 
   }


   void call_send_messages(FactorType& factor, const REAL send_weight)
   {
     meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [&](auto l) {
         // check whether the message supports batch updates. If so, call batch update, else call individual send message
         constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
         if constexpr(l.sends_message_to_adjacent_factor()) {
           auto msg_begin = std::get<n>(msg_).begin();
           auto msg_end = std::get<n>(msg_).end();

           if constexpr(l.CanCallSendMessages()) {

             const INDEX no_messages = std::get<n>(msg_).size(); 

             if(no_messages > 1) { 
               const REAL total_send_weight = no_messages*send_weight;
               l.SendMessages(factor, msg_begin, msg_end, total_send_weight);
             } else {
                l.SendMessage(&factor, *msg_begin, send_weight); 
             } 
           } else {
             for(auto msg_it = msg_begin; msg_it != msg_end; ++msg_it) {
                 l.SendMessage(&factor, *msg_it, send_weight); 
             }
           }
         }
     }); 
   }

   template<typename ITERATOR>
   void CallSendMessages(FactorType& factor, ITERATOR omegaIt) 
   {
     auto omega_begin = omegaIt;
     meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [&](auto l) {
         // check whether the message supports batch updates. If so, call batch update, else call individual send message
         constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
         if constexpr(l.sends_message_to_adjacent_factor()) {
           auto msg_begin = std::get<n>(msg_).begin();
           auto msg_end = std::get<n>(msg_).end();

           if constexpr(l.CanCallSendMessages()) {

             const INDEX no_messages = std::get<n>(msg_).size(); 
             const INDEX no_active_messages = std::count_if(omegaIt, omegaIt+no_messages, [](const REAL x){ return x > 0.0; });

             if(no_active_messages > 1) { 
               const REAL omega_sum = std::accumulate(omegaIt, omegaIt+no_messages, 0.0);
               auto msg_begin_c = conditional_message_iterator_begin(msg_begin, omegaIt);
               auto msg_end_c = conditional_message_iterator_end(msg_begin, msg_end, omegaIt);
               l.SendMessages(factor, msg_begin_c, msg_end_c, omega_sum);
               omegaIt += no_messages;
             } else {
                 for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it, ++omegaIt) {
                     if(*omegaIt > 0) {
                         l.SendMessage(&factor, *msg_it, *omegaIt); 
                     } 
                 }
             } 
           } else {
             for(auto msg_it = msg_begin; msg_it != msg_end; ++msg_it, ++omegaIt) {
               if(*omegaIt != 0.0) {
#ifndef NDEBUG
                 const REAL before_lb = LowerBound() + l.get_adjacent_factor(*msg_it)->LowerBound();
#endif
                 l.SendMessage(&factor, *msg_it, *omegaIt); 
#ifndef NDEBUG
                 const REAL after_lb = LowerBound() + l.get_adjacent_factor(*msg_it)->LowerBound();
                 assert(before_lb <= after_lb + eps);
#endif
               }
             }
           }
         }
     });
     assert(omegaIt - omega_begin == no_send_messages());
   }

#ifdef LP_MP_PARALLEL
   template<typename ITERATOR>
   void CallSendMessagesSynchronized(FactorType& factor, ITERATOR omegaIt) 
   {
     meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [&](auto l) {
         // check whether the message supports batch updates. If so, call batch update.
         // If not, check whether individual updates are supported. If yes, call individual updates. If no, do nothing
         if(l.sends_message_to_adjacent_factor()) {
         if constexpr(l.CanCallSendMessages()) {
             constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
             const REAL omega_sum = std::accumulate(omegaIt, omegaIt + std::get<n>(msg_).size(), 0.0);
             if(omega_sum > 0.0) { 
               f(l).SendMessagesSynchronized(factor, std::get<n>(msg_), omegaIt);
             }
             omegaIt += std::get<n>(msg_).size();
          } else {
                constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
                for(auto it = std::get<n>(msg_).begin(); it != std::get<n>(msg_).end(); ++it, ++omegaIt) {
                  if(*omegaIt != 0.0) {
                    l.SendMessageSynchronized(&factor, *it, *omegaIt); 
                  }
                }
          });
         }
     });
   }
#endif

   void send_messages(const REAL leave_weight)
   {
       assert(leave_weight >= 0.0);
       const auto no_calls = no_send_messages_calls();
       const auto send_weight = 1.0/(leave_weight + no_send_messages());
       if(no_calls == 1) {
           call_send_messages(factor_, send_weight); 
       } else if(no_calls > 1) {
           FactorType tmp_factor(factor_);
           call_send_messages(factor_, send_weight); 
       }
   }

   template<typename WEIGHT_VEC>
   void SendMessages(const WEIGHT_VEC& omega) 
   {
      assert(*std::min_element(omega.begin(), omega.end()) >= 0.0);
      assert(std::accumulate(omega.begin(), omega.end(), 0.0) <= 1.0 + eps);
      assert(std::distance(omega.begin(), omega.end()) == no_send_messages()); 
#ifndef NDEBUG
       const REAL before_lb = LowerBound();
#endif
      // do zrobienia: condition no_send_messages_calls also on omega. whenever omega is zero, we will not send messages
      const INDEX no_calls = no_send_messages_calls();

      if(no_calls == 1) {
        CallSendMessages(factor_, omega.begin());
      } else if( no_calls > 1 ) {
         // make a copy of the current reparametrization. The new messages are computed on it. Messages are updated implicitly and hence possibly the new reparametrization is automatically adjusted, which would interfere with message updates
         FactorType tmp_factor(factor_);

         CallSendMessages(tmp_factor, omega.begin());
      } else {
        assert(omega.size() == 0.0);
      }
#ifndef NDEBUG
       const REAL after_lb = LowerBound();
       assert(before_lb <= after_lb + eps);
#endif

   } 

   template<typename ITERATOR_1, typename ITERATOR_2>
   void adaptive_weight_rescaling_soft_max(ITERATOR_1 dual_improvement_begin, ITERATOR_1 dual_improvement_end, ITERATOR_2 omega_begin, ITERATOR_2 omega_end)
   {
       // soft max
       /*
       const auto scaling_exp = 100.0;
       const auto dual_improvement_sum_help = std::accumulate(dual_improvement.begin(), dual_improvement.end(), 0.0);
       const auto max_dual_improvement = *std::max_element(dual_improvement.begin(), dual_improvement.end());
       const auto offset = scaling_exp*max_dual_improvement/dual_improvement_sum_help;
       const auto dual_improvement_sum = std::accumulate(dual_improvement.begin(), dual_improvement.end(), 0.0, 
               [=](auto sum, auto x) { 
                   return sum + (x > 0 ? std::exp(scaling_exp*x/dual_improvement_sum_help - offset) : 0.0);
               });

       if(dual_improvement_sum_help > eps) {
           for(std::size_t i=0; i<omega.size(); ++i) {
               if(omega[i] > 0) {
                   dual_improvement[i] = 0.1*omega[i] + 0.9*omega_sum*(std::exp(scaling_exp*dual_improvement[i]/dual_improvement_sum_help - offset)/dual_improvement_sum );
               } else {
                   assert(dual_improvement[i] == 0.0);
               }
           }
       } else {
           std::copy(omega.begin(), omega.end(), dual_improvement.begin());
       }
       */
       assert(false);
   }

   template<typename ITERATOR_1, typename ITERATOR_2>
   void adaptive_weight_rescaling(ITERATOR_1 dual_improvement_begin, ITERATOR_1 dual_improvement_end, ITERATOR_2 omega_begin, ITERATOR_2 omega_end)
   {
       assert(std::distance(dual_improvement_begin, dual_improvement_end) == std::distance(omega_begin, omega_end));
       const auto dual_improvement_sum = std::accumulate(dual_improvement_begin, dual_improvement_end, 0.0);
       const auto omega_sum = std::accumulate(omega_begin, omega_end, 0.0);
       if(dual_improvement_sum > 0) {
           auto omega_it = omega_begin;
           for(auto it=dual_improvement_begin; it!=dual_improvement_end; ++it, ++omega_it) {
               *it = 0.5 * (*omega_it) + 0.5 * omega_sum * (*it)/dual_improvement_sum;
           }
       }
   }

   template<typename WEIGHT_VEC>
   void send_messages_with_adaptive_weights(const WEIGHT_VEC& omega)
   {
       assert(*std::min_element(omega.begin(), omega.end()) >= 0.0);
       assert(std::accumulate(omega.begin(), omega.end(), 0.0) <= 1.0 + eps);
       assert(std::distance(omega.begin(), omega.end()) == no_send_messages()); 

       // first go over all messages that will be send (i.e. where omega > 0) and record, how large an improvement will result from sending message with weight 1
       const auto omega_sum = std::accumulate(omega.begin(), omega.end(), 0.0);
       std::vector<REAL> dual_improvement;
       dual_improvement.reserve(omega.size());

       auto omega_it = omega.begin();
       // now record message update
       meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [&](auto l) {
         // check whether the message supports batch updates. If so, call batch update.
         // If not, check whether individual updates are supported. If yes, call individual updates. If no, do nothing
         constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
           if(!std::get<n>(msg_).empty()) {

             auto msg_begin = std::get<n>(msg_).begin();
             auto msg_end = std::get<n>(msg_).end(); 

             if constexpr(l.CanCallSendMessages()) {

               const INDEX no_messages = std::get<n>(msg_).size();
               const INDEX no_active_messages = std::count_if(omega_it, omega_it+no_messages, [](REAL x){ return x > 0.0; });
               REAL impr;
               if(no_active_messages > 0) { 
                 auto msg_begin_c = conditional_message_iterator_begin(msg_begin, omega_it);
                 auto msg_end_c = conditional_message_iterator_end(msg_begin, msg_end, omega_it);
                 impr = l.send_messages_improvement(msg_begin_c, msg_end_c);
                 //impr = l.send_messages_improvement(conditional_message_iterator(msg_begin, omega_it), conditional_message_iterator_end(msg_begin, omega_it));
                 assert(impr >= -eps);
                 impr = std::abs(impr);
               } 

               // record dual improvement
               for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it, ++omega_it) {
                   if(*omega_it > 0) {
                       dual_improvement.push_back(impr/no_active_messages);
                   } else {
                       dual_improvement.push_back(0.0);
                   }
               }

             } else {
                 for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it, ++omega_it) {
                    if(*omega_it > 0.0) {
                        dual_improvement.push_back( l.send_message_improvement(*msg_it) );
                        assert(dual_improvement.back() >= -eps);
                        dual_improvement.back() = std::abs(dual_improvement.back());
                    } else {
                        dual_improvement.push_back(0.0);
                    }
                 }
             }
           }
       });

       // reweight omega according to dual improvement
       assert(dual_improvement.size() == omega.size());
       adaptive_weight_rescaling(dual_improvement.begin(), dual_improvement.end(), omega.begin(), omega.end());

       //assert(std::abs( omega_sum - std::accumulate(dual_improvement.begin(), dual_improvement.end(), 0.0) ) <= eps);

       SendMessages(dual_improvement);
   }

   static constexpr INDEX active_messages_array_size = 16;

   template<typename MSG_ITERATOR, typename ACTIVE_ITERATOR>
   auto get_active_messages_array(MSG_ITERATOR msgs_begin, MSG_ITERATOR msgs_end, ACTIVE_ITERATOR active_begin)
   {
     using message_ptr_type = decltype( &(*msgs_begin) );
     std::array<message_ptr_type,active_messages_array_size> msgs;
     for(auto [no_msgs,msg_it]=std::pair(0,msgs_begin); msg_it!=msgs_end; ++msg_it, ++active_begin) {
       if(*active_begin > 0.0) {
         msgs[no_msgs++] = &(*msg_it);
       }
     }
     return msgs; 
   }

   template<typename MSG_ITERATOR, typename ACTIVE_ITERATOR>
   auto get_active_messages_vector(MSG_ITERATOR msgs_begin, MSG_ITERATOR msgs_end, ACTIVE_ITERATOR active_begin, const INDEX no_active_messages)
   {
     using message_ptr_type = decltype( &(*msgs_begin) );
     std::vector<message_ptr_type> msgs(no_active_messages);
     INDEX no_msgs = 0;
     for(auto msg_it=msgs_begin; msg_it!=msgs_end; ++msg_it, ++active_begin) {
       if(*active_begin > 0.0) {
         msgs[no_msgs++] = &(*msg_it);
       }
     }
     assert(no_active_messages == no_msgs); 
     return msgs; 
   }

   // choose order of messages to be sent and immediately reparametrize after each send message call and increase the remaining weights
   template<typename WEIGHT_VEC>
   void send_messages_residual(const WEIGHT_VEC& omega)
   {
       // first send messages in shared mode
       SendMessages(omega);

       auto omegaIt = omega.begin();
       REAL residual_omega = 0.0;
       meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [&](auto l) {
         // check whether the message supports batch updates. If so, call batch update.
         // If not, check whether individual updates are supported. If yes, call individual updates. If no, do nothing
         constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
         using message_ptr_type = decltype( *(std::get<n>(msg_).begin()) );
         if constexpr(l.sends_message_to_adjacent_factor()) {
           if(!std::get<n>(msg_).empty()) {
             auto& msgs = std::get<n>(msg_);
             auto msg_begin = msgs.begin();
             auto msg_end = msgs.end();

             if constexpr(l.CanCallSendMessages()) {

               const INDEX no_messages = std::get<n>(msg_).size();
               const INDEX no_active_messages = std::count_if(omegaIt, omegaIt+no_messages, [](REAL x){ return x > 0.0; });
               const REAL omega_sum = std::accumulate(omegaIt, omegaIt+no_messages, 0.0);
               residual_omega += omega_sum;
               if(no_active_messages > 0) { 
                 auto msg_begin_c = conditional_message_iterator_begin(msg_begin, omegaIt);
                 auto msg_end_c = conditional_message_iterator_end(msg_begin, msg_end, omegaIt);
                 l.SendMessages(factor_, msg_begin_c, msg_end_c, residual_omega);
                 //l.SendMessages(factor_, conditional_message_iterator(msgs.begin(), omegaIt), conditional_message_iterator_end( msgs.begin(), msgs.end(), omegaIt));
               } 
               omegaIt += no_messages;
  
             } else {
  
               for(auto it = msgs.begin(); it != msgs.end(); ++it, ++omegaIt) {
                 if(*omegaIt != 0.0) {
                   residual_omega += *omegaIt;
                   l.SendMessage(&factor_, *it, residual_omega); 
                 }
               }
             }
           }

         }
         assert(0.0 <= residual_omega && residual_omega <= 1.0 + eps);

       });
   }

#ifdef LP_MP_PARALLEL
   template<typename WEIGHT_VEC>
   void SendMessagesSynchronized(const WEIGHT_VEC& omega) 
   {
      // do zrobienia: condition no_send_messages_calls also on omega. whenever omega is zero, we will not send messages
      const INDEX no_calls = no_send_messages_calls();

      if(no_calls == 1) {
        CallSendMessagesSynchronized(factor_, omega.begin());
      } else if( no_calls > 1 ) {
         // make a copy of the current reparametrization. The new messages are computed on it. Messages are updated implicitly and hence possibly the new reparametrization is automatically adjusted, which would interfere with message updates
         FactorType tmp_factor(factor_);

         CallSendMessagesSynchronized(tmp_factor, omega.begin());
      } else {
        assert(omega.size() == 0.0);
      }
   } 
#endif

   bool ReceivesMessage() const
   {
      bool can_receive = false;
      meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [this,&can_receive](auto l) {
         constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
         if(l.receives_message_from_adjacent_factor() && !std::get<n>(msg_).empty()) {
            can_receive = true;
         }
      });
      return can_receive;
   }

   // check whether actually send messages is called. Can be false, even if CanSendMessages is true, e.g. when no message is present
   bool SendsMessage() const
   {
      bool calls_send = false;
      meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [&](auto l) {
            constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
            if(l.sends_message_to_adjacent_factor() && !std::get<n>(msg_).empty()) {
               calls_send = true;
            }
      });
      return calls_send;
   }


   template<typename ...MESSAGE_DISPATCHER_TYPES_REST>
   bool SendsMessage(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>, const INDEX) const 
   {
      throw std::runtime_error("message index out of bound");
   }
   template<typename MESSAGE_DISPATCHER_TYPE, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   bool SendsMessage(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...>, const INDEX cur_msg_idx) const // to get the current MESSAGE_TYPE
   {
      constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      const INDEX no_msgs = std::get<n>(msg_).size();
      if(cur_msg_idx < no_msgs) {
         MESSAGE_DISPATCHER_TYPE l;
         assert(!std::get<n>(msg_).empty());
         if( l.sends_message_to_adjacent_factor() ) {
            return true;
         } else {
           return false;
         }
      } else {
         return SendsMessage(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{}, cur_msg_idx - no_msgs);
      }
   }

   bool SendsMessage(const INDEX msg_idx) const final
   {
      return SendsMessage(MESSAGE_DISPATCHER_TYPELIST{}, msg_idx);
   }

   template<typename ...MESSAGE_DISPATCHER_TYPES_REST>
   bool ReceivesMessage(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>, const INDEX) const 
   {
      throw std::runtime_error("message index out of bound");
   }
   template<typename MESSAGE_DISPATCHER_TYPE, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   bool ReceivesMessage(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...>, const INDEX cur_msg_idx) const // to get the current MESSAGE_TYPE
   {
      constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      const INDEX no_msgs = std::get<n>(msg_).size();
      if(cur_msg_idx < no_msgs) {
         MESSAGE_DISPATCHER_TYPE l;
         assert(!std::get<n>(msg_).empty());
         if( l.receives_message_from_adjacent_factor() ) {
            return true;
         } else {
           return false;
         }
      } else {
         return ReceivesMessage(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{}, cur_msg_idx - no_msgs);
      }
   }

   bool ReceivesMessage(const INDEX msg_idx) const final
   {
       return ReceivesMessage(MESSAGE_DISPATCHER_TYPELIST{}, msg_idx);
   }

   // check whether actually receive restricted messages is called. Can be false, even if CanReceiveRestrictedMessages is true, e.g. when no message is present
   bool ReceivesRestrictedMessage() const
   {
      bool calls_receive_restricted = false;
      meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [&](auto l) {
            constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
            if(l.CanCallReceiveRestrictedMessage() && !std::get<n>(msg_).empty()) {
               calls_receive_restricted = true;
            }
      });
      return calls_receive_restricted;
   }

   // does factor call {Receive(Restricted)?|Send}Messages or does it compute primal? If not, UpdateFactor need not be called.
   bool FactorUpdated() const final
   {
      if(CanComputePrimal()) {
         return true;
      }
      if(ReceivesMessage()) {
         return true;
      }
      if(SendsMessage()) {
         return true;
      }
      if(ReceivesRestrictedMessage()) {
         return true;
      }
      return false;
   }

   void SetAndPropagatePrimal() const
   {
      assert(false);
     //assert(GetPrimalOffset() + PrimalSize() <= primal.size());
      //for(INDEX i=0; i<PrimalSize(); ++i) {
      //   primal[i + GetPrimalOffset()] = label[i];
      //}
      propagate_primal_through_messages();
   }

   struct apply_subgradient {
      apply_subgradient(double* _w, REAL _sign) : w(_w), sign(_sign) { assert(sign == 1.0 || sign == -1.0); }
      void operator[](const INDEX i) { w[i] = sign; }
      private:
      REAL* const w;
      const REAL sign;
   };
   constexpr static bool can_apply()
   {
      return FunctionExistence::has_apply<FactorType, void, apply_subgradient>(); 
   }
   virtual INDEX subgradient(double* w, const REAL sign) final
   {
      assert(sign == -1.0 || sign == 1.0);
      if constexpr(can_apply()) {
        apply_subgradient a(w,sign);
        factor_.apply(a);
      } else {
        assert(false);
      }
      return 0;
   }

   virtual REAL dot_product(double* w) final
   {
      class apply_dot_product {
      public:
         apply_dot_product(double* _w) : w(_w) {}
         void operator[](const INDEX i) { dp += w[i]; }
         REAL dot_product() const { return dp; }
      private:
         REAL* const w;
         REAL dp = 0;
      };

      apply_dot_product d(w);
      if constexpr(can_apply()) {
          factor_.apply(d);
      } else {
         assert(false);
      }
      return d.dot_product(); 
   }

   virtual void serialize_dual(load_archive& ar) final
   { factor_.serialize_dual(ar); }
   virtual void serialize_primal(load_archive& ar) final
   { factor_.serialize_primal(ar); } 
   virtual void serialize_dual(save_archive& ar) final
   { factor_.serialize_dual(ar); }
   virtual void serialize_primal(save_archive& ar) final
   { factor_.serialize_primal(ar); } 
   virtual void serialize_dual(allocate_archive& ar) final
   { factor_.serialize_dual(ar); }
   virtual void serialize_primal(allocate_archive& ar) final
   { factor_.serialize_primal(ar); } 
   virtual void serialize_dual(addition_archive& ar) final
   { factor_.serialize_dual(ar); }

   // returns size in bytes
   virtual INDEX dual_size() final
   {
      return dual_size_in_bytes()/sizeof(REAL);
   }

   virtual INDEX dual_size_in_bytes() final
   {
      allocate_archive ar;
      factor_.serialize_dual(ar);
      assert(ar.size() % sizeof(REAL) == 0);
      return ar.size();
   }

   virtual void divide(const REAL val) final
   {
      arithmetic_archive<operation::division> ar(val);
      factor_.serialize_dual(ar);
   }

   virtual INDEX primal_size_in_bytes() final
   {
      allocate_archive ar;
      factor_.serialize_primal(ar);
      return ar.size(); 
   }

   REAL LowerBound() const final {
      //return factor_.LowerBound(*this); 
      return factor_.LowerBound(); 
   } 

   REAL EvaluatePrimal() const final
   {
      return factor_.EvaluatePrimal();
   }


   FactorType* GetFactor() const { return &factor_; }
   FactorType* GetFactor() { return &factor_; }

  template<typename MESSAGE_TYPE>
  constexpr static 
  INDEX get_message_number()
  {
     static_assert(MESSAGE_TYPE::leftFactorNumber == FACTOR_NO || MESSAGE_TYPE::rightFactorNumber == FACTOR_NO,"");
     static_assert(MESSAGE_TYPE::leftFactorNumber != MESSAGE_TYPE::rightFactorNumber,""); // otherwise we cannot distinguish

     constexpr bool left = MESSAGE_TYPE::leftFactorNumber == FACTOR_NO;
     using dispatcher_type = typename meta::if_c<left, MessageDispatcher<MESSAGE_TYPE, LeftMessageFuncGetter> , MessageDispatcher<MESSAGE_TYPE, RightMessageFuncGetter>>;
     return  FactorContainerType::FindMessageDispatcherTypeIndex<dispatcher_type>();

     //if(MESSAGE_TYPE::leftFactorNumber == FACTOR_NO) {
     //   return FindMessageDispatcherTypeIndex<MessageDispatcher<MESSAGE_TYPE, LeftMessageFuncGetter>>();
     //} else {
     //   return FindMessageDispatcherTypeIndex<MessageDispatcher<MESSAGE_TYPE, RightMessageFuncGetter>>();
     //} 
  }

  template<typename MESSAGE_TYPE>
  auto get_messages() const 
  {
      std::vector<const MESSAGE_TYPE*> messages;
      constexpr auto n = get_message_number<MESSAGE_TYPE>();
      messages.reserve(std::get<n>(msg_).size());
      auto msg_begin = std::get<n>(msg_).begin();
      auto msg_end = std::get<n>(msg_).end();
      for(auto msg_it=msg_begin; msg_it!=msg_end; ++msg_it) {
          messages.push_back( &*msg_it );
      }
      return messages;
  }
   
protected:
   FactorType factor_; // the factor operation
public:
   INDEX primal_access_ = 0; // counts when primal was accessed last, do zrobienia: make setter and getter for clean interface or make MessageContainer a friend

   virtual void init_primal() final
   {
      factor_.init_primal();
   }
   void conditionally_init_primal(const INDEX timestamp) 
   {
      assert(primal_access_ <= timestamp);
      if(primal_access_ < timestamp) {
         factor_.init_primal();
         primal_access_ = timestamp;
      } 
   }

   INDEX runtime_estimate()
   {
     INDEX runtime = 0;
     assert(false);
     // go over all messages to be received and sum the dual size of connected factors

     // get number of messages to be sent and multiply by dual size of current factor (discount for SendMessages calls?)

     return runtime;
   }

   std::vector<FactorTypeAdapter*> get_adjacent_factors() const
   {
       std::vector<FactorTypeAdapter*> v;
       v.reserve(no_messages());

       meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [&](auto l) {
               constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
               auto msg_begin = std::get<n>(msg_).begin();
               auto msg_end = std::get<n>(msg_).end();
               for(auto it = msg_begin; it != msg_end; ++it) {
                   v.push_back(l.get_adjacent_factor(*it)); 
               }
       });

       return v; 
   }

   std::vector<message_trait> get_messages() const
   {
       std::vector<message_trait> v;
       v.reserve(no_messages());

       meta::for_each(MESSAGE_DISPATCHER_TYPELIST{}, [&](auto l) {
               constexpr INDEX n = FactorContainerType::FindMessageDispatcherTypeIndex<decltype(l)>();
               Chirality c = l.get_chirality();
               auto msg_begin = std::get<n>(msg_).begin();
               auto msg_end = std::get<n>(msg_).end();
               for(auto it = msg_begin; it != msg_end; ++it) {

                   message_trait t;
                   t.adjacent_factor = l.get_adjacent_factor(*it);
                   t.chirality = c;
                   t.sends_to_adjacent_factor = l.sends_message_to_adjacent_factor();
                   t.receives_from_adjacent_factor = l.receives_message_from_adjacent_factor();
                   t.adjacent_factor_sends = l.adjacent_factor_sends_message();
                   t.adjacent_factor_receives = l.adjacent_factor_receives_message();
                   
                   v.push_back(t);
               }
       });

       assert(v.size() == no_messages());
       return v; 
   }

protected:
   // pool memory allocator specific for this factor container
   // note: the below construction is not perfect when more than one solver is run simultaneously: The same allocator is used, yet the optimization problems are different + not thread safe.
   // -> investigate thread_local inline static! inline static however is only supported in C++17
   struct Allocator { // we enclose static allocator in nested class as only there (since C++11) we can access sizeof(FactorContainerType).
      using type = MemoryPool<FactorContainerType,4096*(sizeof(FactorContainerType)+sizeof(void*))>; 
      static type& get() {
         static type allocator;
         return allocator;
      }
   };
   
   // compile time metaprogramming to transform Factor-Message information into lists of which messages this factor must hold
   // first get lists with left and right message types
   struct get_msg_type_list {
      template<typename LIST>
         using invoke = typename LIST::MessageContainerType;
   };
   struct get_left_msg {
      template<class LIST>
         using invoke = typename std::is_same<meta::size_t<LIST::leftFactorNumber>, meta::size_t<FACTOR_NO>>::type;
   };
   struct get_right_msg {
      template<class LIST>
         using invoke = typename std::is_same<meta::size_t<LIST::rightFactorNumber>, meta::size_t<FACTOR_NO> >::type;
   };
   struct get_left_msg_container_type_list {
      template<class MSG_CONTAINER_TYPE> 
         using invoke = typename message_container_selector<MSG_CONTAINER_TYPE, MSG_CONTAINER_TYPE::no_left_factors(), MSG_CONTAINER_TYPE::no_right_factors(), MSG_CONTAINER_TYPE::estimated_no_left_factors(), MSG_CONTAINER_TYPE::estimated_no_right_factors(), Chirality::left>::type;//LIST::LeftMessageContainerStorageType;
   };
   struct get_right_msg_container_type_list {
      template<class MSG_CONTAINER_TYPE>
         using invoke = typename message_container_selector<MSG_CONTAINER_TYPE, MSG_CONTAINER_TYPE::no_left_factors(), MSG_CONTAINER_TYPE::no_right_factors(), MSG_CONTAINER_TYPE::estimated_no_left_factors(), MSG_CONTAINER_TYPE::estimated_no_right_factors(), Chirality::right>::type;//LIST::LeftMessageContainerStorageType;
   };

   using left_msg_list = meta::transform< meta::filter<typename FACTOR_MESSAGE_TRAIT::MessageList, get_left_msg>, get_msg_type_list>;
   using right_msg_list = meta::transform< meta::filter<typename FACTOR_MESSAGE_TRAIT::MessageList, get_right_msg>, get_msg_type_list>;
   using left_msg_container_list = meta::transform< meta::filter<typename FACTOR_MESSAGE_TRAIT::MessageList, get_left_msg>, get_left_msg_container_type_list>;
   using right_msg_container_list = meta::transform< meta::filter<typename FACTOR_MESSAGE_TRAIT::MessageList, get_right_msg>, get_right_msg_container_type_list>;

   // now construct a tuple with left and right dispatcher
   struct left_dispatch {
      template<class LIST>
         using invoke = MessageDispatcher<LIST, LeftMessageFuncGetter>;
   };
   struct right_dispatch {
      template<class LIST>
         using invoke = MessageDispatcher<LIST, RightMessageFuncGetter>;
   };
   using left_dispatcher_list = meta::transform< left_msg_list, left_dispatch >;
   using right_dispatcher_list = meta::transform< right_msg_list, right_dispatch >;

   using MESSAGE_DISPATCHER_TYPELIST = meta::concat<left_dispatcher_list, right_dispatcher_list>;

public:
   // helper function for getting the index in msg_ of given MESSAGE_DISPATCHER_TYPE
   template<typename MESSAGE_DISPATCHER_TYPE>
   static constexpr INDEX FindMessageDispatcherTypeIndex()
   {
      constexpr auto n = meta::find_index<MESSAGE_DISPATCHER_TYPELIST, MESSAGE_DISPATCHER_TYPE>::value;
      static_assert(n < meta::size<MESSAGE_DISPATCHER_TYPELIST>::value,"");
      return n;
   }


   // construct tuple holding messages for left and right dispatch
   // the tuple will hold some container for the message type. The container type is specified in the {Left|Right}MessageContainerStorageType fields of MessageList
   using msg_container_type_list = meta::concat<left_msg_container_list, right_msg_container_list>;

private:

   using msg_storage_type = meta::apply<meta::quote<std::tuple>, msg_container_type_list>;
   msg_storage_type msg_;

#ifdef LP_MP_PARALLEL
   // a recursive mutex is required only for SendMessagesTo{Left|Right}, as multiple messages may be have the same endpoints. Then the corresponding lock is acquired multiple times.
   // if no two messages have the same endpoints, an ordinary mutex is enough.
   std::recursive_mutex mutex_;
#endif

public:

   // functions for interfacing with external solver interface DD_ILP

   template<typename EXTERNAL_SOLVER>
   auto convert_variables_to_external(EXTERNAL_SOLVER& s, REAL x)
   { return s.add_variable(); }

   template<typename EXTERNAL_SOLVER>
   auto convert_variables_to_external(EXTERNAL_SOLVER& s, const vector<REAL>& v)
   { return s.add_vector(v); }

   template<typename EXTERNAL_SOLVER, std::size_t N>
   auto convert_variables_to_external(EXTERNAL_SOLVER& s, const std::array<REAL,N>& v)
   { return s.add_vector(v); }

   template<typename EXTERNAL_SOLVER>
   auto convert_variables_to_external(EXTERNAL_SOLVER& s, const std::vector<REAL>& v)
   { return s.add_vector(v); } 

   template<typename EXTERNAL_SOLVER>
   auto convert_variables_to_external(EXTERNAL_SOLVER& s, const matrix<REAL>& m)
   { return s.add_matrix(m); }

   template<typename EXTERNAL_SOLVER>
   auto convert_variables_to_external(EXTERNAL_SOLVER& s, const tensor3<REAL>& t)
   { return s.add_tensor(t); }

   template<typename EXTERNAL_SOLVER>
   auto load_external_variables(EXTERNAL_SOLVER& s, REAL x)
   { return s.load_variable(); }

   template<typename EXTERNAL_SOLVER>
   auto load_external_variables(EXTERNAL_SOLVER& s, vector<REAL>& x)
   { return s.load_vector(); }

   template<typename EXTERNAL_SOLVER, std::size_t N>
   auto load_external_variables(EXTERNAL_SOLVER& s, std::array<REAL,N>& x)
   { return s.load_vector(); }

   template<typename EXTERNAL_SOLVER>
   auto load_external_variables(EXTERNAL_SOLVER& s, std::vector<REAL>& x)
   { return s.load_vector(); }

   template<typename EXTERNAL_SOLVER>
   auto load_external_variables(EXTERNAL_SOLVER& s, matrix<REAL>& x)
   { return s.load_matrix(); }

   template<typename EXTERNAL_SOLVER>
   auto load_external_variables(EXTERNAL_SOLVER& s, tensor3<REAL>& x)
   { return s.load_tensor(); }

   template<typename EXTERNAL_SOLVER>
   void add_objective(EXTERNAL_SOLVER& s, REAL cost)
   { s.add_variable_objective(cost); }

   template<typename EXTERNAL_SOLVER>
   void add_objective(EXTERNAL_SOLVER& s, const vector<REAL>& cost)
   { s.add_vector_objective(cost); }

   template<typename EXTERNAL_SOLVER, std::size_t N>
   auto add_objective(EXTERNAL_SOLVER& s, std::array<REAL,N>& cost)
   { return s.add_vector_objective(cost); } 

   template<typename EXTERNAL_SOLVER>
   void add_objective(EXTERNAL_SOLVER& s, const std::vector<REAL>& cost)
   { s.add_vector_objective(cost); }

   template<typename EXTERNAL_SOLVER>
   void add_objective(EXTERNAL_SOLVER& s, const matrix<REAL>& cost)
   { s.add_matrix_objective(cost); }

   template<typename EXTERNAL_SOLVER>
   void add_objective(EXTERNAL_SOLVER& s, const tensor3<REAL>& cost)
   { s.add_tensor_objective(cost); }

   // functions for implementing external solver interface
   template<typename EXTERNAL_SOLVER>
   void construct_constraints_impl(EXTERNAL_SOLVER& s)
   {
      // transform exported variables to external solver variables
      auto vars = factor_.export_variables();
      auto external_vars = std::apply([this,&s](auto... x){ return std::make_tuple(this->convert_variables_to_external(s, x)...); }, vars);

      // unpack tuple and call construct_constraints function of factor
      auto construct_constraints_fun = [this,&s](auto... x) { this->factor_.construct_constraints(s, x...); };
      std::apply(construct_constraints_fun, external_vars);
   }


   template<typename EXTERNAL_SOLVER>
   void load_costs_impl(EXTERNAL_SOLVER& s)
   {
      // load external solver variables corresponding to reparametrization ones and add reparametrization as cost
      auto vars = factor_.export_variables();
      std::apply([this,&s](auto... x){ ((this->add_objective(s,x)), ...); },  vars);
      //auto external_vars = std::apply([this,&s](auto... x){ return std::make_tuple(this->leftFactor_->load_external_variables(s, x)...); }, vars);
      // for all variables,
   }


   template<typename SOLVER>
   void convert_primal_impl(SOLVER& s)
   {
      auto vars = factor_.export_variables();
      auto external_vars = std::apply([this,&s](auto... x){ return std::make_tuple(this->load_external_variables(s, x)...); }, vars); 

      auto convert_primal_fun = [this,&s](auto... x) { this->factor_.convert_primal(s, x...); };
      std::apply(convert_primal_fun, external_vars);

      //propagate_primal_through_messages();
   }

   virtual void construct_constraints(DD_ILP::external_solver_interface<DD_ILP::sat_solver>& s) final { construct_constraints_impl(s); }
   virtual void load_costs(DD_ILP::external_solver_interface<DD_ILP::sat_solver>& s) final {  }
   virtual void convert_primal(DD_ILP::external_solver_interface<DD_ILP::sat_solver>& solver) { convert_primal_impl(solver); }

   virtual void construct_constraints(DD_ILP::external_solver_interface<DD_ILP::problem_export>& s) final { construct_constraints_impl(s); }
   virtual void load_costs(DD_ILP::external_solver_interface<DD_ILP::problem_export>& s) final { load_costs_impl(s); } 
   virtual void convert_primal(DD_ILP::external_solver_interface<DD_ILP::problem_export>& solver) { convert_primal_impl(solver); }

#ifdef DD_ILP_WITH_GUROBI
   virtual void construct_constraints(DD_ILP::external_solver_interface<DD_ILP::gurobi_interface>& s) final { construct_constraints_impl(s); }
   virtual void load_costs(DD_ILP::external_solver_interface<DD_ILP::gurobi_interface>& s) final { load_costs_impl(s); } 
   virtual void convert_primal(DD_ILP::external_solver_interface<DD_ILP::gurobi_interface>& solver) { convert_primal_impl(solver); }
#endif
};

} // end namespace LP_MP

#endif // LP_MP_FACTORS_MESSAGES_HXX

