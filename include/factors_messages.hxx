#ifndef LP_MP_FACTORS_MESSAGES_HXX
#define LP_MP_FACTORS_MESSAGES_HXX

#include <vector>
#include <valarray> // do zrobienia: do not use
#include <tuple>
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
#include <assert.h>
#include <cxxabi.h>

#include "template_utilities.hxx"
#include "function_existence.hxx"
#include "meta/meta.hpp"

#include "factors/reparametrization_storage.hxx" 
#include "messages/message_storage.hxx"

#include "LP_MP.h"


// do zrobienia: remove these
#include <fstream>
#include <sstream>

// this file provides message and factor containers. The factors and messages are plugged into the container and then every method call is dispatched correctly with static polymorphism and template tricks.

// do zrobienia: Introduce MessageConstraint and FactorConstraint for templates
// cleanup name inconsistencies: MessageType, MessageDispatcher etc

namespace LP_MP {

// we must check existence of functions in message classes. The necessary test code is concentrated here. 
namespace FunctionExistence {

// Macros to construct help functions for checking existence of member functions of classes
LP_MP_FUNCTION_EXISTENCE_CLASS(HasReceiveMessageFromRight,ReceiveMessageFromRight);
LP_MP_FUNCTION_EXISTENCE_CLASS(HasReceiveMessageFromLeft, ReceiveMessageFromLeft);
   
LP_MP_FUNCTION_EXISTENCE_CLASS(HasReceiveRestrictedMessageFromRight,ReceiveRestrictedMessageFromRight);
LP_MP_FUNCTION_EXISTENCE_CLASS(HasReceiveRestrictedMessageFromLeft, ReceiveRestrictedMessageFromLeft);

LP_MP_FUNCTION_EXISTENCE_CLASS(HasSendMessageToRight,SendMessageToRight);
LP_MP_FUNCTION_EXISTENCE_CLASS(HasSendMessageToLeft, SendMessageToLeft);

LP_MP_FUNCTION_EXISTENCE_CLASS(HasSendMessagesToRight,SendMessagesToRight);
LP_MP_FUNCTION_EXISTENCE_CLASS(HasSendMessagesToLeft, SendMessagesToLeft);

LP_MP_FUNCTION_EXISTENCE_CLASS(HasRepamRight, RepamRight);
LP_MP_FUNCTION_EXISTENCE_CLASS(HasRepamLeft, RepamLeft);

LP_MP_FUNCTION_EXISTENCE_CLASS(HasComputeLeftFromRightPrimal, ComputeLeftFromRightPrimal);
LP_MP_FUNCTION_EXISTENCE_CLASS(HasComputeRightFromLeftPrimal, ComputeRightFromLeftPrimal); 

LP_MP_FUNCTION_EXISTENCE_CLASS(HasCheckPrimalConsistency, CheckPrimalConsistency); 

LP_MP_FUNCTION_EXISTENCE_CLASS(HasPrimalSize,PrimalSize);
LP_MP_FUNCTION_EXISTENCE_CLASS(HasPropagatePrimal, PropagatePrimal);
LP_MP_FUNCTION_EXISTENCE_CLASS(HasMaximizePotential, MaximizePotential);

LP_MP_FUNCTION_EXISTENCE_CLASS(HasCreateConstraints, CreateConstraints);

LP_MP_ASSIGNMENT_FUNCTION_EXISTENCE_CLASS(IsAssignable, operator[]);
}

// function getters for statically dispatching ReceiveMessage and SendMessage to left and right side correctly, used in FactorContainer
template<typename MSG_CONTAINER>
struct LeftMessageFuncGetter
{
   using ConnectedFactorType = typename MSG_CONTAINER::RightFactorContainer;

   constexpr static decltype(&MSG_CONTAINER::GetLeftMessage) GetMessageFunc() { return &MSG_CONTAINER::GetLeftMessage; }

   constexpr static decltype(&MSG_CONTAINER::ReceiveMessageFromRightContainer) GetReceiveFunc() { return &MSG_CONTAINER::ReceiveMessageFromRightContainer; }
   constexpr static decltype(&MSG_CONTAINER::ReceiveRestrictedMessageFromRightContainer) GetReceiveRestrictedFunc() { return &MSG_CONTAINER::ReceiveRestrictedMessageFromRightContainer; }
   template<typename ARRAY>
   constexpr static decltype(&MSG_CONTAINER::template SendMessageToRightContainer<ARRAY>) GetSendFunc() { return &MSG_CONTAINER::template SendMessageToRightContainer<ARRAY>; }

   template<typename LEFT_FACTOR, typename LEFT_REPAM, typename MSG_ARRAY, typename ITERATOR>
   constexpr static decltype(&MSG_CONTAINER::template SendMessagesToRightContainer<LEFT_FACTOR, LEFT_REPAM, MSG_ARRAY, ITERATOR>) GetSendMessagesFunc() 
   { return &MSG_CONTAINER::template SendMessagesToRightContainer<LEFT_FACTOR, LEFT_REPAM, MSG_ARRAY, ITERATOR>; }

   constexpr static bool 
   CanCallReceiveMessage()
   { return MSG_CONTAINER::CanCallReceiveMessageFromRightContainer(); }

   constexpr static bool
   CanCallReceiveRestrictedMessage()
   { return MSG_CONTAINER::CanCallReceiveRestrictedMessageFromRightContainer(); }

   template<typename REPAM_ARRAY>
   constexpr static bool CanCallSendMessage() 
   { return MSG_CONTAINER::template CanCallSendMessageToRightContainer<REPAM_ARRAY>(); }

   template<typename LEFT_FACTOR, typename REPAM_ARRAY, typename MSG_ARRAY, typename ITERATOR>
   constexpr static bool 
   CanCallSendMessages()
   { return MSG_CONTAINER::template CanCallSendMessagesToRightContainer<LEFT_FACTOR, REPAM_ARRAY, MSG_ARRAY, ITERATOR>(); }

   // do zrobienia: rename CanPropagatePrimalThroughMessage
   constexpr static bool CanComputePrimalThroughMessage()
   { return MSG_CONTAINER::CanComputeRightFromLeftPrimal(); }
   constexpr static decltype(&MSG_CONTAINER::ComputeRightFromLeftPrimal) GetComputePrimalThroughMessageFunc()
   { return &MSG_CONTAINER::ComputeRightFromLeftPrimal; }
};

template<typename MSG_CONTAINER>
struct RightMessageFuncGetter
{
   using ConnectedFactorType = typename MSG_CONTAINER::LeftFactorContainer;

   constexpr static decltype(&MSG_CONTAINER::GetRightMessage) GetMessageFunc() { return &MSG_CONTAINER::GetRightMessage; }

   constexpr static decltype(&MSG_CONTAINER::ReceiveMessageFromLeftContainer) GetReceiveFunc() { return &MSG_CONTAINER::ReceiveMessageFromLeftContainer; }
   constexpr static decltype(&MSG_CONTAINER::ReceiveRestrictedMessageFromLeftContainer) GetReceiveRestrictedFunc() { return &MSG_CONTAINER::ReceiveRestrictedMessageFromLeftContainer; }
   template<typename ARRAY>
   constexpr static decltype(&MSG_CONTAINER::template SendMessageToLeftContainer<ARRAY>) GetSendFunc() { return &MSG_CONTAINER::template SendMessageToLeftContainer<ARRAY>; }

   template<typename RIGHT_FACTOR, typename RIGHT_REPAM, typename MSG_ARRAY, typename ITERATOR>
   constexpr static decltype(&MSG_CONTAINER::template SendMessagesToLeftContainer<RIGHT_FACTOR, RIGHT_REPAM, MSG_ARRAY, ITERATOR>) GetSendMessagesFunc() 
   { return &MSG_CONTAINER::template SendMessagesToLeftContainer<RIGHT_FACTOR, RIGHT_REPAM, MSG_ARRAY, ITERATOR>; }

   constexpr static bool CanCallReceiveMessage() 
   { return MSG_CONTAINER::CanCallReceiveMessageFromLeftContainer(); }

   constexpr static bool CanCallReceiveRestrictedMessage() 
   { return MSG_CONTAINER::CanCallReceiveRestrictedMessageFromLeftContainer(); }

   template<typename REPAM_ARRAY>
   constexpr static bool CanCallSendMessage() 
   { return MSG_CONTAINER::template CanCallSendMessageToLeftContainer<REPAM_ARRAY>(); }

   template<typename RIGHT_FACTOR, typename RIGHT_REPAM, typename MSG_ARRAY, typename ITERATOR>
   constexpr static bool
   CanCallSendMessages()
   { return MSG_CONTAINER::template CanCallSendMessagesToLeftContainer<RIGHT_FACTOR, RIGHT_REPAM, MSG_ARRAY, ITERATOR>(); }

   constexpr static bool CanComputePrimalThroughMessage()
   { return MSG_CONTAINER::CanComputeLeftFromRightPrimal(); }
   constexpr static decltype(&MSG_CONTAINER::ComputeLeftFromRightPrimal) GetComputePrimalThroughMessageFunc()
   { return &MSG_CONTAINER::ComputeLeftFromRightPrimal; }
};

template<class MSG_CONTAINER, template<typename> class FuncGetter>
struct MessageDispatcher
{
   using ConnectedFactorType = typename FuncGetter<MSG_CONTAINER>::ConnectedFactorType; // this is the type of factor container to which the message is connected

   constexpr static bool CanCallReceiveMessage() { return FuncGetter<MSG_CONTAINER>::CanCallReceiveMessage(); }
   static void ReceiveMessage(MSG_CONTAINER& t)
   {
      auto staticMemberFunc = FuncGetter<MSG_CONTAINER>::GetReceiveFunc();
      return (t.*staticMemberFunc)();
   }
   constexpr static bool CanCallReceiveRestrictedMessage() { return FuncGetter<MSG_CONTAINER>::CanCallReceiveRestrictedMessage(); }
   static void ReceiveRestrictedMessage(MSG_CONTAINER& t, PrimalSolutionStorage::Element primal)
   {
      auto staticMemberFunc = FuncGetter<MSG_CONTAINER>::GetReceiveRestrictedFunc();
      return (t.*staticMemberFunc)(primal);
   }

   // individual message sending
   template<typename REPAM_ARRAY>
   constexpr static bool CanCallSendMessage() { return FuncGetter<MSG_CONTAINER>::template CanCallSendMessage<REPAM_ARRAY>(); }
   template<typename REPAM_ARRAY>
   static void SendMessage(MSG_CONTAINER& t, const REPAM_ARRAY& repam, const REAL omega)
   {
      auto staticMemberFunc = FuncGetter<MSG_CONTAINER>::template GetSendFunc<REPAM_ARRAY>();
      return (t.*staticMemberFunc)(repam, omega);
   }

   // batch message sending
   template<typename FACTOR, typename REPAM_ARRAY, typename MSG_ARRAY, typename ITERATOR>
   constexpr static bool CanCallSendMessages() { return FuncGetter<MSG_CONTAINER>::template CanCallSendMessages<FACTOR, REPAM_ARRAY, MSG_ARRAY, ITERATOR>(); }

   template<typename FACTOR, typename REPAM_ARRAY, typename MSG_ARRAY, typename ITERATOR>
   static void SendMessages(const FACTOR& f, const REPAM_ARRAY& repam, const MSG_ARRAY& msgs, ITERATOR omegaBegin)
   {
      auto staticMemberFunc = FuncGetter<MSG_CONTAINER>::template GetSendMessagesFunc<FACTOR, REPAM_ARRAY, MSG_ARRAY, ITERATOR>();
      (*staticMemberFunc)(f, repam, msgs, omegaBegin);
   }

   static REAL GetMessage(MSG_CONTAINER& t, const INDEX i)
   {
      auto staticMemberFunc = FuncGetter<MSG_CONTAINER>::GetMessageFunc();
      return (t.*staticMemberFunc)(i);
   }

   constexpr static bool CanComputePrimalThroughMessage() // do zrobienia: return false, if the factor from which this is called computes its own primal already
   {
      return FuncGetter<MSG_CONTAINER>::CanComputePrimalThroughMessage();
   }

   static void ComputePrimalThroughMessage(MSG_CONTAINER& t, typename PrimalSolutionStorage::Element primal) 
   {
      auto staticMemberFunc = FuncGetter<MSG_CONTAINER>::GetComputePrimalThroughMessageFunc();
      return (t.*staticMemberFunc)(primal);
   }
};

// this container does the following: as it can never shrink, we choose the following growth strategy:
// the first sizes are double, when 16 is reached, we increase by blocks of 16. 
// by this approach we can save one variable over std::vector, which needs three pointers
// do zrobienia: use own memory allocator. Possibly use stride 4 instead of 16.
template<typename T>
class VariableSizeMessageContainer {
public:
   VariableSizeMessageContainer() { begin_ = nullptr; end_ = nullptr; }
   ~VariableSizeMessageContainer() { static_assert(std::is_pointer<T>::value, "Message container must hold pointers to messages"); }
   INDEX size() const { return end_ - begin_; }
   const T operator[](const INDEX i) const { assert(i < size()); return begin_[i]; }

   void push_back(T t) {
      INDEX newSize = std::min(size()*2, size() + 16);
      if(size() == 0) {
         begin_ = new T[1];
         end_ = begin_+1;
         begin_[0] = t;
      } else if(size() == 1 || size() == 2 || size() == 4 || size() == 8 || size()%16 == 0) {
         T beginOld = begin_;
         T endOld = end_;
         begin_ = new T[newSize];
         end_ = begin_ + newSize;
         for(INDEX i=0; i<endOld - beginOld; ++i) {
            begin_[i] = beginOld[i];
         }
         begin_[endOld - beginOld] = t;
         delete[] beginOld;
      } else {
         end_[0] = t;
      }
   }
private:
   T* begin_;
   T* end_;
};

template<INDEX NO_ELEMENTS, typename T>
class FixedSizeMessageContainer : public std::array<T,NO_ELEMENTS> {
public: 
   FixedSizeMessageContainer() { this->fill(nullptr); }
   ~FixedSizeMessageContainer() { static_assert(std::is_pointer<T>::value, "Message container must hold pointers to messages"); }
   void push_back(T t) {
      // do zrobienia: possibly use binary search when NO_ELEMENTS is bigger than some threshold
      for(INDEX i=0; i<NO_ELEMENTS; ++i) {
         if(this->operator[](i) == nullptr) {
            this->operator[](i) = t;
            return;
         }
      }
      throw std::range_error("added more messages than can be held");
   }
   constexpr INDEX size() const { return NO_ELEMENTS; }
};

// holds at most NO_ELEMENTS in std::array. Unused entries have nullptr in them
template<INDEX NO_ELEMENTS, typename T>
class UpToFixedSizeMessageContainer : public std::array<T,NO_ELEMENTS> {
public:
   UpToFixedSizeMessageContainer() : size_(0) { this->fill(nullptr); }
   ~UpToFixedSizeMessageContainer() { 
      static_assert(std::is_pointer<T>::value, "Message container must hold pointers to messages"); 
      static_assert(NO_ELEMENTS > 0, "");
   }
   void push_back(T t) {
      assert(size_ < NO_ELEMENTS);
      this->operator[](size_) = t;
      ++size_;
   }
   INDEX size() const { return size_; }

private:
   INDEX size_;
};


// N=0 means variable number of messages, > 0 means compile time fixed number of messages and <0 means at most compile time number of messages
// see config.hxx for shortcuts
template<SIGNED_INDEX N, typename MESSAGE_CONTAINER_TYPE>
struct MessageContainerSelector {
   using type = typename std::conditional<(N > 0), FixedSizeMessageContainer<INDEX(N),MESSAGE_CONTAINER_TYPE*>,
        typename std::conditional<(N < 0), UpToFixedSizeMessageContainer<INDEX(-N),MESSAGE_CONTAINER_TYPE*>, std::vector<MESSAGE_CONTAINER_TYPE*> >::type >::type;
};

// there are two possible choices: (i) knowing message size in advance, (ii) holding message explicitly (e.g. when one factor is reparametrized implicitly)
// provide message storage classes for the four combinations of these cases
// do zrobienia: replace std::vector
// do zrobienia: simdize the two classes below
// do zrobienia: if message knows its size (constexpr), derive automatically from it.
class VariableSizeMessageStorage : public std::vector<REAL>
{
public:
   VariableSizeMessageStorage(const INDEX msg_size) : std::vector<REAL>(msg_size,0.0) {}
}; 

template<INDEX N>
class FixedSizeMessageStorage : public std::array<REAL,N>
{
public:
   FixedSizeMessageStorage(const INDEX msg_size) { assert(msg_size == N); std::array<REAL,N>::fill(0.0); }
   FixedSizeMessageStorage() { std::array<REAL,N>::fill(0.0); }
}; 

template<INDEX N>
class FixedSizeEmptyMessageStorage
{
public:
   FixedSizeEmptyMessageStorage(const INDEX msg_size) { static_assert(msg_size == N,""); }
   FixedSizeEmptyMessageStorage() {}
   constexpr static INDEX size() { return N; }
};
class VariableSizeEmptyMessageStorage
{
public:
   VariableSizeEmptyMessageStorage(const INDEX msg_size) : n_(msg_size) {}
   INDEX size() const { return n_; }
private:
   INDEX n_;
};

template<SIGNED_INDEX SIZE, bool HOLD_MESSAGE>
struct MessageStorageSelector {};

template<SIGNED_INDEX SIZE>
struct MessageStorageSelector<SIZE, false> {
   using type = typename std::conditional<(SIZE >= 0), FixedSizeEmptyMessageStorage<INDEX(SIZE)>, VariableSizeEmptyMessageStorage>::type;
};
template<SIGNED_INDEX SIZE>
struct MessageStorageSelector<SIZE, true> {
   using type = typename std::conditional<(SIZE >= 0), FixedSizeMessageStorage<INDEX(SIZE)>, VariableSizeMessageStorage>::type;
};


// Class holding message and left and right factor
// We simulate holding two messages, one for the left and one for the right factor, each being the negative of the other one. Physically, we just hold one.
// The following sign convention must be followed: 
// {Left|Right}Repam use +
// message computation uses -
// Dispatch of signs to {Left|Right}Repam:
// Left | Right
//  -   |   +
//
template<typename MESSAGE_TYPE, 
         INDEX LEFT_FACTOR_NO, INDEX RIGHT_FACTOR_NO, SIGNED_INDEX NO_OF_LEFT_FACTORS, SIGNED_INDEX NO_OF_RIGHT_FACTORS,
         SIGNED_INDEX MESSAGE_SIZE, 
         typename FACTOR_MESSAGE_TRAIT, 
         INDEX MESSAGE_NO
         >
class MessageContainer : public MessageStorageSelector<MESSAGE_SIZE,true>::type, public MessageTypeAdapter
{
public:
   using leftFactorNumber_t = std::integral_constant<INDEX, LEFT_FACTOR_NO>;
   static constexpr INDEX leftFactorNumber = LEFT_FACTOR_NO;
   static constexpr INDEX rightFactorNumber = RIGHT_FACTOR_NO;

   typedef MessageContainer<MESSAGE_TYPE, LEFT_FACTOR_NO, RIGHT_FACTOR_NO, NO_OF_LEFT_FACTORS, NO_OF_RIGHT_FACTORS, MESSAGE_SIZE, FACTOR_MESSAGE_TRAIT, MESSAGE_NO> MessageContainerType;
   typedef MESSAGE_TYPE MessageType;
   typedef typename MessageStorageSelector<MESSAGE_SIZE,true>::type MessageStorageType; // do zrobienia: true is just for now. In general, message need not hold actual message, except when some factor is reparametrized implicitly

   // structures used in FactorContainer to hold pointers to messages
   using LeftMessageContainerStorageType = typename MessageContainerSelector<NO_OF_LEFT_FACTORS, MessageContainerType>::type;
   using RightMessageContainerStorageType = typename MessageContainerSelector<NO_OF_RIGHT_FACTORS, MessageContainerType>::type;

   // FactorContainer
   using LeftFactorContainer = meta::at_c<typename FACTOR_MESSAGE_TRAIT::FactorList, leftFactorNumber>;
   using RightFactorContainer = meta::at_c<typename FACTOR_MESSAGE_TRAIT::FactorList, rightFactorNumber>;
   // Factor
   using LeftFactorType = typename LeftFactorContainer::FactorType;
   using RightFactorType = typename RightFactorContainer::FactorType;

   // message size known in advance
   MessageContainer(MESSAGE_TYPE msg_op, LeftFactorContainer* const l, RightFactorContainer* const r) 
      :MessageStorageType(),
      msg_op_(msg_op),
      leftFactor_(l),
      rightFactor_(r)
   {
      leftFactor_->template AddMessage<MessageDispatcher<MessageContainerType, LeftMessageFuncGetter>, MessageContainerType>(this);
      rightFactor_->template AddMessage<MessageDispatcher<MessageContainerType, RightMessageFuncGetter>, MessageContainerType>(this);
   }
   // message size not known in advance
   MessageContainer(MESSAGE_TYPE msg_op, LeftFactorContainer* const l, RightFactorContainer* const r, const INDEX msg_size) 
      : MessageStorageType(msg_size),
      msg_op_(msg_op),
      leftFactor_(l), 
      rightFactor_(r) 
   {
      int status;
      //std::cout << "msg holding type = " << abi::__cxa_demangle(typeid(*this).name(),0,0,&status) << "\n";
      //std::cout << FunctionExistence::IsAssignable<RightFactorContainer,REAL,INDEX>() << "\n";
      //std::cout << "msg holding type = " << abi::__cxa_demangle(typeid(msg_op_).name(),0,0,&status) << "\n";
      //std::cout << "left factor number = " << leftFactorNumber << "\n";
      //std::cout << "right factor number = " << rightFactorNumber << "\n";
      //std::cout << "left factor type = " << abi::__cxa_demangle(typeid(LeftFactorContainer).name(),0,0,&status) << "\n";
      //std::cout << "right factor type = " << abi::__cxa_demangle(typeid(RightFactorContainer).name(),0,0,&status) << "\n";
      // register messages in factors
      leftFactor_->template AddMessage<MessageDispatcher<MessageContainerType, LeftMessageFuncGetter>, MessageContainerType>(this);
      rightFactor_->template AddMessage<MessageDispatcher<MessageContainerType, RightMessageFuncGetter>, MessageContainerType>(this);
   }
   ~MessageContainer() {
      static_assert(meta::unique<typename FACTOR_MESSAGE_TRAIT::MessageList>::size() == FACTOR_MESSAGE_TRAIT::MessageList::size(), 
            "Message list must have unique elements");
      static_assert(MESSAGE_NO >= 0 && MESSAGE_NO < FACTOR_MESSAGE_TRAIT::MessageList::size(), "message number must be smaller than length of message list");
      static_assert(leftFactorNumber < FACTOR_MESSAGE_TRAIT::FactorList::size(), "left factor number out of bound");
      static_assert(rightFactorNumber < FACTOR_MESSAGE_TRAIT::FactorList::size(), "right factor number out of bound");
      // do zrobienia: put message constraint here, i.e. which methods MESSAGE_TYPE must minimally implement
   } 

   // we can throw these functions out again
   /*
   bool CanSendMessageToLeft() const final
   {
      if(CanCallSendMessageToLeftContainer<std::vector<REAL>>()) return true;
      if(CanCallSendMessagesToLeftContainer<RightFactorType,RightFactorContainer,std::vector<MessageContainerType*>,std::vector<REAL>::iterator>()) return true;
      else return false;
   }
   bool CanSendMessageToRight() const final
   {
      if(CanCallSendMessageToRightContainer<std::vector<REAL>>() ) return true;//|| CanCallSendMessagesToRightContainer<>()) return true;
      if(CanCallSendMessagesToRightContainer<LeftFactorType,LeftFactorContainer,std::vector<MessageContainerType*>,std::vector<REAL>::iterator>()) return true;
      else return false;
   }
   */

   constexpr static bool
   CanCallReceiveMessageFromRightContainer()
   { 
      return FunctionExistence::HasReceiveMessageFromRight<MessageType, void, 
      decltype(rightFactor_->GetFactor()), decltype(*rightFactor_), MessageContainerType>(); 
   }
   void ReceiveMessageFromRightContainer()
   { msg_op_.ReceiveMessageFromRight(rightFactor_->GetFactor(),*rightFactor_, *static_cast<MessageContainerView<Chirality::right>*>(this) ); }

   // do zrobienia: must use one additional argument for primal storage
   constexpr static bool
   CanCallReceiveRestrictedMessageFromRightContainer()
   { 
      return FunctionExistence::HasReceiveRestrictedMessageFromRight<MessageType, void, 
      decltype(rightFactor_->GetFactor()), decltype(*rightFactor_), MessageContainerType, PrimalSolutionStorage::Element>(); // do zrobienia: signature is slighly different: MessageContainerType is not actually used
   }
   void ReceiveRestrictedMessageFromRightContainer(PrimalSolutionStorage::Element primal)
   {
      msg_op_.ReceiveRestrictedMessageFromRight(rightFactor_->GetFactor(), *rightFactor_, *static_cast<RestrictedMessageContainerView<Chirality::left>*>(this), primal + rightFactor_->GetPrimalOffset());
   }

   constexpr static bool 
   CanCallReceiveMessageFromLeftContainer()
   { 
      return FunctionExistence::HasReceiveMessageFromLeft<MessageType, void, 
      decltype(leftFactor_->GetFactor()), decltype(*leftFactor_), MessageContainerType>(); 
   }
   void ReceiveMessageFromLeftContainer()
   { msg_op_.ReceiveMessageFromLeft(leftFactor_->GetFactor(), *leftFactor_, *static_cast<MessageContainerView<Chirality::left>*>(this) ); }

   // do zrobienia: must use one additional argument for primal storage
   constexpr static bool
   CanCallReceiveRestrictedMessageFromLeftContainer()
   { 
      return FunctionExistence::HasReceiveRestrictedMessageFromLeft<MessageType, void, 
      decltype(leftFactor_->GetFactor()), decltype(*leftFactor_), MessageContainerType, PrimalSolutionStorage::Element>(); 
   }
   void ReceiveRestrictedMessageFromLeftContainer(PrimalSolutionStorage::Element primal)
   {
      msg_op_.ReceiveRestrictedMessageFromLeft(leftFactor_->GetFactor(), *leftFactor_, *static_cast<RestrictedMessageContainerView<Chirality::right>*>(this), primal + leftFactor_->GetPrimalOffset());
   }


   template<typename REPAM_ARRAY>
   constexpr static bool 
   CanCallSendMessageToRightContainer()
   { 
      return FunctionExistence::HasSendMessageToRight<MessageType, void, 
      decltype(leftFactor_->GetFactor()), REPAM_ARRAY, MessageContainerType, REAL>(); 
   }

   template<typename ARRAY>
   void SendMessageToRightContainer(const ARRAY& repam, const REAL omega)
   {
      msg_op_.SendMessageToRight(leftFactor_->GetFactor(), repam, *static_cast<MessageContainerView<Chirality::left>*>(this), omega);
   }

   template<typename REPAM_ARRAY>
   constexpr static bool
   CanCallSendMessageToLeftContainer()
   { 
      return FunctionExistence::HasSendMessageToLeft<MessageType, void, 
      decltype(rightFactor_->GetFactor()), REPAM_ARRAY, MessageContainerType, REAL>(); 
   }

   template<typename ARRAY>
   void SendMessageToLeftContainer(const ARRAY& repam, const REAL omega)
   {
      msg_op_.SendMessageToLeft(rightFactor_->GetFactor(), repam, *static_cast<MessageContainerView<Chirality::right>*>(this), omega);
   }

   template<typename RIGHT_FACTOR, typename RIGHT_REPAM, typename MSG_ARRAY, typename ITERATOR>
   constexpr static bool
   CanCallSendMessagesToLeftContainer()
   { 
      return FunctionExistence::HasSendMessagesToLeft<MessageType, void, RIGHT_FACTOR, RIGHT_REPAM, MSG_ARRAY, ITERATOR>();
   }
   template<typename RIGHT_FACTOR, typename RIGHT_REPAM, typename MSG_ARRAY, typename ITERATOR>
   static void SendMessagesToLeftContainer(const RIGHT_FACTOR& rightFactor, const RIGHT_REPAM& repam, const MSG_ARRAY& msgs, ITERATOR omegaBegin) 
   {
      // this is not nice: heavy static casting!
      // We get msgs an array with pointers to messages. We wrap it so that operator[] gives a reference to the respective message with the correct view
      struct ViewWrapper : public MSG_ARRAY {
         MessageContainerView<Chirality::right>& operator[](const INDEX i) const 
         { 
            return *static_cast<MessageContainerView<Chirality::right>*>( (static_cast<const MSG_ARRAY*>(this)->operator[](i)) ); 
            //return *static_cast<MessageContainerView<Chirality::right>*>( operator[](i) ); 
            //return operator[](i);
         }
      };
      return MessageType::SendMessagesToLeft(rightFactor, repam, *static_cast<const ViewWrapper*>(&msgs), omegaBegin);

      struct ViewWrapper2 {
         ViewWrapper2(const MSG_ARRAY& msgs) : msgs_(msgs) {}
         MessageContainerView<Chirality::right>& operator[](const INDEX i) const {
            return static_cast<MessageContainerView<Chirality::right>&>( *(msgs_[i]) );
         }
         INDEX size() const {
            return msgs_.size();
         }
         private:
         const MSG_ARRAY& msgs_;
      } msgsProxy(msgs);

      return MessageType::SendMessagesToLeft(rightFactor, repam, msgsProxy, omegaBegin);
   }

   template<typename LEFT_FACTOR, typename LEFT_REPAM, typename MSG_ARRAY, typename ITERATOR>
   constexpr static bool
   CanCallSendMessagesToRightContainer()
   { 
      return FunctionExistence::HasSendMessagesToRight<MessageType, void, LEFT_FACTOR, LEFT_REPAM, MSG_ARRAY, ITERATOR>(); 
   }
   template<typename LEFT_FACTOR, typename LEFT_REPAM, typename MSG_ARRAY, typename ITERATOR>
   static void SendMessagesToRightContainer(const LEFT_FACTOR& leftFactor, const LEFT_REPAM& repam, const MSG_ARRAY& msgs, ITERATOR omegaBegin) 
   {
      struct ViewWrapper : public MSG_ARRAY {
         MessageContainerView<Chirality::left>& operator[](const INDEX i) const 
         { 
            return *static_cast<MessageContainerView<Chirality::left>*>( (static_cast<const MSG_ARRAY*>(this)->operator[](i)) ); 
         }
      };
      MessageType::SendMessagesToRight(leftFactor, repam, *static_cast<const ViewWrapper*>(&msgs), omegaBegin);
   }

   constexpr static bool
   CanComputeRightFromLeftPrimal()
   {
      return FunctionExistence::HasComputeRightFromLeftPrimal<MessageType,void,
             PrimalSolutionStorage::Element, decltype(leftFactor_->GetFactor()),
             PrimalSolutionStorage::Element, decltype(rightFactor_->GetFactor())>();
   }
   constexpr static bool
   CanComputeLeftFromRightPrimal()
   {
      return FunctionExistence::HasComputeLeftFromRightPrimal<MessageType,void,
             PrimalSolutionStorage::Element, decltype(leftFactor_->GetFactor()),
             PrimalSolutionStorage::Element, decltype(rightFactor_->GetFactor())>();
   }

   void ComputeRightFromLeftPrimal(typename PrimalSolutionStorage::Element primal) 
   {
      msg_op_.ComputeRightFromLeftPrimal(primal + leftFactor_->GetPrimalOffset(), leftFactor_->GetFactor(), primal + rightFactor_->GetPrimalOffset(), rightFactor_->GetFactor());
      rightFactor_->PropagatePrimal(primal + rightFactor_->GetPrimalOffset());
      rightFactor_->ComputePrimalThroughMessages(primal);

   }

   void ComputeLeftFromRightPrimal(PrimalSolutionStorage::Element primal)
   {
      msg_op_.ComputeLeftFromRightPrimal(primal + leftFactor_->GetPrimalOffset(), leftFactor_->GetFactor(), primal + rightFactor_->GetPrimalOffset(), rightFactor_->GetFactor());
      leftFactor_->PropagatePrimal(primal + leftFactor_->GetPrimalOffset());
      leftFactor_->ComputePrimalThroughMessages(primal);
   }

   constexpr static bool
   CanCheckPrimalConsistency()
   {
      return FunctionExistence::HasCheckPrimalConsistency<MessageType,bool,
          PrimalSolutionStorage::Element, typename LeftFactorContainer::FactorType*,
          PrimalSolutionStorage::Element, typename RightFactorContainer::FactorType*>();
   }
   template<bool ENABLE=CanCheckPrimalConsistency()>
   typename std::enable_if<ENABLE,bool>::type
   CheckPrimalConsistencyImpl(PrimalSolutionStorage::Element primal) const
   {
      static_assert(ENABLE == CanCheckPrimalConsistency(),"");
      return msg_op_.CheckPrimalConsistency(primal + leftFactor_->GetPrimalOffset(), leftFactor_->GetFactor(), primal + rightFactor_->GetPrimalOffset(), rightFactor_->GetFactor());
   }
   template<bool ENABLE=CanCheckPrimalConsistency()>
   typename std::enable_if<!ENABLE,bool>::type
   CheckPrimalConsistencyImpl(PrimalSolutionStorage::Element primal) const
   {
      static_assert(ENABLE == CanCheckPrimalConsistency(),"");
      return true;
   }
   bool CheckPrimalConsistency(PrimalSolutionStorage::Element primal) const final
   { 
      return CheckPrimalConsistencyImpl(primal);
   }

   // do zrobienia: this does not capture write back functions not returning REAL&
   constexpr static bool IsAssignableLeft() {
      return FunctionExistence::IsAssignable<typename LeftFactorContainer::RepamStorageType,REAL,INDEX>();
   }
   constexpr static bool IsAssignableRight() {
      return FunctionExistence::IsAssignable<typename RightFactorContainer::RepamStorageType,REAL,INDEX>();
   }

   template<typename ARRAY, bool IsAssignable = IsAssignableLeft()>
   constexpr static bool CanBatchRepamLeft()
   {
      return FunctionExistence::HasRepamLeft<MessageType,void,LeftFactorContainer,ARRAY>();
   }
   template<typename ARRAY, bool IsAssignable = IsAssignableLeft()>
   typename std::enable_if<CanBatchRepamLeft<ARRAY>() == true && IsAssignable == true>::type
   RepamLeft(const ARRAY& m)
   { 
      assert(false); // no -+ distinguishing
      msg_op_.RepamLeft(*leftFactor_, m);
   }
   template<typename ARRAY, bool IsAssignable = IsAssignableLeft()>
   typename std::enable_if<CanBatchRepamLeft<ARRAY>() == false && IsAssignable == true>::type
   RepamLeft(const ARRAY& m)
   { 
      assert(false); // no -+ distinguishing
      assert(m.size() == this->size());
      for(INDEX i=0; i<m.size(); ++i) {
         msg_op_.RepamLeft(*leftFactor_, m[i], i);
      }
   }
   template<typename ARRAY, bool IsAssignable = IsAssignableLeft()>
   typename std::enable_if<IsAssignable == false>::type
   RepamLeft(const ARRAY& m)
   {}

   template<bool IsAssignable = IsAssignableLeft()>
   typename std::enable_if<IsAssignable == true>::type
   RepamLeft(const REAL diff, const INDEX dim) {
      msg_op_.RepamLeft(*leftFactor_, diff, dim); // note: in right, we reparametrize by +diff, here by -diff
   }
   template<bool IsAssignable = IsAssignableLeft()>
   typename std::enable_if<IsAssignable == false>::type
   RepamLeft(const REAL diff, const INDEX dim)
   {}

   template<typename ARRAY>
   constexpr static bool CanBatchRepamRight()
   {
      return FunctionExistence::HasRepamRight<MessageType,void,RightFactorContainer,ARRAY>();
   }
   template<typename ARRAY, bool IsAssignable = IsAssignableRight()>
   typename std::enable_if<CanBatchRepamRight<ARRAY>() == true && IsAssignable == true>::type
   RepamRight(const ARRAY& m)
   { 
      assert(false); // no -+ distinguishing
      msg_op_.RepamRight(*rightFactor_, m);
   }
   template<typename ARRAY, bool IsAssignable = IsAssignableRight()>
   typename std::enable_if<CanBatchRepamRight<ARRAY>() == false && IsAssignable == true>::type
   RepamRight(const ARRAY& m)
   {
      assert(false); // no -+ distinguishing
      assert(m.size() == this->size());
      for(INDEX i=0; i<m.size(); ++i) {
         msg_op_.RepamRight(*rightFactor_, m[i], i);
      }
   }
   template<typename ARRAY, bool IsAssignable = IsAssignableRight()>
   typename std::enable_if<IsAssignable == false>::type
   RepamRight(const ARRAY& m)
   {}

   template<bool IsAssignable = IsAssignableRight()>
   typename std::enable_if<IsAssignable == true>::type
   RepamRight(const REAL diff, const INDEX dim) {
      msg_op_.RepamRight(*rightFactor_, diff, dim);
   }
   template<bool IsAssignable = IsAssignableRight()>
   typename std::enable_if<IsAssignable == false>::type
   RepamRight(const REAL diff, const INDEX dim)
   {}

   // do zrobienia: better name?
   REAL GetLeftMessage(const INDEX i) const { return msg_op_.GetLeftMessage(i,*this); }
   REAL GetRightMessage(const INDEX i) const { return msg_op_.GetRightMessage(i,*this);  }

   //FactorTypeAdapter* GetLeftFactor() const { return leftFactor_; }
   //FactorTypeAdapter* GetRightFactor() const { return rightFactor_; }
   // do zrobienia: Rename Get{Left|Right}FactorContainer
   LeftFactorContainer* GetLeftFactor() const final { return leftFactor_; }
   RightFactorContainer* GetRightFactor() const final { return rightFactor_; }

   INDEX GetMessageNumber() const final { return MESSAGE_NO; } 
   //REAL GetMessageWeightToRight() const final { return SEND_MESSAGE_TO_RIGHT_WEIGHT::value; }
   //REAL GetMessageWeightToLeft() const final { return SEND_MESSAGE_TO_LEFT_WEIGHT::value;  }
   
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
            static_cast<typename MessageContainerType::MessageStorageType*>(msg_)->operator[](dim_) -= x;
            msg_->RepamLeft( +x, dim_);
            msg_->RepamRight(-x, dim_);
         } else if (CHIRALITY == Chirality::left) { // message is computed by left factor
            static_cast<typename MessageContainerType::MessageStorageType*>(msg_)->operator[](dim_) += x;
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
            static_cast<typename MessageContainerType::MessageStorageType*>(msg_)->operator[](dim_) += x;
            msg_->RepamLeft( x, dim_);
            msg_->RepamRight( x, dim_);
         } else if(CHIRALITY == Chirality::left) { // message is computed by left factor
            static_cast<typename MessageContainerType::MessageStorageType*>(msg_)->operator[](dim_) -= x;
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
      operator REAL() const __attribute__ ((always_inline)) { return static_cast<typename MessageContainerType::MessageStorageType*>(msg_)->operator[](dim_); }
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
   };

   // for primal computation: record message change only in one side and into a special array
   template<Chirality CHIRALITY>
   class RestrictedMsgVal
   {
   public:
      RestrictedMsgVal(MessageContainerType* msg, const INDEX dim) : 
         msg_(msg), 
         dim_(dim)
      {}

      RestrictedMsgVal& operator-=(const REAL x) __attribute__ ((always_inline))
      {
         if(CHIRALITY == Chirality::right) { // message is computed by right factor
            msg_->RepamRight(-x, dim_);
         } else if (CHIRALITY == Chirality::left) { // message is computed by left factor
            msg_->RepamLeft(-x, dim_);
         } else {
            assert(false);
         }
         return *this;
      }

      RestrictedMsgVal& operator+=(const REAL x) __attribute__ ((always_inline))
      {
         assert(false);
         if(CHIRALITY == Chirality::right) {
            msg_->RepamRight(+x, dim_);
         } else if(CHIRALITY == Chirality::left) {
            msg_->RepamLeft(+x, dim_);
         } else {
            assert(false);
         }
         return *this;
      }

   private:
      MessageContainerType* const msg_;
      const INDEX dim_;
   };

   // this view is given to receive restricted message operations. 
   // Reparametrization is recorded only on one side
   template<Chirality CHIRALITY>
   class RestrictedMessageContainerView : public MessageContainerType{
   public:
      //using MessageContainerType;
      RestrictedMsgVal<CHIRALITY> operator[](const INDEX i) 
      {
         return RestrictedMsgVal<CHIRALITY>(this,i);
      }
   };


   // there must be four different implementations of msg updating with SIMD: 
   // (i) If parallel reparametrization is not supported by either left and right factor
   // If (ii) left or (iii) right factor supports reparametrization but not both
   // If (iv) both left and right factor support reparametrization


   template<typename ARRAY>
   MessageContainerType& operator=(const ARRAY& msg) {
      // construct difference to current message and then call +=
      // better do this via expression templates
      // diff = msg - *this;
      MinusExprVec<ARRAY,decltype(*this)> diff(msg, *this);

      //std::vector<REAL> diff(msg.size());
      //for(INDEX i=0; i<diff.size(); ++i) {
      //   diff[i] = msg[i] - operator[](i);
      //}
      operator+=(diff);
      return *this;
   }

   template<typename ARRAY>
   MessageContainerType& operator-=(const ARRAY& diff) {
      assert(false); // update to left right -+
      MinusVec<ARRAY> minus_diff(diff);
      assert(minus_diff.size() == this->size());
      RepamLeft(minus_diff);
      RepamRight(minus_diff);
      return *this;
   }

   template<typename ARRAY>
   MessageContainerType& operator+=(const ARRAY& diff) {
      assert(false); // update to left right -+
      PlusVec<ARRAY> plus_diff(diff); // used to wrap Vc::Memory, otherwise not needed // do zrobienia: change this with better vector architecture
      assert(plus_diff.size() == this->size()); // or entriesCount
      RepamLeft(plus_diff);
      RepamRight(plus_diff);
      return *this;
   }


   /*
   template<typename LAMBDA>
      struct repamOp {
         repamOp(const std::vector<REAL>& msg, LAMBDA& f) : f_(f) {}
         inline REAL operator[](const INDEX i) const { return f_(i); }
         LAMBDA& f_;
      };
      */
   // do zrobienia: change from valarray to std::vector and name SetMessageVal
   // function obsolete in current form
   void SetMessage(const std::valarray<REAL>& m) final
   { 
      assert(m.size() == MessageStorageType::size());
      assert(false); // not supported currently
      
      for(INDEX i=0; i<m.size(); ++i) { this->operator[](i) = m[i]; }

      /*
      // get difference to current message
      // check later if the construction below is faster
      //static std::vector<REAL> d(200);
      std::vector<REAL> d(m.size());
      for(INDEX i=0; i<m.size(); ++i) { d[i] = m[i] - msg_val_[i]; }

      // set message
      msg_val_ = m; 
     
      // update the reparametrization stored by the factor
      auto l = [&](const INDEX i) { return static_cast<MESSAGE_TYPE*>(this)->ApplyLeftRepam(i,d); };
      repamOp<decltype(l)> lOp(m,l);
      auto r = [&](const INDEX i) { return static_cast<MESSAGE_TYPE*>(this)->ApplyRightRepam(i,d); };
      repamOp<decltype(r)> rOp(m,r);

      leftFactor_->UpdateRepam(lOp);
      rightFactor_->UpdateRepam(rOp);
      */
   }
   // do zrobienia: change return type to std::vector, rename to GetMessageVal
   const std::valarray<REAL> GetMessage() const final
   { 
      assert(false); // do zrobienia: function not supported anymore. Only reparametrizations should be of interest
      std::valarray<REAL> m(0.0, MessageStorageType::size());
      for(INDEX i=0; i<MessageStorageType::size(); ++i) {
         m[i] = MessageStorageType::operator[](i);
      }
      return m; 
   }

   // possibly not the best choice: Sometimes msg_op_ needs access to this class
   const MessageType& GetMessageOp() const
   {
      return msg_op_;
   }

   constexpr static bool CanCreateConstraints()
   {
      //return FunctionExistence::HasCreateConstraints<MessageType,LpInterfaceAdapter*, LeftFactorContainer*, RightFactorContainer*>();
      return FunctionExistence::HasCreateConstraints<MessageType,void, LpInterfaceAdapter*, LeftFactorContainer*, RightFactorContainer*>();
   }
   
   template<bool ENABLE = CanCreateConstraints()>
   typename std::enable_if<ENABLE>::type
   CreateConstraintsImpl(LpInterfaceAdapter* l) const
   {
      msg_op_.CreateConstraints(l,leftFactor_->GetFactor(),rightFactor_->GetFactor());
   }

   template<bool ENABLE = CanCreateConstraints()>
   typename std::enable_if<!ENABLE>::type
   CreateConstraintsImpl(LpInterfaceAdapter* l) const
   {
      throw std::runtime_error("create constraints not implemented by message");
   }

   virtual void CreateConstraints(LpInterfaceAdapter* l) final
   {
      CreateConstraintsImpl(l);
   }


protected:
   MessageType msg_op_;
   LeftFactorContainer* const leftFactor_;
   RightFactorContainer* const rightFactor_;
};


// container class for factors. Here we hold the factor, all connected messages, reparametrization storage and perform reparametrization and coordination for sending and receiving messages.
// derives from REPAM_STORAGE_TYPE to mixin a class for storing the reparametrized potential
// implements the interface from FactorTypeAdapter for access from LP_MP
// if COMPUTE_PRIMAL_SOLUTION is true, MaximizePotential is expected to return either an integer of type INDEX or a std::vector<INDEX>
// if WRITE_PRIMAL_SOLUTION is false, WritePrimal will not output anything
// do zrobienia: introduce enum classes for COMPUTE_PRIMAL_SOLUTION and WRITE_PRIMAL_SOLUTION
template<typename FACTOR_TYPE, 
         template<class> class REPAM_STORAGE_TYPE, 
         class FACTOR_MESSAGE_TRAIT,
         INDEX FACTOR_NO,
         bool COMPUTE_PRIMAL_SOLUTION = false,
         bool WRITE_PRIMAL_SOLUTION = false> // do zrobienia: remove this template parameter
class FactorContainer : public REPAM_STORAGE_TYPE<FactorContainer<FACTOR_TYPE, REPAM_STORAGE_TYPE, FACTOR_MESSAGE_TRAIT, FACTOR_NO, COMPUTE_PRIMAL_SOLUTION, WRITE_PRIMAL_SOLUTION> >, public FactorTypeAdapter
{
public:
   using FactorContainerType = FactorContainer<FACTOR_TYPE, REPAM_STORAGE_TYPE, FACTOR_MESSAGE_TRAIT, FACTOR_NO, COMPUTE_PRIMAL_SOLUTION, WRITE_PRIMAL_SOLUTION>;
   using FactorType = FACTOR_TYPE;
   using RepamStorageType = REPAM_STORAGE_TYPE<FactorContainerType>;
   friend class REPAM_STORAGE_TYPE<FactorContainerType>;

   // do zrobienia: templatize cosntructor to allow for more general initialization of reparametrization storage and factor
   template<typename ...ARGS>
   FactorContainer(const FactorType& factor, ARGS... args) : RepamStorageType(factor,args...), factor_(factor) {
      //INDEX status;
      //std::cout << "msg_ type= "  << abi::__cxa_demangle(typeid(msg_).name(),0,0,&status) << "\n";
      //std::cout << "dispatcher list = "  << abi::__cxa_demangle(typeid(MESSAGE_DISPATCHER_TYPELIST).name(),0,0,&status) << "\n";
      //std::cout << "msg_ type= "  << abi::__cxa_demangle(typeid(msg_).name(),0,0,&status) << "\n";
      //std::cout << "left message list = " << abi::__cxa_demangle(typeid(left_message_list_).name(),0,0,&status) << "\n";
      //std::cout << "left message list = " << abi::__cxa_demangle(typeid(left_message_list_1).name(),0,0,&status) << "\n";
   
   }
   ~FactorContainer() { 
      static_assert(meta::unique<MESSAGE_DISPATCHER_TYPELIST>::size() == MESSAGE_DISPATCHER_TYPELIST::size(), 
            "Message dispatcher typelist must have unique elements");
      static_assert(FACTOR_NO >= 0 && FACTOR_NO < FACTOR_MESSAGE_TRAIT::FactorList::size(), "factor number must be smaller than length of factor list");
   }

   template<typename MESSAGE_DISPATCHER_TYPE, typename MESSAGE_TYPE> 
   void AddMessage(MESSAGE_TYPE* m) { 
      constexpr INDEX n = FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      static_assert( n < meta::size<MESSAGE_DISPATCHER_TYPELIST>(), "message dispatcher not supported");
      static_assert( n < std::tuple_size<decltype(msg_)>(), "message dispatcher not supported");
      //INDEX status;
      //std::cout << "msg dispatcher list =\n" << abi::__cxa_demangle(typeid(MESSAGE_DISPATCHER_TYPELIST).name(),0,0,&status) << "\n";
      //std::cout << "dispatcher  type =\n" << abi::__cxa_demangle(typeid(MESSAGE_DISPATCHER_TYPE).name(),0,0,&status) << "\n";
      //std::cout << " number = " << n << "\n" ;
      //std::cout << "message type = " << abi::__cxa_demangle(typeid(MESSAGE_TYPE).name(),0,0,&status) << "\n";

      std::get<n>(msg_).push_back(m);
   }

   // get sum of all messages in dimension i
   const REAL GetMessageSum(const INDEX i) const {
      return GetMessageSum(MESSAGE_DISPATCHER_TYPELIST{},i);
   }
   template<typename... MESSAGE_DISPATCHER_TYPES_REST>
   const REAL GetMessageSum(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>, const INDEX i) const { return 0.0; }
   template<typename MESSAGE_DISPATCHER_TYPE, typename... MESSAGE_DISPATCHER_TYPES_REST>
   const REAL GetMessageSum(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...>, const INDEX i) const {
      //INDEX status;
      //std::cout << "return message for " << abi::__cxa_demangle(typeid(MESSAGE_DISPATCHER_TYPE).name(),0,0,&status) << "\n";
      // current number of MESSAGE_DISPATCHER_TYPE
      constexpr INDEX n = FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      REAL msg_val = 0.0;
      for(auto it=std::get<n>(msg_).begin(); it!=std::get<n>(msg_).end(); ++it) {
         msg_val += MESSAGE_DISPATCHER_TYPE::GetMessage(*(*it),i);
      }
      // receive messages for subsequent MESSAGE_DISPATCER_TYPES
      return msg_val + GetMessageSum(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{},i);
   }

   void UpdateFactor(const std::vector<REAL>& omega) final
   {
      ReceiveMessages(omega);
      MaximizePotential();
      SendMessages(omega);
   }

   // do zrobienia: possibly also check if method present
   constexpr static bool
   CanComputePrimal()
   {
      return COMPUTE_PRIMAL_SOLUTION;
   }

   constexpr static bool
   CanPropagatePrimal()
   {
      return FunctionExistence::HasPropagatePrimal<FactorType,void,PrimalSolutionStorage::Element>();
   }

   template<bool ENABLE=CanPropagatePrimal()>
   typename std::enable_if<!ENABLE,void>::type
   PropagatePrimalImpl(PrimalSolutionStorage::Element primal) 
   {}
   template<bool ENABLE=CanPropagatePrimal()>
   typename std::enable_if<ENABLE,void>::type
   PropagatePrimalImpl(PrimalSolutionStorage::Element primal) 
   {
      factor_.PropagatePrimal(primal);
   }

   void PropagatePrimal(PrimalSolutionStorage::Element primal) 
   {
      PropagatePrimalImpl(primal);
   }

   constexpr static bool
   CanMaximizePotential()
   {
      return FunctionExistence::HasMaximizePotential<FactorType,void,FactorContainerType>();
   }

   // do zrobienia: remove this function, obsolete: New primal computation mode
   void UpdateFactor(const std::vector<REAL>& omega, typename PrimalSolutionStorage::Element primal) final
   {

      if(CanComputePrimal()) { // do zrobienia: for now
         if(CanReceiveRestrictedMessages()) {
            std::vector<REAL> tmpRepam(this->size()); // temporary structure where repam is stored before it is reverted back.
            for(INDEX i=0; i<tmpRepam.size(); ++i) {
               tmpRepam[i] = this->operator[](i);
            }
            // first we compute restricted incoming messages, on which to compute the primal
            ReceiveRestrictedMessages(primal);
            // now we compute primal
            MaximizePotentialAndComputePrimal(primal);
            // restore original reparametrization
            for(INDEX i=0; i<tmpRepam.size(); ++i) {
               this->operator[](i) = tmpRepam[i];
            }
            MaximizePotential();
         } else {
            MaximizePotentialAndComputePrimal(primal);
         }
         // now prapagate primal to adjacent factors
         ComputePrimalThroughMessages(primal);
      } else {
         MaximizePotential();
      } 

      ReceiveMessages(omega);
      SendMessages(omega);
   }

   template<bool ENABLE=CanMaximizePotential()>
   typename std::enable_if<ENABLE,void>::type
   MaximizePotential() 
   {
      factor_.MaximizePotential(*this);
   }
   template<bool ENABLE=CanMaximizePotential()>
   typename std::enable_if<!ENABLE,void>::type
   MaximizePotential() {}

   // if primal solution is to be computed by this factor, then we must take primal solution from some other factor and derive the solution through the messages
   template<bool COMPUTE_PRIMAL_SOLUTION_TMP = COMPUTE_PRIMAL_SOLUTION>
   typename std::enable_if<COMPUTE_PRIMAL_SOLUTION_TMP == false>::type 
   MaximizePotentialAndComputePrimal(typename PrimalSolutionStorage::Element primal)
   {
      assert(false); // this should not occur at all
   }

   // if primal solution is to be computed by this factor, we store the solution in primal
   template<bool COMPUTE_PRIMAL_SOLUTION_TMP = COMPUTE_PRIMAL_SOLUTION>
   typename std::enable_if<COMPUTE_PRIMAL_SOLUTION_TMP == true>::type 
   MaximizePotentialAndComputePrimal(typename PrimalSolutionStorage::Element primal)
   {
      static_assert(COMPUTE_PRIMAL_SOLUTION_TMP == COMPUTE_PRIMAL_SOLUTION,"");
      factor_.MaximizePotentialAndComputePrimal(*this, primal + primalOffset_);
   }


   template<typename MESSAGE_DISPATCHER_TYPE>
   typename std::enable_if<MESSAGE_DISPATCHER_TYPE::CanComputePrimalThroughMessage() == true>::type 
   ComputePrimalThroughMessagesImpl(MESSAGE_DISPATCHER_TYPE, typename PrimalSolutionStorage::Element primal) const
   {
      constexpr INDEX n = FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      for(INDEX i=0; i<std::get<n>(msg_).size(); ++i) {
         MESSAGE_DISPATCHER_TYPE::ComputePrimalThroughMessage(*(std::get<n>(msg_)[i]), primal);
      }
      //for(auto it=std::get<n>(msg_).cbegin(); it!=std::get<n>(msg_).cend(); ++it) {
      //   MESSAGE_DISPATCHER_TYPE::ComputePrimalThroughMessage(*(*it), primal);
      //}
   }
   template<typename MESSAGE_DISPATCHER_TYPE>
   typename std::enable_if<MESSAGE_DISPATCHER_TYPE::CanComputePrimalThroughMessage() == false>::type 
   ComputePrimalThroughMessagesImpl(MESSAGE_DISPATCHER_TYPE, typename PrimalSolutionStorage::Element primal) const {}

   template<typename... MESSAGE_DISPATCHER_TYPES_REST>
   void ComputePrimalThroughMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>, typename PrimalSolutionStorage::Element primal) const {}
   template<typename MESSAGE_DISPATCHER_TYPE, typename... MESSAGE_DISPATCHER_TYPES_REST>
   void ComputePrimalThroughMessages(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...>, typename PrimalSolutionStorage::Element primal) const 
   {
      ComputePrimalThroughMessagesImpl(MESSAGE_DISPATCHER_TYPE{}, primal);
      ComputePrimalThroughMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{}, primal);
   }
   // do zrobienia: rename PropagatePrimalThroughMessages
   void ComputePrimalThroughMessages(typename PrimalSolutionStorage::Element primal) const
   {
      ComputePrimalThroughMessages(MESSAGE_DISPATCHER_TYPELIST{}, primal);
   }
   
   

   // SFINAE-based selection whether we will perform message updates for receiving   
   template<typename MESSAGE_DISPATCHER_TYPE, typename MSG_ARRAY, typename ITERATOR>
   typename std::enable_if<MESSAGE_DISPATCHER_TYPE::CanCallReceiveMessage() == true>::type 
   ReceiveMessagesImpl(MESSAGE_DISPATCHER_TYPE msg_dispatcher, const MSG_ARRAY& msgs, ITERATOR omegaIt)
   {
      // receive messages for current MESSAGE_DISPATCER_TYPE
      constexpr INDEX n = FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      // do zrobienia: note that msgs array is not used!
      for(INDEX i=0; i<std::get<n>(msg_).size(); ++i) {
         // this is not valid. Instead, use a vector of bools which indicates whether to receive messages and let it be computed by Compute...Weights
         //if(*omegaIt == 0.0) { // makes large difference for cosegmentation_bins, why?
         MESSAGE_DISPATCHER_TYPE::ReceiveMessage(*(std::get<n>(msg_)[i]));
         //}
      }
      //for(auto it=std::get<n>(msg_).cbegin(); it!=std::get<n>(msg_).cend(); ++it, ++omegaIt) {
      //   // this is not valid. Instead, use a vector of bools which indicates whether to receive messages and let it be computed by Compute...Weights
      //   //if(*omegaIt == 0.0) { // makes large difference for cosegmentation_bins, why?
      //      MESSAGE_DISPATCHER_TYPE::ReceiveMessage(*(*it));
      //   //}
      //}
   }
   // or do not perform receiving message updates (if no receive message is implemented)
   template<typename MESSAGE_DISPATCHER_TYPE, typename MSG_ARRAY, typename ITERATOR>
   typename std::enable_if<MESSAGE_DISPATCHER_TYPE::CanCallReceiveMessage() == false>::type 
   ReceiveMessagesImpl(MESSAGE_DISPATCHER_TYPE msg_dispatcher, const MSG_ARRAY& msgs, ITERATOR omegaIt)
   {}

   template<typename ITERATOR, typename... MESSAGE_DISPATCHER_TYPES_REST>
   void ReceiveMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>, ITERATOR omegaIt) {}
   template<typename ITERATOR, typename MESSAGE_DISPATCHER_TYPE, typename... MESSAGE_DISPATCHER_TYPES_REST>
   void ReceiveMessages(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...>, ITERATOR omegaIt) 
   {
      // receive messages for current MESSAGE_DISPATCER_TYPE
      constexpr INDEX n = FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      ReceiveMessagesImpl(MESSAGE_DISPATCHER_TYPE{}, std::get<n>(msg_), omegaIt);
      //omegaIt += std::get<n>(msg_).size(); 
      // receive messages for subsequent MESSAGE_DISPATCHER_TYPES
      ReceiveMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{}, omegaIt);
   }
   void ReceiveMessages(const std::vector<REAL>& omega) 
   {
      // note: currently all messages are received, even if not needed. Change this again.
      //assert(omega.size() == GetNoMessages());
      ReceiveMessages(MESSAGE_DISPATCHER_TYPELIST{}, omega.cbegin());
   }



   // SFINAE-based restricted message updates
   template<typename MESSAGE_DISPATCHER_TYPE>
   typename std::enable_if<MESSAGE_DISPATCHER_TYPE::CanCallReceiveRestrictedMessage() == true>::type 
   ReceiveRestrictedMessagesImpl(MESSAGE_DISPATCHER_TYPE msg_dispatcher, PrimalSolutionStorage::Element primal)
   {
      // receive restricted messages for current MESSAGE_DISPATCHER_TYPE
      constexpr INDEX n = FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      // do zrobienia: note that msgs array is not used!
      for(INDEX i=0; i<std::get<n>(msg_).size(); ++i) {
         MESSAGE_DISPATCHER_TYPE::ReceiveRestrictedMessage(*(std::get<n>(msg_)[i]), primal); // do zrobienia: only receive messages from sensible ones
      }
      //for(auto it=std::get<n>(msg_).cbegin(); it!=std::get<n>(msg_).cend(); ++it) { // do zrobienia: only receive messages from sensible ones
      //   MESSAGE_DISPATCHER_TYPE::ReceiveRestrictedMessage(*(*it),primal);
      //}
   }
   // or do not perform receiving message updates (if no receive restricted message is implemented)
   template<typename MESSAGE_DISPATCHER_TYPE>
   typename std::enable_if<MESSAGE_DISPATCHER_TYPE::CanCallReceiveRestrictedMessage() == false>::type 
   ReceiveRestrictedMessagesImpl(MESSAGE_DISPATCHER_TYPE msg_dispatcher, PrimalSolutionStorage::Element)
   {}

   // we write message change not into original reparametrization, but into temporary one named pot
   template<typename... MESSAGE_DISPATCHER_TYPES_REST>
   void ReceiveRestrictedMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>, PrimalSolutionStorage::Element) {}
   template<typename MESSAGE_DISPATCHER_TYPE, typename... MESSAGE_DISPATCHER_TYPES_REST>
   void ReceiveRestrictedMessages(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...>, PrimalSolutionStorage::Element primal) 
   {
      // receive messages for current MESSAGE_DISPATCER_TYPE
      ReceiveRestrictedMessagesImpl(MESSAGE_DISPATCHER_TYPE{}, primal);
      ReceiveRestrictedMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{}, primal);
   }

   void ReceiveRestrictedMessages(PrimalSolutionStorage::Element primal) 
   {
      ReceiveRestrictedMessages(MESSAGE_DISPATCHER_TYPELIST{}, primal);
   }


   void SendMessages(const std::vector<REAL>& omega) 
   {
      //assert(omega.size() == GetNoMessages()); // this is not true: omega.size() is the number of messages that implement a send function
      static constexpr INDEX n = NumberOfSendMessagesCalls<std::vector<REAL>, decltype(omega.begin())>(MESSAGE_DISPATCHER_TYPELIST{});
      // do zrobienia: also do not construct currentRepam, if exactly one message update call will be issued. 
      // Check if there is one message dispatcher such that its size can be called via a constexpr function and is 1 -> complicated!
      // also possible: check whether omega has only one nonnegative entry
      if( n > 0 ) { // no need to construct currentRepam, if it will not be used at all
         // make a copy of the current reparametrization. The new messages are computed on it. Messages are updated implicitly and hence possibly the new reparametrization is automatically adjusted, which would interfere with message updates
         // do zrobienia: use static memory or custom memory allocator for this, do not always allocate new memory via system call
         std::vector<REAL> currentRepam(RepamStorageType::size());
         for(INDEX i=0; i<currentRepam.size(); ++i) {
            currentRepam[i] = RepamStorageType::operator[](i); 
         }
         SendMessages(MESSAGE_DISPATCHER_TYPELIST{}, currentRepam, omega.cbegin());
      }
   }

   template<typename REPAM_ARRAY, typename ITERATOR, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   constexpr static INDEX NumberOfSendMessagesCalls(meta::list<MESSAGE_DISPATCHER_TYPES_REST...> t) 
   { return 0; }
   template<typename REPAM_ARRAY, typename ITERATOR, typename MESSAGE_DISPATCHER_TYPE, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   constexpr static INDEX NumberOfSendMessagesCalls(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...> t) 
   { 
      constexpr INDEX n = FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      constexpr INDEX no = MESSAGE_DISPATCHER_TYPE::template CanCallSendMessages<FactorType, REPAM_ARRAY, decltype(std::get<n>(msg_)), ITERATOR>()
         || MESSAGE_DISPATCHER_TYPE::template CanCallSendMessage<REPAM_ARRAY>()
         ? 1 : 0;
      return no + NumberOfSendMessagesCalls<REPAM_ARRAY,ITERATOR>(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{});
   }

   template<typename ...MESSAGE_DISPATCHER_TYPES_REST>
   constexpr static bool CanReceiveRestrictedMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...> t) 
   { return false; }
   template<typename MESSAGE_DISPATCHER_TYPE, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   constexpr static bool CanReceiveRestrictedMessages(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...> t) 
   {
      constexpr bool canReceive = MESSAGE_DISPATCHER_TYPE::CanCallReceiveRestrictedMessage();
      if(canReceive == true) {
         return true;
      } else {
         return CanReceiveRestrictedMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{});
      }
   }
   constexpr static bool CanReceiveRestrictedMessages() 
   {
      return CanReceiveRestrictedMessages(MESSAGE_DISPATCHER_TYPELIST{});
   }

   // SFINAE-based seletion whether we will do batch or individual message updates for sending
   // batch message update
   // do zrobienia: note that CanCallSendMessages depends on more tempalte arguments than CanCalSendMessage. Reduce template usage of first one by passign tempalte arguments later.
   template<typename MESSAGE_DISPATCHER_TYPE, typename ITERATOR, typename REPAM_ARRAY, typename MSG_ARRAY>
   typename std::enable_if<MESSAGE_DISPATCHER_TYPE::template CanCallSendMessages<FactorType, REPAM_ARRAY, MSG_ARRAY, ITERATOR>() == true,INDEX>::type 
   SendMessagesImpl(MESSAGE_DISPATCHER_TYPE msg_dispatcher, const MSG_ARRAY& msgs, const REPAM_ARRAY& repam, ITERATOR omegaIt)
   {
      constexpr INDEX n = FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      const REAL omega_sum = std::accumulate(omegaIt, omegaIt + std::get<n>(msg_).size(), 0.0);
      if(omega_sum > 0.0) { // do zrobienia: possibly not allowed with MessageReplicatorFactor
         // do zrobienia: construct proxy object for msgs, so that it directly points to &(msgs[i]->msg_op_), make msg_op_ protected in MessageContainer again
                  //MsgProxy<MSG_ARRAY> msgProxy({msgs});
               
         /*
         struct MsgProxy {
            MsgProxy(const MSG_ARRAY& msgs) : msgs_(msgs) {}
            decltype(*(std::declval<const MSG_ARRAY&>().operator[](0)))& operator[](const INDEX i) const { return *(msgs_[i]); }
            INDEX size() const { return msgs_.size(); }
            private:
            const MSG_ARRAY& msgs_;
         } msgProxy(msgs);
         MESSAGE_DISPATCHER_TYPE::SendMessages(factor_, repam, msgProxy, omegaIt);
         */

         MESSAGE_DISPATCHER_TYPE::SendMessages(factor_, repam, msgs, omegaIt);
      }
      return msgs.size();
   }
   // individual message update
   template<typename MESSAGE_DISPATCHER_TYPE, typename ITERATOR, typename REPAM_ARRAY, typename MSG_ARRAY>
   typename std::enable_if<
   MESSAGE_DISPATCHER_TYPE::template CanCallSendMessages<FactorType,REPAM_ARRAY,MSG_ARRAY,ITERATOR>() == false && 
   MESSAGE_DISPATCHER_TYPE::template CanCallSendMessage <REPAM_ARRAY>() == true
   ,INDEX>::type 
   SendMessagesImpl(MESSAGE_DISPATCHER_TYPE msg_dispatcher, const MSG_ARRAY& msgs, const REPAM_ARRAY& repam, ITERATOR omegaIt)
   {
      // call individual message updates
      for(INDEX i=0; i<msgs.size(); ++i, ++omegaIt) {
         if(*omegaIt != 0.0) {
            MESSAGE_DISPATCHER_TYPE::SendMessage(*msgs[i], repam, *omegaIt); // do zrobienia: only receive messages from sensible ones
         }
      }
      //for(auto it=msgs.cbegin(); it!=msgs.cend(); ++it, ++omegaIt) {
      //   if(*omegaIt != 0.0) {
      //      MESSAGE_DISPATCHER_TYPE::SendMessage(*(*it), repam, *omegaIt);
      //   }
      //}
      return msgs.size(); 
   }
   // no updates if they are not implemented 
   template<typename MESSAGE_DISPATCHER_TYPE, typename ITERATOR, typename REPAM_ARRAY, typename MSG_ARRAY>
   typename std::enable_if<
   MESSAGE_DISPATCHER_TYPE::template CanCallSendMessages<FactorType,REPAM_ARRAY,MSG_ARRAY,ITERATOR>() == false &&
   MESSAGE_DISPATCHER_TYPE::template CanCallSendMessage <REPAM_ARRAY>() == false
   ,INDEX>::type 
   SendMessagesImpl(MESSAGE_DISPATCHER_TYPE msg_dispatcher, const MSG_ARRAY& msgs, const REPAM_ARRAY& repam, ITERATOR omegaIt)
   { return 0;}

   // note that messages must be iterated over in the same order as done by MessageIterator
   template<typename ITERATOR, typename ARRAY, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   void SendMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...> t, const ARRAY& repam, ITERATOR omegaIt) {}
   template<typename ITERATOR, typename ARRAY, typename MESSAGE_DISPATCHER_TYPE, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   void SendMessages(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...> t, const ARRAY& repam, ITERATOR omegaIt) // to get the current MESSAGE_TYPE
   { 
      // receive messages for current MESSAGE_DISPATCHER_TYPE
      constexpr INDEX n = FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();

      // check whether the message supports batch updates. If so, call batch update. If not, check whether individual updates are supported. If yes, call individual updates. If no, do nothing
      INDEX noCalls = SendMessagesImpl(MESSAGE_DISPATCHER_TYPE{}, std::get<n>(msg_), repam, omegaIt);
      omegaIt += noCalls;

      // receive messages for subsequent MESSAGE_DISPATCHER_TYPES
      SendMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{}, repam, omegaIt);
   }

   // methods used by MessageIterator
   template<typename ...MESSAGE_DISPATCHER_TYPES_REST>
   const INDEX GetNoMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...> t) const
   {
      return 0;
   }
   template<typename MESSAGE_DISPATCHER_TYPE, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   const INDEX GetNoMessages(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...> t) const 
   {
      constexpr INDEX n = FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      const INDEX no_msgs = std::get<n>(msg_).size(); 
      return no_msgs + GetNoMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{});
   }
   const INDEX GetNoMessages() const final
   {
      return GetNoMessages(MESSAGE_DISPATCHER_TYPELIST{});
   }

   template<typename ...MESSAGE_DISPATCHER_TYPES_REST>
   MessageTypeAdapter* GetMessage(meta::list<MESSAGE_DISPATCHER_TYPES_REST...> t, const INDEX msgNo) const 
   {
      assert(false);
      throw std::runtime_error("index out of bound");
      return nullptr;
   }
   template<typename MESSAGE_DISPATCHER_TYPE, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   MessageTypeAdapter* GetMessage(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...> t, const INDEX msgNo) const 
   {
      constexpr INDEX n = FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      if(msgNo < std::get<n>(msg_).size()) { return  std::get<n>(msg_)[msgNo]; }
      else return GetMessage(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{}, msgNo - std::get<n>(msg_).size());
   }
   MessageTypeAdapter* GetMessage(const INDEX n) const final
   {
      assert(n<GetNoMessages());
      return GetMessage(MESSAGE_DISPATCHER_TYPELIST{}, n);
   }

   template<typename ...MESSAGE_DISPATCHER_TYPES_REST>
   FactorTypeAdapter* GetConnectedFactor(meta::list<MESSAGE_DISPATCHER_TYPES_REST...> t, const INDEX cur_msg_idx) const 
   {
      throw std::runtime_error("message index out of bound");
   }
   template<typename MESSAGE_DISPATCHER_TYPE, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   FactorTypeAdapter* GetConnectedFactor(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...> t, const INDEX cur_msg_idx) const // to get the current message_type
   {
      constexpr INDEX n = FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      const INDEX no_msgs = std::get<n>(msg_).size();
      if(cur_msg_idx < no_msgs) {
         auto msg = std::get<n>(msg_)[cur_msg_idx];
         assert(msg != nullptr);
         if(msg->GetLeftFactor() == static_cast<const FactorTypeAdapter*>(this)) { return msg->GetRightFactor(); }
         else { return msg->GetLeftFactor(); }
      } else {
         return GetConnectedFactor(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{}, cur_msg_idx - no_msgs);
      }
   }
   FactorTypeAdapter* GetConnectedFactor (const INDEX msg_idx) const final
   { 
      auto f = GetConnectedFactor(MESSAGE_DISPATCHER_TYPELIST{}, msg_idx);
      assert(f != this);
      return f;
   }

   template<typename ...MESSAGE_DISPATCHER_TYPES_REST>
   constexpr static bool CanReceiveMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...> t) 
   {
      return false;
   }
   template<typename MESSAGE_DISPATCHER_TYPE, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   constexpr static bool CanReceiveMessages(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...> t)
   {
      constexpr INDEX n = FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      if( MESSAGE_DISPATCHER_TYPE::CanCallReceiveMessage()) {
         return true;
      } else {
         return CanReceiveMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{});
      }
   }
   constexpr static bool CanReceiveMessages() 
   {
      return CanReceiveMessages(MESSAGE_DISPATCHER_TYPELIST{});
   }
   // check whether actually receive restricted messages is called. Can be false, even if CanReceiveRestrictedMessages is true, e.g. when no message is present
   template<typename ...MESSAGE_DISPATCHER_TYPES_REST>
   bool CallsReceiveMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...> t) const
   { return false; }
   template<typename MESSAGE_DISPATCHER_TYPE, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   bool CallsReceiveMessages(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...> t)  const
   {
      constexpr bool canReceive = MESSAGE_DISPATCHER_TYPE::CanCallReceiveMessage();
      constexpr INDEX n = FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      if(canReceive == true && std::get<n>(msg_).size() > 0) {
         return true;
      } else {
         return CallsReceiveMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{});
      }
   }
   bool CallsReceiveMessages() const
   {
      return CallsReceiveMessages(MESSAGE_DISPATCHER_TYPELIST{});
   }


   template<typename ...MESSAGE_DISPATCHER_TYPES_REST>
   constexpr static bool CanSendMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...> t) 
   {
      return false;
   }
   template<typename MESSAGE_DISPATCHER_TYPE, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   constexpr static bool CanSendMessages(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...> t)
   {
      constexpr INDEX n = FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      if( MESSAGE_DISPATCHER_TYPE::template CanCallSendMessage<std::vector<REAL>>() ||
             MESSAGE_DISPATCHER_TYPE::template CanCallSendMessages<FactorContainerType, std::vector<REAL>, decltype(std::get<n>(msg_)), std::vector<REAL>::iterator>() ) {
         return true;
      } else {
         return CanSendMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{});
      }
   }
   constexpr static bool CanSendMessages() 
   {
      return CanSendMessages(MESSAGE_DISPATCHER_TYPELIST{});
   }

   // check whether actually send messages is called. Can be false, even if CanSendMessages is true, e.g. when no message is present
   template<typename ...MESSAGE_DISPATCHER_TYPES_REST>
   bool CallsSendMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...> t) const
   { return false; }
   template<typename MESSAGE_DISPATCHER_TYPE, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   bool CallsSendMessages(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...> t)  const
   {
      constexpr INDEX n = FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      constexpr bool canSend = MESSAGE_DISPATCHER_TYPE::template CanCallSendMessage<std::vector<REAL>>() || 
             MESSAGE_DISPATCHER_TYPE::template CanCallSendMessages<decltype(*this), std::vector<REAL>, decltype(std::get<n>(msg_)), std::vector<REAL>::iterator>();
      if(canSend == true && std::get<n>(msg_).size() > 0) {
         return true;
      } else {
         return CallsSendMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{});
      }
   }
   bool CallsSendMessages() const
   {
      return CallsSendMessages(MESSAGE_DISPATCHER_TYPELIST{});
   }


   template<typename ...MESSAGE_DISPATCHER_TYPES_REST>
   bool CanSendMessage(meta::list<MESSAGE_DISPATCHER_TYPES_REST...> t, const INDEX cur_msg_idx) const 
   {
      throw std::runtime_error("message index out of bound");
   }
   template<typename MESSAGE_DISPATCHER_TYPE, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   bool CanSendMessage(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...> t, const INDEX cur_msg_idx) const // to get the current MESSAGE_TYPE
   {
      constexpr INDEX n = FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      const INDEX no_msgs = std::get<n>(msg_).size();
      if(cur_msg_idx < no_msgs) {
         if( MESSAGE_DISPATCHER_TYPE::template CanCallSendMessage<std::vector<REAL>>() || 
             MESSAGE_DISPATCHER_TYPE::template CanCallSendMessages<decltype(*this), std::vector<REAL>, decltype(std::get<n>(msg_)), std::vector<REAL>::iterator>() )
            return true;
         else return false;
      } else {
         return CanSendMessage(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{}, cur_msg_idx - no_msgs);
      }
   }

   bool CanSendMessage(const INDEX msg_idx) const final
   {
      return CanSendMessage(MESSAGE_DISPATCHER_TYPELIST{}, msg_idx);
   }


   // check whether actually receive restricted messages is called. Can be false, even if CanReceiveRestrictedMessages is true, e.g. when no message is present
   template<typename ...MESSAGE_DISPATCHER_TYPES_REST>
   bool CallsReceiveRestrictedMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...> t) const
   { return false; }
   template<typename MESSAGE_DISPATCHER_TYPE, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   bool CallsReceiveRestrictedMessages(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...> t)  const
   {
      constexpr bool canReceive = MESSAGE_DISPATCHER_TYPE::CanCallReceiveRestrictedMessage();
      constexpr INDEX n = FindMessageDispatcherTypeIndex<MESSAGE_DISPATCHER_TYPE>();
      if(canReceive == true && std::get<n>(msg_).size() > 0) {
         return true;
      } else {
         return CallsReceiveRestrictedMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{});
      }
   }
   bool CallsReceiveRestrictedMessages() const
   {
      return CallsReceiveRestrictedMessages(MESSAGE_DISPATCHER_TYPELIST{});
   }


   // does factor call {Receive|Send}Messages? If not, it need, UpdateFactor need not be called
   bool FactorUpdated() const final
   {
      if(CanComputePrimal()) {
         return true;
      }
      if(CanReceiveMessages() && CallsReceiveMessages()) {
         return true;
      }
      if(CanSendMessages() && CallsSendMessages()) {
         return true;
      }
      if(CanReceiveRestrictedMessages() && CallsReceiveRestrictedMessages()) {
         return true;
      }
      return false;
   }

   template<typename ITERATOR>
   void SetAndPropagatePrimal(PrimalSolutionStorage::Element primal, ITERATOR label) const
   {
      for(INDEX i=0; i<PrimalSize(); ++i) {
         primal[i + GetPrimalOffset()] = label[i];
      }
      ComputePrimalThroughMessages(primal);
   }

   // do zrobienia: possibly do it with std::result_of
   //auto begin() -> decltype(std::declval<RepamStorageType>().begin()) { return RepamStorageType::begin(); }
   //auto end()   -> decltype(std::declval<RepamStorageType>().end()) { return RepamStorageType::end(); }
   //auto cbegin() -> decltype(std::declval<RepamStorageType>().cbegin()) const { return RepamStorageType::cbegin(); } // do zrobienia: somehow discards const qualifiers
   //auto cend()   -> decltype(std::declval<RepamStorageType>().cend()) const { return RepamStorageType::cend(); } 

   /*
   const INDEX size() const { return factor_.size(); }
   // get reparametrized cost
   const REAL operator[](const INDEX i) const { return RepamStorageType::operator[](i); }
   //auto operator[](const INDEX i) -> decltype(std::declval<RepamStorageType>().operator[](0)) { return RepamStorageType::operator[](i); }
   REAL& operator[](const INDEX i) { return RepamStorageType::operator[](i); }
   */

   std::vector<REAL> GetReparametrizedPotential() const final
   {
      std::vector<REAL> repam(size());
      for(INDEX i=0; i<repam.size(); ++i) {
         repam[i] = RepamStorageType::operator[](i);
      }
      return repam;
   }

   INDEX size() const final { return RepamStorageType::size(); }

   constexpr static bool CanComputePrimalSize()
   {
      return FunctionExistence::HasPrimalSize<FactorType,INDEX>();
   }
   template<bool ENABLE = CanComputePrimalSize()>
   typename std::enable_if<!ENABLE,INDEX>::type
   PrimalSizeImpl() const
   {
      return this->size();
   }
   template<bool ENABLE = CanComputePrimalSize()>
   typename std::enable_if<ENABLE,INDEX>::type
   PrimalSizeImpl() const
   {
      return factor_.PrimalSize();
   }
   // return size for primal storage
   INDEX PrimalSize() const final { return PrimalSizeImpl(); }

   REAL LowerBound() const final { return factor_.LowerBound(*this); } 

   FactorType* GetFactor() const { return &factor_; }
   FactorType* GetFactor() { return &factor_; }
   void SetPrimalOffset(const INDEX n) final { primalOffset_ = n; } // this function is used in AddFactor in LP class
   INDEX GetPrimalOffset() const final { return primalOffset_; }

  void SetAuxOffset(const INDEX n) final { auxOffset_ = n; }
  INDEX GetAuxOffset() const final { return auxOffset_; }
   
protected:
   FactorType factor_; // the factor operation
   INDEX primalOffset_;
   INDEX auxOffset_;
   

   // we do this via templates only (no values), as types are incomplete yet. This is more cumbersome than a direct value computation as can be done usually.
   /*
   template<typename MESSAGE_LIST>
   constexpr static auto MetaComputeMessageTuple()
   {
      //auto message_tuple = hana::transform(msg_list, [](auto m) { return std::declval<typename decltype(m)::type>(); }); // go from type_c to actual type
      //auto left_message_type_tuple = hana::filter(msg_list, [](auto m) { return std::declval<typename decltype(m)::type>().LeftFactorNumber() == FACTOR_NO; });
      auto msg_list_t = hana::type_c<MESSAGE_LIST>; // MESSAGE_LIST cannot be instantiated yet (incomplete type) -> hold it as type_c
      // transform message list to tuple_t holding all the tuples
      constexpr auto
      auto has_same_left_factor = hana::is_valid([](auto&& x) -> decltype((void)x.leftFactorNumber) { });
      //MESSAGE_LIST* msg_list;
      //auto left_message_type_tuple = hana::filter(*msg_list, [](auto&& m) { return decltype(m)::leftFactorNumber == FACTOR_NO; });
      //auto left_message_dispatcher_tuple = hana::transform(left_message_type_tuple, [](auto m) { return hana::type_c<MessageDispatcher<typename decltype(m)::type, LeftMessageFuncGetter>>; });
      //auto right_message_tuple = hana::filter(msg_list, [](auto m) { return decltype(m)::type::rightFactorNumber == FACTOR_NO; });
      //auto right_message_dispatcher_tuple = hana::transform(right_message_tuple, [](auto m) { return hana::type_c<MessageDispatcher<typename decltype(m)::type, RightMessageFuncGetter>>; });
      //return hana::concat<left_message_dispatcher_tuple, right_message_dispatcher_tuple>;

      //return left_message_dispatcher_tuple;
      return hana::type_c<char>;
   }
   */

   //decltype(MetaComputeMessageTuple<typename FACTOR_MESSAGE_TRAIT::MessageListHana>()) test;

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
      template<class LIST>
         using invoke = typename LIST::LeftMessageContainerStorageType;
   };
   struct get_right_msg_container_type_list {
      template<class LIST>
         using invoke = typename LIST::RightMessageContainerStorageType;
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

   // helper function for getting the index in msg_ of given MESSAGE_DISPATCHER_TYPE
   template<typename MESSAGE_DISPATCHER_TYPE>
   static constexpr INDEX FindMessageDispatcherTypeIndex()
   {
      constexpr INDEX n = meta::find_index<MESSAGE_DISPATCHER_TYPELIST, MESSAGE_DISPATCHER_TYPE>::value;
      static_assert(n < meta::size<MESSAGE_DISPATCHER_TYPELIST>::value,"");
      return n;
   }


   // construct tuple holding messages for left and right dispatch
   // the tuple will hold some container for the message type. The container type is specified in the {Left|Right}MessageContainerStorageType fields of MessageList
   using msg_container_type_list = meta::concat<left_msg_container_list, right_msg_container_list>;

   tuple_from_list<msg_container_type_list> msg_;

public:
   REAL EvaluatePrimal(typename PrimalSolutionStorage::Element primalIt) const final
   {
      return factor_.EvaluatePrimal(*this,primalIt + primalOffset_);
   }

   constexpr static bool CanCreateConstraints()
   {
      //return FunctionExistence::HasCreateConstraints<FactorType,LpInterfaceAdapter*>();
      return FunctionExistence::HasCreateConstraints<FactorType,void,LpInterfaceAdapter*>();
   }
   
   template<bool ENABLE = CanCreateConstraints()>
   typename std::enable_if<ENABLE>::type
   CreateConstraintsImpl(LpInterfaceAdapter* l) const
   {
      factor_.CreateConstraints(l);
   }
   
   INDEX GetNumberOfAuxVariables() const 
   { 
    return factor_.GetNumberOfAuxVariables();
   }

   template<bool ENABLE = CanCreateConstraints()>
   typename std::enable_if<!ENABLE>::type
   CreateConstraintsImpl(LpInterfaceAdapter* l) const
   {
      throw std::runtime_error("create constraints not implemented by factor");
   }


   void CreateConstraints(LpInterfaceAdapter* l) const final
   {
      CreateConstraintsImpl(l);
   }
};


} // end namespace LP_MP

#endif // LP_MP_FACTORS_MESSAGES_HXX

