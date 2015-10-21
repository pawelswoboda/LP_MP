#ifndef LP_MP_FACTORS_MESSAGES_HXX
#define LP_MP_FACTORS_MESSAGES_HXX

#include <vector>
#include <valarray>
#include <map>
#include <list>
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

#include "reparametrization_storage.hxx"
#include "messages/message_storage.hxx"

#include "LP_MP.h"

// this file provides message and factor containers. The factors and messages are plugged into the container and then every method call is dispatched correctly with static polymorphism and template tricks.


// do zrobienia: make those functions const which are logically so
// do zrobienia: Introduce MessageConstraint and FactorConstraint for templates
// cleanup name inconsistencies: MessageType, MessageDispatcher etc
// rename MultiplexMargMessage into MultiplexMsg, EqualityMessage into EqualityMsg
namespace LP_MP {

// we must check existence of functions in message classes. The necessary test code is concentrated here. 
namespace FunctionExistence {

// Macros to construct help functions for checking member functions of classes
LP_MP_FUNCTION_EXISTENCE_CLASS(HasReceiveMessageFromRight,ReceiveMessageFromRight);
LP_MP_FUNCTION_EXISTENCE_CLASS(HasReceiveMessageFromLeft, ReceiveMessageFromLeft);
   
LP_MP_FUNCTION_EXISTENCE_CLASS(HasSendMessageToRight,SendMessageToRight);
LP_MP_FUNCTION_EXISTENCE_CLASS(HasSendMessageToLeft, SendMessageToLeft);

LP_MP_FUNCTION_EXISTENCE_CLASS(HasSendMessagesToRight,SendMessagesToRight);
LP_MP_FUNCTION_EXISTENCE_CLASS(HasSendMessagesToLeft, SendMessagesToLeft);

LP_MP_FUNCTION_EXISTENCE_CLASS(HasRepamRight, RepamRight);
LP_MP_FUNCTION_EXISTENCE_CLASS(HasRepamLeft, RepamLeft);

LP_MP_FUNCTION_EXISTENCE_CLASS(IsAssignable, operator[]);
}

// function getters for statically dispatching ReceiveMessage and SendMessage to left and right side correctly, used in FactorContainer
// do zrobienia: are left and right correctly dispatched?
template<typename MSG_CONTAINER>
struct LeftMessageFuncGetter
{
   constexpr static decltype(&MSG_CONTAINER::GetLeftMessageContainer) GetMessageFunc() { return &MSG_CONTAINER::GetLeftMessageContainer; }

   constexpr static decltype(&MSG_CONTAINER::ReceiveMessageFromRightContainer) GetReceiveFunc() { return &MSG_CONTAINER::ReceiveMessageFromRightContainer; }
   template<typename ARRAY>
   constexpr static decltype(&MSG_CONTAINER::template SendMessageToRightContainer<ARRAY>) GetSendFunc() { return &MSG_CONTAINER::template SendMessageToRightContainer<ARRAY>; }

   template<typename MSG_ARRAY, typename REPAM_ARRAY, typename ITERATOR>
   constexpr static decltype(&MSG_CONTAINER::template SendMessagesToRightContainer<MSG_ARRAY, REPAM_ARRAY, ITERATOR>) GetSendMessagesFunc() 
   { return &MSG_CONTAINER::template SendMessagesToRightContainer<MSG_ARRAY,REPAM_ARRAY,ITERATOR>; }

   constexpr static bool 
   CanCallReceiveMessage()
   { return MSG_CONTAINER::CanCallReceiveMessageFromRightContainer(); }

   template<typename REPAM_ARRAY>
   constexpr static bool CanCallSendMessage() 
   { return MSG_CONTAINER::template CanCallSendMessageToRightContainer<REPAM_ARRAY>(); }

   template<typename MSG_ARRAY, typename REPAM_ARRAY, typename ITERATOR>
   constexpr static bool 
   CanCallSendMessages()
   { return MSG_CONTAINER::template CanCallSendMessagesToRightContainer<MSG_ARRAY, REPAM_ARRAY, ITERATOR>(); }
};

template<typename MSG_CONTAINER>
struct RightMessageFuncGetter
{
   constexpr static decltype(&MSG_CONTAINER::GetRightMessageContainer) GetMessageFunc() { return &MSG_CONTAINER::GetRightMessageContainer; }

   constexpr static decltype(&MSG_CONTAINER::ReceiveMessageFromLeftContainer) GetReceiveFunc() { return &MSG_CONTAINER::ReceiveMessageFromLeftContainer; }
   template<typename ARRAY>
   constexpr static decltype(&MSG_CONTAINER::template SendMessageToLeftContainer<ARRAY>) GetSendFunc() { return &MSG_CONTAINER::template SendMessageToLeftContainer<ARRAY>; }

   template<typename MSG_ARRAY, typename REPAM_ARRAY, typename ITERATOR>
   constexpr static decltype(&MSG_CONTAINER::template SendMessagesToLeftContainer<MSG_ARRAY, REPAM_ARRAY, ITERATOR>) GetSendMessagesFunc() 
      { return &MSG_CONTAINER::template SendMessagesToLeftContainer<MSG_ARRAY,REPAM_ARRAY,ITERATOR>; }

   constexpr static bool CanCallReceiveMessage() 
   { return MSG_CONTAINER::CanCallReceiveMessageFromLeftContainer(); }

   template<typename REPAM_ARRAY>
   constexpr static bool CanCallSendMessage() 
   { return MSG_CONTAINER::template CanCallSendMessageToLeftContainer<REPAM_ARRAY>(); }

   template<typename MSG_ARRAY, typename REPAM_ARRAY, typename ITERATOR>
   constexpr static bool
   CanCallSendMessages()
   { return MSG_CONTAINER::template CanCallSendMessagesToLeftContainer<MSG_ARRAY, REPAM_ARRAY, ITERATOR>(); }
};

template<class MSG_CONTAINER, template<typename> class FuncGetter>
struct MessageDispatcher
{
   constexpr static bool CanCallReceiveMessage() { return FuncGetter<MSG_CONTAINER>::CanCallReceiveMessage(); }
   static void ReceiveMessage(MSG_CONTAINER& t)
   {
      auto staticMemberFunc = FuncGetter<MSG_CONTAINER>::GetReceiveFunc();
      return (t.*staticMemberFunc)();
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

   // batch based message sending
   template<typename MSG_ARRAY, typename REPAM_ARRAY, typename ITERATOR>
   constexpr static bool CanCallSendMessages() { return FuncGetter<MSG_CONTAINER>::template CanCallSendMessages<MSG_ARRAY, REPAM_ARRAY, ITERATOR>(); }

   template<typename MSG_ARRAY, typename REPAM_ARRAY, typename ITERATOR>
   static void SendMessages(const MSG_ARRAY& msgs, const REPAM_ARRAY& repam, ITERATOR omegaBegin)
   {
      auto staticMemberFunc = FuncGetter<MSG_CONTAINER>::template GetSendMessagesFunc<MSG_ARRAY, REPAM_ARRAY, ITERATOR>();
      (*staticMemberFunc)(msgs, repam, omegaBegin);
   }

   static REAL GetMessage(MSG_CONTAINER& t, const INDEX i)
   {
      auto staticMemberFunc = FuncGetter<MSG_CONTAINER>::GetMessageFunc();
      return (t.*staticMemberFunc)(i);
   }
};

// message container class
template<typename MESSAGE_TYPE, typename MESSAGE_STORAGE, typename FACTOR_MESSAGE_TRAIT, INDEX MESSAGE_NO>
class MessageContainer : public MESSAGE_STORAGE, public MessageTypeAdapter
{
   constexpr static INDEX left_factor_number = meta::at_c< meta::at_c<typename FACTOR_MESSAGE_TRAIT::msg_list, MESSAGE_NO>, 1>::value;
   constexpr static INDEX right_factor_number = meta::at_c< meta::at_c<typename FACTOR_MESSAGE_TRAIT::msg_list, MESSAGE_NO>, 2 >::value;

public:
   typedef MessageContainer<MESSAGE_TYPE, MESSAGE_STORAGE, FACTOR_MESSAGE_TRAIT, MESSAGE_NO> Container;
   typedef MESSAGE_TYPE MessageType;
   typedef MESSAGE_STORAGE MessageStorageType;

   // FactorContainer
   using LEFT_FACTOR_TYPE = meta::at_c<typename FACTOR_MESSAGE_TRAIT::factor_list, left_factor_number>;
   using RIGHT_FACTOR_TYPE = meta::at_c<typename FACTOR_MESSAGE_TRAIT::factor_list, right_factor_number>;
   // FactorType of factor held by FactorContainer
   // do zrobienia: clean up naming mess: WRITING_STYLE like this should only refer to templates. NormalClassName should refer to class names otherwise
   using LEFT_FACTOR_TYPE_OP = typename LEFT_FACTOR_TYPE::FactorType;
   using RIGHT_FACTOR_TYPE_OP = typename RIGHT_FACTOR_TYPE::FactorType;

   MessageContainer(MESSAGE_TYPE msg_op, LEFT_FACTOR_TYPE* const l, RIGHT_FACTOR_TYPE* const r, const INDEX msg_size) 
      : MESSAGE_STORAGE(msg_size),
      msg_op_(msg_op),
      leftFactor_(l), 
      rightFactor_(r) 
   {
      //INDEX status;
      //std::cout << "msg holding type = " << abi::__cxa_demangle(typeid(msg_op_).name(),0,0,&status) << "\n";
      //std::cout << "left factor number = " << left_factor_number << "\n";
      //std::cout << "right factor number = " << right_factor_number << "\n";
      //std::cout << "left factor type = " << abi::__cxa_demangle(typeid(LEFT_FACTOR_TYPE).name(),0,0,&status) << "\n";
      //std::cout << "right factor type = " << abi::__cxa_demangle(typeid(RIGHT_FACTOR_TYPE).name(),0,0,&status) << "\n";
      // register messages in factors
      leftFactor_->template AddMessage<MessageDispatcher<Container, LeftMessageFuncGetter>, Container>(this);
      rightFactor_->template AddMessage<MessageDispatcher<Container, RightMessageFuncGetter>, Container>(this);
   }
   ~MessageContainer() {
      static_assert(meta::unique<typename FACTOR_MESSAGE_TRAIT::msg_list>::size() == FACTOR_MESSAGE_TRAIT::msg_list::size(), 
            "Message list must have unique elements");
      static_assert(MESSAGE_NO >= 0 && MESSAGE_NO < FACTOR_MESSAGE_TRAIT::msg_list::size(), "message number must be smaller than length of message list");
      static_assert(left_factor_number < FACTOR_MESSAGE_TRAIT::factor_list::size(), "left factor number out of bound");
      static_assert(right_factor_number < FACTOR_MESSAGE_TRAIT::factor_list::size(), "right factor number out of bound");
   } // put message constraINDEX here, i.e. which methods MESSAGE_TYPE must minimally implement


   constexpr static bool 
   CanCallReceiveMessageFromRightContainer()
   { 
      return FunctionExistence::HasReceiveMessageFromRight<MessageType, void, 
      decltype(rightFactor_->GetFactor()), decltype(*rightFactor_), Container>(); 
   }
   void ReceiveMessageFromRightContainer()
   { msg_op_.ReceiveMessageFromRight(rightFactor_->GetFactor(),*rightFactor_, *this); }

   constexpr static bool 
   CanCallReceiveMessageFromLeftContainer()
   { 
      return FunctionExistence::HasReceiveMessageFromLeft<MessageType, void, 
      decltype(leftFactor_->GetFactor()), decltype(*leftFactor_), Container>(); 
   }
   void ReceiveMessageFromLeftContainer()
   { msg_op_.ReceiveMessageFromLeft(leftFactor_->GetFactor(), *leftFactor_, *this); }


   template<typename REPAM_ARRAY>
   constexpr static bool 
   CanCallSendMessageToRightContainer()
   { 
      return FunctionExistence::HasSendMessageToRight<MessageType, void, 
      decltype(leftFactor_->GetFactor()), decltype(rightFactor_->GetFactor()), REPAM_ARRAY, decltype(*rightFactor_), std::vector<REAL>, REAL>(); 
   }

   template<typename ARRAY>
   void SendMessageToRightContainer(const ARRAY& repam, const REAL omega)
   {
      msg_op_.SendMessageToRight(leftFactor_->GetFactor(), rightFactor_->GetFactor(), repam, *rightFactor_, *this, omega);
   }

   template<typename REPAM_ARRAY>
   constexpr static bool
   CanCallSendMessageToLeftContainer()
   { 
      return FunctionExistence::HasSendMessageToLeft<MessageType, void, 
      decltype(leftFactor_->GetFactor()), decltype(rightFactor_->GetFactor()), decltype(*leftFactor_), REPAM_ARRAY, std::vector<REAL>, REAL>(); 
   }

   template<typename ARRAY>
   void SendMessageToLeftContainer(const ARRAY& repam, const REAL omega)
   {
      msg_op_.SendMessageToLeft(leftFactor_->GetFactor(), rightFactor_->GetFactor(), *leftFactor_, repam, *this, omega);
   }

   template<typename MSG_ARRAY, typename REPAM_ARRAY, typename ITERATOR>
   constexpr static bool
   CanCallSendMessagesToLeftContainer()
   { 
      return FunctionExistence::HasSendMessagesToLeft<MessageType, void, MSG_ARRAY, REPAM_ARRAY, ITERATOR>();
   }
   template<typename MSG_ARRAY, typename REPAM_ARRAY, typename ITERATOR>
   static void SendMessagesToLeftContainer(const MSG_ARRAY& msgs, const REPAM_ARRAY& repam, ITERATOR omegaBegin) {
   MessageType::SendMessagesToLeft(msgs, repam, omegaBegin);
   }

   template<typename... ARGS>
   constexpr static bool
   CanCallSendMessagesToRightContainer()
   { 
      return FunctionExistence::HasSendMessagesToRight<MessageType,void, ARGS...>();
   }
   template<typename MSG_ARRAY, typename REPAM_ARRAY, typename ITERATOR>
   static void SendMessagesToRightContainer(const MSG_ARRAY& msgs, const REPAM_ARRAY& repam, ITERATOR omegaBegin) {
   MessageType::SendMessagesToRight(msgs, repam, omegaBegin);
   }

   constexpr static bool IsAssignableLeft() {
      return FunctionExistence::IsAssignable<typename LEFT_FACTOR_TYPE::RepamStorageType,REAL&,INDEX>();
   }
   constexpr static bool IsAssignableRight() {
      return FunctionExistence::IsAssignable<typename RIGHT_FACTOR_TYPE::RepamStorageType,REAL&,INDEX>();
   }

   template<typename ARRAY, bool IsAssignable = IsAssignableLeft()>
   constexpr static bool CanBatchRepamLeft()
   {
      return FunctionExistence::HasRepamLeft<MessageType,void,LEFT_FACTOR_TYPE,ARRAY>();
   }
   template<typename ARRAY, bool IsAssignable = IsAssignableLeft()>
   typename std::enable_if<CanBatchRepamLeft<ARRAY>() == true && IsAssignable == true>::type
   RepamLeft(const ARRAY& m)
   { 
      msg_op_.RepamLeft(*leftFactor_, m);
   }
   template<typename ARRAY, bool IsAssignable = IsAssignableLeft()>
   typename std::enable_if<CanBatchRepamLeft<ARRAY>() == false && IsAssignable == true>::type
   RepamLeft(const ARRAY& m)
   { 
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
      msg_op_.RepamLeft(*leftFactor_, diff, dim);
   }
   template<bool IsAssignable = IsAssignableLeft()>
   typename std::enable_if<IsAssignable == false>::type
   RepamLeft(const REAL diff, const INDEX dim)
   {}


   template<typename ARRAY>
   constexpr static bool CanBatchRepamRight()
   {
      return FunctionExistence::HasRepamRight<MessageType,void,RIGHT_FACTOR_TYPE,ARRAY>();
   }
   template<typename ARRAY, bool IsAssignable = IsAssignableRight()>
   typename std::enable_if<CanBatchRepamRight<ARRAY>() == true && IsAssignable == true>::type
   RepamRight(const ARRAY& m)
   { 
      msg_op_.RepamRight(*rightFactor_, m);
   }
   template<typename ARRAY, bool IsAssignable = IsAssignableRight()>
   typename std::enable_if<CanBatchRepamRight<ARRAY>() == false && IsAssignable == true>::type
   RepamRight(const ARRAY& m)
   {
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


   REAL GetLeftMessageContainer(const INDEX i) const { return msg_op_.GetLeftMessage(i,*this); }
   REAL GetRightMessageContainer(const INDEX i) const { return msg_op_.GetRightMessage(i,*this);  }

   FactorTypeAdapter* GetLeftFactor() const { return leftFactor_; }
   FactorTypeAdapter* GetRightFactor() const { return rightFactor_; }
   //LEFT_FACTOR_TYPE* GetLeftFactor() { return leftFactor_; }
   //RIGHT_FACTOR_TYPE* GetRightFactor() { return rightFactor_; }
   
   //const REAL operator[](const INDEX i) const { return MessageStorageType::operator[](i); }

   // class for storing a callback upon new assignment of message: update left and right factors
   class MsgVal {
   public:
      MsgVal(Container* msg, const INDEX dim) : 
         msg_(msg), 
         dim_(dim)
      {}
      MsgVal& operator=(const REAL x) __attribute__ ((always_inline))
      {
         const REAL diff = x - msg_->operator[](dim_);
         // set new message
         static_cast<Container::MessageStorageType*>(msg_)->operator[](dim_) = x;
         // propagate difference to left and right factor
         msg_->RepamLeft( diff, dim_);
         msg_->RepamRight( diff, dim_);
         return *this;
      }
      MsgVal& operator-=(const REAL x) __attribute__ ((always_inline))
      {
         static_cast<Container::MessageStorageType*>(msg_)->operator[](dim_) -= x;
         msg_->RepamLeft( -x, dim_);
         msg_->RepamRight( -x, dim_);
         return *this;
      }
      MsgVal& operator+=(const REAL x) __attribute__ ((always_inline))
      {
         static_cast<Container::MessageStorageType*>(msg_)->operator[](dim_) += x;
         msg_->RepamLeft( x, dim_);
         msg_->RepamRight( x, dim_);
         return *this;
      }
      operator REAL() const __attribute__ ((always_inline)) { return static_cast<Container::MessageStorageType*>(msg_)->operator[](dim_); }
   private:
      Container* const msg_;
      const INDEX dim_;
   };
   MsgVal operator[](const INDEX i) {
      return MsgVal(this,i);
   }
   const REAL operator[](const INDEX i) const {
      return MessageStorageType::operator[](i); // do zrobienia: needed?
   }


   // there must be four different implementations of msg updating with SIMD: 
   // (i) If parallel reparametrization is not supported by either left and right factor
   // If (ii) left or (iii) right factor supports reparametrization
   // If (iv) left and right factor support reparametrization


   template<typename ARRAY>
   Container& operator=(const ARRAY& diff) {
      throw std::runtime_error("not implemented yet");
      return *this;
   }

   template<typename ARRAY>
   Container& operator-=(const ARRAY& diff) {
      MinusVec<ARRAY> minus_diff(diff);
      assert(minus_diff.size() == this->size());
      RepamLeft(minus_diff);
      RepamRight(minus_diff);
      return *this;
   }

   template<typename ARRAY>
   Container& operator+=(const ARRAY& diff) {
      PlusVec<ARRAY> plus_diff(diff); // used to wrap Vc::Memory, otherwise not needed
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
   void SetMessage(const std::valarray<REAL>& m) { 
      assert(m.size() == MessageStorageType::size());
      
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
   const std::valarray<REAL> GetMessage() const { 
      std::valarray<REAL> m(0.0, MessageStorageType::size());
      for(INDEX i=0; i<MessageStorageType::size(); ++i) {
         m[i] = MessageStorageType::operator[](i);
      }
      return m; 
   }
   // possibly not the best choice: Sometimes msg_op_ needs access to this class
   MessageType msg_op_;

protected:
   LEFT_FACTOR_TYPE* const leftFactor_;
   RIGHT_FACTOR_TYPE* const rightFactor_;

};


// container class for factors. Here we hold the factor, all connected messages, reparametrization storage and perform reparametrization and coordination for sending and receiving messages.
// derives from REPAM_STORAGE_TYPE to mixin a class for storing the reparametrized potential
// implements the interface from FactorTypeAdapter for acces from LP_MP
template<typename FACTOR_TYPE, 
         template<class> class REPAM_STORAGE_TYPE, 
         class FACTOR_MESSAGE_TRAIT,
         INDEX FACTOR_NO>
class FactorContainer : public REPAM_STORAGE_TYPE<FactorContainer<FACTOR_TYPE, REPAM_STORAGE_TYPE, FACTOR_MESSAGE_TRAIT, FACTOR_NO> >, public FactorTypeAdapter
{
public:
   //typedef FactorContainer<FACTOR_TYPE,REPAM_STORAGE_TYPE, MESSAGE_TYPELIST, FACTOR_SPEC> ContainerType;
   typedef FACTOR_TYPE FactorType;
   typedef REPAM_STORAGE_TYPE<FactorContainer<FACTOR_TYPE, REPAM_STORAGE_TYPE, FACTOR_MESSAGE_TRAIT, FACTOR_NO> > RepamStorageType;
   friend class REPAM_STORAGE_TYPE<FactorContainer<FACTOR_TYPE, REPAM_STORAGE_TYPE, FACTOR_MESSAGE_TRAIT, FACTOR_NO> >;
   
   FactorContainer(FactorType factor, const std::vector<double>& cost) : RepamStorageType(factor,cost), factor_(factor) {
      //INDEX status;
      //std::cout << "msg_ type= "  << abi::__cxa_demangle(typeid(msg_).name(),0,0,&status) << "\n";
      //std::cout << "dispatcher list = "  << abi::__cxa_demangle(typeid(MESSAGE_DISPATCHER_TYPELIST).name(),0,0,&status) << "\n";
      //std::cout << "msg_ type= "  << abi::__cxa_demangle(typeid(msg_).name(),0,0,&status) << "\n";
      //std::cout << "left message list = " << abi::__cxa_demangle(typeid(left_message_list_).name(),0,0,&status) << "\n";
      //std::cout << "left message list = " << abi::__cxa_demangle(typeid(left_message_list_1).name(),0,0,&status) << "\n";
   
   }
   /*
   FactorContainer(const std::vector<REAL>& pot) : RepamStorageType(pot) {
      INDEX status;
      std::cout << "msg holding type = " << abi::__cxa_demangle(typeid(msg_).name(),0,0,&status) << "\n";
      std::cout << "msg type list = " << abi::__cxa_demangle(typeid(msg_type_list).name(),0,0,&status) << "\n";
      std::cout << "msg dispatcher list = " << abi::__cxa_demangle(typeid(MESSAGE_DISPATCHER_TYPELIST).name(),0,0,&status) << "\n";
      //std::get<0>(msg_);
      std::cout << "Constructor for FactorContainer " << std::endl; }
      */
   virtual ~FactorContainer() { 
      static_assert(meta::unique<MESSAGE_DISPATCHER_TYPELIST>::size() == MESSAGE_DISPATCHER_TYPELIST::size(), 
            "Message dispatcher typelist must have unique elements");
      static_assert(FACTOR_NO >= 0 && FACTOR_NO < FACTOR_MESSAGE_TRAIT::factor_list::size(), "factor number must be smaller than length of factor list");
   }

   template<typename MESSAGE_DISPATCHER_TYPE, typename MESSAGE_TYPE> 
      void AddMessage(MESSAGE_TYPE* m) { 
         constexpr INDEX n = meta::find_index<MESSAGE_DISPATCHER_TYPELIST, MESSAGE_DISPATCHER_TYPE>::value;
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
      constexpr INDEX n = meta::find_index<MESSAGE_DISPATCHER_TYPELIST, MESSAGE_DISPATCHER_TYPE>::value;
      REAL msg_val = 0.0;
      for(auto it=std::get<n>(msg_).begin(); it!=std::get<n>(msg_).end(); ++it) {
         msg_val += MESSAGE_DISPATCHER_TYPE::GetMessage(*(*it),i);
      }
      // receive messages for subsequent MESSAGE_DISPATCER_TYPES
      return msg_val + GetMessageSum(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{},i);
   }


   void UpdateFactor(const std::vector<REAL>& omega) {
      ReceiveMessages(omega);
      factor_.MaximizePotential(); // do zrobienia: pass the reparametrized potential here. Currently, MaximizePotential makes no sense, as we first need to have a "self-message" to modify the reparametrization. This is not implemented now
      SendMessages(omega);
   }

   // SFINAE-based seletion whether we will perform message updates for receiving   
   template<typename MESSAGE_DISPATCHER_TYPE, typename MSG_ARRAY, typename ITERATOR>
   typename std::enable_if<MESSAGE_DISPATCHER_TYPE::CanCallReceiveMessage() == true>::type 
   ReceiveMessagesImpl(MESSAGE_DISPATCHER_TYPE msg_dispatcher, const MSG_ARRAY& msgs, ITERATOR omegaIt)
   {
      // receive messages for current MESSAGE_DISPATCER_TYPE
      constexpr INDEX n = meta::find_index<MESSAGE_DISPATCHER_TYPELIST, MESSAGE_DISPATCHER_TYPE>::value;
      for(auto it=std::get<n>(msg_).cbegin(); it!=std::get<n>(msg_).cend(); ++it, ++omegaIt) {
         //if(*omegaIt != 0.0) {
         MESSAGE_DISPATCHER_TYPE::ReceiveMessage(*(*it));
         //}
      }
   }
   // or do not perform receiving message updates (if no receive message is implemented)
   template<typename MESSAGE_DISPATCHER_TYPE, typename MSG_ARRAY, typename ITERATOR>
   typename std::enable_if<MESSAGE_DISPATCHER_TYPE::CanCallReceiveMessage() == false>::type 
   ReceiveMessagesImpl(MESSAGE_DISPATCHER_TYPE msg_dispatcher, const MSG_ARRAY& msgs, ITERATOR omegaIt)
   {}

   void ReceiveMessages(const std::vector<REAL>& omega) {
      assert(omega.size() == GetNoMessages());
      ReceiveMessages(MESSAGE_DISPATCHER_TYPELIST{}, omega.cbegin());
   }
   template<typename ITERATOR, typename... MESSAGE_DISPATCHER_TYPES_REST>
   void ReceiveMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>, ITERATOR omegaIt) {}
   template<typename ITERATOR, typename MESSAGE_DISPATCHER_TYPE, typename... MESSAGE_DISPATCHER_TYPES_REST>
   void ReceiveMessages(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...>, ITERATOR omegaIt) {
      // receive messages for current MESSAGE_DISPATCER_TYPE
      constexpr INDEX n = meta::find_index<MESSAGE_DISPATCHER_TYPELIST, MESSAGE_DISPATCHER_TYPE>::value;
      ReceiveMessagesImpl(MESSAGE_DISPATCHER_TYPE{}, std::get<n>(msg_), omegaIt);
      omegaIt += std::get<n>(msg_).size(); 
      // receive messages for subsequent MESSAGE_DISPATCER_TYPES
      ReceiveMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{}, omegaIt);
   }


   void SendMessages(const std::vector<REAL>& omega) {
      assert(omega.size() == GetNoMessages());
      static constexpr INDEX n = NumberOfSendMessagesCalls<std::valarray<REAL>, decltype(omega.begin())>(MESSAGE_DISPATCHER_TYPELIST{});
      // do zrobienia: also do not construct currentRepam, if exactly one message update call will be issued. 
      // Check if there is one message dispatcher such that its size can be called via a constexpr function and is 1 -> complicated!
      if( n > 0 ) { // no need to construct currentRepam, if it will not be used at all
         // make a copy of the current reparametrization. The new messages are computed on it. Messages are updated implicitly and hence possibly the new reparametrization is automatically adjusted, which would INDEXerfere with message updates
         // do zrobienia: use static memory for this, do not always allocate new memory
         std::valarray<REAL> currentRepam(RepamStorageType::size());
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
   constexpr static INDEX NumberOfSendMessagesCalls(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...> t) { 
      constexpr INDEX n = meta::find_index<MESSAGE_DISPATCHER_TYPELIST, MESSAGE_DISPATCHER_TYPE>::value;
      constexpr INDEX no = MESSAGE_DISPATCHER_TYPE::template CanCallSendMessages<decltype(std::get<n>(msg_)), REPAM_ARRAY, ITERATOR>()
         || MESSAGE_DISPATCHER_TYPE::template CanCallSendMessage<REPAM_ARRAY>()
         ? 1 : 0;
      return no + NumberOfSendMessagesCalls<REPAM_ARRAY,ITERATOR>(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{});
   }

   // SFINAE-based seletion whether we will do batch or individual message updates for sending
   // batch message update
   template<typename MESSAGE_DISPATCHER_TYPE, typename ITERATOR, typename REPAM_ARRAY, typename MSG_ARRAY>
   typename std::enable_if<MESSAGE_DISPATCHER_TYPE::template CanCallSendMessages<MSG_ARRAY,REPAM_ARRAY,ITERATOR>() == true>::type 
   SendMessagesImpl(MESSAGE_DISPATCHER_TYPE msg_dispatcher, const MSG_ARRAY& msgs, const REPAM_ARRAY& repam, ITERATOR omegaIt)
   {
      constexpr INDEX n = meta::find_index<MESSAGE_DISPATCHER_TYPELIST, MESSAGE_DISPATCHER_TYPE>::value;
      const REAL omega_sum = std::accumulate(omegaIt, omegaIt + std::get<n>(msg_).size(), 0.0);
      if(omega_sum > 0.0) {
         // do zrobienia: construct proxy object for msgs, so that it directly poINDEXs to &(msgs[i]->msg_op_), make msg_op_ protected in MessageContainer again
         MESSAGE_DISPATCHER_TYPE::SendMessages(msgs, repam, omegaIt);
      }
   }
   // individual message update
   template<typename MESSAGE_DISPATCHER_TYPE, typename ITERATOR, typename REPAM_ARRAY, typename MSG_ARRAY>
   typename std::enable_if<
   MESSAGE_DISPATCHER_TYPE::template CanCallSendMessages<MSG_ARRAY,REPAM_ARRAY,ITERATOR>() == false && 
   MESSAGE_DISPATCHER_TYPE::template CanCallSendMessage <REPAM_ARRAY>() == true
   >::type 
   SendMessagesImpl(MESSAGE_DISPATCHER_TYPE msg_dispatcher, const MSG_ARRAY& msgs, const REPAM_ARRAY& repam, ITERATOR omegaIt)
   {
      // call individual message updates
      for(auto it=msgs.cbegin(); it!=msgs.cend(); ++it, ++omegaIt) {
         if(*omegaIt != 0.0) {
            MESSAGE_DISPATCHER_TYPE::SendMessage(*(*it), repam, *omegaIt);
         }
      }
   }
   // no updates if they are not implemented 
   template<typename MESSAGE_DISPATCHER_TYPE, typename ITERATOR, typename REPAM_ARRAY, typename MSG_ARRAY>
   typename std::enable_if<
   MESSAGE_DISPATCHER_TYPE::template CanCallSendMessages<MSG_ARRAY,REPAM_ARRAY,ITERATOR>() == false &&
   MESSAGE_DISPATCHER_TYPE::template CanCallSendMessage <REPAM_ARRAY>() == false
   >::type 
   SendMessagesImpl(MESSAGE_DISPATCHER_TYPE msg_dispatcher, const MSG_ARRAY& msgs, const REPAM_ARRAY& repam, ITERATOR omegaIt)
   {}



   // note that messages must be iterated over in the same order as done by MessageIterator
   template<typename ITERATOR, typename ARRAY, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   void SendMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...> t, const ARRAY& repam, ITERATOR omegaIt) {}
   template<typename ITERATOR, typename ARRAY, typename MESSAGE_DISPATCHER_TYPE, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   void SendMessages(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...> t, const ARRAY& repam, ITERATOR omegaIt) { // to get the current MESSAGE_TYPE

      // receive messages for current MESSAGE_DISPATCHER_TYPE
      constexpr INDEX n = meta::find_index<MESSAGE_DISPATCHER_TYPELIST, MESSAGE_DISPATCHER_TYPE>::value;

      // check whether the message supports batch updates. If so, call batch update. If not, check whether individual updates are supported. If yes, call individual updates. If no, do nothing
      SendMessagesImpl(MESSAGE_DISPATCHER_TYPE{}, std::get<n>(msg_), repam, omegaIt);
      omegaIt += std::get<n>(msg_).size();

      // receive messages for subsequent MESSAGE_DISPATCHER_TYPES
      SendMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{}, repam, omegaIt);
   }

   // methods used by MessageIterator
   template<typename ...MESSAGE_DISPATCHER_TYPES_REST>
   const INDEX GetNoMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...> t) const {
      return 0;
   }
   template<typename MESSAGE_DISPATCHER_TYPE, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   const INDEX GetNoMessages(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...> t) const {
      constexpr INDEX n = meta::find_index<MESSAGE_DISPATCHER_TYPELIST, MESSAGE_DISPATCHER_TYPE>::value;
      const INDEX no_msgs = std::get<n>(msg_).size(); 
      return no_msgs + GetNoMessages(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{});
      return 0;
   }
   const INDEX GetNoMessages() const {
      return GetNoMessages(MESSAGE_DISPATCHER_TYPELIST{});
   }

   template<typename ...MESSAGE_DISPATCHER_TYPES_REST>
   FactorTypeAdapter* GetConnectedFactor(meta::list<MESSAGE_DISPATCHER_TYPES_REST...> t, const INDEX cur_msg_idx) const {
      throw std::runtime_error("message index out of bound");
   }
   template<typename MESSAGE_DISPATCHER_TYPE, typename ...MESSAGE_DISPATCHER_TYPES_REST>
   FactorTypeAdapter* GetConnectedFactor(meta::list<MESSAGE_DISPATCHER_TYPE, MESSAGE_DISPATCHER_TYPES_REST...> t, const INDEX cur_msg_idx) const { // to get the current MESSAGE_TYPE
      constexpr INDEX n = meta::find_index<MESSAGE_DISPATCHER_TYPELIST, MESSAGE_DISPATCHER_TYPE>::value;
      const INDEX no_msgs = std::get<n>(msg_).size();
      if(cur_msg_idx < no_msgs) {
         auto msg = std::get<n>(msg_)[cur_msg_idx];
         if(msg->GetLeftFactor() == this) { return msg->GetRightFactor(); }
         else { return msg->GetLeftFactor(); }
      } else {
         return GetConnectedFactor(meta::list<MESSAGE_DISPATCHER_TYPES_REST...>{}, cur_msg_idx - no_msgs);
      }
   }
   FactorTypeAdapter* GetConnectedFactor (const INDEX msg_idx) const
   { 
      auto f = GetConnectedFactor(MESSAGE_DISPATCHER_TYPELIST{}, msg_idx);
      assert(f != this);
      return f;
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

   std::vector<REAL> GetReparametrizedPotential() const {
      std::vector<REAL> repam(factor_.size());
      for(INDEX i=0; i<repam.size(); ++i) {
         repam[i] = RepamStorageType::operator[](i);
      }
      return repam;
   }

   /*
   // get original cost. do zrobienia: not really needed and not always supported, better remove
   REAL GetOrigCost(const INDEX i) const { return factor_[i]; } 
   const std::vector<REAL> GetOrigPotential() const { return factor_.GetPotential(); } 
   */

   REAL LowerBound() const { return factor_.LowerBound(*this); }

   const FactorType* GetFactor() const { return &factor_; }

protected:
   FactorType factor_; // the factor operation

   // compile time metaprogramming to transform Factor-Message information into lists of which messages this factor must hold
   // first get lists with left and right message types
   struct get_msg_type_list {
      template<typename... E>
         using apply = meta::front<E...>;
   };
   struct get_left_msg {
      template<class... E>
         using apply = typename std::is_same< meta::at_c<E..., 1>, meta::size_t<FACTOR_NO> >::type;
   };
   struct get_right_msg {
      template<class... E>
         using apply = typename std::is_same< meta::at_c<E..., 2>, meta::size_t<FACTOR_NO> >::type;
   };
   struct get_left_msg_container_type_list {
      template<class... E>
         using apply = meta::at_c<E...,3>;
   };
   struct get_right_msg_container_type_list {
      template<class... E>
         using apply = meta::at_c<E...,4>;
   };

   using left_msg_list = meta::transform< meta::filter<typename FACTOR_MESSAGE_TRAIT::msg_list, get_left_msg>, get_msg_type_list>;
   using right_msg_list = meta::transform< meta::filter<typename FACTOR_MESSAGE_TRAIT::msg_list, get_right_msg>, get_msg_type_list>;
   using left_msg_container_list = meta::transform< meta::filter<typename FACTOR_MESSAGE_TRAIT::msg_list, get_left_msg>, get_left_msg_container_type_list>;
   using right_msg_container_list = meta::transform< meta::filter<typename FACTOR_MESSAGE_TRAIT::msg_list, get_right_msg>, get_right_msg_container_type_list>;

   // now construct a tuple with left and right dispatcher
   struct left_dispatch {
      template<class... E>
         using apply = MessageDispatcher<E..., LeftMessageFuncGetter>;
   };
   struct right_dispatch {
      template<class... E>
         using apply = MessageDispatcher<E..., RightMessageFuncGetter>;
   };
   using left_dispatcher_list = meta::transform< left_msg_list, left_dispatch >;
   using right_dispatcher_list = meta::transform< right_msg_list, right_dispatch >;

   using MESSAGE_DISPATCHER_TYPELIST = meta::concat<left_dispatcher_list, right_dispatcher_list>;

   // construct tuple holding messages for left and right dispatch
   template<class MSG> using get_vector_of_pointers = std::vector<MSG*>; // templatize this: e.g. standard pairwise factors hold exactly two messages.
   // the tuple will hold some container for the message type. The container type is specified in the fourth and fifth component of the message list
   using msg_list = meta::concat<left_msg_list, right_msg_list>;
   using msg_type_list = meta::transform< msg_list, meta::quote<get_vector_of_pointers> >;
   using msg_container_type_list = meta::concat<left_msg_container_list, right_msg_container_list>;

   tuple_from_list<msg_container_type_list> msg_;
};


// various containers: must implement push_back and indexing and iterators
// used for holding messages

// FixedSizeContainer must hold exactly N messages
// do zrobienia: use partial specialization as on page 732 in Stroustrup to eliminate code bloat
// do zrobienia: possibly overload operator[] to check for accessing nullptr
template<class T, INDEX N>
class FixedSizeContainer : public std::array<T,N>
{
public: 
   FixedSizeContainer() { this->fill(nullptr); }
   ~FixedSizeContainer() { static_assert(std::is_pointer<T>::value, "Message container must hold pointers to messages"); }
   void push_back(T t) {
      for(INDEX i=0; i<N; ++i) {
         if(this->operator[](i) == nullptr) {
            this->operator[](i) = t;
            return;
         }
      }
      throw std::range_error("added more messages than can be held");
   }
   constexpr size_t size() const { return N; }
};

} // end namespace LP_MP

#endif // LP_MP_FACTORS_MESSAGES_HXX

