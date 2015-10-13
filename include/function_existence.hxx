#ifndef LP_MP_FUNCTION_EXISTENCE_HXX
#define LP_MP_FUNCTION_EXISTENCE_HXX

#include <type_traits>


// generates helper classes and constexpr function with which it is possible to detect  existence and callability of member functions

#define LP_MP_FUNCTION_EXISTENCE_CLASS(TESTER_NAME, MEMBER) \
template<typename, typename T> \
struct struct_##TESTER_NAME { \
   static_assert( \
         std::integral_constant<T, false>::value, \
         "Second template parameter needs to be of function type."); \
};  \
\
template<typename C, typename Ret, typename... Args> \
struct struct_##TESTER_NAME <C,Ret(Args...)> { \
private: \
   template<typename T> \
   static constexpr auto check(T*) \
   -> typename \
      std::is_same< \
      decltype( std::declval<T>() .MEMBER ( std::declval<Args>()...) ), \
      Ret \
         >::type; \
\
   template<typename> \
   static constexpr std::false_type check(...); \
\
   typedef decltype(check<C>(0)) type; \
\
public: \
   static constexpr bool value = type::value; \
}; \
template<class C, typename RET, typename... ARGS> \
constexpr static bool TESTER_NAME () { \
  return struct_##TESTER_NAME<C,RET(ARGS&...)>::value; \
} \





#endif // LP_MP_FUNCTION_EXISTENCE_HXX
