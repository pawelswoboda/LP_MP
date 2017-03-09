#ifndef LP_MP_TEMPLATE_UTILITIES
#define LP_MP_TEMPLATE_UTILITIES

#include <type_traits>
#include <utility>
#include <tuple>
#include "config.hxx"
#include <cstddef>
#include <utility>
#include "meta/meta.hpp"

// mostly obsolete due to the meta-library
// do zrobienia: remove altogether

namespace LP_MP {

   // tuple iteration


   template <typename Tuple, typename F, std::size_t ...Indices>
      void for_each_tuple_impl(Tuple&& tuple, F&& f, std::index_sequence<Indices...>) {
         using swallow = int[];
         (void)swallow{1,
            (f(std::get<Indices>(std::forward<Tuple>(tuple))), void(), int{})...
         };
      }

   template <typename Tuple, typename F>
      void for_each_tuple(Tuple&& tuple, F&& f) {
         constexpr std::size_t N = std::tuple_size<std::remove_reference_t<Tuple>>::value;
         for_each_tuple_impl(std::forward<Tuple>(tuple), std::forward<F>(f),
               std::make_index_sequence<N>{});
      }

   template<class LIST> using tuple_from_list = meta::apply<meta::quote<std::tuple>, LIST>;

   /* remove
      template<SIGNED_INDEX NUMERATOR, SIGNED_INDEX DENOMINATOR>
      struct RationalNumberTemplate
      {
      static constexpr REAL value = REAL(NUMERATOR)/REAL(DENOMINATOR);

   };
   */


/*
   template<class LIST> struct tuple_from_list_impl {};
   template<template<class...> class LIST, class...T> 
      struct tuple_from_list_impl<LIST<T...> > 
      { using type = std::tuple<T...>; };
   template<class LIST> using tuple_from_list = typename tuple_from_list_impl<LIST>::type;
*/

   // expression templates for various array transformations
   // introduce traits for expression templates such that Vc::Memory compatible expression templates are taken care of recursively
   // do zrobienia: put this all in expression_template.hxx
   // do zrobienia: expression templates not working recursively yet, as only references are stored, but references to expressions might vanish. They need to be stored by copy.

   // do zrobienia: detect for const expressions and let them not have assignment operators

   template<typename T>
   struct VecSlice {
      VecSlice(T& a, const INDEX begin, const INDEX end)
         : a_(a), begin_(begin), end_(end)
      {
         assert(end_ > begin_);
         assert(begin_ >= 0);
         assert(end_ <= a_.size());
      }
      const REAL operator[](const INDEX i) const {
         return a_[i + begin_];
      }
      REAL& operator[](const INDEX i) {
         return a_[i + begin_];
      }
      const INDEX size() const {
         return end_ - begin_;
      }
      private:
      T& a_;
      const INDEX begin_;
      const INDEX end_;
   };

   template<typename T>
   struct ScaledVec {
      ScaledVec(const REAL& omega, const T& a) : omega_(omega), a_(a) {}
      const REAL operator[](const INDEX i) const {
         return omega_*a_[i];
      }
      const INDEX size() const {
         return a_.size();
      }
      private:
      const T& a_;
      const REAL& omega_;
   };

   //template<>
   //struct ScaledVec<Vc::Memory<REAL_SIMD>> {
   //
   //};

   template<typename T>
      struct PlusVec {
         PlusVec(const T& a) : a_(a) {}
         const REAL operator[](const INDEX i) const {
            return a_[i];
         }
         const INDEX size() const {
            return a_.size();
         }
         private:
         const T& a_;
      };

   /*
   template<>
      struct PlusVec<Vc::Memory<REAL_SIMD>> {
         PlusVec(const Vc::Memory<REAL_SIMD>& a) : a_(a) {}
         //auto vector(const INDEX i) const -> decltype(- std::declval<const Vc::Memory<REAL_SIMD>>().vector(0)) {
         const REAL_SIMD vector(const INDEX i) const {
            REAL_SIMD tmp = a_.vector(i);
            return tmp;
         }
         const REAL operator[](const INDEX i) const {
            return a_[i];
         }
         const INDEX size() const {
            return a_.entriesCount();
         }
         const INDEX vectorsCount() const {
            return a_.vectorsCount();
         }
         private:
         const Vc::Memory<REAL_SIMD> a_;
      };  
      */


   template<typename T>
      struct MinusVec {
         MinusVec(const T& a) : a_(a) {}
         const REAL operator[](const INDEX i) const {
            return -a_[i];
         }
         const INDEX size() const {
            return a_.size();
         }
         private:
         const T& a_;
      };

/*
   template<>
      struct MinusVec<Vc::Memory<REAL_SIMD>> {
         MinusVec(const Vc::Memory<REAL_SIMD>& a) : a_(a) {}
         //auto vector(const INDEX i) const -> decltype(- std::declval<const Vc::Memory<REAL_SIMD>>().vector(0)) {
         const REAL_SIMD vector(const INDEX i) const {
            REAL_SIMD tmp = a_.vector(i);
            return -tmp;
         }
         const REAL operator[](const INDEX i) const {
            return -a_[i];
         }
         const INDEX size() const {
            return a_.entriesCount();
         }
         const INDEX vectorsCount() const {
            return a_.vectorsCount();
         }
         private:
         const Vc::Memory<REAL_SIMD> a_;
      };  
*/

   template<typename T1, typename T2>
   struct PlusExprVec {
      PlusExprVec(const T1& a, const T2& b) : a_(a), b_(b) {
         assert(a_.size() == b_.size());
      }
      const REAL operator[](const INDEX i) const {
         return a_[i] + b_[i];
      }
      const INDEX size() const {
         return a_.size();
      }
      private:
      const T1& a_;
      const T2& b_;
   };


   template<typename T1, typename T2>
   struct MinusExprVec {
      MinusExprVec(const T1& a, const T2& b) : a_(a), b_(b) {
         assert(a_.size() == b_.size());
      }
      const REAL operator[](const INDEX i) const {
         return a_[i] - b_[i];
      }
      const INDEX size() const {
         return a_.size();
      }
      private:
      const T1& a_;
      const T2& b_;
   };




//===============================================================================
// Determine function return type, number of arguments and argument types
template<class F>
struct function_traits;
template<class R, class... Args>
struct function_traits<R(*)(Args...)> : public function_traits<R(Args...)>
{};
// function pointer
template<class R, class... Args>
struct function_traits<R(Args...)>
{
   using return_type = R;
   static constexpr std::size_t arity = sizeof...(Args);
   template<std::size_t N>
      struct argument
      {
         static_assert(N < arity, "error: invalid parameter index.");
         using type = typename std::tuple_element<N,std::tuple<Args...> >::type;
      };
};
// member function pointer
template<class C, class R, class... Args>
struct function_traits<R(C::*)(Args...)> : public function_traits<R(C&,Args...)>
{};
// const member function pointer
template<class C, class R, class... Args>
struct function_traits<R(C::*)(Args...) const> : public function_traits<R(C&,Args...)>
{};
// member object pointer
template<class C, class R>
struct function_traits<R(C::*)> : public function_traits<R(C&)>
{};



//===============================================================================

namespace integer_sequence {

template <typename T, T... N>
struct integer_sequence
{
  typedef T value_type;
  static_assert(
    std::is_integral<T>::value,
    "std::integer_sequence can only be instantiated with an integral type" );

  static inline
  std::size_t size(
  ) {
    return (sizeof...(N));
  }
};

template <std::size_t... N>
using index_sequence = integer_sequence<std::size_t, N...>;

namespace integer_sequence_detail {

template <typename T, std::size_t ..._Extra>
struct repeat;

template <typename T, T ...N, std::size_t ..._Extra>
struct repeat<integer_sequence<T, N...>, _Extra...>
{
  typedef integer_sequence<T, N...,
    1 * sizeof...(N) + N...,
    2 * sizeof...(N) + N...,
    3 * sizeof...(N) + N...,
    4 * sizeof...(N) + N...,
    5 * sizeof...(N) + N...,
    6 * sizeof...(N) + N...,
    7 * sizeof...(N) + N...,
    _Extra...> type;
};

template <std::size_t N> struct parity;
template <std::size_t N> struct make:parity<N % 8>::template pmake<N> {};

template <> struct make<0> { typedef integer_sequence<std::size_t> type; };
template <> struct make<1> { typedef integer_sequence<std::size_t, 0> type; };
template <> struct make<2> { typedef integer_sequence<std::size_t, 0, 1> type; };
template <> struct make<3> { typedef integer_sequence<std::size_t, 0, 1, 2> type; };
template <> struct make<4> { typedef integer_sequence<std::size_t, 0, 1, 2, 3> type; };
template <> struct make<5> { typedef integer_sequence<std::size_t, 0, 1, 2, 3, 4> type; };
template <> struct make<6> { typedef integer_sequence<std::size_t, 0, 1, 2, 3, 4, 5> type; };
template <> struct make<7> { typedef integer_sequence<std::size_t, 0, 1, 2, 3, 4, 5, 6> type; };

template <> struct parity<0> { template <std::size_t N> struct pmake:repeat<typename make<N / 8>::type> {}; };
template <> struct parity<1> { template <std::size_t N> struct pmake:repeat<typename make<N / 8>::type, N - 1> {}; };
template <> struct parity<2> { template <std::size_t N> struct pmake:repeat<typename make<N / 8>::type, N - 2, N - 1> {}; };
template <> struct parity<3> { template <std::size_t N> struct pmake:repeat<typename make<N / 8>::type, N - 3, N - 2, N - 1> {}; };
template <> struct parity<4> { template <std::size_t N> struct pmake:repeat<typename make<N / 8>::type, N - 4, N - 3, N - 2, N - 1> {}; };
template <> struct parity<5> { template <std::size_t N> struct pmake:repeat<typename make<N / 8>::type, N - 5, N - 4, N - 3, N - 2, N - 1> {}; };
template <> struct parity<6> { template <std::size_t N> struct pmake:repeat<typename make<N / 8>::type, N - 6, N - 5, N - 4, N - 3, N - 2, N - 1> {}; };
template <> struct parity<7> { template <std::size_t N> struct pmake:repeat<typename make<N / 8>::type, N - 7, N - 6, N - 5, N - 4, N - 3, N - 2, N - 1> {}; };

template <typename T, typename U>
struct convert
{
  template <typename>
  struct result;

  template <T ...N>
  struct result<integer_sequence<T, N...> >
  {
    typedef integer_sequence<U, N...> type;
  };
};

template <typename T>
struct convert<T, T>
{
  template <typename U>
  struct result
  {
    typedef U type;
  };
};

template <typename T, T N>
using make_integer_sequence_unchecked =
typename convert<std::size_t, T>::template result<typename make<N>::type>::type;

template <typename T, T N>
struct make_integer_sequence
{
  static_assert(std::is_integral<T>::value,
    "std::make_integer_sequence can only be instantiated with an integral type");
  static_assert(0 <= N,"std::make_integer_sequence input shall not be negative");

  typedef make_integer_sequence_unchecked<T, N> type;
};

} // namespace integer_sequence_detail


template <typename T, T N>
using make_integer_sequence = typename integer_sequence_detail::make_integer_sequence<T, N>::type;

template <std::size_t N>
using make_index_sequence = make_integer_sequence<std::size_t, N>;

template <typename... T>
using index_sequence_for = make_index_sequence<sizeof...(T)>;

} // end namespace integer_sequence


// get tuple element by type
template<class T, std::size_t N, class ... Args>
struct get_number_of_element_from_tuple_by_type_impl
{
   static constexpr auto value = N;
};

template<class T, std::size_t N, class ... Args>
struct get_number_of_element_from_tuple_by_type_impl<T,N,T,Args...>
{
   static constexpr auto value = N;
};

template<class T, std::size_t N, class U, class ... Args>
struct get_number_of_element_from_tuple_by_type_impl<T,N,U,Args...>
{
   static constexpr auto value = get_number_of_element_from_tuple_by_type_impl<T, N+1, Args...>::value;
};

// S = std::tuple<Args...>
template<class S>
constexpr S* get_null_pointer() { return nullptr; }
template<class T, class S>
constexpr int get_number_of_element_by_type()
{
   // static const S* t{nullptr}; // one can directly declare this in C++14 and omit get_null_pointer
   return get_number_of_element_by_type<T>(get_null_pointer<S>()); // this is ugly!
}

template<class T, class... Args>
constexpr int get_number_of_element_by_type(std::tuple<Args...>* t)
{
   return get_number_of_element_from_tuple_by_type_impl<T, 0, Args...>::value;
}

template<class T, class... Args>
T& get_element_by_type(std::tuple<Args...>& t)
{
   return std::get<get_number_of_element_from_tuple_by_type_impl<T, 0, Args...>::value>(t);
}


/*
template<typename S, typename T, typename ...REST>
struct tuple_no<std::tuple<T,REST...> { 
   //std::enable_if<std::is_same<S,T>>
   static constexpr int value = 0; 
};
*/
/*
template<typename S, typename T, typename ...REST>
struct tuple_no<std::tuple<T,REST...> { 
   int value = 1 + tuple_no<S, REST...>;
};

template<typename S, typename ...REST>
struct tuple_no<S,REST...> {
   int value = -1000000000;
};
*/





//===============================================================================
// META-FUNCTIONS FOR EXTRACTING THE n-th TYPE OF A PARAMETER PACK

// Declare primary template
template<int I, typename... Ts>
struct nth_type_of
{
};

// Base step
template<typename T, typename... Ts>
struct nth_type_of<0, T, Ts...>
{
    using type = T;
};

// Induction step
template<int I, typename T, typename... Ts>
struct nth_type_of<I, T, Ts...>
{
    using type = typename nth_type_of<I - 1, Ts...>::type;
};

// Helper meta-function for retrieving the first type in a parameter pack
template<typename... Ts>
struct first_type_of
{
    using type = typename nth_type_of<0, Ts...>::type;
};

// Helper meta-function for retrieving the last type in a parameter pack
template<typename... Ts>
struct last_type_of
{
    using type = typename nth_type_of<sizeof...(Ts) - 1, Ts...>::type;
};

//===============================================================================
// FUNCTIONS FOR EXTRACTING THE n-th VALUE OF AN ARGUMENT PACK

// Base step
template<int I, typename T, typename... Ts>
auto nth_value_of(T&& t, Ts&&... args) ->
    typename std::enable_if<(I == 0), decltype(std::forward<T>(t))>::type
{
    return std::forward<T>(t);
}

// Induction step
template<int I, typename T, typename... Ts>
auto nth_value_of(T&& t, Ts&&... args) ->
    typename std::enable_if<(I > 0), decltype(
        std::forward<typename nth_type_of<I, T, Ts...>::type>(
            std::declval<typename nth_type_of<I, T, Ts...>::type>()
            )
        )>::type
{
    using return_type = typename nth_type_of<I, T, Ts...>::type;
    return std::forward<return_type>(nth_value_of<I - 1>((std::forward<Ts>(args))...));
}

// Helper function for retrieving the first value of an argument pack
template<typename... Ts>
auto first_value_of(Ts&&... args) ->
    decltype(
        std::forward<typename first_type_of<Ts...>::type>(
            std::declval<typename first_type_of<Ts...>::type>()
            )
        )
{
    using return_type = typename first_type_of<Ts...>::type;
    return std::forward<return_type>(nth_value_of<0>((std::forward<Ts>(args))...));
}

// Helper function for retrieving the last value of an argument pack
template<typename... Ts>
auto last_value_of(Ts&&... args) ->
    decltype(
        std::forward<typename last_type_of<Ts...>::type>(
            std::declval<typename last_type_of<Ts...>::type>()
            )
        )
{
    using return_type = typename last_type_of<Ts...>::type;
    return std::forward<return_type>(nth_value_of<sizeof...(Ts) - 1>((std::forward<Ts>(args))...));
}

//===============================================================================
// METAFUNCTION FOR COMPUTING THE UNDERLYING TYPE OF HOMOGENEOUS PARAMETER PACKS

// Used as the underlying type of non-homogeneous parameter packs
struct null_type
{
};

// Declare primary template
template<typename... Ts>
struct homogeneous_type;

// Base step
template<typename T>
struct homogeneous_type<T>
{
    using type = T;
    static const bool isHomogeneous = true;
};

// Induction step
template<typename T, typename... Ts>
struct homogeneous_type<T, Ts...>
{
    // The underlying type of the tail of the parameter pack
    using type_of_remaining_parameters = typename homogeneous_type<Ts...>::type;

    // True if each parameter in the pack has the same type
    static const bool isHomogeneous = std::is_same<T, type_of_remaining_parameters>::value;

    // If isHomogeneous is "false", the underlying type is the fictitious null_type
    using type = typename std::conditional<isHomogeneous, T, null_type>::type;
};

// Meta-function to determine if a parameter pack is homogeneous
template<typename... Ts>
struct is_homogeneous_pack
{
    static const bool value = homogeneous_type<Ts...>::isHomogeneous;
};

//===============================================================================
// META-FUNCTIONS FOR CREATING INDEX LISTS

// The structure that encapsulates index lists
template <unsigned... Is>
struct index_list
{
};

// Collects internal details for generating index ranges [MIN, MAX)
namespace detail
{
    // Declare primary template for index range builder
    template <unsigned MIN, unsigned N, unsigned... Is>
    struct range_builder;

    // Base step
    template <unsigned MIN, unsigned... Is>
    struct range_builder<MIN, MIN, Is...>
    {
        typedef index_list<Is...> type;
    };

    // Induction step
    template <unsigned MIN, unsigned N, unsigned... Is>
    struct range_builder : public range_builder<MIN, N - 1, N - 1, Is...>
    {
    };
}

// Meta-function that returns a [MIN, MAX) index range
template<unsigned MIN, unsigned MAX>
using index_range = typename detail::range_builder<MIN, MAX>::type;

//===============================================================================
// CLASSES AND FUNCTIONS FOR REALIZING LOOPS ON ARGUMENT PACKS

// Implementation inspired by @jogojapan's answer to this question:
// http://stackoverflow.com/questions/14089637/return-several-arguments-for-another-function-by-a-single-function

// Collects internal details for implementing functor invocation
namespace detail
{
    // Functor invocation is realized through variadic inheritance.
    // The constructor of each base class invokes an input functor.
    // An functor invoker for an argument pack has one base class
    // for each argument in the pack

    // Realizes the invocation of the functor for one parameter
    template<unsigned I, typename T>
    struct invoker_base
    {
        template<typename F, typename U>
        invoker_base(F&& f, U&& u) { f(u); }
    };

    // Necessary because a class cannot inherit the same class twice
    template<unsigned I, typename T>
    struct indexed_type
    {
        static const unsigned int index = I;
        using type = T;
    };

    // The functor invoker: inherits from a list of base classes.
    // The constructor of each of these classes invokes the input
    // functor with one of the arguments in the pack.
    template<typename... Ts>
    struct invoker : public invoker_base<Ts::index, typename Ts::type>...
    {
        template<typename F, typename... Us>
        invoker(F&& f, Us&&... args)
            :
            invoker_base<Ts::index, typename Ts::type>(std::forward<F>(f), std::forward<Us>(args))...
        {
        }
    };
}

// The functor provided in the first argument is invoked for each
// argument in the pack whose index is contained in the index list
// specified in the second argument
template<typename F, unsigned... Is, typename... Ts>
void for_each_in_arg_pack_subset(F&& f, index_list<Is...> const& i, Ts&&... args)
{
    // Constructors of invoker's sub-objects will invoke the functor.
    // Note that argument types must be paired with numbers because the
    // implementation is based on inheritance, and one class cannot
    // inherit the same base class twice.
    detail::invoker<detail::indexed_type<Is, typename nth_type_of<Is, Ts...>::type>...> invoker(
        f,
        (nth_value_of<Is>(std::forward<Ts>(args)...))...
        );
}

// The functor provided in the first argument is invoked for each
// argument in the pack
template<typename F, typename... Ts>
void for_each_in_arg_pack(F&& f, Ts&&... args)
{
    for_each_in_arg_pack_subset(f, index_range<0, sizeof...(Ts)>(), std::forward<Ts>(args)...);
}

// The functor provided in the first argument is given in input the
// arguments in whose index is contained in the index list specified
// as the second argument.
template<typename F, unsigned... Is, typename... Ts>
void forward_subpack(F&& f, index_list<Is...> const& i, Ts&&... args)
{
    f((nth_value_of<Is>(std::forward<Ts>(args)...))...);
}

// The functor provided in the first argument is given in input all the
// arguments in the pack.
template<typename F, typename... Ts>
void forward_pack(F&& f, Ts&&... args)
{
    f(std::forward<Ts>(args)...);
}


}

#endif // LP_MP_TEMPLATE_UTILITIES
