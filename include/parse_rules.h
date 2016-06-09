#ifndef LP_MP_PARSE_RULES_HXX
#define LP_MP_PARSE_RULES_HXX

#include "pegtl.hh"

// elementary grammar rules for PEGTL

namespace LP_MP {

      struct mand_whitespace : pegtl::plus< pegtl::blank > {}; 
      struct opt_whitespace : pegtl::star< pegtl::blank > {}; 
      struct positive_integer : pegtl::plus< pegtl::digit > {};

      struct real_number_standard : pegtl::sor<
	pegtl::seq< pegtl::opt< pegtl::one<'+','-'> >, pegtl::plus<pegtl::digit>, pegtl::opt< pegtl::seq< pegtl::string<'.'>, pegtl::star<pegtl::digit> > > >,
	pegtl::string<'I','n','f'>,
	pegtl::string<'i','n','f'>
	> {}; 
      struct real_number_smaller1 : pegtl::seq< pegtl::opt< pegtl::one<'+','-'> >, pegtl::string<'.'>, pegtl::plus< pegtl::digit > > {};
      struct real_number_exponential : pegtl::seq< pegtl::opt< pegtl::one<'+','-'> >, pegtl::star< pegtl::digit >, pegtl::opt<pegtl::seq< pegtl::string<'.'>, pegtl::star< pegtl::digit>>>, pegtl::string<'e'>, pegtl::opt< pegtl::one<'+','-'> >, pegtl::plus< pegtl::digit > > {};
      struct real_number : pegtl::sor<real_number_exponential, real_number_standard, real_number_smaller1> {};

      struct vector : pegtl::seq< pegtl::string<'['>, opt_whitespace, real_number, pegtl::star< pegtl::seq< mand_whitespace, real_number > >, opt_whitespace, pegtl::string<']'> > {};


   template<typename COMMENT_PREFIX>
      struct comment_line : pegtl::seq< opt_whitespace, COMMENT_PREFIX, pegtl::until< pegtl::eol >> {};
   // this does not work yet
   //template<char COMMENT_PREFIX> 
   //   struct comment_lines : pegtl::seq< opt_whitespace, pegtl::string<COMMENT_PREFIX>, pegtl::until< pegtl::eol >> {};
}

// do zrobienia: test the above definitions with unit testing whether they accept and reject what they are supposed to


#endif// LP_MP_PARSE_RULES_HXX
