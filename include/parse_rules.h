#ifndef LP_MP_PARSE_RULES_HXX
#define LP_MP_PARSE_RULES_HXX

#include "pegtl.hh"
#include <sstream>

// elementary grammar rules for PEGTL

namespace LP_MP {
namespace Parsing {

      struct mand_whitespace : pegtl::plus< pegtl::blank > {}; 
      struct opt_whitespace : pegtl::star< pegtl::blank > {}; 
      struct opt_invisible : pegtl::star< pegtl::sor< pegtl::blank, pegtl::eol > > {};
      struct mand_invisible : pegtl::plus< pegtl::sor< pegtl::blank, pegtl::eol > > {};
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



   // variable length sequence of numbers, where first number denotes length (integer), followed by real numbers to be read in
   /*
   namespace variable_length_number_seq { 

      template< typename Rule >
         struct action
         : pegtl::nothing< Rule > {};

      template<> struct action< positive_integer >
      {
         static void apply(const pegtl::action_input & in, INDEX& length, std::vector<REAL>& real_seq )
         {
            length = std::stoul(in.string());
            real_seq.reserve(length);
            real_seq.resize(0);
         }
      };

      template<> struct action< real_number >
      {
         static void apply(const pegtl::action_input & in, INDEX& length, std::vector<REAL>& real_seq )
         {
            real_seq..push_back(std::stod(in.string()));
         }
      };

      struct sequence_end
      {
         template< pegtl::apply_mode A, template< typename ... > class Action, template< typename ... > class Control, typename Input >
            static bool match( Input & in, const INDEX& length, const std::vector<REAL>& real_seq )
            {
               return length <= real_seq.size();
            }
      };

      struct grammar
         : pegtl::seq< positive_integer, pegtl::until< sequence_end, pegtl::sor< invisible, real_number > > > {};

      template< typename Rule >
         struct action
         : pegtl::nothing< Rule > {};

   } // namespace variable_length_number_seq
   */

} // end namespace Parsing
}

// do zrobienia: test the above definitions with unit testing whether they accept and reject what they are supposed to


#endif// LP_MP_PARSE_RULES_HXX
