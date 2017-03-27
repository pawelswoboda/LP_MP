//#include "graph_matching.h"
//#include "visitors/standard_visitor.hxx"
//using FMC_INST = FMC_MP_T<PairwiseConstruction::Left>;
//using BaseSolverType = Solver<FMC_INST,LP,StandardTighteningVisitor>;
//LP_MP_CONSTRUCT_SOLVER_WITH_INPUT_AND_VISITOR_MP_ROUNDING(FMC_INST, TorresaniEtAlInput::ParseProblemMP<BaseSolverType>, StandardTighteningVisitor);

#include <iostream>
#include <array>
#include <type_traits>
#include <utility>

template<size_t... LABELS>
struct labeling {

   template<size_t LABEL_NO, size_t LABEL, size_t... LABELS_REST>
   constexpr 
   typename std::enable_if<LABEL_NO == 0,size_t>::type get_label()
   {
      return LABEL;
   }

   template<size_t LABEL_NO, size_t LABEL, size_t... LABELS_REST>
   constexpr 
   typename std::enable_if<(LABEL_NO > 0),size_t>::type get_label()
   {
      return get_label<LABEL_NO-1, LABELS_REST...>();
   }

   template<size_t LABEL_NO>
   constexpr size_t label()
   { 
      static_assert(LABEL_NO < sizeof...(LABELS), "label number must be smaller than number of labels");
      return get_label<LABEL_NO, LABELS...>();
   } 
};

template<typename... LABELINGS> // all labels must be instances of labeling
struct labelings
{
   template<size_t LABELING_NO, size_t LABEL_NO, typename LABELING, typename... LABELINGS_REST>
   constexpr 
   typename std::enable_if<LABELING_NO == 0,size_t>::type get_label()
   {
      return LABELING{}.template label<LABEL_NO>();
   }

   template<size_t LABELING_NO, size_t LABEL_NO, typename LABELING, typename... LABELINGS_REST>
   constexpr 
   typename std::enable_if<(LABELING_NO > 0),size_t>::type get_label()
   {
      return get_label<LABELING_NO-1, LABEL_NO, LABELINGS_REST...>();
   }

   template<size_t LABELING_NO, size_t LABEL_NO>
   constexpr size_t label()
   {
      static_assert(LABELING_NO < sizeof...(LABELINGS), "labeling number must be smaller than number of labelings");
      return get_label<LABELING_NO,LABEL_NO,LABELINGS...>();
   }
};



int main()
{
   using l = labelings<
      labeling<0,1,1,1,0>,
      labeling<1,0,0,1,1>,
      labeling<0,0,0,0,1>
         >;


   std::cout << l{}.label<0,0>() << "," << l{}.label<0,1>() << "," << l{}.label<0,2>() << "," << l{}.label<0,3>() << "," << l{}.label<0,4>() << "\n";
   std::cout << l{}.label<1,0>() << "," << l{}.label<1,1>() << "," << l{}.label<1,2>() << "," << l{}.label<1,3>() << "," << l{}.label<1,4>() << "\n";
   std::cout << l{}.label<2,0>() << "," << l{}.label<2,1>() << "," << l{}.label<2,2>() << "," << l{}.label<2,3>() << "," << l{}.label<2,4>() << "\n";
   return 0;
}
