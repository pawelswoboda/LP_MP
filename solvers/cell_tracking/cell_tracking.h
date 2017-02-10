#ifndef LP_MP_CELL_TRACKING__HXX
#define LP_MP_CELL_TRACKING__HXX

#include "factors_messages.hxx"
#include "detection_factor.hxx"
#include "cell_tracking_constructor.hxx"
#include "conservation_tracking_constructor.hxx"
#include "LP_MP.h"

namespace LP_MP {
  /*
struct FMC_CELL_TRACKING {
  constexpr static char* name = "Cell tracking";

  using detection_factor_container = FactorContainer<detection_factor, FMC_CELL_TRACKING, 0, false>;
  using at_most_one_hypothesis_container = FactorContainer<at_most_one_cell_factor, FMC_CELL_TRACKING, 1, false>;
  using transition_factor_container = FactorContainer<transition_factor, FMC_CELL_TRACKING, 2, false>;
  using split_factor_container = FactorContainer<split_factor, FMC_CELL_TRACKING, 3, false>;

  using at_most_one_cell_message_container = MessageContainer<at_most_one_cell_message, 0, 1, variableMessageNumber, variableMessageNumber, FMC_CELL_TRACKING, 0>;
  using outgoing_edge_message_container = MessageContainer<outgoing_edge_message, 0, 2, variableMessageNumber, atMostOneMessage, FMC_CELL_TRACKING, 1>;
  using incoming_edge_message_container = MessageContainer<incoming_edge_message, 2, 0, atMostOneMessage, variableMessageNumber, FMC_CELL_TRACKING, 2>;
  using transition_split_message_out_container = MessageContainer<transition_split_message_out , 2, 3, atMostOneMessage, 1, FMC_CELL_TRACKING, 3>;
  using transition_split_message_in_container = MessageContainer<transition_split_message_in , 3, 2, 2, atMostOneMessage, FMC_CELL_TRACKING, 4>;

  using FactorList = meta::list< detection_factor_container, at_most_one_hypothesis_container, transition_factor_container, split_factor_container >;
  using MessageList = meta::list< 
    at_most_one_cell_message_container,
    outgoing_edge_message_container,
    incoming_edge_message_container,
    transition_split_message_out_container,
    transition_split_message_in_container
      >;

  using constructor = cell_tracking_constructor<
    detection_factor_container, at_most_one_hypothesis_container, transition_factor_container, split_factor_container,
    at_most_one_cell_message_container, outgoing_edge_message_container, incoming_edge_message_container, transition_split_message_out_container, transition_split_message_in_container
      >;

  using ProblemDecompositionList = meta::list<constructor>;
};
*/

struct FMC_CELL_TRACKING_MOTHER_MACHINE {
  constexpr static char* name = "Cell tracking in the mother machine";

  using detection_factor_container = FactorContainer<detection_factor, FMC_CELL_TRACKING_MOTHER_MACHINE, 0, true>;
  using at_most_one_hypothesis_container = FactorContainer<at_most_one_cell_factor, FMC_CELL_TRACKING_MOTHER_MACHINE, 1, false>;
  using exit_constraint = FactorContainer<exit_constraint_factor, FMC_CELL_TRACKING_MOTHER_MACHINE, 2, false>;

  using transition_message_container = MessageContainer<transition_message, 0, 0, variableMessageNumber, variableMessageNumber, FMC_CELL_TRACKING_MOTHER_MACHINE, 0>;
  using at_most_one_cell_message_container = MessageContainer<at_most_one_cell_message, 0, 1, variableMessageNumber, variableMessageNumber, FMC_CELL_TRACKING_MOTHER_MACHINE, 1>;
  using exit_constraint_lower_message = MessageContainer<exit_constraint_message<exit_constraint_position::lower>, 0, 2, variableMessageNumber, 1, FMC_CELL_TRACKING_MOTHER_MACHINE, 2>;
  using exit_constraint_upper_message = MessageContainer<exit_constraint_message<exit_constraint_position::upper>, 0, 2, variableMessageNumber, 1, FMC_CELL_TRACKING_MOTHER_MACHINE, 3>;

  using FactorList = meta::list< detection_factor_container, at_most_one_hypothesis_container, exit_constraint >;
  using MessageList = meta::list< transition_message_container, at_most_one_cell_message_container, exit_constraint_lower_message, exit_constraint_upper_message >;

  using constructor = cell_tracking_constructor<
    detection_factor_container, at_most_one_hypothesis_container, 
    transition_message_container, at_most_one_cell_message_container
      >;

  using constructor_mother_machine = cell_tracking_mother_machine_constructor<
    constructor,
    exit_constraint, exit_constraint_lower_message, exit_constraint_upper_message
      >;

  using ProblemDecompositionList = meta::list<constructor_mother_machine>;
};

struct FMC_CONSERVATION_TRACKING {
  constexpr static char* name = "Conservation tracking";

  using detection_factor_container = FactorContainer<multiple_detection_factor, FMC_CONSERVATION_TRACKING, 0>;

  using transition_message_container = MessageContainer<transition_message_multiple, 0, 0, variableMessageNumber, variableMessageNumber, FMC_CONSERVATION_TRACKING, 0>;

  using FactorList = meta::list< detection_factor_container >;

  using MessageList = meta::list< transition_message_container >;

  using constructor = conservation_tracking_constructor<detection_factor_container, transition_message_container>;

  using ProblemDecompositionList = meta::list<constructor>;
};

} // end namespace LP_MP
#endif // LP_MP_CELL_TRACKING__HXX
