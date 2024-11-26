#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>

#include "boundary_conditions/reflect_functions.hpp"
#include "debug/debug.hpp"
#include "error/error.hpp"
#include "io/io.hpp"
#include "macro.hpp"
#include "math/constants.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#ifdef mpi_
#include "mpi.h"
#endif
#include "mpi/mpi_function.hpp"
#include "parser/parser_functions.hpp"
#include "reactions/association/association.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "split.cpp"
#include "system_setup/system_setup.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

void check_perform_zeroth_first_order_reactions(
    unsigned simItr, Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, SimulVolume& simulVolume,
    std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays,
    Membrane& membraneObject, std::vector<double>& IL2DbindingVec,
    std::vector<double>& IL2DUnbindingVec, std::vector<double>& ILTableIDs,
    MpiContext& mpiContext) {
  if (VERBOSE) printf("Checking 1st order reaction...\n");
  char fnameProXYZ[256];
  snprintf(fnameProXYZ, sizeof(fnameProXYZ), "DATA/assoc_dissoc_time.dat");
  std::ofstream assocDissocFile(fnameProXYZ);
  check_for_unimolecular_reactions(
      simItr, params, moleculeList, complexList, simulVolume, forwardRxns,
      backRxns, createDestructRxns, molTemplateList, observablesList,
      counterArrays, membraneObject, IL2DbindingVec, IL2DUnbindingVec,
      ILTableIDs, assocDissocFile);

  if (DEBUG) {
    DEBUG_FIND_MOL("4.41");
  }

  // Zeroth order reactions (creation)
  if (VERBOSE) printf("Checking 0th order reaction...\n");
  check_for_zeroth_order_creation(simItr, params, simulVolume, forwardRxns,
                                  createDestructRxns, moleculeList, complexList,
                                  molTemplateList, observablesList,
                                  counterArrays, membraneObject);
  if (DEBUG) {
    DEBUG_FIND_MOL("4.42");
  }

  // Update member lists after creation and destruction
  // if (VERBOSE)
  //   printf(
  //       "Updating memberMolLists for simulVolume after creation and "
  //       "destruction...\n");
  // simulVolume.update_memberMolLists(params, moleculeList, complexList,
  //                                   molTemplateList, membraneObject, simItr,
  //                                   mpiContext);

  if (params.hasUniMolStateChange == true) {
    if (VERBOSE) printf("Checking 1st state change reactions...\n");
    check_for_unimolstatechange_reactions(
        simItr, params, moleculeList, complexList, simulVolume, forwardRxns,
        backRxns, createDestructRxns, molTemplateList, observablesList,
        counterArrays, membraneObject);
  }
  if (DEBUG) {
    DEBUG_FIND_MOL("4.43");
  }

  // Skip this entire loop if the system has no implicit lipids.
  if (params.implicitLipid == true) {
    // check dissociation (implicit)
    if (VERBOSE) printf("Checking dissociation with the implicit lipid...\n");
    for (unsigned molItr{0}; molItr < moleculeList.size(); ++molItr) {
      // only do checks if the Molecule exists
      if (moleculeList[molItr].isEmpty ||
          moleculeList[molItr].isImplicitLipid == true ||
          complexList[moleculeList[molItr].myComIndex].OnSurface == false ||
          params.implicitLipid == false ||
          moleculeList[molItr].isGhosted == true)
        continue;
      check_dissociation_implicitlipid(
          simItr, params, simulVolume, molTemplateList, observablesList, molItr,
          moleculeList, complexList, backRxns, forwardRxns, createDestructRxns,
          counterArrays, membraneObject, IL2DbindingVec, IL2DUnbindingVec,
          ILTableIDs, assocDissocFile);
    }
  }
}