#include <unistd.h>

#include <algorithm>
#include <chrono>
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
#include "mpi.h"
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

void remove_empty_slots(
    unsigned simItr, Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, SimulVolume& simulVolume,
    std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays,
    Membrane& membraneObject, MpiContext& mpiContext) {
  //------------------------------------------------------------------------------------
  // Remove empty complexes:
  // Put the last non-empty complex in the list to the first empty slot:
  sort(Complex::emptyComList.begin(), Complex::emptyComList.end());
  int lastNonEmptyIndex = complexList.size() - 1;

  for (auto& firstEmptyIndex : Complex::emptyComList) {
    // Find the last non-empty complex:
    while (complexList[lastNonEmptyIndex].isEmpty) lastNonEmptyIndex--;
    if (lastNonEmptyIndex <= firstEmptyIndex) {
      break;
    }
    // Move last non-empty complex to the first empty position:
    // cout << "  Moving a complex from " << lastNonEmptyIndex << " to " <<
    // firstEmptyIndex << "... "
    //         << "complexList[firstEmptyIndex].isEmpty=" <<
    //         complexList[firstEmptyIndex].isEmpty << endl;
    complexList[firstEmptyIndex] = complexList[lastNonEmptyIndex];

    // Update complex index to match new position:
    complexList[firstEmptyIndex].index = firstEmptyIndex;

    // Update myComIndex for all complex members to match firstEmptyIndex:
    for (auto& molIndex : complexList[firstEmptyIndex].memberList) {
      moleculeList[molIndex].myComIndex = firstEmptyIndex;
    }
    lastNonEmptyIndex--;
  }

  // Remove last elements in complexList
  for (int k = 0; k < Complex::emptyComList.size(); k++) complexList.pop_back();

  // Empty Complex::emptyComList:
  Complex::emptyComList.clear();

  // if(DEBUG) {debug_bndpartner_interface(mpiContext, "35: (before remove
  // empty molecules)");}

  //------------------------------------------------------------------------------------
  // Remove empty molecules:
  sort(Molecule::emptyMolList.begin(), Molecule::emptyMolList.end());
  // Always copy last occupied element from moleculeList to the first empty
  // index:
  lastNonEmptyIndex = moleculeList.size() - 1;
  // While there is a molecule to delete:
  for (auto& firstEmptyIndex : Molecule::emptyMolList) {
    // cout << "firstEmptyIndex=" << firstEmptyIndex << endl;
    if (DEBUG && (firstEmptyIndex >= moleculeList.size()))
      error("10: firstEmptyIndex(" + to_string(firstEmptyIndex) +
            ") >= moleculeList.size() (" + to_string(moleculeList.size()) +
            ")");
    // Find the last non-empty molecule:
    while (moleculeList[lastNonEmptyIndex].myComIndex == -1)
      lastNonEmptyIndex--;
    if (lastNonEmptyIndex <= firstEmptyIndex) break;
    // Move last non-empty molecule to the first empty position:
    moleculeList[firstEmptyIndex] = moleculeList[lastNonEmptyIndex];
    // Update molecule index to match new position:
    moleculeList[firstEmptyIndex].index = firstEmptyIndex;
    // Update complex member list elements to match new firstEmptyIndex:
    int tmpComIndex{moleculeList[firstEmptyIndex].myComIndex};
    for (auto& tmpMember : complexList[tmpComIndex].memberList) {
      if (tmpMember == lastNonEmptyIndex) tmpMember = firstEmptyIndex;
    }
    // For all partners,
    // update their interface.interaction.partnerIndex and mol.bndpartner
    // to match new firstEmptyIndex:
    for (auto& tmpPartner : moleculeList[firstEmptyIndex].bndpartner) {
      if (DEBUG && (tmpPartner >= moleculeList.size()) && tmpPartner != -1)
        error("13: tmpPartner(" + to_string(tmpPartner) +
              ") >= moleculeList.size()(" + to_string(moleculeList.size()) +
              ")");
      for (auto& partner : moleculeList[tmpPartner].bndpartner) {
        if (DEBUG && (partner >= moleculeList.size()))
          error("14: partner >= moleculeList.size()");
        if (partner == lastNonEmptyIndex) partner = firstEmptyIndex;
      }
      for (auto& tmpIface : moleculeList[tmpPartner].interfaceList) {
        if (DEBUG &&
            (tmpIface.interaction.partnerIndex >= moleculeList.size()) &&
            (tmpIface.interaction.partnerIndex != -1))
          error("15: tmpIface.interaction.partnerIndex >= moleculeList.size()");
        if (tmpIface.interaction.partnerIndex == lastNonEmptyIndex)
          tmpIface.interaction.partnerIndex = firstEmptyIndex;
      }
    }

    // update Molecule::emptyMolList
    // Molecule::emptyMolList.pop_back();
    lastNonEmptyIndex--;
  }

  // Remove last elements in moleculeList
  for (int k = 0; k < Molecule::emptyMolList.size(); k++)
    moleculeList.pop_back();
  Molecule::emptyMolList.clear();
}