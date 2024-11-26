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
#ifdef mpi_
#include "mpi.h"
#endif
#include "mpi/mpi_function.hpp"
#include "parser/parser_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "split.cpp"
#include "system_setup/system_setup.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

/**
 * Parses add parameters file for a restart simulation.
 *
 * @param addFileNameInput Name of the add.inp file.
 * @param params Simulation parameters.
 * @param observablesList Map of observables to track during simulation.
 * @param forwardRxns List of forward reactions.
 * @param backRxns List of backward reactions.
 * @param createDestructRxns List of create/destruct reactions.
 * @param molTemplateList List of molecular templates.
 * @param membraneObject Membrane object.
 * @param moleculeList List of individual molecules.
 * @param complexList List of complexes.
 * @param numMolTemplateBeforeAdd Number of template before adding.
 * @param numDoubleBeforeAdd Number of bimolecular species before adding.
 */
void parse_input_for_add_file(
    std::string addFileNameInput, Parameters& params,
    std::map<std::string, int>& observablesList,
    std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<MolTemplate>& molTemplateList, Membrane& membraneObject,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    int& numMolTemplateBeforeAdd, int& numDoubleBeforeAdd) {
  std::cout << "This is a restart simulation with add file: "
            << addFileNameInput << std::endl;

  int numForwardRxnBdeforeAdd{0};
  int numBackRxnBdeforeAdd{0};
  int numCreatDestructRxnBdeforeAdd{0};
  int tempLastStateIndexBeforeAdd{0};
  int tempLastStateIndexAfterAdd{0};
  int numStateAdd{0};  ///< num of added states

  numMolTemplateBeforeAdd = molTemplateList.size();
  numForwardRxnBdeforeAdd = static_cast<int>(forwardRxns.size());
  numBackRxnBdeforeAdd = static_cast<int>(backRxns.size());
  numCreatDestructRxnBdeforeAdd = static_cast<int>(createDestructRxns.size());

  if (numMolTemplateBeforeAdd > 0) {
    for (int forwardRxnIndex{0}; forwardRxnIndex < numForwardRxnBdeforeAdd;
         forwardRxnIndex++) {
      ForwardRxn oneRxn{};
      oneRxn = forwardRxns[forwardRxnIndex];
      if (oneRxn.rxnType == ReactionType::bimolecular) {
        numDoubleBeforeAdd++;
      }
    }
  }

  tempLastStateIndexBeforeAdd =
      molTemplateList.back().interfaceList.back().stateList.back().index;
  std::vector<TransmissionRxn> transmissionRxns{};
  parse_input_for_add(addFileNameInput, params, observablesList, forwardRxns,
                      backRxns, createDestructRxns, transmissionRxns, molTemplateList,
                      membraneObject, numDoubleBeforeAdd);
  // Check the implicit lipid is the first, and unpdate mol.molTypeIndex
  for (auto& tempMolTemplate : molTemplateList) {
    if (tempMolTemplate.isImplicitLipid == true &&
        tempMolTemplate.molTypeIndex != 0) {
      error("Implicit Lipid must be the first molecule type!");
    }
  }

  tempLastStateIndexAfterAdd =
      molTemplateList.back().interfaceList.back().stateList.back().index;

  numStateAdd = tempLastStateIndexAfterAdd - tempLastStateIndexBeforeAdd;

  MolTemplate::numMolTypes = molTemplateList.size();

  unsigned long reservation{};
  for (auto& molTemp : molTemplateList) reservation += molTemp.copies;
  moleculeList.reserve(reservation);
  complexList.reserve(reservation);

  // Update the iface.index for all molecules that the index is larger than
  // tempLastStateIndexBeforeAdd
  for (auto& tempMol : moleculeList) {
    for (auto& tempIface : tempMol.interfaceList) {
      if (tempIface.index > tempLastStateIndexBeforeAdd)
        tempIface.index += numStateAdd;
    }
  }

  // Update the index of reacts and products of forward reactions
  for (size_t forwardReactIndex{0}; forwardReactIndex < forwardRxns.size();
       forwardReactIndex++) {
    if (forwardReactIndex < numForwardRxnBdeforeAdd) {
      if (forwardRxns[forwardReactIndex].rxnType == ReactionType::bimolecular) {
        // This is a old forward biomolecule react, update the products'
        // index by + numStateAdd
        for (auto& tempProduct :
             forwardRxns[forwardReactIndex].productListNew) {
          tempProduct.absIfaceIndex += numStateAdd;
        }
      }
    }
    if (forwardReactIndex >= numForwardRxnBdeforeAdd) {
      if (forwardRxns[forwardReactIndex].rxnType == ReactionType::bimolecular) {
        // This is a new forward biomolecule react, update the products'
        // index by + numDoubleBeforeAdd
        for (auto& tempProduct :
             forwardRxns[forwardReactIndex].productListNew) {
          tempProduct.absIfaceIndex += numDoubleBeforeAdd;
        }
      }
    }
  }

  // Update the index of reacts and products of back reactions
  for (size_t backReactIndex{0}; backReactIndex < backRxns.size();
       backReactIndex++) {
    if (backReactIndex < numBackRxnBdeforeAdd) {
      // This is a old back react, update the reactants' index by +
      // numStateAdd
      for (auto& tempReactant : backRxns[backReactIndex].reactantListNew) {
        tempReactant.absIfaceIndex += numStateAdd;
      }
    } else {
      // This is a new back react, update the reactants' index by +
      // numDoubleBeforeAdd
      for (auto& tempReactant : backRxns[backReactIndex].reactantListNew) {
        tempReactant.absIfaceIndex += numDoubleBeforeAdd;
      }
    }
  }

  // Generate the coordinates, write out coordinate and topology files for
  // added molecules
  generate_coordinates_for_restart(
      params, moleculeList, complexList, molTemplateList, forwardRxns,
      membraneObject, numMolTemplateBeforeAdd, numForwardRxnBdeforeAdd);
  for (auto& tmpComplex : complexList) {
    tmpComplex.numEachMol.clear();
    tmpComplex.numEachMol.resize(molTemplateList.size());
    for (auto& memMol : tmpComplex.memberList)
      ++tmpComplex.numEachMol[moleculeList[memMol].molTypeIndex];

    tmpComplex.lastNumberUpdateItrEachMol.resize(molTemplateList.size());
  }

  write_psf(params, moleculeList, molTemplateList);
}