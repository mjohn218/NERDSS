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

void perform_bimolecular_reactions(
    unsigned simItr, Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, SimulVolume& simulVolume,
    std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays,
    Membrane& membraneObject, std::vector<int>& region, bool isLeft, std::ofstream& debugFile, int& associationCount) {
  for (auto& molItr : region) {
    // Only continue if the molecule actually exists, and isn't implicit-lipid, and can react
    if (moleculeList[molItr].isEmpty || moleculeList[molItr].isImplicitLipid || moleculeList[molItr].isDissociated || moleculeList[molItr].isAssociated)
      continue;

    if (moleculeList[molItr].isRightEdge || moleculeList[molItr].isRightGhost) {
      continue;
    }

    if (moleculeList[molItr].crossbase.size() > 0) {
      // Evaluate whether to perform a reaction with protein, and with whom.
      // Flag=1 means reaction is performed. Returns correct ci1 and ci2 for
      // this rxn. Loop over all reactions individually, instead of summing
      // probabilities
      //
      // these are indices in crossbase/mycrossint/crossrxn of the reaction
      // for molecules 1 and 2, should it occur
      int crossIndex1{0};
      int crossIndex2{0};
      bool willReact{determine_if_reaction_occurs(
          crossIndex1, crossIndex2, Constants::iRandMax, moleculeList[molItr],
          moleculeList, forwardRxns)};

      if (willReact) {
        // This molecule will perform a bimolecular reaction
        // either physically associate two molecules into a complex (A+B->AB)
        // or change state of one (or both) reactants (A+B->A+B')
        //
        int molItr2{moleculeList[molItr].crossbase[crossIndex1]};
        int ifaceIndex1{moleculeList[molItr].mycrossint[crossIndex1]};
        int ifaceIndex2;
        if (moleculeList[molItr2].isImplicitLipid == false) {
          ifaceIndex2 = moleculeList[molItr2].mycrossint[crossIndex2];
        } else {
          ifaceIndex2 = 0;
        }
        std::array<int, 3> rxnIndex =
            moleculeList[molItr].crossrxn[crossIndex1];

        if (moleculeList[molItr].isDissociated || moleculeList[molItr].isAssociated || moleculeList[molItr2].isDissociated || moleculeList[molItr2].isAssociated) {
          continue;
        }

        // First if statement is to determine if reactants are physically
        // associating
        if (forwardRxns[rxnIndex[0]].rxnType == ReactionType::bimolecular) {
          if (moleculeList[molItr2].isImplicitLipid) {
            if (molTemplateList
                    [forwardRxns[rxnIndex[0]].reactantListNew[0].molTypeIndex]
                        .isImplicitLipid ==
                false) {  // IL is listed second as the reactant.
              associate_implicitlipid(
                  ifaceIndex1, ifaceIndex2, moleculeList[molItr],
                  moleculeList[molItr2],
                  complexList[moleculeList[molItr].myComIndex],
                  complexList[moleculeList[molItr2].myComIndex], params,
                  forwardRxns[rxnIndex[0]], moleculeList, molTemplateList,
                  observablesList, counterArrays, complexList, membraneObject,
                  forwardRxns, backRxns);
            } else {  // IL is listed first as the reactant.
              associate_implicitlipid(
                  ifaceIndex2, ifaceIndex1, moleculeList[molItr2],
                  moleculeList[molItr],
                  complexList[moleculeList[molItr2].myComIndex],
                  complexList[moleculeList[molItr].myComIndex], params,
                  forwardRxns[rxnIndex[0]], moleculeList, molTemplateList,
                  observablesList, counterArrays, complexList, membraneObject,
                  forwardRxns, backRxns);
            }
          } else {
            if (moleculeList[molItr].interfaceList[ifaceIndex1].index ==
                forwardRxns[rxnIndex[0]].reactantListNew[0].absIfaceIndex) {
              associate(simItr, ifaceIndex1, ifaceIndex2, moleculeList[molItr],
                        moleculeList[molItr2],
                        complexList[moleculeList[molItr].myComIndex],
                        complexList[moleculeList[molItr2].myComIndex], params,
                        forwardRxns[rxnIndex[0]], moleculeList, molTemplateList,
                        observablesList, counterArrays, complexList,
                        membraneObject, forwardRxns, backRxns, debugFile, associationCount);
            } else {
              associate(simItr, ifaceIndex2, ifaceIndex1, moleculeList[molItr2],
                        moleculeList[molItr],
                        complexList[moleculeList[molItr2].myComIndex],
                        complexList[moleculeList[molItr].myComIndex], params,
                        forwardRxns[rxnIndex[0]], moleculeList, molTemplateList,
                        observablesList, counterArrays, complexList,
                        membraneObject, forwardRxns, backRxns, debugFile, associationCount);
            }
          }
        } else if (forwardRxns[rxnIndex[0]].rxnType ==
                   ReactionType::biMolStateChange) {
          size_t facilMolIndex{static_cast<size_t>(molItr)};
          int facilIfaceIndex{ifaceIndex1};
          int facilComIndex{moleculeList[molItr].myComIndex};
          int stateMolIndex{molItr2};
          int stateIfaceIndex{ifaceIndex2};
          int stateComIndex{moleculeList[molItr2].myComIndex};

          // figure out if the current molecule is the facilitator molecule or
          // the molecule which has its interface state changed
          if (moleculeList[molItr].molTypeIndex !=
              forwardRxns[rxnIndex[0]].reactantListNew[0].molTypeIndex) {
            facilMolIndex = molItr2;
            facilIfaceIndex = ifaceIndex2;
            facilComIndex = moleculeList[molItr2].myComIndex;
            stateMolIndex = molItr;
            stateComIndex = moleculeList[molItr].myComIndex;
            stateIfaceIndex = ifaceIndex1;
          }

          if (moleculeList[molItr2].isImplicitLipid) {
            // In this case, one implicit lipid changes state.
            perform_implicitlipid_state_change(
                stateIfaceIndex, facilIfaceIndex, rxnIndex,
                moleculeList[stateMolIndex], moleculeList[facilMolIndex],
                complexList[stateComIndex], complexList[facilComIndex],
                counterArrays, params, forwardRxns, backRxns, moleculeList,
                complexList, molTemplateList, observablesList, membraneObject);
          } else {
            perform_bimolecular_state_change(
                stateIfaceIndex, facilIfaceIndex, rxnIndex,
                moleculeList[stateMolIndex], moleculeList[facilMolIndex],
                complexList[stateComIndex], complexList[facilComIndex],
                counterArrays, params, forwardRxns, backRxns, moleculeList,
                complexList, molTemplateList, observablesList, membraneObject);
          }
        } else {
          error(
              "Attemping bimolecular reaction which has no reaction type. "
              "Exiting...");
        }
      } else {
        // No reaction was chosen for this molecule
        if (moleculeList[molItr].trajStatus ==
            TrajStatus::none) {  // diffusing a complex, while all its members
                                 // diffused
          // Propagate complex with random translational and rotational motion
          create_complex_propagation_vectors(
              params, complexList[moleculeList[molItr].myComIndex],
              moleculeList, complexList, molTemplateList, membraneObject);
          for (auto& memMol :
               complexList[moleculeList[molItr].myComIndex].memberList)
            moleculeList[memMol].trajStatus = TrajStatus::canBeResampled;
        }
        // Set probability of this protein to zero in all reactions so it
        // doesn't try to react again but the partners still will avoid
        // overlapping.
        for (unsigned crossItr{0};
             crossItr < moleculeList[molItr].crossbase.size(); ++crossItr) {
          int skipMol{moleculeList[molItr].crossbase[crossItr]};
          for (unsigned crossItr2{0};
               crossItr2 < moleculeList[skipMol].crossbase.size();
               ++crossItr2) {
            if (moleculeList[skipMol].crossbase[crossItr2] ==
                moleculeList[molItr].index)
              moleculeList[skipMol].probvec[crossItr2] = 0;
          }
        }
      }

    } else {  // else, moleculeList[molItr].crossbase.size() == 0
      // This molecule shoud diffuse
      // this protein has ncross=0,
      // meaning it neither dissociated nor tried to associate.
      // however, it could have movestat=2 if it is part of a multi-protein
      // complex that already displaced.
      //
      if (moleculeList[molItr].trajStatus ==
          TrajStatus::none) {  // No actions on this molecule so far, so it
                               // can diffuse
        // Move the whole complex this molecule belongs to:
        create_complex_propagation_vectors(
            params, complexList[moleculeList[molItr].myComIndex], moleculeList,
            complexList, molTemplateList, membraneObject);
        for (auto& memMol :
             complexList[moleculeList[molItr].myComIndex].memberList) {
          // Need to be checked against overlapping with others:
          moleculeList[memMol].trajStatus = TrajStatus::canBeResampled;
        }
      }
    }
  }  // done testing all molecules for bimolecular reactions
}
