#include "parser/parser_functions.hpp"

void determine_iface_indices(int specieIndex, int& totSpecies, ParsedMol& targMol, ParsedRxn& parsedRxn,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<MolTemplate>& molTemplateList)
{
    for (auto tempItr = molTemplateList.begin(); tempItr != molTemplateList.end(); ++tempItr) {
        if (tempItr->molName == targMol.molName) {
            /* -now since we know its the same protein, let's find the indices
             * TODO: find a better way of doing this
             */

            // iterate over the target ifaces
            // TODO: change to for-range
            for (auto targIfaceItr = targMol.interfaceList.begin(); targIfaceItr != targMol.interfaceList.end();
                 ++targIfaceItr) {
                // iterate over the template ifaces
                for (auto tempIfaceItr = tempItr->interfaceList.begin(); tempIfaceItr != tempItr->interfaceList.end();
                     ++tempIfaceItr) {
                    if (tempIfaceItr->name == targIfaceItr->ifaceName) {
                        targIfaceItr->molTypeIndex = targMol.molTypeIndex;
                        targIfaceItr->speciesIndex = specieIndex;
                        targIfaceItr->relIndex = tempIfaceItr->index;

                        // only look for the index if it isn't yet known
                        // TODO: can split this up into several if statements to stop it from being so nested
                        if (targIfaceItr->absIndex == -1) {
                            // if the Interface::states vector is larger than 1, the states are named
                            if (tempIfaceItr->stateList.size() == 1) {
                                if (targIfaceItr->ifaceName == tempIfaceItr->name) {
                                    // set the index and tell the parser that the index was found
                                    targIfaceItr->absIndex = tempIfaceItr->stateList[0].index;
                                    targIfaceItr->relIndex = tempIfaceItr->index;
                                }
                            } else if (targIfaceItr->state == '\0') {
                                if (parsedRxn.rxnType != ReactionType::zerothOrderCreation
                                    && parsedRxn.rxnType != ReactionType::uniMolCreation
                                    && parsedRxn.rxnType != ReactionType::destruction) {
                                    parsedRxn.willBeMultipleRxns = true;
                                } else {
                                    std::cerr
                                        << "WARNING: Creation and Destruction reactions must have states explicitly "
                                           "stated. Using default states.\n";
                                }
                            } else {
                                for (auto& tempState : tempIfaceItr->stateList) {
                                    if (targIfaceItr->state == tempState.iden) {
                                        // set the index and tell the parser that the index was found
                                        targIfaceItr->absIndex = tempState.index;
                                        targIfaceItr->relIndex = tempIfaceItr->index;
                                        break;
                                    }
                                }
                            }

                            // if the iface is bound (and not a wildcard), we know it has to be an iface that has
                            // changed interaction
                            if (targIfaceItr->isBound && targIfaceItr->bondIndex != 0) {
                                determine_bound_iface_index(
                                    totSpecies, *targIfaceItr, parsedRxn, forwardRxns, molTemplateList);
                            }
                        }
                    }
                }
            }
        }
    }
}
