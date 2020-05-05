#include "parser/parser_functions.hpp"

void populate_reaction_lists(const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    const std::vector<CreateDestructRxn>& createDestructRxns, std::vector<MolTemplate>& molTemplateList)
{
    for (const auto& forwardRxn : forwardRxns) {
        if (forwardRxn.rxnType == ReactionType::bimolecular) {
            // Set first interface index
            Interface& reactIface1 = molTemplateList[forwardRxn.reactantListNew[0].molTypeIndex]
                                         .interfaceList[forwardRxn.reactantListNew[0].relIfaceIndex];

            unsigned stateIndex1 { 0 };
            if (reactIface1.stateList.size() > 1) {
                // if there are more than one possible states of the interface, find which one the reaction uses
                for (unsigned stateIndex { 0 }; stateIndex < reactIface1.stateList.size(); ++stateIndex) {
                    if (reactIface1.stateList[stateIndex].index == forwardRxn.reactantListNew[0].absIfaceIndex) {
                        reactIface1.stateList[stateIndex].myForwardRxns.push_back(forwardRxn.relRxnIndex);
                        stateIndex1 = stateIndex;
                        break;
                    }
                }
            } else {
                // if not, just add it to the list on the first, default, state
                reactIface1.stateList[0].myForwardRxns.push_back(forwardRxn.relRxnIndex);
            }

            // If the reaction is non-symmetric, set the second interface index as well
            if (!forwardRxn.isSymmetric) {
                Interface& reactIface2 = molTemplateList[forwardRxn.reactantListNew[1].molTypeIndex]
                                             .interfaceList[forwardRxn.reactantListNew[1].relIfaceIndex];

                if (reactIface2.stateList.size() > 1) {
                    // if there are more than one possible states of the interface, find which one the reaction uses
                    for (unsigned stateIndex { 0 }; stateIndex < reactIface2.stateList.size(); ++stateIndex) {
                        if (reactIface2.stateList[stateIndex].index == forwardRxn.reactantListNew[1].absIfaceIndex) {
                            reactIface2.stateList[stateIndex].myForwardRxns.push_back(forwardRxn.relRxnIndex);
                            reactIface2.stateList[stateIndex].rxnPartners.push_back(
                                reactIface1.stateList[stateIndex1].index);
                            reactIface1.stateList[stateIndex1].rxnPartners.push_back(
                                reactIface2.stateList[stateIndex].index);
                            break;
                        }
                    }
                } else {
                    // if not, just add it to the list on the first, default, state
                    reactIface2.stateList[0].myForwardRxns.push_back(forwardRxn.relRxnIndex);

                    reactIface1.stateList[stateIndex1].rxnPartners.push_back(reactIface2.stateList[0].index);
                    reactIface2.stateList[0].rxnPartners.push_back(reactIface1.stateList[stateIndex1].index);
                }
            } else {
                reactIface1.stateList[stateIndex1].rxnPartners.push_back(forwardRxn.reactantListNew[0].absIfaceIndex);
            }
        } else if (forwardRxn.rxnType == ReactionType::biMolStateChange) {
            /* BIMOLECULAR STATE CHANGE REACTIONS */
            Interface& reactIface0 = molTemplateList[forwardRxn.reactantListNew[0].molTypeIndex]
                                         .interfaceList[forwardRxn.reactantListNew[0].relIfaceIndex];
            Interface& reactIface1 = molTemplateList[forwardRxn.reactantListNew[1].molTypeIndex]
                                         .interfaceList[forwardRxn.reactantListNew[1].relIfaceIndex];

            for (auto& oneState : reactIface1.stateList) {
                if (oneState.index == forwardRxn.reactantListNew[1].absIfaceIndex) {
                    oneState.stateChangeRxns.emplace_back(forwardRxn.relRxnIndex, forwardRxn.conjBackRxnIndex);
                    oneState.myForwardRxns.emplace_back(forwardRxn.relRxnIndex);
                    oneState.rxnPartners.emplace_back(forwardRxn.reactantListNew[0].absIfaceIndex);
                } else if (forwardRxn.isReversible && oneState.index == forwardRxn.productListNew[1].absIfaceIndex) {
                    oneState.stateChangeRxns.emplace_back(forwardRxn.conjBackRxnIndex, forwardRxn.relRxnIndex);
                    oneState.myForwardRxns.emplace_back(forwardRxn.relRxnIndex);
                    oneState.rxnPartners.emplace_back(forwardRxn.productListNew[0].absIfaceIndex);
                }
            }

            for (auto& oneState : reactIface0.stateList) {
                if (oneState.index == forwardRxn.reactantListNew[0].absIfaceIndex) {
                    oneState.stateChangeRxns.emplace_back(forwardRxn.relRxnIndex, forwardRxn.conjBackRxnIndex);
                    oneState.myForwardRxns.emplace_back(forwardRxn.relRxnIndex);
                    oneState.rxnPartners.emplace_back(forwardRxn.reactantListNew[1].absIfaceIndex);
                } else if (forwardRxn.isReversible && oneState.index == forwardRxn.productListNew[0].absIfaceIndex) {
                    oneState.stateChangeRxns.emplace_back(forwardRxn.conjBackRxnIndex, forwardRxn.relRxnIndex);
                    oneState.myForwardRxns.emplace_back(forwardRxn.relRxnIndex);
                    oneState.rxnPartners.emplace_back(forwardRxn.productListNew[1].absIfaceIndex);
                }
            }

            // don't add it to the facilitator's list of state change reactions
            //            for (auto& oneState : reactIface0.stateList) {
            //                if (oneState.index == forwardRxn.reactantListNew[0].absIfaceIndex) {
            //                    oneState.rxnPartners.emplace_back(forwardRxn.reactantListNew[1].absIfaceIndex);
            //                } else if (forwardRxn.isReversible && oneState.index == forwardRxn.productListNew[0].absIfaceIndex) {
            //                    oneState.rxnPartners.emplace_back(forwardRxn.reactantListNew[1].absIfaceIndex);
            //                }
            //            }
            //
            //            for (auto& oneState : reactIface1.stateList) {
            //                if (oneState.index == forwardRxn.reactantListNew[1].absIfaceIndex) {
            //                    oneState.rxnPartners.emplace_back(forwardRxn.reactantListNew[0].absIfaceIndex);
            //                } else if (forwardRxn.isReversible && oneState.index == forwardRxn.productListNew[1].absIfaceIndex) {
            //                    oneState.rxnPartners.emplace_back(forwardRxn.reactantListNew[0].absIfaceIndex);
            //                }
            //            }
        } else if (forwardRxn.rxnType == ReactionType::uniMolStateChange) {
            /* UNIMOLECULAR STATE CHANGE REACTIONS */
            Interface& reactIface = molTemplateList[forwardRxn.reactantListNew[0].molTypeIndex]
                                        .interfaceList[forwardRxn.reactantListNew[0].relIfaceIndex];

            // add state change reactions only to the stateChangeRxns list, not myForwardRxns
            for (auto& oneState : reactIface.stateList) {
                if (oneState.index == forwardRxn.reactantListNew[0].absIfaceIndex)
                    oneState.stateChangeRxns.emplace_back(forwardRxn.relRxnIndex, forwardRxn.conjBackRxnIndex);
                else if (forwardRxn.isReversible && oneState.index == forwardRxn.productListNew[0].absIfaceIndex)
                    oneState.stateChangeRxns.emplace_back(forwardRxn.conjBackRxnIndex, forwardRxn.relRxnIndex);
            }
        }
    }
}

void populate_reaction_lists_for_add(const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    const std::vector<CreateDestructRxn>& createDestructRxns, std::vector<MolTemplate>& molTemplateList, int addForwardRxnNum, int addBackRxnNum, int addCreateDestructRxnNum)
{
    for (int forwardRxnsIndex { addForwardRxnNum }; forwardRxnsIndex < forwardRxns.size(); forwardRxnsIndex++) {
        ForwardRxn forwardRxn {};
        forwardRxn = forwardRxns[forwardRxnsIndex];
        if (forwardRxn.rxnType == ReactionType::bimolecular) {
            // Set first interface index

            Interface& reactIface1 = molTemplateList[forwardRxn.reactantListNew[0].molTypeIndex]
                                         .interfaceList[forwardRxn.reactantListNew[0].relIfaceIndex];

            unsigned stateIndex1 { 0 };

            for (unsigned index { 0 }; index < reactIface1.stateList[0].myForwardRxns.size(); index++) {
                std::cout << forwardRxns[reactIface1.stateList[0].myForwardRxns[index]].rateList[0].rate << "\t";
            }

            if (reactIface1.stateList.size() > 1) {
                // if there are more than one possible states of the interface, find which one the reaction uses
                for (unsigned stateIndex { 0 }; stateIndex < reactIface1.stateList.size(); ++stateIndex) {
                    if (reactIface1.stateList[stateIndex].index == forwardRxn.reactantListNew[0].absIfaceIndex) {
                        reactIface1.stateList[stateIndex].myForwardRxns.push_back(forwardRxn.relRxnIndex);
                        stateIndex1 = stateIndex;
                        break;
                    }
                }
            } else {
                // if not, just add it to the list on the first, default, state
                reactIface1.stateList[0].myForwardRxns.push_back(forwardRxn.relRxnIndex);
            }

            // If the reaction is non-symmetric, set the second interface index as well
            if (!forwardRxn.isSymmetric) {
                Interface& reactIface2 = molTemplateList[forwardRxn.reactantListNew[1].molTypeIndex]
                                             .interfaceList[forwardRxn.reactantListNew[1].relIfaceIndex];

                if (reactIface2.stateList.size() > 1) {
                    // if there are more than one possible states of the interface, find which one the reaction uses
                    for (unsigned stateIndex { 0 }; stateIndex < reactIface2.stateList.size(); ++stateIndex) {
                        if (reactIface2.stateList[stateIndex].index == forwardRxn.reactantListNew[1].absIfaceIndex) {
                            reactIface2.stateList[stateIndex].myForwardRxns.push_back(forwardRxn.relRxnIndex);
                            reactIface2.stateList[stateIndex].rxnPartners.push_back(
                                reactIface1.stateList[stateIndex1].index);
                            reactIface1.stateList[stateIndex1].rxnPartners.push_back(
                                reactIface2.stateList[stateIndex].index);
                            break;
                        }
                    }
                } else {
                    // if not, just add it to the list on the first, default, state
                    reactIface2.stateList[0].myForwardRxns.push_back(forwardRxn.relRxnIndex);

                    reactIface1.stateList[stateIndex1].rxnPartners.push_back(reactIface2.stateList[0].index);
                    reactIface2.stateList[0].rxnPartners.push_back(reactIface1.stateList[stateIndex1].index);
                }
            } else {
                reactIface1.stateList[stateIndex1].rxnPartners.push_back(forwardRxn.reactantListNew[0].absIfaceIndex);
            }

        } else if (forwardRxn.rxnType == ReactionType::biMolStateChange) {
            /* BIMOLECULAR STATE CHANGE REACTIONS */
            Interface& reactIface0 = molTemplateList[forwardRxn.reactantListNew[0].molTypeIndex]
                                         .interfaceList[forwardRxn.reactantListNew[0].relIfaceIndex];
            Interface& reactIface1 = molTemplateList[forwardRxn.reactantListNew[1].molTypeIndex]
                                         .interfaceList[forwardRxn.reactantListNew[1].relIfaceIndex];

            for (auto& oneState : reactIface1.stateList) {
                if (oneState.index == forwardRxn.reactantListNew[1].absIfaceIndex) {
                    oneState.stateChangeRxns.emplace_back(forwardRxn.relRxnIndex, forwardRxn.conjBackRxnIndex);
                    oneState.myForwardRxns.emplace_back(forwardRxn.relRxnIndex);
                    oneState.rxnPartners.emplace_back(forwardRxn.reactantListNew[0].absIfaceIndex);
                } else if (forwardRxn.isReversible && oneState.index == forwardRxn.productListNew[1].absIfaceIndex) {
                    oneState.stateChangeRxns.emplace_back(forwardRxn.conjBackRxnIndex, forwardRxn.relRxnIndex);
                    oneState.myForwardRxns.emplace_back(forwardRxn.relRxnIndex);
                    oneState.rxnPartners.emplace_back(forwardRxn.productListNew[0].absIfaceIndex);
                }
            }

            for (auto& oneState : reactIface0.stateList) {
                if (oneState.index == forwardRxn.reactantListNew[0].absIfaceIndex) {
                    oneState.stateChangeRxns.emplace_back(forwardRxn.relRxnIndex, forwardRxn.conjBackRxnIndex);
                    oneState.myForwardRxns.emplace_back(forwardRxn.relRxnIndex);
                    oneState.rxnPartners.emplace_back(forwardRxn.reactantListNew[1].absIfaceIndex);
                } else if (forwardRxn.isReversible && oneState.index == forwardRxn.productListNew[0].absIfaceIndex) {
                    oneState.stateChangeRxns.emplace_back(forwardRxn.conjBackRxnIndex, forwardRxn.relRxnIndex);
                    oneState.myForwardRxns.emplace_back(forwardRxn.relRxnIndex);
                    oneState.rxnPartners.emplace_back(forwardRxn.productListNew[1].absIfaceIndex);
                }
            }

            // don't add it to the facilitator's list of state change reactions
            //            for (auto& oneState : reactIface0.stateList) {
            //                if (oneState.index == forwardRxn.reactantListNew[0].absIfaceIndex) {
            //                    oneState.rxnPartners.emplace_back(forwardRxn.reactantListNew[1].absIfaceIndex);
            //                } else if (forwardRxn.isReversible && oneState.index == forwardRxn.productListNew[0].absIfaceIndex) {
            //                    oneState.rxnPartners.emplace_back(forwardRxn.reactantListNew[1].absIfaceIndex);
            //                }
            //            }
            //
            //            for (auto& oneState : reactIface1.stateList) {
            //                if (oneState.index == forwardRxn.reactantListNew[1].absIfaceIndex) {
            //                    oneState.rxnPartners.emplace_back(forwardRxn.reactantListNew[0].absIfaceIndex);
            //                } else if (forwardRxn.isReversible && oneState.index == forwardRxn.productListNew[1].absIfaceIndex) {
            //                    oneState.rxnPartners.emplace_back(forwardRxn.reactantListNew[0].absIfaceIndex);
            //                }
            //            }
        } else if (forwardRxn.rxnType == ReactionType::uniMolStateChange) {
            /* UNIMOLECULAR STATE CHANGE REACTIONS */
            Interface& reactIface = molTemplateList[forwardRxn.reactantListNew[0].molTypeIndex]
                                        .interfaceList[forwardRxn.reactantListNew[0].relIfaceIndex];

            // add state change reactions only to the stateChangeRxns list, not myForwardRxns
            for (auto& oneState : reactIface.stateList) {
                if (oneState.index == forwardRxn.reactantListNew[0].absIfaceIndex)
                    oneState.stateChangeRxns.emplace_back(forwardRxn.relRxnIndex, forwardRxn.conjBackRxnIndex);
                else if (forwardRxn.isReversible && oneState.index == forwardRxn.productListNew[0].absIfaceIndex)
                    oneState.stateChangeRxns.emplace_back(forwardRxn.conjBackRxnIndex, forwardRxn.relRxnIndex);
            }
        }
    }
}
