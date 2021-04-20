#include "math/constants.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "tracing.hpp"

void check_for_unimolstatechange_reactions(unsigned simItr, Parameters& params, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    SimulVolume& simulVolume, const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    const std::vector<CreateDestructRxn>& createDestructRxns, std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays, Membrane& membraneObject)
{
    for (unsigned molItr { 0 }; molItr < moleculeList.size(); ++molItr) {
        // only do checks if the Molecule exists
        if (moleculeList[molItr].isEmpty)
            continue;

        //check for state changes
        if (params.hasUniMolStateChange == true) {
            if (molTemplateList[moleculeList[molItr].molTypeIndex].ifacesWithStates.empty() == false) {
                for (int relIndex : molTemplateList[moleculeList[molItr].molTypeIndex].ifacesWithStates) {
                    const auto& stateList = molTemplateList[moleculeList[molItr].molTypeIndex].interfaceList[relIndex].stateList;
                    if (moleculeList[molItr].isImplicitLipid == true) {
                        if (moleculeList[molItr].interfaceList[relIndex].isBound) {
                            // TODO
                        } else {
                            for (int tmpImplicitLipidStateIndex = 0; tmpImplicitLipidStateIndex < membraneObject.nStates; tmpImplicitLipidStateIndex++) {
                                int tmpOneStateLipid = membraneObject.numberOfFreeLipidsEachState[tmpImplicitLipidStateIndex];
                                for (int tmpItr = 0; tmpItr < tmpOneStateLipid; tmpItr++) {
                                    moleculeList[molItr].interfaceList[relIndex].index += tmpImplicitLipidStateIndex;
                                    moleculeList[molItr].interfaceList[relIndex].relIndex += tmpImplicitLipidStateIndex;
                                    moleculeList[molItr].interfaceList[relIndex].stateIndex += tmpImplicitLipidStateIndex;

                                    int rxnIndex { -1 };
                                    int rateIndex { -1 };
                                    bool isStateChangeBackRxn { false };
                                    find_which_state_change_reaction(relIndex, rxnIndex, rateIndex, isStateChangeBackRxn,
                                        moleculeList[molItr], molTemplateList[moleculeList[molItr].molTypeIndex].interfaceList[0].stateList[tmpImplicitLipidStateIndex], forwardRxns, backRxns);

                                    if (rxnIndex != -1 && rateIndex != -1) {
                                        double rate { (!isStateChangeBackRxn) ? forwardRxns[rxnIndex].rateList[rateIndex].rate
                                                                              : backRxns[rxnIndex].rateList[rateIndex].rate };
                                        double prob { 1 - exp(-rate * params.timeStep * Constants::usToSeconds) };

                                        // make sure that the time step is resonable according to the prob of reaction
                                        if (prob > 1.000001) {
                                            std::cerr << "Error: prob of reaction > 1. Avoid this using a smaller time step." << std::endl;
                                            exit(1);
                                        }
                                        if (prob > 0.5) {
                                            // std::cout << "WARNING: prob of reaction > 0.5. If this is a reaction for a bimolecular binding with multiple binding sites, please use a smaller time step." << std::endl;
                                        }

                                        if (prob > rand_gsl()) {
                                            // to get higher resolution random numbers
                                            // double rnum2 { rnum + rand_gsl() * Constants::iRandMax };

                                            // if (prob > rnum2) {
                                            // std::cout << "State change at iteration " << simItr << ". Molecule "
                                            //           << moleculeList[molItr].index << " changing state from " << molTemplateList[moleculeList[molItr].molTypeIndex].interfaceList[0].stateList[tmpImplicitLipidStateIndex].iden
                                            //           << '\n';

                                            const auto& newState = (!isStateChangeBackRxn) ? forwardRxns[rxnIndex].productListNew[0]
                                                                                           : backRxns[rxnIndex].productListNew[0];
                                            int relStateIndex { -1 };
                                            for (auto& state : stateList) {
                                                if (state.index == newState.absIfaceIndex) {
                                                    relStateIndex = static_cast<int>(&state - &stateList[0]);
                                                    break;
                                                }
                                            }

                                            // check observables
                                            bool isObserved { false };
                                            std::string observeLabel {};
                                            if (!isStateChangeBackRxn) {
                                                isObserved = forwardRxns[rxnIndex].isObserved;
                                                observeLabel = forwardRxns[rxnIndex].observeLabel;
                                            } else {
                                                isObserved = backRxns[rxnIndex].isObserved;
                                                observeLabel = backRxns[rxnIndex].observeLabel;
                                            }
                                            if (isObserved) {
                                                auto observeItr = observablesList.find(observeLabel);
                                                if (observeItr == observablesList.end()) {
                                                    // std::cerr << "WARNING: Observable " << observeLabel << " not defined.\n";
                                                } else {
                                                    ++observeItr->second;
                                                }
                                            }

                                            --counterArrays.copyNumSpecies[tmpImplicitLipidStateIndex];
                                            ++counterArrays.copyNumSpecies[newState.absIfaceIndex];
                                            --membraneObject.numberOfFreeLipidsEachState[tmpImplicitLipidStateIndex];
                                            ++membraneObject.numberOfFreeLipidsEachState[relStateIndex];
                                            // }
                                        }
                                    }
                                    moleculeList[molItr].interfaceList[relIndex].index -= tmpImplicitLipidStateIndex;
                                    moleculeList[molItr].interfaceList[relIndex].relIndex -= tmpImplicitLipidStateIndex;
                                    moleculeList[molItr].interfaceList[relIndex].stateIndex -= tmpImplicitLipidStateIndex;
                                }
                            }
                        }
                    } else {
                        if (moleculeList[molItr].trajStatus == TrajStatus::none) {
                            auto stateItr = std::find_if(stateList.begin(), stateList.end(), [&](const Interface::State& oneState) {
                                return oneState.iden == moleculeList[molItr].interfaceList[relIndex].stateIden;
                            });
                            if (stateItr == stateList.end()) {
                                std::cerr << "Error, state is invalid. Exiting...\n";
                                exit(1);
                            }

                            if (moleculeList[molItr].interfaceList[relIndex].isBound) {
                                // TODO
                            } else {
                                int rxnIndex { -1 };
                                int rateIndex { -1 };
                                bool isStateChangeBackRxn { false };
                                find_which_state_change_reaction(relIndex, rxnIndex, rateIndex, isStateChangeBackRxn,
                                    moleculeList[molItr], *stateItr, forwardRxns, backRxns);

                                if (rxnIndex != -1 && rateIndex != -1) {
                                    double rate { (!isStateChangeBackRxn) ? forwardRxns[rxnIndex].rateList[rateIndex].rate
                                                                          : backRxns[rxnIndex].rateList[rateIndex].rate };
                                    double prob { 1 - exp(-rate * params.timeStep * Constants::usToSeconds) };

                                    // make sure that the time step is resonable according to the prob of reaction
                                    if (prob > 1.000001) {
                                        std::cerr << "Error: prob of reaction > 1. Avoid this using a smaller time step." << std::endl;
                                        exit(1);
                                    }
                                    if (prob > 0.5) {
                                        // std::cout << "WARNING: prob of reaction > 0.5. If this is a reaction for a bimolecular binding with multiple binding sites, please use a smaller time step." << std::endl;
                                    }

                                    if (prob > rand_gsl()) {
                                        // to get higher resolution random numbers
                                        // double rnum2 { rnum + rand_gsl() * Constants::iRandMax };

                                        // if (prob > rnum2) {
                                        // std::cout << "State change at iteration " << simItr << ". Molecule "
                                        //           << moleculeList[molItr].index << " changing state from " << stateItr->iden
                                        //           << '\n';

                                        const auto& newState = (!isStateChangeBackRxn) ? forwardRxns[rxnIndex].productListNew[0]
                                                                                       : backRxns[rxnIndex].productListNew[0];
                                        int relStateIndex { -1 };
                                        for (auto& state : stateList) {
                                            if (state.index == newState.absIfaceIndex) {
                                                relStateIndex = static_cast<int>(&state - &stateList[0]);
                                                break;
                                            }
                                        }

                                        // check observables
                                        bool isObserved { false };
                                        std::string observeLabel {};
                                        if (!isStateChangeBackRxn) {
                                            isObserved = forwardRxns[rxnIndex].isObserved;
                                            observeLabel = forwardRxns[rxnIndex].observeLabel;
                                        } else {
                                            isObserved = backRxns[rxnIndex].isObserved;
                                            observeLabel = backRxns[rxnIndex].observeLabel;
                                        }
                                        if (isObserved) {
                                            auto observeItr = observablesList.find(observeLabel);
                                            if (observeItr == observablesList.end()) {
                                                // std::cerr << "WARNING: Observable " << observeLabel << " not defined.\n";
                                            } else {
                                                ++observeItr->second;
                                            }
                                        }

                                        --counterArrays.copyNumSpecies[stateItr->index];
                                        ++counterArrays.copyNumSpecies[newState.absIfaceIndex];

                                        /*State change happens here.*/
                                        moleculeList[molItr].interfaceList[relIndex].change_state(
                                            relStateIndex, newState.absIfaceIndex, newState.requiresState);

                                        moleculeList[molItr].trajStatus = TrajStatus::propagated;
                                        complexList[moleculeList[molItr].myComIndex].trajStatus = TrajStatus::propagated;
                                        break; // break out of interface loop
                                        // }
                                    }
                                }
                            }
                        }
                    } // only try each pair dissociating once
                } // finished trying out stateChange reactions for this proteins
            }
        }
    }
}
