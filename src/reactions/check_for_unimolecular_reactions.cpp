#include "math/constants.hpp"
//#include "math/math_functions.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"

void check_for_unimolecular_reactions(unsigned simItr, Parameters& params, std::vector<int>& emptyMolList,
                                      std::vector<int>& emptyComList, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
                                      SimulVolume& simulVolume, const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
                                      const std::vector<CreateDestructRxn>& createDestructRxns, const std::vector<MolTemplate>& molTemplateList,
                                      std::map<std::string, int>& observablesList, copyCounters& counterArrays, const Membrane &membraneObject, std::vector<double> & IL2DbindingVec, std::vector<double> & IL2DUnbindingVec, std::vector<double> &ILTableIDs)
{
    // Note: cannot use a for-range or vector iterator loop here, since we are changing the moleculeList vector
    // while iterating over it
    for (unsigned molItr { 0 }; molItr < moleculeList.size(); ++molItr) {
        // only do checks if the Molecule exists
        if (moleculeList[molItr].isEmpty || moleculeList[molItr].isImplicitLipid)
            continue;
        // first check for state changes
        for (int relIndex : molTemplateList[moleculeList[molItr].molTypeIndex].ifacesWithStates) {
            const auto& stateList = molTemplateList[moleculeList[molItr].molTypeIndex].interfaceList[relIndex].stateList;
            
            if (moleculeList[molItr].trajStatus == TrajStatus::none) {
                auto stateItr = std::find_if(stateList.begin(), stateList.end(), [&](const Interface::State& oneState) {
                    return oneState.iden == moleculeList[molItr].interfaceList[relIndex].stateIden;
                });
                if (stateItr == stateList.end()) {
                    std::cerr << "Error, state is invalid. Exiting...\n";
                    exit(1);
                }

                if (moleculeList[molItr].interfaceList[0].isBound) {
                    // TODO
                }

                int rxnIndex { -1 };
                int rateIndex { -1 };
                bool isStateChangeBackRxn { false };
                find_which_state_change_reaction(relIndex, rxnIndex, rateIndex, isStateChangeBackRxn,
                    moleculeList[molItr], *stateItr, forwardRxns, backRxns);

                if (rxnIndex != -1 && rateIndex != -1) {
                    double rate { (!isStateChangeBackRxn) ? forwardRxns[rxnIndex].rateList[rateIndex].rate
                                                          : backRxns[rxnIndex].rateList[rateIndex].rate };
                    double prob { 1 - exp(-rate * params.timeStep * Constants::usToSeconds) };
                    double rnum { 1.0 * rand_gsl() };

                    if (prob > rnum) {
                        // to get higher resolution random numbers
                        double rnum2 { rnum + rand_gsl() * Constants::iRandMax };

                        if (prob > rnum2) {
                            std::cout << "State change at iteration " << simItr << ". Molecule "
                                      << moleculeList[molItr].index << " changing state from " << stateItr->iden
                                      << '\n';

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
                                    std::cerr << "WARNING: Observable " << observeLabel << " not defined.\n";
                                } else {
                                    ++observeItr->second;
                                }
                            }

                            --counterArrays.copyNumSpecies[stateItr->index];
                            ++counterArrays.copyNumSpecies[newState.absIfaceIndex];
                            moleculeList[molItr].interfaceList[relIndex].change_state(
                                relStateIndex, newState.absIfaceIndex, newState.requiresState);

                            moleculeList[molItr].trajStatus = TrajStatus::propagated;
                            complexList[moleculeList[molItr].myComIndex].trajStatus = TrajStatus::propagated;
                            break; // break out of interface loop
                        }
                    }
                }
            } // only try each pair dissociating once
        } // finished trying out stateChange reactions for this proteins

        // now check for unimolecular creation and destruction
        // TODO: this is currently constrained to only consider a Molecule and only have one rate
        for (const auto& oneRxn : createDestructRxns) {
            if (oneRxn.rxnType == ReactionType::destruction && moleculeList[molItr].trajStatus == TrajStatus::none
                && isReactant(moleculeList[molItr], complexList[moleculeList[molItr].myComIndex], oneRxn, moleculeList) ) {
                double prob { 1 - exp(-oneRxn.rateList[0].rate * params.timeStep * Constants::usToSeconds) };
                double rNum { 1.0 * rand_gsl() };
                if (prob > rNum) {
                    double rNum2 { rNum + rand_gsl() * Constants::iRandMax }; // to get higher resolution
                    if (prob > rNum2) {
                        std::cout << "Destroying molecule of type "
                                  << molTemplateList[moleculeList[molItr].molTypeIndex].molName << " with reaction "
                                  << oneRxn.absRxnIndex << " at iteration " << simItr << '\n';

                        // check if observable, before destroying
                        //                        for (auto& observable : observablesList) {
                        //                            if (observable == moleculeList[molItr])
                        //                                --observable.currNum;
                        //                        }

                        // decrement the copy number array for everything in complex
                        for (auto& memMol : complexList[moleculeList[molItr].myComIndex].memberList) {
                            for (auto& iface : moleculeList[memMol].interfaceList) {
                                --counterArrays.copyNumSpecies[iface.index];
                            }
                        }

                        complexList[moleculeList[molItr].myComIndex].destroy(moleculeList, emptyMolList, emptyComList);
                        // remove the molecule from the SimulVolume subsCellList
                        // have this here to avoid circular header calls with SimulVolume and Molecule_Complex
                        for (auto itr
                             = simulVolume.subCellList[moleculeList[molItr].mySubVolIndex].memberMolList.begin();
                             itr != simulVolume.subCellList[moleculeList[molItr].mySubVolIndex].memberMolList.end();
                             ++itr) {
                            if (*itr == moleculeList[molItr].index) {
                                simulVolume.subCellList[moleculeList[molItr].mySubVolIndex].memberMolList.erase(itr);
                                break;
                            }
                        }
                        moleculeList[molItr].mySubVolIndex = -1;

                        if (oneRxn.isObserved) {
                            auto observeItr = observablesList.find(oneRxn.observeLabel);
                            if (observeItr == observablesList.end()) {
                                std::cerr << "WARNING: Observable " << oneRxn.observeLabel << " not defined.\n";
                            } else {
                                --observeItr->second;
                            }
                        }

                        break; // don't do anything else with the molecule this timestep
                    }
                }
            } else if (oneRxn.rxnType == ReactionType::uniMolCreation && moleculeList[molItr].trajStatus == TrajStatus::none
                       && isReactant(moleculeList[molItr], complexList[moleculeList[molItr].myComIndex], oneRxn, moleculeList))
            {   
                double rNum { 1.0 * rand_gsl() };
                long double lambda { oneRxn.rateList.at(0).rate * params.timeStep * Constants::usToSeconds };
                long double prob { 1 - exp(-lambda) };
                if (prob > rNum) {
                    int newMolIndex { 0 };
                    int newComIndex { 0 };
                    std::cout << "Creating molecule of type " << oneRxn.productMolList.back().molName
                              << " from reaction " << oneRxn.absRxnIndex << '\n';
                    create_molecule_and_complex_from_rxn(molItr, newMolIndex, newComIndex, true,
                        molTemplateList[oneRxn.productMolList.back().molTypeIndex], params, emptyMolList, emptyComList,
							 oneRxn, simulVolume, moleculeList, complexList, molTemplateList, forwardRxns, membraneObject);
                    moleculeList[molItr].trajStatus = TrajStatus::propagated;
                    complexList[moleculeList[molItr].myComIndex].trajStatus = TrajStatus::propagated;

                    // update the copy number array
                    for (const auto& iface : moleculeList[newMolIndex].interfaceList)
                        ++counterArrays.copyNumSpecies[iface.index];

                    if (oneRxn.isObserved) {
                        auto observeItr = observablesList.find(oneRxn.observeLabel);
                        if (observeItr == observablesList.end()) {
                            std::cerr << "WARNING: Observable " << oneRxn.observeLabel << " not defined.\n";
                        } else {
                            ++observeItr->second;
                        }
                    }
                    break; // don't do anything else with the molecule this timestep
                }
            }
        }

        // now check for dissociation
        check_dissociation(simItr, params, simulVolume, molTemplateList, observablesList, molItr, emptyMolList,
			   emptyComList, moleculeList, complexList, backRxns, forwardRxns, createDestructRxns, counterArrays, membraneObject);
	      
	check_dissociation_implicitlipid(simItr, params, simulVolume, molTemplateList, observablesList, molItr, emptyMolList,
					 emptyComList, moleculeList, complexList, backRxns, forwardRxns, createDestructRxns, counterArrays, membraneObject, IL2DbindingVec, IL2DUnbindingVec, ILTableIDs);
    }
}
