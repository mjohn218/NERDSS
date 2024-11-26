#include "boundary_conditions/reflect_functions.hpp"
#include "math/constants.hpp"
#include "math/math_functions.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "tracing.hpp"
#include <iostream>
#include <vector>

void check_for_unimolecular_reactions_population(long long int simItr, Parameters& params, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    SimulVolume& simulVolume, const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    const std::vector<CreateDestructRxn>& createDestructRxns, std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays, Membrane& membraneObject, std::ofstream& assocDissocFile)
{
    // In this routine, we treat the unimolecular reactions based on population by looping over reaction list, which is more acurate for large k*dt;
    // we need to track the molecules that can be destroyed or dissociated, hope this can be faster than the loop over molecule list

    for (auto& oneRxn : createDestructRxns) {
        // check destruction, e.g. A -> null
        if (oneRxn.rxnType == ReactionType::destruction) {
            MolTemplate& oneTemp = molTemplateList[oneRxn.reactantMolList.at(0).molTypeIndex];
            if (oneTemp.isImplicitLipid == false) { // A is not implicit lipid

                // Decide how many reactions occur based on lambda=k*dt*NA, where p(m_events)=exp(-lambda)*lamba^m /m!
                // determine the number of the A that can be destroyed
                int NA { 0 };
                NA = static_cast<int>(oneTemp.monomerList.size());

                // determine the reaction number in this time step, numEvents
                long double lambda { oneRxn.rateList.at(0).rate * params.timeStep * Constants::usToSeconds * NA };
                double explam=exp(-lambda);
				long double prob { explam };
				int numEvents { 0 };
				double rNum { 1.0 * rand_gsl() };
				long double lampow=lambda;
				long double factor = 1;
								
                while (rNum > prob) {
                    ++numEvents;
                    //prob += (exp(-lambda) * pow(lambda, numEvents)) / MathFuncs::gammFactorial(numEvents);
					/*this method is much faster, same accuracy*/
					prob += explam *lampow /factor;
					lampow = lampow*lambda;// accumulate the powers
					factor = factor * (numEvents + 1);//accumulate the factorial
				
                }
                if (numEvents > NA) {
                    // std::cout << " WARNING: CANNOT perform N destruction Events of : " << numEvents << " Not enough species exist, perform: " << NA << std::endl;
                    numEvents = NA;
                }

                // whether to use the strong destruction
                if ( oneRxn.rateList[0].rate == -1 ){ // && membraneObject.hasCompartment == true ){
                    numEvents =  oneTemp.monomerList.size() - oneTemp.copies;
                    if ( numEvents < 0 ) // for the monomers fewer than the target number, then no need desctruction but need creation(titration)
                        numEvents = 0;
                }

                // destory numEvents A molecule in this time step, the cooresponding numEvents A are chosen randomly from the pool stored in oneTemp.monomerList
                if (numEvents > 0) {
                    // choose numEvents A from the poolAList randomly
                    // std::cout << "Destroying " << numEvents << " molecule(s) of type " << oneRxn.reactantMolList.at(0).molName
                    //           << " from reaction " << oneRxn.absRxnIndex << " at iteration " << simItr << std::endl;

                    // record the index of the destroyed A
                    std::vector<int> destoryMolIndex {};
                    while (numEvents > 0) {
                        // chose one A whose TrajStatus should be None
                        int randIntNum { static_cast<int>(1.0 * NA * rand_gsl()) };
                        if (randIntNum == NA) {
                            randIntNum = NA - 1;
                        }

                        int molItr { oneTemp.monomerList[randIntNum] };

                        // avoid to destroy the already-destroyed particle
                        std::vector<int>::iterator result { std::find(std::begin(destoryMolIndex), std::end(destoryMolIndex), oneTemp.monomerList[randIntNum]) }; //check whether this A has been destroyed in this step
                        if (result == std::end(destoryMolIndex) && moleculeList[molItr].mySubVolIndex != -1)
                        {
                            // decrement the copy number array for everything in complex
                            for (auto& memMol : complexList[moleculeList[oneTemp.monomerList[randIntNum]].myComIndex].memberList) {
                                for (auto& iface : moleculeList[memMol].interfaceList) {
                                    --counterArrays.copyNumSpecies[iface.index];
                                }
                            }

                            complexList[moleculeList[oneTemp.monomerList[randIntNum]].myComIndex].destroy(moleculeList, complexList);

                            // remove the molecule from the SimulVolume subsCellList
                            // have this here to avoid circular header calls with SimulVolume and Molecule_Complex
                            for (auto itr = simulVolume.subCellList[moleculeList[molItr].mySubVolIndex].memberMolList.begin(); itr != simulVolume.subCellList[moleculeList[molItr].mySubVolIndex].memberMolList.end(); ++itr) {
                                if (*itr == moleculeList[molItr].index) {
                                    simulVolume.subCellList[moleculeList[molItr].mySubVolIndex].memberMolList.erase(itr);
                                    break;
                                }
                            }
                            moleculeList[molItr].mySubVolIndex = -1;
                            destoryMolIndex.emplace_back(oneTemp.monomerList[randIntNum]);
                            --numEvents;
                            if (oneRxn.isObserved) {
                                auto observeItr = observablesList.find(oneRxn.observeLabel);
                                if (observeItr == observablesList.end()) {
                                    // std::cerr << "WARNING: Observable " << oneRxn.observeLabel << " not defined.\n";
                                } else {
                                    --observeItr->second;
                                }
                            }
                        }
                    }
                    // remove the As from the monomerList because they are destroyed
                    for (auto oneDestroy : destoryMolIndex) {
                        oneTemp.monomerList.erase(std::find_if(oneTemp.monomerList.begin(), oneTemp.monomerList.end(), [&](const size_t& mol) { return mol == oneDestroy; }));
                    }
                }
            } else { // A is implicit lipid, do not consider the case with multiple implicit lipid states
                // Decide how many reactions occur based on lambda=k*dt*NA, where p(m_events)=exp(-lambda)*lamba^m /m!
                // determine the number of the A that can be destroyed
                int NA { 0 };
                NA = membraneObject.numberOfFreeLipidsEachState[0];

                // determine the reaction number in this time step, numEvents
                long double lambda { oneRxn.rateList.at(0).rate * params.timeStep * Constants::usToSeconds * NA };
                double explam=exp(-lambda);
				long double prob { explam };
				unsigned numEvents { 0 };
				double rNum { 1.0 * rand_gsl() };
				long double lampow=lambda;
				long double factor = 1;
				
                while (rNum > prob) {
                    ++numEvents;
                    //prob += (exp(-lambda) * pow(lambda, numEvents)) / MathFuncs::gammFactorial(numEvents);
					/*this method is much faster, same accuracy*/
					prob += explam *lampow /factor;
					lampow = lampow*lambda;// accumulate the powers
					factor = factor * (numEvents + 1);//accumulate the factorial
				
                }
                if (numEvents > NA) {
                    // std::cout << " WARNING: CANNOT perform N IL Events of : " << numEvents << " Not enough species exist, perform: " << NA << std::endl;
                    numEvents = NA;
                }
                // destory numEvents A molecule in this time step, just need to update the membraneObject.numberOfFreeLipidsEachState because A is implicit lipid
                if (numEvents > 0) {
                    // std::cout << "Destroying " << numEvents << " molecule(s) of type " << oneRxn.reactantMolList.at(0).molName
                    //           << " from reaction " << oneRxn.absRxnIndex << " at iteration " << simItr << std::endl;

                    while (numEvents > 0) {
                        if (counterArrays.copyNumSpecies[0] > 0) {
                            --counterArrays.copyNumSpecies[0];
                        }
                        if (membraneObject.numberOfFreeLipidsEachState[0] > 0) {
                            --membraneObject.numberOfFreeLipidsEachState[0];
                        }
                        --numEvents;
                        if (oneRxn.isObserved) {
                            auto observeItr = observablesList.find(oneRxn.observeLabel);
                            if (observeItr == observablesList.end()) {
                                // std::cerr << "WARNING: Observable " << oneRxn.observeLabel << " not defined.\n";
                            } else {
                                --observeItr->second;
                            }
                        }
                    }
                }
            }
        }
    }

    for (auto& oneRxn : createDestructRxns) {
        // check uniMolCreation, e.g. A -> A+B
        if (oneRxn.rxnType == ReactionType::uniMolCreation) {
            // Decide how many reactions occur based on lambda=k*dt*NA, where p(m_events)=exp(-lambda)*lamba^m /m!
            // determine the number of the template molecule A
            int NA { 0 };
            std::vector<int> poolAList {}; // this is all the A molecule index with TrajStatus::None and is the reaction's reactant
            for (auto& oneMol : moleculeList) {
                if (oneMol.trajStatus == TrajStatus::none && isReactant(oneMol, complexList[oneMol.myComIndex], oneRxn, moleculeList)) {
                    poolAList.push_back(oneMol.index);
                }
            }
            NA = static_cast<int>(poolAList.size());

            // determine the reaction number in this time step, numEvents
            long double lambda { oneRxn.rateList.at(0).rate * params.timeStep * Constants::usToSeconds * NA };
            double explam=exp(-lambda);
			long double prob { explam };
			unsigned numEvents { 0 };
			double rNum { 1.0 * rand_gsl() };
			long double lampow=lambda;
			long double factor = 1;
			
            while (rNum > prob) {
                ++numEvents;
				// prob += (exp(-lambda) * pow(lambda, numEvents)) / MathFuncs::gammFactorial(numEvents);
				/*this method is much faster, same accuracy*/
				prob += explam *lampow /factor;
				lampow = lampow*lambda;// accumulate the powers
				factor = factor * (numEvents + 1);//accumulate the factorial
					
            }
            //Here we are still limited in creation, because once one molecule creates, it cannot create again.
            //Allow one to create more than one in this case, since it would occur with a smaller time-step.

          //if(numEvents>NA){
	          //std::cout <<" WARNING: CANNOT perform N Creation Events of : "<<numEvents<<" Not enough species exist, perform: "<<NA<<std::endl;
	          //numEvents=NA;
	      //}

            // create numEvents B molecule in this time step, the cooresponding numEvents A are chosen randomly from the pool
            if (numEvents > 0) {
                // choose numEvents A from the poolAList randomly
                // std::cout << "Creating " << numEvents << " molecule(s) of type " << oneRxn.productMolList.back().molName
                //           << " from reaction " << oneRxn.absRxnIndex << " at iteration " << simItr << std::endl;
                while (numEvents > 0) {
                    // chose one A whose TrajStatus should be None

                    int randIntNum { static_cast<int>(1.0 * NA * rand_gsl()) };
                    if (randIntNum == NA) {
                        randIntNum = NA - 1;
                    }
                    if (moleculeList[poolAList[randIntNum]].trajStatus == TrajStatus::none) {
                        int newMolIndex { 0 };
                        int newComIndex { 0 };
                        create_molecule_and_complex_from_rxn(poolAList[randIntNum], newMolIndex, newComIndex, true,
                            molTemplateList[oneRxn.productMolList.back().molTypeIndex], params,
                            oneRxn, simulVolume, moleculeList, complexList, molTemplateList, forwardRxns, membraneObject);
                        if (numEvents <= NA) {
                            //need to allow at least one molecule to create more than one if numEvents>NA
                            //so do not update TrajStatus in that case.
                            moleculeList[poolAList[randIntNum]].trajStatus = TrajStatus::propagated;
                            complexList[moleculeList[poolAList[randIntNum]].myComIndex].trajStatus = TrajStatus::propagated;
                        }

                        // update the copy number array
                        for (const auto& iface : moleculeList[newMolIndex].interfaceList)
                            ++counterArrays.copyNumSpecies[iface.index];

                        --numEvents;

                        if (oneRxn.isObserved) {
                            auto observeItr = observablesList.find(oneRxn.observeLabel);
                            if (observeItr == observablesList.end()) {
                                // std::cerr << "WARNING: Observable " << oneRxn.observeLabel << " not defined.\n";
                            } else {
                                ++observeItr->second;
                            }
                        }
                    }
                }
            }
        }
    }
	bool breakLinkComplex;
	bool cancelDissociation;
	//std::vector<int> dissociateMolIndex {};
	//dissociateMolIndex.push_back(1);//allocate before clearing.

    // std::cout << "Start check dissociation on population..." << std::endl;
    for (auto& oneRxn : backRxns) {
        // In order to treat the dissociation based on population, we set the specie's canDissociate = true if it can dissociate
        // and track the specie's mol index in specie's bindPairList
        if (forwardRxns[oneRxn.conjForwardRxnIndex].rxnType == ReactionType::bimolecular) {
            bool isImplicit { false };
            if (molTemplateList[oneRxn.productListNew[0].molTypeIndex].isImplicitLipid || molTemplateList[oneRxn.productListNew[1].molTypeIndex].isImplicitLipid)
                isImplicit = true;

            if (isImplicit == false) {
                // check explicit dissociation, e.g. A-B bind pair -> A + B
                // Decide how many reactions occur based on lambda=k*dt*NAB, where p(m_events)=exp(-lambda)*lamba^m /m!
                // determine the number of the bind A-B
                int NAB { 0 };
                NAB = static_cast<int>(counterArrays.bindPairList[oneRxn.reactantListNew[0].absIfaceIndex].size());
                // if ( counterArrays.copyNumSpecies[oneRxn.reactantListNew[0].absIfaceIndex] != NAB ){
                //     std::cout<<"Wrong! in dissociation, copynum = "<<counterArrays.copyNumSpecies[oneRxn.reactantListNew[0].absIfaceIndex]<<", but NAB="<<NAB<<std::endl;
                //     exit(0);
                // }
                if (NAB > 0) { // determine the reaction rate, now only work for all offrates of one backRxn are the same; TODO: extend this to the general case is tricky, for we need to track the boundPairs including required ancillary interfaces
                    double rate { 0.0 };
                    rate = oneRxn.rateList.at(0).rate;
                    // std::cout << rate << std::endl;
                    for (auto oneRate : oneRxn.rateList) {
                        if (std::abs(oneRate.rate - rate) / std::abs(rate) > 1E-2) {
                          std::cerr
                              << "The current version only works for the case "
                                 "that dissociation reaction has one off rate "
                                 "value. i.e. the required ancillary "
                                 "interfaces does not influence the off rate."
                              << std::endl;
                          exit(1);
                        }
                    }

                    // determine the reaction number in this time step, numEvents
                    long double lambda { rate * params.timeStep * Constants::usToSeconds * NAB };
                    double explam = exp(-lambda);
                    long double prob { explam };
                    // std::cout << rate << "," << lambda << "," << prob << std::endl;
                    unsigned numEvents { 0 };
                    double rNum { 1.0 * rand_gsl() };
					long double lampow=lambda;
					long double factor = 1;
					
                    while (rNum > prob) {
                        ++numEvents;
						// prob += (exp(-lambda) * pow(lambda, numEvents)) / MathFuncs::gammFactorial(numEvents);
						/*this method is much faster, same accuracy*/
						prob += explam *lampow /factor;
						lampow = lampow*lambda;// accumulate the powers
						factor = factor * (numEvents + 1);//accumulate the factorial
					}

                    //if (params.debugParams.forceDissoc)
					//   numEvents = NAB;
                    if (numEvents > NAB) {
                        std::cout << " WARNING: CANNOT perform N Events of : " << numEvents << " Not enough species exist, perform: " << NAB << std::endl;
                        numEvents = NAB;
                    }

                    // choose numEvents A-B from the poolABList (stored in counterArrays.bindPairList) randomly, then break the interaction between the pair
                    if (numEvents > 0) {
                        // std::cout << "Before dissociation: " << std::endl
                        //           << oneRxn.productListNew[0].ifaceName << "-" << oneRxn.productListNew[1].ifaceName << " has " << static_cast<int>(counterArrays.bindPairList[oneRxn.reactantListNew[0].absIfaceIndex].size())
                        //           << "pairs. They are listed following..." << std::endl;
                        // for (auto one : counterArrays.bindPairList[oneRxn.reactantListNew[0].absIfaceIndex]) {
                        //     std::cout << one << "\t";
                        // }
                        // std::cout << std::endl;
					  
                        // choose numEvents A-B from the poolABList randomly
					  //std::cout << "Dissociating " << numEvents <<" Out of possible N: "<<NAB << " pair(s) " << oneRxn.productListNew[0].ifaceName << "-" << oneRxn.productListNew[1].ifaceName << " from reaction " << oneRxn.absRxnIndex << " at iteration " << simItr << std::endl;

					  // record the index of the selected mol
					  //dissociateMolIndex.clear();
					  std::vector<int> dissociateMolIndex {};
					  while (numEvents > 0) {
                            int randIntNum { static_cast<int>(1.0 * NAB * rand_gsl()) };
                            if (randIntNum == NAB) {
                                randIntNum = NAB - 1;
                            }

                            // int molIndexA { counterArrays.bindPairList[oneRxn.reactantListNew[0].absIfaceIndex][randIntNum] };
                            // figure out the iface rel index and mol index of B, molIndexA must be the first productant
                            // int ifaceIndexA { oneRxn.productListNew[0].relIfaceIndex };
                            // int molIndexB  { moleculeList[molIndexA].interfaceList[ifaceIndexA].interaction.partnerIndex };
                            // int ifaceIndexB { moleculeList[molIndexA].interfaceList[ifaceIndexA].interaction.partnerIfaceIndex };

                            // check whether this A-B has been selected in this step
                            // if (moleculeList[molIndexA].myComIndex != moleculeList[molIndexB].myComIndex) {
                            //     std::cerr << "ERROR: Molecules in different complexes are attempting to dissociate." << std::endl;
                            //     exit(1);
                            // }

                            std::vector<int>::iterator result { std::find(std::begin(dissociateMolIndex), std::end(dissociateMolIndex), counterArrays.bindPairList[oneRxn.reactantListNew[0].absIfaceIndex][randIntNum]) }; //check whether this A-B has been selected in this step
							if (result == std::end(dissociateMolIndex)) {
                                // dissociate
                                int molIndexA { counterArrays.bindPairList[oneRxn.reactantListNew[0].absIfaceIndex][randIntNum] };
                                // figure out the iface rel index and mol index of B, molIndexA must be the first productant
                                int ifaceIndexA { oneRxn.productListNew[0].relIfaceIndex };
                                int molIndexB { -1 };
                                int ifaceIndexB { -1 };
                                molIndexB = moleculeList[molIndexA].interfaceList[ifaceIndexA].interaction.partnerIndex;
                                ifaceIndexB = moleculeList[molIndexA].interfaceList[ifaceIndexA].interaction.partnerIfaceIndex;

                                if(moleculeList[molIndexA].isDissociated || moleculeList[molIndexB].isDissociated){
                                    continue;
                                }

								//std::cout << "Dissociation at iteration: " << simItr << " protein: " << molIndexA
                                //           << " partner: " << molIndexB << std::endl;
								//std::cout << "Complex " << moleculeList[molIndexA].myComIndex << ", composed of "
								//         << complexList[moleculeList[molIndexA].myComIndex].memberList.size() << " molecules" << std::endl;

                                if (moleculeList[molIndexA].myComIndex != moleculeList[molIndexB].myComIndex) {
                                    std::cerr << "ERROR: Molecules in different complexes are attempting to dissociate." << std::endl;
                                    exit(1);
                                }
								/*If a loop is broken but the complex remains closed, the rebinding of the loop is a unimolecular reaction, "evaluate_binding_within_complex". In case the forward rate and time-step, k*dt, is large (i.e. ~10^5 or 10^6), then there is a growing error in the actual vs theoretical binding, as poisson mean is k*dt, but we will sample only a mean of 1-exp(-k*dt) when we evaluate this bond one-at-a-time.
								 So we here will introduce a correction term, the purpose of the correction term is to ensure the correct equilibrium if the rebinding probability is too large, by lowering the unbinding probability. Theoretically tested, this can recover the exact eequilibrium although the kinetics will be shifted because both rates are a bit lower than expected.
								*/

                                // std::cout << "positions before dissociation: "
                                //           << moleculeList[molIndexA].comCoord
                                //           << "," << moleculeList[molIndexB].comCoord
                                //           << std::endl;
                                // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                                // This is commented out when merging
                                // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                                // bool breakLinkComplex
                                //     = break_interaction(simItr, ifaceIndexA, ifaceIndexB, moleculeList[molIndexA], moleculeList[molIndexB],
                                //         oneRxn, forwardRxns, moleculeList, complexList, molTemplateList, membraneObject.implicitlipidIndex);
                                // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                // std::cout << "positions AFTER dissociation: "
                                //           << moleculeList[molIndexA].comCoord
                                //           << ","
                                //           << moleculeList[molIndexB].comCoord
                                //           << std::endl;
								cancelDissociation
                                    = break_interaction(simItr, ifaceIndexA, ifaceIndexB, moleculeList[molIndexA], moleculeList[molIndexB],
														oneRxn, moleculeList, complexList, molTemplateList, membraneObject.implicitlipidIndex, 
                                                        forwardRxns[oneRxn.conjForwardRxnIndex], breakLinkComplex, params.timeStep, assocDissocFile);

								//if(cancelDissociation) std::cout <<" Canceled dissociation! "<<std::endl;
                                if(cancelDissociation == false){
                                    //std::cout <<" Returned to check_for_unimolecular_reactions_population after break_interactino "<<std::endl;
                                    //Continue on with dissociation event
                                    if (breakLinkComplex)
                                        counterArrays.nLoops--;

                                    // Change the number of bound pairs in the system
                                    update_Nboundpairs(moleculeList[molIndexA].molTypeIndex, moleculeList[molIndexB].molTypeIndex, -1, params, counterArrays);

                                    counterArrays.copyNumSpecies[oneRxn.reactantListNew[0].absIfaceIndex]--;
                                    counterArrays.copyNumSpecies[oneRxn.productListNew[0].absIfaceIndex]++;
                                    counterArrays.copyNumSpecies[oneRxn.productListNew[1].absIfaceIndex]++;

                                    //**Commented out since molecules undergoing dissociation do not move**//
                                    // consider the reflecting-surface movement
                                    // reflect_complex_rad_rot(membraneObject, complexList[moleculeList[molIndexA].myComIndex], moleculeList, 0.0);
                                    // reflect_complex_rad_rot(membraneObject, complexList[moleculeList[molIndexB].myComIndex], moleculeList, 0.0);

                                    // std::cout << "Coords of p1 (COM): " << moleculeList[molIndexA].comCoord << std::endl;
                                    // std::cout << "Coords of p2 (COM): " << moleculeList[molIndexB].comCoord << std::endl;

                                    for (auto memMol : complexList[moleculeList[molIndexA].myComIndex].memberList){
                                        moleculeList[memMol].trajStatus = TrajStatus::propagated;
                                    }
                                    for (auto memMol : complexList[moleculeList[molIndexB].myComIndex].memberList){
                                        moleculeList[memMol].trajStatus = TrajStatus::propagated;
                                    }
                                    complexList[moleculeList[molIndexA].myComIndex].trajStatus = TrajStatus::propagated;
                                    complexList[moleculeList[molIndexB].myComIndex].trajStatus = TrajStatus::propagated;
                                    moleculeList[molIndexA].isDissociated = true;
                                    moleculeList[molIndexB].isDissociated = true;

                                  if (oneRxn.isCoupled) {
                                    if (oneRxn.coupledRxn.probCoupled > rand_gsl()) {
									  // std::cout << "PERFORMING Coupled RXN after Dissociation! probability: " << oneRxn.coupledRxn.probCoupled << std::endl;
									  if (oneRxn.coupledRxn.rxnType == ReactionType::destruction) {
										int destroyProIndex { -1 };
										const CreateDestructRxn& coupledRxn
										  = createDestructRxns[oneRxn.coupledRxn.relRxnIndex]; // which reaction is being performed.
										if (moleculeList[molIndexA].molTypeIndex == coupledRxn.reactantMolList[0].molTypeIndex) {
										  destroyProIndex = molIndexA; // Is it A or B?
										} else {
										  destroyProIndex = molIndexB;
										}
										MolTemplate& oneTemp { molTemplateList[moleculeList[destroyProIndex].molTypeIndex] };

										// std::cout << "Performing coupled destruction reaction." << std::endl;
										// decrement the copy number array for everything in complex
										for (auto& memMol : complexList[moleculeList[destroyProIndex].myComIndex].memberList) {
										  for (auto& iface : moleculeList[memMol].interfaceList) {
											--counterArrays.copyNumSpecies[iface.index];
										  }
										}
										complexList[moleculeList[destroyProIndex].myComIndex].destroy(moleculeList, complexList); // destroying the entire complex that this molecule is a part of.

										// remove the molecule from the SimulVolume subsCellList, have this here to avoid circular header calls with SimulVolume and Molecule_Complex
										for (auto itr = simulVolume.subCellList[moleculeList[destroyProIndex].mySubVolIndex].memberMolList.begin();
											 itr != simulVolume.subCellList[moleculeList[destroyProIndex].mySubVolIndex].memberMolList.end();
											 ++itr) {
										  if (*itr == moleculeList[destroyProIndex].index) {
											simulVolume.subCellList[moleculeList[destroyProIndex].mySubVolIndex].memberMolList.erase(itr);
											break;
										  }
										}
										moleculeList[destroyProIndex].mySubVolIndex = -1; // reinitialize index

										// std::cout << "destroyProIndex: " << destroyProIndex;
										// std::cout << "Before coupled destruction, the monomerList: ";
										// for (auto one : oneTemp.monomerList) {
										//     std::cout << one << "\t";
										// }
										// std::cout << std::endl;
										// exit(1);

										oneTemp.monomerList.erase(std::find_if(oneTemp.monomerList.begin(), oneTemp.monomerList.end(), [&](const size_t& mol) { return mol == destroyProIndex; }));

										if (coupledRxn.isObserved) {
										  auto observeItr = observablesList.find(coupledRxn.observeLabel);
										  if (observeItr == observablesList.end()) {
											// std::cerr << "WARNING: Observable " << coupledRxn.observeLabel << " not defined.\n";
										  } else {
											--observeItr->second;
										  }
										}
									  }

									  if (oneRxn.coupledRxn.rxnType == ReactionType::uniMolStateChange) {
										int stateChangeProIndex { -1 };
										int relIndex { -1 };
										const ForwardRxn& coupledRxn = forwardRxns[oneRxn.coupledRxn.relRxnIndex]; // which reaction is being performed.
										// make sure the molecule iface have the same state with the uniMolStateChange reaction's reactant
										for (auto& tmpIface : moleculeList[molIndexA].interfaceList) {
										  if (tmpIface.index == coupledRxn.reactantListNew[0].absIfaceIndex) {
											stateChangeProIndex = molIndexA; // Is it A or B?
											relIndex = tmpIface.relIndex;
										  }
										}
										for (auto& tmpIface : moleculeList[molIndexB].interfaceList) {
										  //std::cout <<" mol: "<<molIndexB<<" ifaceIndex: "<<tmpIface.index<<"\t";
										  if (tmpIface.index == coupledRxn.reactantListNew[0].absIfaceIndex) {
											stateChangeProIndex = molIndexB; // Is it A or B?
											relIndex = tmpIface.relIndex;
										  }
										}

										if (stateChangeProIndex == -1) {
										  std::cerr << "The products of the disscociation do not match the corresponding uniMolStateChange reactant." << std::endl;
										  exit(1);
										}

										// std::cout << "Performing coupled uniMolStateChange reaction on protein: " << stateChangeProIndex << std::endl;

										const auto& stateList = molTemplateList[moleculeList[stateChangeProIndex].molTypeIndex].interfaceList[relIndex].stateList;
										const auto& newState = coupledRxn.productListNew[0];
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

										isObserved = coupledRxn.isObserved;
										observeLabel = coupledRxn.observeLabel;

										if (isObserved) {
										  auto observeItr = observablesList.find(observeLabel);
										  if (observeItr == observablesList.end()) {
											// std::cerr << "WARNING: Observable " << observeLabel << " not defined.\n";
										  } else {
											++observeItr->second;
										  }
										}

										--counterArrays.copyNumSpecies[coupledRxn.reactantListNew[0].absIfaceIndex];
										++counterArrays.copyNumSpecies[coupledRxn.productListNew[0].absIfaceIndex];

										// Change the state happens here
										moleculeList[stateChangeProIndex].interfaceList[relIndex].change_state(relStateIndex, newState.absIfaceIndex, newState.requiresState);

										moleculeList[stateChangeProIndex].trajStatus = TrajStatus::propagated;
										complexList[moleculeList[stateChangeProIndex].myComIndex].trajStatus = TrajStatus::propagated;
									  }
                                    } //isCoupled happened: prob>rnum
								  } // finished with IsCoupled?

								  dissociateMolIndex.emplace_back(counterArrays.bindPairList[oneRxn.reactantListNew[0].absIfaceIndex][randIntNum]);


								  if (oneRxn.isObserved) {
									auto observeItr = observablesList.find(oneRxn.observeLabel);
									if (observeItr == observablesList.end()) {
									  // std::cerr << "WARNING: Observable " << oneRxn.observeLabel << " not defined.\n";
                                    } else {
									  --observeItr->second;
                                    }
								  }
								}//IF CANCELLED DISSOCIATION, do not update anything
								--numEvents;//even cancelled dissociation counts as an event
							}
					  }
					  // remove the A-B from the bindPairList because they are seperated
					  for (auto oneDissociate : dissociateMolIndex) {
						counterArrays.bindPairList[oneRxn.reactantListNew[0].absIfaceIndex].erase(std::find_if(counterArrays.bindPairList[oneRxn.reactantListNew[0].absIfaceIndex].begin(), counterArrays.bindPairList[oneRxn.reactantListNew[0].absIfaceIndex].end(), [&](const size_t& mol) { return mol == oneDissociate; }));
					  }
					  // std::cout << "After dissociation: " << std::endl
					  //           << oneRxn.productListNew[0].ifaceName << "-" << oneRxn.productListNew[1].ifaceName << " has " << static_cast<int>(counterArrays.bindPairList[oneRxn.reactantListNew[0].absIfaceIndex].size())
					  //           << "pairs. They are listed following..." << std::endl;
					  // for (auto one : counterArrays.bindPairList[oneRxn.reactantListNew[0].absIfaceIndex]) {
					  //     std::cout << one << "\t";
					  // }
					  // std::cout << std::endl;
                    }//If we will perform any dissociation events
                }//If any bonds exist to break
            } // This is for the case when it is not implicit lipid. For implicit lipids it is using the unimolecular_reaction method (on at a time)
        }
    }

}
