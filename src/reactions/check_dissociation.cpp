#include "boundary_conditions/reflect_functions.hpp"
#include "io/io.hpp"
#include "math/constants.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"

void check_dissociation(unsigned int simItr, const Parameters& params, SimulVolume& simulVolume,
    const std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, unsigned int molItr,
    std::vector<int>& emptyMolList, std::vector<int>& emptyComList, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, const std::vector<BackRxn>& backRxns, const std::vector<ForwardRxn>& forwardRxns,
			const std::vector<CreateDestructRxn>& createDestructRxns, copyCounters& counterArrays, const Membrane &membraneObject)
{
    for (int relIface1Itr { 0 }; relIface1Itr < moleculeList[molItr].bndlist.size(); ++relIface1Itr) {
        int relIface1 { moleculeList[molItr].bndlist[relIface1Itr] };
        if (moleculeList[molItr].interfaceList[relIface1].isBound) { // make sure it's actually bound
            int pro2Index { moleculeList[molItr].interfaceList[relIface1].interaction.partnerIndex };
            int relIface2 { moleculeList[molItr].interfaceList[relIface1].interaction.partnerIfaceIndex };

            if (molItr > pro2Index || moleculeList[pro2Index].isImplicitLipid)
                continue;

	    
            // Make sure the other protein knows it's bound to this protein
            bool candissociate = false;
            if (moleculeList[pro2Index].interfaceList[relIface2].interaction.partnerIndex == molItr
                && moleculeList[pro2Index].interfaceList[relIface2].interaction.partnerIfaceIndex == relIface1)
                    candissociate = true;
            
            if ( candissociate ) {
                int mu { moleculeList[molItr].interfaceList[relIface1].interaction.conjBackRxn };

                // if the reaction is irreversible, the conjBackRxn index will be -1, so continue
                if (mu == -1)
                    continue;

                int rateItr { find_reaction_rate_state(simItr, relIface1, relIface2, moleculeList[molItr],
                    moleculeList[pro2Index], backRxns[mu], molTemplateList) };
                if (rateItr == -1)
                    continue;

                /*This if statement is redundant to the continue if molItr>pro2Index above*/
                // if (molItr < pro2Index) {
                /*then this is the first protein in the reaction, so try it.
                This if statement ensures we do not try to dissociate the same complex twice
                */

                double kb { backRxns[mu].rateList[rateItr].rate }; // <- kr[mu]
                double prob = 1 - exp(-kb * params.timeStep * Constants::usToSeconds);
                double rnum { 1.0 * rand_gsl() };

                if (params.debugParams.forceDissoc)
                    prob = 1.0;

                if (prob > rnum) {
                    double rnum2 { rnum + rand_gsl() * Constants::iRandMax }; // to get higher resolution

                    if (prob > rnum2) {
                        std::cout << "Dissociation at iteration: " << simItr << " protein: " << molItr
                                  << " partner: " << pro2Index << '\n';
                        std::cout << "Complex " << moleculeList[molItr].myComIndex << ", composed of "
                                  << complexList[moleculeList[molItr].myComIndex].memberList.size() << " molecules\n";
                        /*Perform this dissociation reaction.
                        Sometimes it is a bond broken, not a full dissociation to two complexes if
                        the two interfaces are part of the same complex
                        */

                        /*std::cout << "Complex members:";
                        for (const auto& memMol : complexList[moleculeList[molItr].myComIndex].memberList)
                            std::cout << ' ' << memMol;
                        std::cout << '\n';
			*/
                        if (moleculeList[molItr].myComIndex != moleculeList[pro2Index].myComIndex) {
                            std::cerr << "ERROR: Molecules in different complexes are attempting to dissociate.\n";
                            exit(1);
                        }
                        bool breakLinkComplex
                            = break_interaction(relIface1, relIface2, moleculeList[molItr], moleculeList[pro2Index],
						backRxns[mu], emptyComList, moleculeList, complexList, molTemplateList, membraneObject.implicitlipidIndex);

                        if (breakLinkComplex)
                            counterArrays.nLoops--;
                        --relIface1Itr; // replaced this reaction, so stay on this one

                        /*Change the number of bound pairs in the system.*/
                        update_Nboundpairs(moleculeList[molItr].molTypeIndex, moleculeList[pro2Index].molTypeIndex, -1,
                            params, counterArrays);
                        /*Update species copy numbers*/
                        // decrement bound state
                        counterArrays.copyNumSpecies[backRxns[mu].reactantListNew[0].absIfaceIndex]--;
                        // increment free species
                        counterArrays.copyNumSpecies[backRxns[mu].productListNew[0].absIfaceIndex]++;
                        counterArrays.copyNumSpecies[backRxns[mu].productListNew[1].absIfaceIndex]++;

                        /*If dissociated products are removed from overlap lists, use ncross=-1.
                        If they remain in list to avoid overlap, use movestat=2 and also
                        ensure that they are not allowed to diffuse again, by, for example,
                        temporarily setting D=0.
                        */
                        // consider the reflecting-surface movement
                        reflect_complex_rad_rot(membraneObject, complexList[moleculeList[molItr].myComIndex], moleculeList);

                        std::cout << "Coords of p1 (COM): " << moleculeList[molItr].comCoord << '\n';
                        std::cout << "Coords of p2 (COM): " << moleculeList[pro2Index].comCoord << '\n';

                        for (auto memMol : complexList[moleculeList[molItr].myComIndex].memberList)
                            moleculeList[memMol].trajStatus = TrajStatus::propagated;
                        for (auto memMol : complexList[moleculeList[pro2Index].myComIndex].memberList)
                                moleculeList[memMol].trajStatus = TrajStatus::propagated;

                        // change complex traj Status
                        complexList[moleculeList[molItr].myComIndex].trajStatus = TrajStatus::propagated;
                        complexList[moleculeList[pro2Index].myComIndex].trajStatus = TrajStatus::propagated;

                        // TODO: Temporary implementation for destruction coupled to dissociation
                        /*Must be a C->A+B
                          Find out if A or B is being destroyed.
                         */
                        if (backRxns[mu].isCoupled) {
                            int destroyProIndex { -1 };
                            const CreateDestructRxn& coupledRxn
                                = createDestructRxns[backRxns[mu].coupledRxn.relRxnIndex]; // which reaction is being
                                                                                           // performed.
                            if (moleculeList[molItr].molTypeIndex == coupledRxn.reactantMolList[0].molTypeIndex) {
                                destroyProIndex = molItr; // Is it A or B?
                            } else{
                                destroyProIndex = pro2Index;
                            }

                            std::cout << "Performing coupled reaction.\n";
                            // decrement the copy number array for everything in complex
                            for (auto& memMol : complexList[moleculeList[destroyProIndex].myComIndex].memberList) {
                                for (auto& iface : moleculeList[memMol].interfaceList) {
                                    --counterArrays.copyNumSpecies[iface.index];
                                }
                            }
                            complexList[moleculeList[destroyProIndex].myComIndex].destroy(moleculeList, emptyMolList,
                                emptyComList); // destroying the entire complex that this molecule is a part of.

                            // remove the molecule from the SimulVolume subsCellList
                            // have this here to avoid circular header calls with SimulVolume and
                            // Molecule_Complex
                            for (auto itr = simulVolume.subCellList[moleculeList[destroyProIndex].mySubVolIndex].memberMolList.begin();
                                      itr!= simulVolume.subCellList[moleculeList[destroyProIndex].mySubVolIndex].memberMolList.end();
                                      ++itr) {
                                if (*itr == moleculeList[destroyProIndex].index) {
                                    simulVolume.subCellList[moleculeList[destroyProIndex].mySubVolIndex].memberMolList.erase(itr);
                                    break;
                                }
                            }
                            moleculeList[destroyProIndex].mySubVolIndex = -1; // reinitialize index

                            if (coupledRxn.isObserved) {
                                auto observeItr = observablesList.find(coupledRxn.observeLabel);
                                if (observeItr == observablesList.end()) {
                                    std::cerr << "WARNING: Observable " << coupledRxn.observeLabel << " not defined.\n";
                                } else {
                                    --observeItr->second;
                                }
                            }
                        } // finished with IsCoupled?

                        if (backRxns[mu].isObserved) {
                            auto obsItr = observablesList.find(backRxns[mu].observeLabel);
                            if (obsItr != observablesList.end())
                                --obsItr->second;
                        }
                    }
                }
                //} // only try each pair dissociating once

            } // finished trying out stateChange reactions for this proteins
        }
    }
}
