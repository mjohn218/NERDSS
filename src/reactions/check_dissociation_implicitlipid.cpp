#include "boundary_conditions/reflect_functions.hpp"
#include "io/io.hpp"
#include "math/constants.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "tracing.hpp"

void check_dissociation_implicitlipid(unsigned int simItr, const Parameters& params, SimulVolume& simulVolume,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, unsigned int molItr,
    std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, const std::vector<BackRxn>& backRxns, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<CreateDestructRxn>& createDestructRxns, copyCounters& counterArrays, Membrane& membraneObject, std::vector<double>& IL2DbindingVec, std::vector<double>& IL2DUnbindingVec, std::vector<double>& ILTableIDs)
{
    // TRACE();
    double Dtot;
    for (int relIface1Itr { 0 }; relIface1Itr < moleculeList[molItr].bndlist.size(); ++relIface1Itr) {
        int relIface1 { moleculeList[molItr].bndlist[relIface1Itr] };
        if (moleculeList[molItr].interfaceList[relIface1].isBound) { // make sure it's actually bound
            int pro2Index { moleculeList[molItr].interfaceList[relIface1].interaction.partnerIndex };
            int relIface2 { moleculeList[molItr].interfaceList[relIface1].interaction.partnerIfaceIndex };

            //int com1Index { moleculeList[molItr].myComIndex };
            //int com2Index { moleculeList[pro2Index].myComIndex };
            // only consider when pro2 is implicit-lipid
            if (complexList[moleculeList[molItr].myComIndex].OnSurface == false || moleculeList[pro2Index].isImplicitLipid == false)
                continue;

            int mu { moleculeList[molItr].interfaceList[relIface1].interaction.conjBackRxn };
            // if the reaction is irreversible, the conjBackRxn index will be -1, so continue
            if (mu == -1)
                continue;

            int rateItr { find_reaction_rate_state(simItr, relIface1, relIface2, moleculeList[molItr],
                moleculeList[pro2Index], backRxns[mu], molTemplateList) };
            if (rateItr == -1)
                continue;

            double kb { backRxns[mu].rateList[rateItr].rate }; // <- kr[mu]
            int kfIndex = backRxns[mu].conjForwardRxnIndex; //!< the index of this reaction's ForwardRxn counterpart
            double ka = forwardRxns[kfIndex].rateList[rateItr].rate;
            double prob;

            /*Calculate Dtot for the binding reaction that gave rise to this complex. Assume Dcomplex=(1/D1+1/D2)-1 (two components)
	      D1=molTemplate[implicitLipidIndex].D, and Dcomplex is current D, solve for D2, and Dtot=D1+D2.
	     */
            //std::cout <<" In check implicit lipid Unbinding Complex D: "<<Dcomplex_x<<' '<<Dcomplex_y<<' '<<Dcomplex_z<<" And DIL:" <<DIL_x<<' '<<DIL_y<<' '<<DIL_z<<std::endl;
            //std::cout <<" In check implicit lipid Unbinding, calculated Dtot of protein bound to IL: "<<D2_x<<' '<<D2_y<<' '<<D2_z<<" And Dtot:" <<Dtot<<std::endl;
            // the expression of unbinding probability varies as implicit-lipid model (2D/3D) and explicit-lipid model
            // if the complex has only one protein-implicit.lipid bond, then the unbinding will make the complex into solution, 2D->3D,
            // then the reflect-surface RS3D needs to be introduced; otherwise, it is just 2D->2D, we donnot need use RS3D.
            if (complexList[moleculeList[molItr].myComIndex].linksToSurface == 1 && membraneObject.TwoD == false) {
                // 2D->3D
                double Dcomplex_x = complexList[moleculeList[molItr].myComIndex].D.x;
                double DIL_x = molTemplateList[moleculeList[pro2Index].molTypeIndex].D.x;
                double D2_x = 1.0 / (1.0 / Dcomplex_x - 1.0 / DIL_x);
                //If it is a protein bound to a lipid only, then DIL might be== Dcomplex due to forced rounding
                if (std::abs(Dcomplex_x - DIL_x) < 1E-10)
                    D2_x = molTemplateList[moleculeList[molItr].molTypeIndex].D.x;

                double Dcomplex_y = complexList[moleculeList[molItr].myComIndex].D.y;
                double DIL_y = molTemplateList[moleculeList[pro2Index].molTypeIndex].D.y;
                double D2_y = 1.0 / (1.0 / Dcomplex_y - 1.0 / DIL_y);
                if (std::abs(Dcomplex_y - DIL_y) < 1E-10)
                    D2_y = molTemplateList[moleculeList[molItr].molTypeIndex].D.y;
                //For z, both are currently on the membrane, so Dcomplex_z and DIL_z==0
                double Dcomplex_z = complexList[moleculeList[molItr].myComIndex].D.z;
                double DIL_z = molTemplateList[moleculeList[pro2Index].molTypeIndex].D.z;
                double D2_z = molTemplateList[moleculeList[molItr].molTypeIndex].D.z;

                Dtot = 1.0 / 3.0 * (D2_x + DIL_x) + 1.0 / 3.0 * (D2_y + DIL_y) + 1.0 / 3.0 * (D2_z + DIL_z);

                double sigma = forwardRxns[kfIndex].bindRadius;
                double h = params.timeStep;
                prob = dissociate3D(h, Dtot, sigma, 2.0 * ka, kb); //membraneObject.Punbinding3D;

                // make sure that the time step is resonable according to the prob of reaction
                if (prob > 1.000001) {
                    std::cerr << "Error: prob of reaction > 1. Avoid this using a smaller time step." << std::endl;
                    exit(1);
                }
                if (prob > 0.5) {
                    // std::cout << "WARNING: prob of reaction > 0.5. If this is a reaction for a bimolecular binding with multiple binding sites, please use a smaller time step." << std::endl;
                }
            } else {
                // 2D->2D
                //get the free number of the protein in the system
                int numberProtein { 0 };
                int proteinIndex { -1 };
                if (molTemplateList[forwardRxns[kfIndex].reactantListNew[0].molTypeIndex].isImplicitLipid == false) {
                    proteinIndex = forwardRxns[kfIndex].reactantListNew[0].absIfaceIndex;
                } else {
                    proteinIndex = forwardRxns[kfIndex].reactantListNew[1].absIfaceIndex;
                }
                //for (auto& tempMol : moleculeList) {
                //    for (auto& tempItf : tempMol.interfaceList) {
                //        if (tempItf.index == proteinIndex) {
                //            numberProtein++;
                //        }
                //    }
                //}
                numberProtein = counterArrays.copyNumSpecies[proteinIndex];

                long long int probMatrixIndex { 0 };
                double Dcomplex_x = complexList[moleculeList[molItr].myComIndex].D.x;
                double DIL_x = molTemplateList[moleculeList[pro2Index].molTypeIndex].D.x;
                double D2_x = 1.0 / (1.0 / Dcomplex_x - 1.0 / DIL_x);
                //If both are on membrane, D2 should be close to DIL
                if (std::abs(Dcomplex_x - DIL_x) < 1E-10)
                    D2_x = DIL_x; //molTemplateList[moleculeList[molItr].molTypeIndex].D.x;

                double Dcomplex_y = complexList[moleculeList[molItr].myComIndex].D.y;
                double DIL_y = molTemplateList[moleculeList[pro2Index].molTypeIndex].D.y;
                double D2_y = 1.0 / (1.0 / Dcomplex_y - 1.0 / DIL_y);
                if (std::abs(Dcomplex_y - DIL_y) < 1E-10)
                    D2_y = DIL_y; //molTemplateList[moleculeList[molItr].molTypeIndex].D.y;
                //For z, both are currently on the membrane, so Dcomplex_z and DIL_z==0
                double Dcomplex_z = complexList[moleculeList[molItr].myComIndex].D.z;
                double DIL_z = molTemplateList[moleculeList[pro2Index].molTypeIndex].D.z;
                double D2_z = 0;
                //If Links>1, then it is stuck in 2D even if this link breaks, so Dz=0.//otherwise, use the bound protein's Dz

                Dtot = 1.0 / 2.0 * (D2_x + DIL_x) + 1.0 / 2.0 * (D2_y + DIL_y);

                // declare intrinsic binding rate of 2D->2D case.
                double ktemp { ka / forwardRxns[kfIndex].length3Dto2D };

                paramsIL params2D {};
                params2D.R2D = 0.0;
                params2D.sigma = forwardRxns[kfIndex].bindRadius;
                params2D.Dtot = Dtot;
                params2D.ka = ktemp;
                int backIndex = forwardRxns[kfIndex].conjBackRxnIndex;
                params2D.kb = kb; //backRxns[backIndex].rateList[rateIndex].rate;
                params2D.area = membraneObject.totalSA;
                params2D.dt = params.timeStep;
                params2D.Na = numberProtein; //free protein

                const ForwardRxn& currRxn = forwardRxns[kfIndex];
                RxnIface implicitLipidState {};
                const auto& implicitLipidStateList = molTemplateList[moleculeList[membraneObject.implicitlipidIndex].molTypeIndex].interfaceList[0].stateList;
                if (molTemplateList[currRxn.reactantListNew[1].molTypeIndex].isImplicitLipid == true) {
                    implicitLipidState = currRxn.reactantListNew[1];
                } else {
                    implicitLipidState = currRxn.reactantListNew[0];
                }
                int relStateIndex { -1 };
                for (auto& state : implicitLipidStateList) {
                    if (state.index == implicitLipidState.absIfaceIndex) {
                        relStateIndex = static_cast<int>(&state - &implicitLipidStateList[0]);
                        break;
                    }
                }
                params2D.Nlipid = membraneObject.numberOfFreeLipidsEachState[relStateIndex]; // free lipids

                prob = dissociate2D(params2D); //2D-2D unbinding

                // make sure that the time step is resonable according to the prob of reaction
                if (prob > 1.000001) {
                    std::cerr << "Error: prob of reaction > 1. Avoid this using a smaller time step." << std::endl;
                    exit(1);
                }
                if (prob > 0.5) {
                    // std::cout << "WARNING: prob of reaction > 0.5. If this is a reaction for a bimolecular binding with multiple binding sites, please use a smaller time step." << std::endl;
                }
            }
            if (params.debugParams.forceDissoc)
                prob = 1.0;

            const ForwardRxn& currRxn = forwardRxns[kfIndex];
            double RS3D { -1.0 };
            for (int RS3Di = 0; RS3Di < 100; RS3Di++) {
                if ((std::abs(membraneObject.RS3Dvect[RS3Di] - currRxn.bindRadius) < 1E-15) && (std::abs(membraneObject.RS3Dvect[RS3Di + 100] - currRxn.rateList[0].rate) < 1E-15) && std::abs(membraneObject.RS3Dvect[RS3Di + 200] - (1.0 / 3.0 * (molTemplateList[currRxn.reactantListNew[0].molTypeIndex].D.x + molTemplateList[currRxn.reactantListNew[1].molTypeIndex].D.x) + 1.0 / 3.0 * (molTemplateList[currRxn.reactantListNew[0].molTypeIndex].D.y + molTemplateList[currRxn.reactantListNew[1].molTypeIndex].D.y) + 1.0 / 3.0 * (molTemplateList[currRxn.reactantListNew[0].molTypeIndex].D.z + molTemplateList[currRxn.reactantListNew[1].molTypeIndex].D.z))) < 1E-15) {
                    RS3D = membraneObject.RS3Dvect[RS3Di + 300];
                    break;
                }
            }

            double rnum { gsl_rng_uniform(r) };
            if (prob > rnum) {
                double rnum2 { rnum + gsl_rng_uniform(r) * 1.0 / (gsl_rng_max(r) + 1.0) }; // to get higher resolution
                if (prob > rnum2) {
                    // std::cout << "Dissociation at iteration: " << simItr << " protein: " << molItr
                    //           << " from the membrane, with probability: " << prob << '\n';
                    // std::cout << "Complex " << moleculeList[molItr].myComIndex << ", composed of "
                    //           << complexList[moleculeList[molItr].myComIndex].memberList.size() << " molecules\n";
                    // std::cout << " protein coords for : " << molItr << std::endl;
                    // moleculeList[molItr].display_my_coords("proteinonsurface");
                    /*std::cout << "Complex members:";
                        for (const auto& memMol : complexList[moleculeList[molItr].myComIndex].memberList)
                            std::cout << ' ' << memMol;
                        std::cout << '\n';
			*/

                    // 'break_interaction' frees certain protein and its interface
                    break_interaction_implicitlipid(relIface1, relIface2, moleculeList[molItr], moleculeList[pro2Index],
                        backRxns[mu], moleculeList, complexList, molTemplateList);

                    // Change the number of bound pairs in the system.
                    update_Nboundpairs(moleculeList[molItr].molTypeIndex, moleculeList[pro2Index].molTypeIndex, -1,
                        params, counterArrays);
                    //Update species copy numbers
                    //std::cout <<"Add new to index:" <<backRxns[mu].productListNew[0].absIfaceIndex<<' '<<backRxns[mu].productListNew[1].absIfaceIndex;
                    //std::cout <<"Subtract from index: "<<backRxns[mu].reactantListNew[0].absIfaceIndex<<std::endl;
                    counterArrays.copyNumSpecies[backRxns[mu].reactantListNew[0].absIfaceIndex] -= 1;
                    counterArrays.copyNumSpecies[backRxns[mu].productListNew[0].absIfaceIndex] += 1;
                    counterArrays.copyNumSpecies[backRxns[mu].productListNew[1].absIfaceIndex] += 1;
                    //counterArrays.copyNumSpecies[moleculeList[molItr].interfaceList[relIface1].index]  += 1;
                    //counterArrays.copyNumSpecies[moleculeList[pro2Index].interfaceList[relIface2].index]  += 1;

                    //update No_free_lipids according to state of IL
                    const BackRxn& currRxn = backRxns[mu];
                    RxnIface implicitLipidState {};
                    const auto& implicitLipidStateList = molTemplateList[moleculeList[membraneObject.implicitlipidIndex].molTypeIndex].interfaceList[0].stateList;
                    if (molTemplateList[currRxn.productListNew[1].molTypeIndex].isImplicitLipid == true) {
                        implicitLipidState = currRxn.productListNew[1];
                    } else {
                        implicitLipidState = currRxn.productListNew[0];
                    }
                    int relStateIndex { -1 };
                    for (auto& state : implicitLipidStateList) {
                        if (state.index == implicitLipidState.absIfaceIndex) {
                            relStateIndex = static_cast<int>(&state - &implicitLipidStateList[0]);
                            break;
                        }
                    }
                    membraneObject.numberOfFreeLipidsEachState[relStateIndex] += 1;

                    // update the number of bonds that this complex has connected to the membrane surface.
                    //this also needs to be done for the individual proteins.
                    complexList[moleculeList[molItr].myComIndex].linksToSurface -= 1;
                    moleculeList[molItr].linksToSurface--;
                    // consider the reflecting-surface movement
                    Vector transVec1 {};
                    if (complexList[moleculeList[molItr].myComIndex].linksToSurface < 1 && membraneObject.TwoD == false) {
                        if (membraneObject.isSphere) {
                            Coord coord = moleculeList[molItr].interfaceList[relIface1].coord;
                            double rtmp = membraneObject.sphereR - RS3D;
                            Coord coordnew = rtmp / coord.get_magnitude() * coord;
                            transVec1.x = coordnew.x - coord.x;
                            transVec1.y = coordnew.y - coord.y;
                            transVec1.z = coordnew.z - coord.z;
                        } else {
                            transVec1.x = 0;
                            transVec1.y = 0;
                            //transVec1.z = 0;
                            transVec1.z = (-membraneObject.waterBox.z / 2.0 + RS3D) - moleculeList[molItr].interfaceList[relIface1].coord.z; //move this interface to the RS3D position in z.
                            // std::cout << " unbinding POSITION UPDATE IN Z: " << transVec1.z << std::endl;
                        }
                        complexList[moleculeList[molItr].myComIndex].OnSurface = false; // the complex is not bound to the membrane anymore.
                    } else {
                        transVec1.x = 0;
                        transVec1.y = 0;
                        transVec1.z = 0;
                        complexList[moleculeList[molItr].myComIndex].OnSurface = true;
                    }
                    //BINDING CURRENT PUTS PROTEINS AT RS3D, SO NO DISPLACEMENT AFTER DISSOCIATIONS
                    // update the temporary coordinates. If binding places it at RS, do not move here!!
                    for (auto& mp : complexList[moleculeList[molItr].myComIndex].memberList)
                        moleculeList[mp].update_association_coords(transVec1);

                    for (auto memMol : complexList[moleculeList[molItr].myComIndex].memberList) {
                        moleculeList[memMol].comCoord = moleculeList[memMol].tmpComCoord;
                        for (unsigned int i { 0 }; i < moleculeList[memMol].interfaceList.size(); ++i)
                            moleculeList[memMol].interfaceList[i].coord = moleculeList[memMol].tmpICoords[i];
                        moleculeList[memMol].clear_tmp_association_coords();
                    }
                    //*/
                    complexList[moleculeList[molItr].myComIndex].update_properties(moleculeList, molTemplateList); // recalculate the properties of the first complex
                    reflect_complex_rad_rot(membraneObject, complexList[moleculeList[molItr].myComIndex], moleculeList, RS3D);

                    // std::cout << "Coords of p1 (COM) after dissociation: \n"; // << moleculeList[molItr].comCoord << '\n';
                    // moleculeList[molItr].display_my_coords("proteinatRS");
                    // std::cout << "The partner was an implicit lipid" << '\n';
                    // complexList[moleculeList[molItr].myComIndex].display();
                    // change the traj Status
                    //complexList[moleculeList[molItr].myComIndex].trajStatus = TrajStatus::propagated;
                    for (auto memMol : complexList[moleculeList[molItr].myComIndex].memberList)
                        moleculeList[memMol].trajStatus = TrajStatus::propagated;
                    complexList[moleculeList[molItr].myComIndex].trajStatus = TrajStatus::propagated;

                    // TODO: Temporary implementation for destruction coupled to dissociation
                    /*Must be a C->A+B
                          Find out if A or B is being destroyed.
                         */
                    if (backRxns[mu].isCoupled) {
                        if (backRxns[mu].coupledRxn.rxnType == ReactionType::destruction) {
                            int destroyProIndex { -1 };
                            const CreateDestructRxn& coupledRxn
                                = createDestructRxns[backRxns[mu].coupledRxn.relRxnIndex]; // which reaction is being
                            // performed.
                            if (moleculeList[molItr].molTypeIndex == coupledRxn.reactantMolList[0].molTypeIndex) {
                                destroyProIndex = molItr; // Is it A or B?
                            } else {
                                destroyProIndex = pro2Index;
                            }
                            if (moleculeList[destroyProIndex].isImplicitLipid == true) {
                                int indexIlState = coupledRxn.reactantMolList.at(0).interfaceList.at(0).absIfaceIndex;
                                // std::cout << "Performing coupled IL destruction reaction.\n";
                                --counterArrays.copyNumSpecies[indexIlState];
                                --membraneObject.numberOfFreeLipidsEachState[indexIlState];
                            } else {
                                // std::cout << "Performing coupled destruction reaction.\n";
                                // decrement the copy number array for everything in complex
                                for (auto& memMol : complexList[moleculeList[destroyProIndex].myComIndex].memberList) {
                                    for (auto& iface : moleculeList[memMol].interfaceList) {
                                        --counterArrays.copyNumSpecies[iface.index];
                                    }
                                }
                                complexList[moleculeList[destroyProIndex].myComIndex].destroy(moleculeList,
                                    complexList); // destroying the entire complex that this molecule is a part of.

                                // remove the molecule from the SimulVolume subsCellList
                                // have this here to avoid circular header calls with SimulVolume and
                                // Molecule_Complex
                                for (auto itr = simulVolume.subCellList[moleculeList[destroyProIndex].mySubVolIndex].memberMolList.begin();
                                     itr != simulVolume.subCellList[moleculeList[destroyProIndex].mySubVolIndex].memberMolList.end();
                                     ++itr) {
                                    if (*itr == moleculeList[destroyProIndex].index) {
                                        simulVolume.subCellList[moleculeList[destroyProIndex].mySubVolIndex].memberMolList.erase(itr);
                                        break;
                                    }
                                }
                                moleculeList[destroyProIndex].mySubVolIndex = -1; // reinitialize index

                                MolTemplate& oneTemp { molTemplateList[moleculeList[destroyProIndex].molTypeIndex] };
                                oneTemp.monomerList.erase(std::find_if(oneTemp.monomerList.begin(), oneTemp.monomerList.end(), [&](const size_t& mol) { return mol == destroyProIndex; }));
                            }
                            if (coupledRxn.isObserved) {
                                auto observeItr = observablesList.find(coupledRxn.observeLabel);
                                if (observeItr == observablesList.end()) {
                                    // std::cerr << "WARNING: Observable " << coupledRxn.observeLabel << " not defined.\n";
                                } else {
                                    --observeItr->second;
                                }
                            }
                        }

                        if (backRxns[mu].coupledRxn.rxnType == ReactionType::uniMolStateChange) {
                            int stateChangeProIndex { -1 };
                            int relIndex { -1 };
                            const ForwardRxn& coupledRxn = forwardRxns[backRxns[mu].coupledRxn.relRxnIndex]; // which reaction is being performed.
                            if (molTemplateList[coupledRxn.reactantListNew[0].molTypeIndex].isImplicitLipid == true) {
                                int indexIlState = coupledRxn.reactantListNew[0].absIfaceIndex;
                                int indexIlStateNew = coupledRxn.productListNew[0].absIfaceIndex;

                                // std::cout << "Performing coupled IL uniMolStateChange reaction.\n";

                                --counterArrays.copyNumSpecies[indexIlState];
                                ++counterArrays.copyNumSpecies[indexIlStateNew];
                                --membraneObject.numberOfFreeLipidsEachState[indexIlState];
                                ++membraneObject.numberOfFreeLipidsEachState[indexIlStateNew];

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
                            } else {
                                // make sure the molecule iface have the same state with the uniMolStateChange reaction's reactant
                                for (auto& tmpIface : moleculeList[molItr].interfaceList) {
                                    if (tmpIface.index == coupledRxn.reactantListNew[0].absIfaceIndex) {
                                        stateChangeProIndex = molItr; // Is it A or B?
                                        relIndex = tmpIface.relIndex;
                                    }
                                }
                                for (auto& tmpIface : moleculeList[pro2Index].interfaceList) {
                                    if (tmpIface.index == coupledRxn.reactantListNew[0].absIfaceIndex) {
                                        stateChangeProIndex = pro2Index; // Is it A or B?
                                        relIndex = tmpIface.relIndex;
                                    }
                                }

                                if (stateChangeProIndex == -1) {
                                    std::cerr << "The products of the disscociation do not match the corresponding uniMolStateChange reactant." << std::endl;
                                    exit(1);
                                }

                                // std::cout << "Performing coupled uniMolStateChange reaction.\n";

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
                                moleculeList[stateChangeProIndex].interfaceList[relIndex].change_state(
                                    relStateIndex, newState.absIfaceIndex, newState.requiresState);

                                moleculeList[stateChangeProIndex].trajStatus = TrajStatus::propagated;
                                complexList[moleculeList[stateChangeProIndex].myComIndex].trajStatus = TrajStatus::propagated;
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
        }
    }
}
