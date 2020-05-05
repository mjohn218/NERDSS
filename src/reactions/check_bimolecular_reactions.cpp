#include "math/rand_gsl.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"
#include <algorithm>

void check_bimolecular_reactions(int pro1Index, int pro2Index, int simItr, double* tableIDs, unsigned& DDTableIndex,
    const Parameters& params, std::vector<gsl_matrix*>& normMatrices, std::vector<gsl_matrix*>& survMatrices,
    std::vector<gsl_matrix*>& pirMatrices, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, copyCounters& counterArrays, Membrane& membraneObject)
{
    TRACE();
    int pro1MolType = moleculeList[pro1Index].molTypeIndex;
    /* if(pro1Index== 21)
 	std::cout <<"In check bimolecular Reaction for protein 21 ! to  "<<pro2Index<<std::endl;
     if(pro1Index== 27)
 	std::cout <<"In check bimolecular Reaction for protein 27 ! to "<<pro2Index<<std::endl;
    */
    bool canInteract { false };
    for (auto partner : molTemplateList[pro1MolType].rxnPartners) {
        if (partner == moleculeList[pro2Index].molTypeIndex) {
            canInteract = true;
            break;
        }
    }

    // only consider when pro2 is NOT implicit-lipid
    if (canInteract) {
        if (moleculeList[pro2Index].isImplicitLipid)
            canInteract = false;
    }

    // If this pair of proteins are already bound together, don't test for binding OR overlap
    if (canInteract) {
        /*If this pair of proteins are already bound together, don't test for binding OR overlap, set canInteract=0*/
        if (moleculeList[pro1Index].bndlist.size() > 0 && moleculeList[pro2Index].bndlist.size() > 0) {
            bool boundPro1 { false };
            bool boundPro2 { false };
            for (auto& partner : moleculeList[pro1Index].bndpartner) {
                if (partner == pro2Index) {
                    boundPro1 = true;
                    break;
                }
            }
            for (auto& partner : moleculeList[pro2Index].bndpartner) {
                if (partner == pro1Index) {
                    boundPro2 = true;
                    break;
                }
            }
            // TODO: add other case
            if (boundPro1 && boundPro2)
                canInteract = false;
        }
    }

    if (canInteract) {
        /* CALCULATE ASSOCIATION PROBABILITIES */
        /*if(pro1Index== 21)
 	    std::cout <<" calculate Association prob to ! "<<pro2Index<<std::endl;
 	if(pro1Index== 27)
 	    std::cout <<" calculate Association prob to ! "<<pro2Index<<std::endl;
	*/
        for (int relIface1Itr { 0 }; relIface1Itr < moleculeList[pro1Index].freelist.size(); ++relIface1Itr) {
            /*test all of i1's binding partners to see whether they are on protein pro2 */
            int relIface1 { moleculeList[pro1Index].freelist[relIface1Itr] };
            int absIface1 { moleculeList[pro1Index].interfaceList[relIface1].index };
            int stateIndex1 { moleculeList[pro1Index].interfaceList[relIface1].stateIndex };
            const Interface::State& state {
                molTemplateList[pro1MolType].interfaceList[relIface1].stateList[stateIndex1]
            };
            //if(pro1Index == 21 || pro1Index ==27)
            //std::cout <<"reliface, abs, stateindex, npartners: "<<relIface1<<' '<<absIface1<<' '<<stateIndex1<<' '<<state.rxnPartners.size()<<std::endl;
            for (auto statePartner : state.rxnPartners) {
                //if(pro1Index == 21 || pro1Index ==27)
                //  std::cout <<"statePertner: "<<statePartner<<std::endl;
                for (int relIface2Idx = 0; relIface2Idx < moleculeList[pro2Index].freelist.size(); ++relIface2Idx) {
                    int relIface2 { moleculeList[pro2Index].freelist[relIface2Idx] };
                    int absIface2 { moleculeList[pro2Index].interfaceList[relIface2].index };

                    if (absIface2 == statePartner) { // both binding interfaces are available!
                        /*if(pro1Index== 21)
 			    std::cout <<" Both Ifaces available! "<<pro2Index<<std::endl;
 			if(pro1Index== 27)
 			    std::cout <<" Both Ifaces available!! "<<pro2Index<<std::endl;
			*/
                        /*Here now we evaluate the probability of binding*/
                        /*Different interfaces can have different diffusion constants if they
                        also rotate. <theta^2>=6Drparams.timeStep.
                        In that case, <d>=sin(sqrt(6Drparams.timeStep)/2)*2R
                        so add in <d>^2=4R^2sin^2(sqrt(6Drparams.timeStep)/2)
                        for a single clathrin, R=arm length,
                        otherwise R will be from pivot point of rotation, COM,
                        calculate distance from interface to the complex COM
                        */
                        int rxnIndex { -1 };
                        int rateIndex { -1 };
                        bool isStateChangeBackRxn { false };

                        find_which_reaction(relIface1, relIface2, rxnIndex, rateIndex, isStateChangeBackRxn, state,
                            moleculeList[pro1Index], moleculeList[pro2Index], forwardRxns, backRxns, molTemplateList);
                        /*if(pro1Index== 21)
 			    std::cout <<" found which reaction! "<<rxnIndex<<" "<<rateIndex<<std::endl;
 			if(pro1Index== 27)
 			    std::cout <<" found which reaction! "<<rxnIndex<<" "<<rateIndex<<std::endl;
			*/
                        if (rxnIndex != -1 && rateIndex != -1) {
                            //int com1Index { moleculeList[pro1Index].myComIndex };
                            //int com2Index { moleculeList[pro2Index].myComIndex };

                            if (moleculeList[pro1Index].myComIndex == moleculeList[pro2Index].myComIndex) {
                                evaluate_binding_within_complex(pro1Index, pro2Index, relIface1, relIface2, rxnIndex,
                                    rateIndex, isStateChangeBackRxn, params, moleculeList,
                                    complexList, molTemplateList,
                                    forwardRxns[rxnIndex], backRxns, membraneObject, counterArrays);
                            } else {
                                Vector ifaceVec { moleculeList[pro1Index].interfaceList[relIface1].coord
                                    - complexList[moleculeList[pro1Index].myComIndex].comCoord };
                                Vector ifaceVec2 { moleculeList[pro2Index].interfaceList[relIface2].coord
                                    - complexList[moleculeList[pro2Index].myComIndex].comCoord };
                                double magMol1 { ifaceVec.x * ifaceVec.x + ifaceVec.y * ifaceVec.y
                                    + ifaceVec.z * ifaceVec.z };
                                double magMol2 { ifaceVec2.x * ifaceVec2.x + ifaceVec2.y * ifaceVec2.y
                                    + ifaceVec2.z * ifaceVec2.z };

                                // binding with explicit-lipid model.
                                //write_rng_state();
                                if (std::abs(complexList[moleculeList[pro1Index].myComIndex].D.z) < 1E-10 && std::abs(complexList[moleculeList[pro2Index].myComIndex].D.z) < 1E-10) {
                                    // both Complexes are on the membrane, evaluate as 2D reaction
                                    double Dtot = 1.0 / 2.0 * (complexList[moleculeList[pro1Index].myComIndex].D.x + complexList[moleculeList[pro2Index].myComIndex].D.x)
                                        + 1.0 / 2.0 * (complexList[moleculeList[pro1Index].myComIndex].D.y + complexList[moleculeList[pro2Index].myComIndex].D.y);

                                    BiMolData biMolData { pro1Index, pro2Index, moleculeList[pro1Index].myComIndex, moleculeList[pro2Index].myComIndex, relIface1, relIface2,
                                        absIface1, absIface2, Dtot, magMol1, magMol2 };

                                    determine_2D_bimolecular_reaction_probability(simItr, rxnIndex, rateIndex,
                                        isStateChangeBackRxn, DDTableIndex, tableIDs, biMolData, params, moleculeList,
                                        complexList, forwardRxns, backRxns, membraneObject, normMatrices, survMatrices, pirMatrices);
                                } else {
                                    //3D reaction
                                    double Dtot = 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.x + complexList[moleculeList[pro2Index].myComIndex].D.x)
                                        + 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.y + complexList[moleculeList[pro2Index].myComIndex].D.y)
                                        + 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.z + complexList[moleculeList[pro2Index].myComIndex].D.z);

                                    BiMolData biMolData { pro1Index, pro2Index, moleculeList[pro1Index].myComIndex, moleculeList[pro2Index].myComIndex, relIface1, relIface2,
                                        absIface1, absIface2, Dtot, magMol1, magMol2 };

                                    determine_3D_bimolecular_reaction_probability(simItr, rxnIndex, rateIndex,
                                        isStateChangeBackRxn, DDTableIndex, tableIDs, biMolData, params, moleculeList,
                                        complexList, forwardRxns, backRxns, membraneObject, normMatrices, survMatrices, pirMatrices);
                                } //end else 3D
                                //read_rng_state();
                            }
                        }
                    }
                }
            }
        } // These protein partners interact
    }
}
