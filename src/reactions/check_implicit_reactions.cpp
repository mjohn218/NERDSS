#include "classes/class_Membrane.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"
#include <algorithm>

// pro2Index is a lipid's index.

void check_implicit_reactions(int pro1Index, int pro2Index, int simItr,
    const Parameters& params, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, copyCounters& counterArrays, Membrane& membraneObject, std::vector<double>& IL2DbindingVec, std::vector<double>& IL2DUnbindingVec, std::vector<double>& ILTableIDs)
{
    // TRACE();
    // only consider when pro2 IS implicit-lipid
    if (moleculeList[pro2Index].isImplicitLipid == false)
        return;

    int pro1MolType = moleculeList[pro1Index].molTypeIndex;

    bool canInteract { true };
    int index = moleculeList[pro2Index].interfaceList[0].index; // protein2 must be implicit-lipid, and has only one interface.
    int number_of_lipids = 0; //sum of all states of IL
    for (int i = 0; i < membraneObject.numberOfFreeLipidsEachState.size(); i++) {
        number_of_lipids += membraneObject.numberOfFreeLipidsEachState[i];
    }
    if (number_of_lipids <= 0) {
        canInteract = false;
    }

    //allow multiple links to the surface.
    /*      if (canInteract) {
        // If this protein is already bound to surface, don't test for binding OR overlap, set canInteract= false/
        if (moleculeList[pro1Index].bndlist.size() > 0 ) {
            bool boundPro1 { false };
            for (auto& partner : moleculeList[pro1Index].bndpartner) {
                if (partner == pro2Index) {
                    boundPro1 = true;
                    break;
                }
            }
            if (boundPro1)
                canInteract = false;
        }
   }
   */
    if (canInteract) {
        /* CALCULATE ASSOCIATION PROBABILITIES */
        for (int relIface1Itr { 0 }; relIface1Itr < moleculeList[pro1Index].freelist.size(); ++relIface1Itr) {
            //test all of i1's binding partners to see whether they are on protein pro2
            int relIface1 { moleculeList[pro1Index].freelist[relIface1Itr] };
            int absIface1 { moleculeList[pro1Index].interfaceList[relIface1].index };
            int stateIndex1 { moleculeList[pro1Index].interfaceList[relIface1].stateIndex };
            const Interface::State& state { molTemplateList[pro1MolType].interfaceList[relIface1].stateList[stateIndex1] };

            for (auto statePartner : state.rxnPartners) {
                for (int relIface2Idx = 0; relIface2Idx < moleculeList[pro2Index].freelist.size(); ++relIface2Idx) {
                    int relIface2 { moleculeList[pro2Index].freelist[relIface2Idx] };
                    int absIface2 { moleculeList[pro2Index].interfaceList[relIface2].index };
                    for (int tmpImplicitLipidStateIndex = 0; tmpImplicitLipidStateIndex < membraneObject.nStates; tmpImplicitLipidStateIndex++) {
                        if ((absIface2 == statePartner - tmpImplicitLipidStateIndex) && membraneObject.numberOfFreeLipidsEachState[tmpImplicitLipidStateIndex] > 0) { // both binding interfaces are available!
                            int rxnIndex { -1 };
                            int rateIndex { -1 };
                            bool isStateChangeBackRxn { false };

                            relIface2 += tmpImplicitLipidStateIndex;
                            absIface2 += tmpImplicitLipidStateIndex;

                            find_which_reaction(relIface1, relIface2, rxnIndex, rateIndex, isStateChangeBackRxn, state,
                                moleculeList[pro1Index], moleculeList[pro2Index], forwardRxns, backRxns, molTemplateList);

                            relIface2 -= tmpImplicitLipidStateIndex;
                            absIface2 -= tmpImplicitLipidStateIndex;

                            if (rxnIndex != -1 && rateIndex != -1) {
                                //int com1Index { moleculeList[pro1Index].myComIndex };
                                //int com2Index { moleculeList[pro2Index].myComIndex };

                                {
                                    Vector ifaceVec { moleculeList[pro1Index].interfaceList[relIface1].coord
                                        - complexList[moleculeList[pro1Index].myComIndex].comCoord };
                                    Vector ifaceVec2 { moleculeList[pro2Index].interfaceList[relIface2].coord
                                        - complexList[moleculeList[pro2Index].myComIndex].comCoord };
                                    double magMol1 { ifaceVec.x * ifaceVec.x + ifaceVec.y * ifaceVec.y
                                        + ifaceVec.z * ifaceVec.z };
                                    double magMol2 { ifaceVec2.x * ifaceVec2.x + ifaceVec2.y * ifaceVec2.y
                                        + ifaceVec2.z * ifaceVec2.z };

                                    //find the state of IL
                                    int forwardRxnIndex { rxnIndex };
                                    const ForwardRxn& currRxn = forwardRxns[forwardRxnIndex];
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

                                    // binding with implicit-lipid model
                                    if (std::abs(complexList[moleculeList[pro1Index].myComIndex].D.z) < 1E-10) {
                                        // both Complexes are on the membrane, evaluate as 2D reaction
                                        double Dtot = 1.0 / 2.0 * (complexList[moleculeList[pro1Index].myComIndex].D.x + complexList[moleculeList[pro2Index].myComIndex].D.x)
                                            + 1.0 / 2.0 * (complexList[moleculeList[pro1Index].myComIndex].D.y + complexList[moleculeList[pro2Index].myComIndex].D.y);

                                        BiMolData biMolData { pro1Index, pro2Index, moleculeList[pro1Index].myComIndex, moleculeList[pro2Index].myComIndex, relIface1, relIface2,
                                            absIface1, absIface2, Dtot, magMol1, magMol2 };
                                        determine_2D_implicitlipid_reaction_probability(simItr, rxnIndex, rateIndex, isStateChangeBackRxn,
                                            ILTableIDs, biMolData, params, moleculeList, complexList, forwardRxns, backRxns,
                                            IL2DbindingVec, IL2DUnbindingVec, membraneObject, relStateIndex);
                                    } else {
                                        //3D->2D reaction
                                        double Dtot = 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.x + complexList[moleculeList[pro2Index].myComIndex].D.x)
                                            + 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.y + complexList[moleculeList[pro2Index].myComIndex].D.y)
                                            + 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.z + complexList[moleculeList[pro2Index].myComIndex].D.z);

                                        BiMolData biMolData { pro1Index, pro2Index, moleculeList[pro1Index].myComIndex, moleculeList[pro2Index].myComIndex, relIface1, relIface2,
                                            absIface1, absIface2, Dtot, magMol1, magMol2 };

                                        determine_3D_implicitlipid_reaction_probability(simItr, rxnIndex, rateIndex, isStateChangeBackRxn,
                                            biMolData, params, moleculeList, complexList, forwardRxns, backRxns,
                                            membraneObject, relStateIndex);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } // These protein partners interact
    }
}
