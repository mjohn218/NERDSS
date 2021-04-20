#include "classes/class_Rxns.hpp"
#include "reactions/association/association.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"

/*
  Loop over all interfaces of base1 and then base2.
  If they are within bindrad of each other (just set here to 1!!!) then set flagCancel=1, so the association will be
  Cancelled. base1 is not associating in this step, it is part of the system, base2 is performing association this step.
  For baseTmp, access its Tmp coords, as it is testing its new orientation!
 */
void measure_overlap_free_protein_interfaces(Molecule base1, Molecule baseTmp, bool& flagCancel,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns)
{
    // TRACE();
    int pro1MolType = base1.molTypeIndex;

    bool canInteract { false };
    for (auto partner : molTemplateList[pro1MolType].rxnPartners) {
        if (partner == baseTmp.molTypeIndex) {
            canInteract = true;
            break;
        }
    }

    // only consider when pro2 is NOT implicit-lipid
    if (canInteract) {
        if (baseTmp.isImplicitLipid)
            canInteract = false;
    }

    /*We now here proteins are in separate complexes. So not yet bound together. */

    if (canInteract) {
        /* CALCULATE ASSOCIATION PROBABILITIES */
        for (int relIface1Itr { 0 }; relIface1Itr < base1.freelist.size(); ++relIface1Itr) {
            /*test all of i1's binding partners to see whether they are on protein pro2 */
            int relIface1 { base1.freelist[relIface1Itr] };
            int absIface1 { base1.interfaceList[relIface1].index };
            int stateIndex1 { base1.interfaceList[relIface1].stateIndex };
            const Interface::State& state {
                molTemplateList[pro1MolType].interfaceList[relIface1].stateList[stateIndex1]
            };
            for (auto statePartner : state.rxnPartners) {
                for (int relIface2Idx = 0; relIface2Idx < baseTmp.freelist.size(); ++relIface2Idx) {
                    int relIface2 { baseTmp.freelist[relIface2Idx] };
                    int absIface2 { baseTmp.interfaceList[relIface2].index };
                    if (absIface2 == statePartner) { // both binding interfaces are available!
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
                            base1, baseTmp, forwardRxns, backRxns, molTemplateList);

                        if (rxnIndex != -1 && rateIndex != -1) {

                            /*distance between these two interfaces*/
                            double d2;
                            double bindrad2 = forwardRxns[rxnIndex].bindRadius * forwardRxns[rxnIndex].bindRadius;
                            // THIS SHOULD BE AN ACTUAL BINDING RADIUS squared BETWEEN INTERFACES WITHIN THIS PROTEIN
                            double dx = base1.interfaceList[relIface1].coord.x - baseTmp.tmpICoords[relIface2].x;
                            double dy = base1.interfaceList[relIface1].coord.y - baseTmp.tmpICoords[relIface2].y;
                            double dz = base1.interfaceList[relIface1].coord.z - baseTmp.tmpICoords[relIface2].z;

                            d2 = dx * dx + dy * dy + dz * dz;
                            if (d2 < bindrad2) {
                                // std::cout << " WARNING: Cancel Association, overlap with INTERFACES in SYSTEM ! " << relIface1 << ' ' << relIface2
                                //           << " separation : " << sqrt(d2) << " from proteins: " << '\t';
                                // cout <<" UNBOUND PROTEINSs "<<mp<<' '<<mp2 <<" Interfaces: "<<myIfaceIndex<<' '<<m<<  " WITHIN A
                                // COMPLEX ARE CLOSE TOGETHER: "<<d<<endl;
                                flagCancel = true;
                            }
                        }
                    }
                } //free interfaces on pro2 baseTmp
            }
        } //free interfaces on pro1 base1

    } //proteins are interactors.
}
