#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "classes/class_Membrane.hpp"
#include <algorithm>
// pro2Index is a lipid's index.

void check_implicit_reactions(int pro1Index, int pro2Index, int simItr, double* tableIDs, unsigned& DDTableIndex,
    const Parameters& params, std::vector<gsl_matrix*>& normMatrices, std::vector<gsl_matrix*>& survMatrices,
    std::vector<gsl_matrix*>& pirMatrices, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns,
			      const std::vector<BackRxn>& backRxns, copyCounters& counterArrays, Membrane &membraneObject, std::vector<double> &IL2DbindingVec, std::vector<double> &IL2DUnbindingVec, std::vector<double> &ILTableIDs)
{
	 // only consider when pro2 IS implicit-lipid
   if(moleculeList[pro2Index].isImplicitLipid == false)
    	 return; 
   
   int pro1MolType = moleculeList[pro1Index].molTypeIndex;
      
   bool canInteract { false };
   int index = moleculeList[pro2Index].interfaceList[0].index; // protein2 must be implicit-lipid, and has only one interface.
   int number_of_lipids = counterArrays.copyNumSpecies[index]; //membraneObject.nSites;
   if (number_of_lipids > 0)
       canInteract = true;    

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
            const Interface::State& state {molTemplateList[pro1MolType].interfaceList[relIface1].stateList[stateIndex1]};
            
            for (auto statePartner : state.rxnPartners) {
                for (int relIface2Idx = 0; relIface2Idx < moleculeList[pro2Index].freelist.size(); ++relIface2Idx) {
                    int relIface2 { moleculeList[pro2Index].freelist[relIface2Idx] };
                    int absIface2 { moleculeList[pro2Index].interfaceList[relIface2].index };
                    if (absIface2 == statePartner) { // both binding interfaces are available!
                        int rxnIndex { -1 };
                        int rateIndex { -1 };
                        bool isStateChangeBackRxn { false };
			
			find_which_reaction(relIface1, relIface2, rxnIndex, rateIndex, isStateChangeBackRxn, state,
					    moleculeList[pro1Index], moleculeList[pro2Index], forwardRxns, backRxns, molTemplateList);
						
                        if (rxnIndex != -1 && rateIndex != -1) {
                            int com1Index { moleculeList[pro1Index].myComIndex };
                            int com2Index { moleculeList[pro2Index].myComIndex };

                            {
                                Vector ifaceVec { moleculeList[pro1Index].interfaceList[relIface1].coord
                                    - complexList[com1Index].comCoord };
                                Vector ifaceVec2 { moleculeList[pro2Index].interfaceList[relIface2].coord
                                    - complexList[com2Index].comCoord };
                                double magMol1 { ifaceVec.x * ifaceVec.x + ifaceVec.y * ifaceVec.y
                                    + ifaceVec.z * ifaceVec.z };
                                double magMol2 { ifaceVec2.x * ifaceVec2.x + ifaceVec2.y * ifaceVec2.y
                                    + ifaceVec2.z * ifaceVec2.z };

                                double Dtot = 1.0 / 3.0 * (complexList[com1Index].D.x + complexList[com2Index].D.x)
                                    + 1.0 / 3.0 * (complexList[com1Index].D.y + complexList[com2Index].D.y)
                                    + 1.0 / 3.0 * (complexList[com1Index].D.z + complexList[com2Index].D.z);

                                BiMolData biMolData { pro1Index, pro2Index, com1Index, com2Index, relIface1, relIface2,
                                    absIface1, absIface2, Dtot, magMol1, magMol2 };
                                    
				// binding with implicit-lipid model
				if (complexList[com1Index].D.z == 0 ){
				    // both Complexes are on the membrane, evaluate as 2D reaction
				    determine_2D_implicitlipid_reaction_probability(simItr, rxnIndex, rateIndex, isStateChangeBackRxn, 
										    ILTableIDs, biMolData, params, moleculeList, complexList, forwardRxns, backRxns, 
										    IL2DbindingVec, IL2DUnbindingVec, membraneObject);
				} else {
				    //3D->2D reaction
				    determine_3D_implicitlipid_reaction_probability(simItr, rxnIndex, rateIndex, isStateChangeBackRxn,
										    biMolData, params, moleculeList, complexList, forwardRxns, backRxns, 
										    membraneObject);
				}
                            }
                        }
                    }
                }
            }
        } // These protein partners interact
    }
}
