#include "math/rand_gsl.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"
#include <algorithm>
#include <vector>

void check_bimolecular_reactions(int pro1Index, int pro2Index, int simItr, double* tableIDs, unsigned& DDTableIndex,
    const Parameters& params, std::vector<gsl_matrix*>& normMatrices, std::vector<gsl_matrix*>& survMatrices,
    std::vector<gsl_matrix*>& pirMatrices, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, copyCounters& counterArrays, Membrane& membraneObject)
{
    // TRACE();
    //  int track1 = 143;
    //int track2= 182;
    int pro1MolType = moleculeList[pro1Index].molTypeIndex;
    /* if(pro1Index== track1)
       std::cout <<"In check bimolecular Reaction for protein track1 !"<<track1<<" to  "<<pro2Index<<"mytype: "<<pro1MolType<<" pro2type: "<<moleculeList[pro2Index].molTypeIndex<<std::endl;
     if(pro1Index== track2)
       std::cout <<"In check bimolecular Reaction for protein track2 !" <<track2<< " to "<<pro2Index<<" mytype: "<<pro1MolType<<" pro2type: "<<moleculeList[pro2Index].molTypeIndex<<std::endl;
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

    bool canExclude { false };
    if ((moleculeList[pro1Index].bndlist.size() > 0 && molTemplateList[moleculeList[pro1Index].molTypeIndex].excludeVolumeBound == true)
        || (moleculeList[pro2Index].bndlist.size() > 0 && molTemplateList[moleculeList[pro2Index].molTypeIndex].excludeVolumeBound == true)) {
        canExclude = true;
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
        /* if(pro1Index== track1)
 	    std::cout <<" calculate Association prob to ! "<<pro2Index<<std::endl;
 	if(pro1Index== track2)
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
            //if(pro1Index == track1 || pro1Index ==track2)
            //std::cout <<"reliface, abs, stateindex, npartners: "<<relIface1<<' '<<absIface1<<' '<<stateIndex1<<' '<<state.rxnPartners.size()<<std::endl;
            for (auto statePartner : state.rxnPartners) {
                //if(pro1Index == track1 || pro1Index ==track2)
                //  std::cout <<"statePertner: "<<statePartner<<std::endl;
                for (int relIface2Idx = 0; relIface2Idx < moleculeList[pro2Index].freelist.size(); ++relIface2Idx) {
                    int relIface2 { moleculeList[pro2Index].freelist[relIface2Idx] };
                    int absIface2 { moleculeList[pro2Index].interfaceList[relIface2].index };

                    if (absIface2 == statePartner) { // both binding interfaces are available!
                        /*if(pro1Index== track1)
			std::cout <<" Both Ifaces available! "<<pro2Index<<" abs1: "<<absIface1<<" abs2: "<<absIface2<<std::endl;
		      if(pro1Index== track2)
			std::cout <<" Both Ifaces available!! "<<pro2Index<<" abs1: "<<absIface1<<" abs2: "<<absIface2<<std::endl;
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
                                if (std::abs(complexList[moleculeList[pro1Index].myComIndex].D.z) < 1E-16 && std::abs(complexList[moleculeList[pro2Index].myComIndex].D.z) < 1E-16) {
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
        }
    } // These protein partners interact

    if (canExclude == true) {
        // pro1 is bound
        if (moleculeList[pro1Index].bndlist.size() > 0 && molTemplateList[moleculeList[pro1Index].molTypeIndex].excludeVolumeBound == true) {
            // loop the bindlist to find the exclude partner of each interface
            for (int relIface1Itr { 0 }; relIface1Itr < moleculeList[pro1Index].bndlist.size(); ++relIface1Itr) {
                int relIface1 { moleculeList[pro1Index].bndlist[relIface1Itr] };
                int absIface1 { moleculeList[pro1Index].interfaceList[relIface1].index };
                int molTypeIndex1 { moleculeList[pro1Index].molTypeIndex };
                if (moleculeList[pro1Index].interfaceList[relIface1].isBound && molTemplateList[moleculeList[pro1Index].molTypeIndex].interfaceList[relIface1].excludeVolumeBoundList.empty() == false) { // make sure it's actually bound and need excludeVolumeBound
                    // loop the interfaceList of pro2
                    for (int relIface2Idx = 0; relIface2Idx < moleculeList[pro2Index].interfaceList.size(); ++relIface2Idx) {
                        int relIface2 { relIface2Idx };
                        int absIface2 { moleculeList[pro2Index].interfaceList[relIface2].index };
                        int molTypeIndex2 { moleculeList[pro2Index].molTypeIndex };
                        for (int indexItr { 0 }; indexItr < molTemplateList[moleculeList[pro1Index].molTypeIndex].interfaceList[relIface1].excludeVolumeBoundList.size(); ++indexItr) {
                            if (molTemplateList[moleculeList[pro1Index].molTypeIndex].interfaceList[relIface1].excludeVolumeBoundList[indexItr] == molTypeIndex2) {
                                double bindRadius { molTemplateList[moleculeList[pro1Index].molTypeIndex].interfaceList[relIface1].excludeRadiusList[indexItr] };
                                int rxnIndex { molTemplateList[moleculeList[pro1Index].molTypeIndex].interfaceList[relIface1].excludeVolumeBoundReactList[indexItr] };
                                if (molTemplateList[moleculeList[pro1Index].molTypeIndex].interfaceList[relIface1].excludeVolumeBoundIfaceList[indexItr] == relIface2) {
                                    // calculate the distance between the two infaces
                                    Vector ifaceVec { moleculeList[pro1Index].interfaceList[relIface1].coord
                                        - complexList[moleculeList[pro1Index].myComIndex].comCoord };
                                    Vector ifaceVec2 { moleculeList[pro2Index].interfaceList[relIface2].coord
                                        - complexList[moleculeList[pro2Index].myComIndex].comCoord };
                                    double magMol1 { ifaceVec.x * ifaceVec.x + ifaceVec.y * ifaceVec.y
                                        + ifaceVec.z * ifaceVec.z };
                                    double magMol2 { ifaceVec2.x * ifaceVec2.x + ifaceVec2.y * ifaceVec2.y
                                        + ifaceVec2.z * ifaceVec2.z };
                                    if (std::abs(complexList[moleculeList[pro1Index].myComIndex].D.z) < 1E-16 && std::abs(complexList[moleculeList[pro2Index].myComIndex].D.z) < 1E-16) {
                                        // both Complexes are on the membrane, evaluate as 2D reaction
                                        double Dtot = 1.0 / 2.0 * (complexList[moleculeList[pro1Index].myComIndex].D.x + complexList[moleculeList[pro2Index].myComIndex].D.x)
                                            + 1.0 / 2.0 * (complexList[moleculeList[pro1Index].myComIndex].D.y + complexList[moleculeList[pro2Index].myComIndex].D.y);

                                        BiMolData biMolData { pro1Index, pro2Index, moleculeList[pro1Index].myComIndex, moleculeList[pro2Index].myComIndex, relIface1, relIface2,
                                            absIface1, absIface2, Dtot, magMol1, magMol2 };

                                        double Dr1 {};
                                        {
                                            double cf { cos(sqrt(2.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep)) };
                                            Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
                                        }

                                        double Dr2 {};
                                        {
                                            double cf = cos(sqrt(2.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
                                            Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
                                        }

                                        biMolData.Dtot += (Dr1 + Dr2) / (4.0 * params.timeStep); // add in contributions from rotation

                                        {
                                            // Only allow 2D. diffusion at certain intervals, to avoid generating too many 2D. Tables
                                            // Keep only one sig fig for <0.1, 2 for 0.1<d<10, 3 for 10<d<100, etc
                                            double dtmp;
                                            if (biMolData.Dtot < 0.0001)
                                                dtmp = biMolData.Dtot * 100000;
                                            else if (biMolData.Dtot < 0.001)
                                                dtmp = biMolData.Dtot * 10000;
                                            else if (biMolData.Dtot < 0.01)
                                                dtmp = biMolData.Dtot * 1000;
                                            else if (biMolData.Dtot < 0.1)
                                                dtmp = biMolData.Dtot * 100;
                                            else
                                                dtmp = biMolData.Dtot * 100;

                                            int d_ones = int(round(dtmp));

                                            if (biMolData.Dtot < 0.0001)
                                                biMolData.Dtot = d_ones * 0.00001;
                                            else if (biMolData.Dtot < 0.001)
                                                biMolData.Dtot = d_ones * 0.0001;
                                            else if (biMolData.Dtot < 0.01)
                                                biMolData.Dtot = d_ones * 0.001;
                                            else if (biMolData.Dtot < 0.1)
                                                biMolData.Dtot = d_ones * 0.01;
                                            else
                                                biMolData.Dtot = d_ones * 0.01;

                                            if (biMolData.Dtot < 1E-50)
                                                biMolData.Dtot = 0;
                                        }

                                        double RMax { 3.5 * sqrt(4.0 * biMolData.Dtot * params.timeStep) + bindRadius };
                                        double R1 { 0.0 };
                                        if (membraneObject.isSphere == true) {
                                            Coord iface11 = moleculeList[pro1Index].interfaceList[relIface1].coord;
                                            Coord iface22 = moleculeList[pro2Index].interfaceList[relIface2].coord;
                                            double r1 = iface11.get_magnitude();
                                            double r2 = iface22.get_magnitude();
                                            double r = (r1 + r2) / 2.0; //membraneObject.sphereR; //
                                            double theta = acos((iface11.x * iface22.x + iface11.y * iface22.y + iface11.z * iface22.z) / r1 / r2);
                                            R1 = r * theta;
                                        } else {
                                            double dx = moleculeList[pro1Index].interfaceList[relIface1].coord.x - moleculeList[pro2Index].interfaceList[relIface2].coord.x;
                                            double dy = moleculeList[pro1Index].interfaceList[relIface1].coord.y - moleculeList[pro2Index].interfaceList[relIface2].coord.y;
                                            double dz { (std::abs(complexList[moleculeList[pro1Index].myComIndex].D.z - 0) < 1E-10
                                                            && std::abs(complexList[moleculeList[pro2Index].myComIndex].D.z - 0) < 1E-10)
                                                    ? 0
                                                    : moleculeList[pro1Index].interfaceList[relIface1].coord.z - moleculeList[pro2Index].interfaceList[relIface2].coord.z };
                                            R1 = sqrt((dx * dx) + (dy * dy) + (dz * dz));
                                        }
                                        if (R1 < RMax * 10.0) {
                                            moleculeList[pro1Index].crossbase.push_back(pro2Index);
                                            moleculeList[pro2Index].crossbase.push_back(pro1Index);
                                            moleculeList[pro1Index].mycrossint.push_back(relIface1);
                                            moleculeList[pro2Index].mycrossint.push_back(relIface2);
                                            moleculeList[pro1Index].crossrxn.push_back(
                                                std::array<int, 3> { rxnIndex, 0, false });
                                            moleculeList[pro2Index].crossrxn.push_back(
                                                std::array<int, 3> { rxnIndex, 0, false });
                                            ++complexList[moleculeList[pro1Index].myComIndex].ncross;
                                            ++complexList[moleculeList[pro2Index].myComIndex].ncross;
                                            moleculeList[pro1Index].probvec.push_back(0);
                                            moleculeList[pro2Index].probvec.push_back(0);
                                        }
                                    } else {
                                        //3D reaction
                                        double Dtot = 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.x + complexList[moleculeList[pro2Index].myComIndex].D.x)
                                            + 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.y + complexList[moleculeList[pro2Index].myComIndex].D.y)
                                            + 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.z + complexList[moleculeList[pro2Index].myComIndex].D.z);

                                        BiMolData biMolData { pro1Index, pro2Index, moleculeList[pro1Index].myComIndex, moleculeList[pro2Index].myComIndex, relIface1, relIface2,
                                            absIface1, absIface2, Dtot, magMol1, magMol2 };

                                        double Dr1 {};
                                        if (std::abs(complexList[biMolData.com1Index].D.z - 0) < 1E-10) {
                                            double cf = cos(sqrt(2.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep));
                                            Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
                                            biMolData.Dtot += Dr1 / (4.0 * params.timeStep);
                                        } else {
                                            double cf = cos(sqrt(4.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep));
                                            Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
                                            biMolData.Dtot += Dr1 / (6.0 * params.timeStep);
                                        }

                                        double Dr2;
                                        if (std::abs(complexList[biMolData.com2Index].D.z - 0) < 1E-10) {
                                            double cf = cos(sqrt(2.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
                                            Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
                                            biMolData.Dtot += Dr2 / (4.0 * params.timeStep);
                                        } else {
                                            double cf = cos(sqrt(4.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
                                            Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
                                            biMolData.Dtot += Dr2 / (6.0 * params.timeStep);
                                        }

                                        double RMax { 3.0 * sqrt(6.0 * biMolData.Dtot * params.timeStep) + bindRadius };
                                        double R1 { 0.0 };
                                        double dx = moleculeList[pro1Index].interfaceList[relIface1].coord.x - moleculeList[pro2Index].interfaceList[relIface2].coord.x;
                                        double dy = moleculeList[pro1Index].interfaceList[relIface1].coord.y - moleculeList[pro2Index].interfaceList[relIface2].coord.y;
                                        double dz { (std::abs(complexList[moleculeList[pro1Index].myComIndex].D.z - 0) < 1E-10
                                                        && std::abs(complexList[moleculeList[pro2Index].myComIndex].D.z - 0) < 1E-10)
                                                ? 0
                                                : moleculeList[pro1Index].interfaceList[relIface1].coord.z - moleculeList[pro2Index].interfaceList[relIface2].coord.z };
                                        R1 = sqrt((dx * dx) + (dy * dy) + (dz * dz));
                                        if (R1 < RMax) {
                                            moleculeList[pro1Index].crossbase.push_back(pro2Index);
                                            moleculeList[pro2Index].crossbase.push_back(pro1Index);
                                            moleculeList[pro1Index].mycrossint.push_back(relIface1);
                                            moleculeList[pro2Index].mycrossint.push_back(relIface2);
                                            moleculeList[pro1Index].crossrxn.push_back(
                                                std::array<int, 3> { rxnIndex, 0, false });
                                            moleculeList[pro2Index].crossrxn.push_back(
                                                std::array<int, 3> { rxnIndex, 0, false });
                                            ++complexList[moleculeList[pro1Index].myComIndex].ncross;
                                            ++complexList[moleculeList[pro2Index].myComIndex].ncross;
                                            moleculeList[pro1Index].probvec.push_back(0);
                                            moleculeList[pro2Index].probvec.push_back(0);
                                        }
                                    } //end else 3D
                                }
                            }
                        }
                    }
                }
            }
        }

        // pro2 is bound
        if (moleculeList[pro2Index].bndlist.size() > 0 && molTemplateList[moleculeList[pro2Index].molTypeIndex].excludeVolumeBound == true) {
            // loop the bindlist to find the exclude partner of each interface
            for (int relIface1Itr { 0 }; relIface1Itr < moleculeList[pro2Index].bndlist.size(); ++relIface1Itr) {
                int relIface2 { moleculeList[pro2Index].bndlist[relIface1Itr] };
                int absIface2 { moleculeList[pro2Index].interfaceList[relIface2].index };
                int molTypeIndex2 { moleculeList[pro2Index].molTypeIndex };
                if (moleculeList[pro2Index].interfaceList[relIface2].isBound && molTemplateList[moleculeList[pro2Index].molTypeIndex].interfaceList[relIface2].excludeVolumeBoundList.empty() == false) { // make sure it's actually bound and need excludeVolumeBound
                    // loop the interfaceList of pro2
                    for (int relIface2Idx = 0; relIface2Idx < moleculeList[pro1Index].interfaceList.size(); ++relIface2Idx) {
                        int relIface1 { relIface2Idx };
                        int absIface1 { moleculeList[pro1Index].interfaceList[relIface1].index };
                        int molTypeIndex1 { moleculeList[pro1Index].molTypeIndex };
                        for (int indexItr { 0 }; indexItr < molTemplateList[moleculeList[pro2Index].molTypeIndex].interfaceList[relIface2].excludeVolumeBoundList.size(); ++indexItr) {
                            if (molTemplateList[moleculeList[pro2Index].molTypeIndex].interfaceList[relIface2].excludeVolumeBoundList[indexItr] == molTypeIndex1) {
                                double bindRadius { molTemplateList[moleculeList[pro2Index].molTypeIndex].interfaceList[relIface2].excludeRadiusList[indexItr] };
                                int rxnIndex { molTemplateList[moleculeList[pro2Index].molTypeIndex].interfaceList[relIface2].excludeVolumeBoundReactList[indexItr] };
                                if (molTemplateList[moleculeList[pro2Index].molTypeIndex].interfaceList[relIface2].excludeVolumeBoundIfaceList[indexItr] == relIface1) {
                                    // calculate the distance between the two infaces
                                    Vector ifaceVec { moleculeList[pro1Index].interfaceList[relIface1].coord
                                        - complexList[moleculeList[pro1Index].myComIndex].comCoord };
                                    Vector ifaceVec2 { moleculeList[pro2Index].interfaceList[relIface2].coord
                                        - complexList[moleculeList[pro2Index].myComIndex].comCoord };
                                    double magMol1 { ifaceVec.x * ifaceVec.x + ifaceVec.y * ifaceVec.y
                                        + ifaceVec.z * ifaceVec.z };
                                    double magMol2 { ifaceVec2.x * ifaceVec2.x + ifaceVec2.y * ifaceVec2.y
                                        + ifaceVec2.z * ifaceVec2.z };
                                    if (std::abs(complexList[moleculeList[pro1Index].myComIndex].D.z) < 1E-16 && std::abs(complexList[moleculeList[pro2Index].myComIndex].D.z) < 1E-16) {
                                        // both Complexes are on the membrane, evaluate as 2D reaction
                                        double Dtot = 1.0 / 2.0 * (complexList[moleculeList[pro1Index].myComIndex].D.x + complexList[moleculeList[pro2Index].myComIndex].D.x)
                                            + 1.0 / 2.0 * (complexList[moleculeList[pro1Index].myComIndex].D.y + complexList[moleculeList[pro2Index].myComIndex].D.y);

                                        BiMolData biMolData { pro1Index, pro2Index, moleculeList[pro1Index].myComIndex, moleculeList[pro2Index].myComIndex, relIface1, relIface2,
                                            absIface1, absIface2, Dtot, magMol1, magMol2 };

                                        double Dr1 {};
                                        {
                                            double cf { cos(sqrt(2.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep)) };
                                            Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
                                        }

                                        double Dr2 {};
                                        {
                                            double cf = cos(sqrt(2.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
                                            Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
                                        }

                                        biMolData.Dtot += (Dr1 + Dr2) / (4.0 * params.timeStep); // add in contributions from rotation

                                        {
                                            // Only allow 2D. diffusion at certain intervals, to avoid generating too many 2D. Tables
                                            // Keep only one sig fig for <0.1, 2 for 0.1<d<10, 3 for 10<d<100, etc
                                            double dtmp;
                                            if (biMolData.Dtot < 0.0001)
                                                dtmp = biMolData.Dtot * 100000;
                                            else if (biMolData.Dtot < 0.001)
                                                dtmp = biMolData.Dtot * 10000;
                                            else if (biMolData.Dtot < 0.01)
                                                dtmp = biMolData.Dtot * 1000;
                                            else if (biMolData.Dtot < 0.1)
                                                dtmp = biMolData.Dtot * 100;
                                            else
                                                dtmp = biMolData.Dtot * 100;

                                            int d_ones = int(round(dtmp));

                                            if (biMolData.Dtot < 0.0001)
                                                biMolData.Dtot = d_ones * 0.00001;
                                            else if (biMolData.Dtot < 0.001)
                                                biMolData.Dtot = d_ones * 0.0001;
                                            else if (biMolData.Dtot < 0.01)
                                                biMolData.Dtot = d_ones * 0.001;
                                            else if (biMolData.Dtot < 0.1)
                                                biMolData.Dtot = d_ones * 0.01;
                                            else
                                                biMolData.Dtot = d_ones * 0.01;

                                            if (biMolData.Dtot < 1E-50)
                                                biMolData.Dtot = 0;
                                        }

                                        double RMax { 3.5 * sqrt(4.0 * biMolData.Dtot * params.timeStep) + bindRadius };
                                        double R1 { 0.0 };
                                        if (membraneObject.isSphere == true) {
                                            Coord iface11 = moleculeList[pro1Index].interfaceList[relIface1].coord;
                                            Coord iface22 = moleculeList[pro2Index].interfaceList[relIface2].coord;
                                            double r1 = iface11.get_magnitude();
                                            double r2 = iface22.get_magnitude();
                                            double r = (r1 + r2) / 2.0; //membraneObject.sphereR; //
                                            double theta = acos((iface11.x * iface22.x + iface11.y * iface22.y + iface11.z * iface22.z) / r1 / r2);
                                            R1 = r * theta;
                                        } else {
                                            double dx = moleculeList[pro1Index].interfaceList[relIface1].coord.x - moleculeList[pro2Index].interfaceList[relIface2].coord.x;
                                            double dy = moleculeList[pro1Index].interfaceList[relIface1].coord.y - moleculeList[pro2Index].interfaceList[relIface2].coord.y;
                                            double dz { (std::abs(complexList[moleculeList[pro1Index].myComIndex].D.z - 0) < 1E-10
                                                            && std::abs(complexList[moleculeList[pro2Index].myComIndex].D.z - 0) < 1E-10)
                                                    ? 0
                                                    : moleculeList[pro1Index].interfaceList[relIface1].coord.z - moleculeList[pro2Index].interfaceList[relIface2].coord.z };
                                            R1 = sqrt((dx * dx) + (dy * dy) + (dz * dz));
                                        }
                                        if (R1 < RMax * 10.0) {
                                            moleculeList[pro1Index].crossbase.push_back(pro2Index);
                                            moleculeList[pro2Index].crossbase.push_back(pro1Index);
                                            moleculeList[pro1Index].mycrossint.push_back(relIface1);
                                            moleculeList[pro2Index].mycrossint.push_back(relIface2);
                                            moleculeList[pro1Index].crossrxn.push_back(
                                                std::array<int, 3> { rxnIndex, 0, false });
                                            moleculeList[pro2Index].crossrxn.push_back(
                                                std::array<int, 3> { rxnIndex, 0, false });
                                            ++complexList[moleculeList[pro1Index].myComIndex].ncross;
                                            ++complexList[moleculeList[pro2Index].myComIndex].ncross;
                                            moleculeList[pro1Index].probvec.push_back(0);
                                            moleculeList[pro2Index].probvec.push_back(0);
                                        }
                                    } else {
                                        //3D reaction
                                        double Dtot = 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.x + complexList[moleculeList[pro2Index].myComIndex].D.x)
                                            + 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.y + complexList[moleculeList[pro2Index].myComIndex].D.y)
                                            + 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.z + complexList[moleculeList[pro2Index].myComIndex].D.z);

                                        BiMolData biMolData { pro1Index, pro2Index, moleculeList[pro1Index].myComIndex, moleculeList[pro2Index].myComIndex, relIface1, relIface2,
                                            absIface1, absIface2, Dtot, magMol1, magMol2 };

                                        double Dr1 {};
                                        if (std::abs(complexList[biMolData.com1Index].D.z - 0) < 1E-10) {
                                            double cf = cos(sqrt(2.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep));
                                            Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
                                            biMolData.Dtot += Dr1 / (4.0 * params.timeStep);
                                        } else {
                                            double cf = cos(sqrt(4.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep));
                                            Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
                                            biMolData.Dtot += Dr1 / (6.0 * params.timeStep);
                                        }

                                        double Dr2;
                                        if (std::abs(complexList[biMolData.com2Index].D.z - 0) < 1E-10) {
                                            double cf = cos(sqrt(2.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
                                            Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
                                            biMolData.Dtot += Dr2 / (4.0 * params.timeStep);
                                        } else {
                                            double cf = cos(sqrt(4.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
                                            Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
                                            biMolData.Dtot += Dr2 / (6.0 * params.timeStep);
                                        }

                                        double RMax { 3.0 * sqrt(6.0 * biMolData.Dtot * params.timeStep) + bindRadius };
                                        double R1 { 0.0 };
                                        double dx = moleculeList[pro1Index].interfaceList[relIface1].coord.x - moleculeList[pro2Index].interfaceList[relIface2].coord.x;
                                        double dy = moleculeList[pro1Index].interfaceList[relIface1].coord.y - moleculeList[pro2Index].interfaceList[relIface2].coord.y;
                                        double dz { (std::abs(complexList[moleculeList[pro1Index].myComIndex].D.z - 0) < 1E-10
                                                        && std::abs(complexList[moleculeList[pro2Index].myComIndex].D.z - 0) < 1E-10)
                                                ? 0
                                                : moleculeList[pro1Index].interfaceList[relIface1].coord.z - moleculeList[pro2Index].interfaceList[relIface2].coord.z };
                                        R1 = sqrt((dx * dx) + (dy * dy) + (dz * dz));
                                        if (R1 < RMax) {
                                            moleculeList[pro1Index].crossbase.push_back(pro2Index);
                                            moleculeList[pro2Index].crossbase.push_back(pro1Index);
                                            moleculeList[pro1Index].mycrossint.push_back(relIface1);
                                            moleculeList[pro2Index].mycrossint.push_back(relIface2);
                                            moleculeList[pro1Index].crossrxn.push_back(
                                                std::array<int, 3> { rxnIndex, 0, false });
                                            moleculeList[pro2Index].crossrxn.push_back(
                                                std::array<int, 3> { rxnIndex, 0, false });
                                            ++complexList[moleculeList[pro1Index].myComIndex].ncross;
                                            ++complexList[moleculeList[pro2Index].myComIndex].ncross;
                                            moleculeList[pro1Index].probvec.push_back(0);
                                            moleculeList[pro2Index].probvec.push_back(0);
                                        }
                                    } //end else 3D
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // // now loop the bndlist of mol1 to determine excludeVolumeBound or not
    // if (molTemplateList[moleculeList[pro1Index].molTypeIndex].excludeVolumeBound == true) {
    //     for (int relIface1Itr { 0 }; relIface1Itr < moleculeList[pro1Index].bndlist.size(); ++relIface1Itr) {
    //         int relIface1 { moleculeList[pro1Index].bndlist[relIface1Itr] };
    //         int absIface1 { moleculeList[pro1Index].interfaceList[relIface1].index };
    //         if (moleculeList[pro1Index].interfaceList[relIface1].isBound && molTemplateList[moleculeList[pro1Index].molTypeIndex].interfaceList[relIface1].excludeVolumeBoundList.empty() == false) { // make sure it's actually bound and need excludeVolumeBound
    //             // loop the interfaceList of pro2
    //             for (int relIface2Idx = 0; relIface2Idx < moleculeList[pro2Index].interfaceList.size(); ++relIface2Idx) {
    //                 int relIface2 { relIface2Idx };
    //                 int absIface2 { moleculeList[pro2Index].interfaceList[relIface2].index };
    //                 for (auto one : molTemplateList[moleculeList[pro1Index].molTypeIndex].interfaceList[relIface1].excludeVolumeBoundList) {
    //                     if (one == absIface2) {
    //                         // calculate the distance between the two infaces
    //                         Vector ifaceVec { moleculeList[pro1Index].interfaceList[relIface1].coord
    //                             - complexList[moleculeList[pro1Index].myComIndex].comCoord };
    //                         Vector ifaceVec2 { moleculeList[pro2Index].interfaceList[relIface2].coord
    //                             - complexList[moleculeList[pro2Index].myComIndex].comCoord };
    //                         double magMol1 { ifaceVec.x * ifaceVec.x + ifaceVec.y * ifaceVec.y
    //                             + ifaceVec.z * ifaceVec.z };
    //                         double magMol2 { ifaceVec2.x * ifaceVec2.x + ifaceVec2.y * ifaceVec2.y
    //                             + ifaceVec2.z * ifaceVec2.z };
    //                         if (std::abs(complexList[moleculeList[pro1Index].myComIndex].D.z) < 1E-16 && std::abs(complexList[moleculeList[pro2Index].myComIndex].D.z) < 1E-16) {
    //                             // both Complexes are on the membrane, evaluate as 2D reaction
    //                             double Dtot = 1.0 / 2.0 * (complexList[moleculeList[pro1Index].myComIndex].D.x + complexList[moleculeList[pro2Index].myComIndex].D.x)
    //                                 + 1.0 / 2.0 * (complexList[moleculeList[pro1Index].myComIndex].D.y + complexList[moleculeList[pro2Index].myComIndex].D.y);

    //                             BiMolData biMolData { pro1Index, pro2Index, moleculeList[pro1Index].myComIndex, moleculeList[pro2Index].myComIndex, relIface1, relIface2,
    //                                 absIface1, absIface2, Dtot, magMol1, magMol2 };

    //                             double Dr1 {};
    //                             {
    //                                 double cf { cos(sqrt(2.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep)) };
    //                                 Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
    //                             }

    //                             double Dr2 {};
    //                             {
    //                                 double cf = cos(sqrt(2.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
    //                                 Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
    //                             }

    //                             biMolData.Dtot += (Dr1 + Dr2) / (4.0 * params.timeStep); // add in contributions from rotation

    //                             {
    //                                 // Only allow 2D. diffusion at certain intervals, to avoid generating too many 2D. Tables
    //                                 // Keep only one sig fig for <0.1, 2 for 0.1<d<10, 3 for 10<d<100, etc
    //                                 double dtmp;
    //                                 if (biMolData.Dtot < 0.0001)
    //                                     dtmp = biMolData.Dtot * 100000;
    //                                 else if (biMolData.Dtot < 0.001)
    //                                     dtmp = biMolData.Dtot * 10000;
    //                                 else if (biMolData.Dtot < 0.01)
    //                                     dtmp = biMolData.Dtot * 1000;
    //                                 else if (biMolData.Dtot < 0.1)
    //                                     dtmp = biMolData.Dtot * 100;
    //                                 else
    //                                     dtmp = biMolData.Dtot * 100;

    //                                 int d_ones = int(round(dtmp));

    //                                 if (biMolData.Dtot < 0.0001)
    //                                     biMolData.Dtot = d_ones * 0.00001;
    //                                 else if (biMolData.Dtot < 0.001)
    //                                     biMolData.Dtot = d_ones * 0.0001;
    //                                 else if (biMolData.Dtot < 0.01)
    //                                     biMolData.Dtot = d_ones * 0.001;
    //                                 else if (biMolData.Dtot < 0.1)
    //                                     biMolData.Dtot = d_ones * 0.01;
    //                                 else
    //                                     biMolData.Dtot = d_ones * 0.01;

    //                                 if (biMolData.Dtot < 1E-50)
    //                                     biMolData.Dtot = 0;
    //                             }

    //                             double RMax { 3.5 * sqrt(4.0 * biMolData.Dtot * params.timeStep) + 1.0 }; // TODO: here use bindRadius = 1.0
    //                             double R1 { 0.0 };
    //                             if (membraneObject.isSphere == true) {
    //                                 Coord iface11 = moleculeList[pro1Index].interfaceList[relIface1].coord;
    //                                 Coord iface22 = moleculeList[pro2Index].interfaceList[relIface2].coord;
    //                                 double r1 = iface11.get_magnitude();
    //                                 double r2 = iface22.get_magnitude();
    //                                 double r = (r1 + r2) / 2.0; //membraneObject.sphereR; //
    //                                 double theta = acos((iface11.x * iface22.x + iface11.y * iface22.y + iface11.z * iface22.z) / r1 / r2);
    //                                 R1 = r * theta;
    //                             } else {
    //                                 double dx = moleculeList[pro1Index].interfaceList[relIface1].coord.x - moleculeList[pro2Index].interfaceList[relIface2].coord.x;
    //                                 double dy = moleculeList[pro1Index].interfaceList[relIface1].coord.y - moleculeList[pro2Index].interfaceList[relIface2].coord.y;
    //                                 double dz { (std::abs(complexList[moleculeList[pro1Index].myComIndex].D.z - 0) < 1E-10
    //                                                 && std::abs(complexList[moleculeList[pro2Index].myComIndex].D.z - 0) < 1E-10)
    //                                         ? 0
    //                                         : moleculeList[pro1Index].interfaceList[relIface1].coord.z - moleculeList[pro2Index].interfaceList[relIface2].coord.z };
    //                                 R1 = sqrt((dx * dx) + (dy * dy) + (dz * dz));
    //                             }
    //                             if (R1 < RMax) {
    //                                 moleculeList[pro1Index].interfaceList[relIface1].excludeVolume = true;
    //                                 moleculeList[pro1Index].interfaceList[relIface1].excludeMolList.push_back(pro2Index);
    //                                 moleculeList[pro1Index].interfaceList[relIface1].excludeInfList.push_back(relIface2);
    //                             }
    //                         } else {
    //                             //3D reaction
    //                             double Dtot = 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.x + complexList[moleculeList[pro2Index].myComIndex].D.x)
    //                                 + 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.y + complexList[moleculeList[pro2Index].myComIndex].D.y)
    //                                 + 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.z + complexList[moleculeList[pro2Index].myComIndex].D.z);

    //                             BiMolData biMolData { pro1Index, pro2Index, moleculeList[pro1Index].myComIndex, moleculeList[pro2Index].myComIndex, relIface1, relIface2,
    //                                 absIface1, absIface2, Dtot, magMol1, magMol2 };

    //                             double Dr1 {};
    //                             if (std::abs(complexList[biMolData.com1Index].D.z - 0) < 1E-10) {
    //                                 double cf = cos(sqrt(2.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep));
    //                                 Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
    //                                 biMolData.Dtot += Dr1 / (4.0 * params.timeStep);
    //                             } else {
    //                                 double cf = cos(sqrt(4.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep));
    //                                 Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
    //                                 biMolData.Dtot += Dr1 / (6.0 * params.timeStep);
    //                             }

    //                             double Dr2;
    //                             if (std::abs(complexList[biMolData.com2Index].D.z - 0) < 1E-10) {
    //                                 double cf = cos(sqrt(2.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
    //                                 Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
    //                                 biMolData.Dtot += Dr2 / (4.0 * params.timeStep);
    //                             } else {
    //                                 double cf = cos(sqrt(4.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
    //                                 Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
    //                                 biMolData.Dtot += Dr2 / (6.0 * params.timeStep);
    //                             }

    //                             double RMax { 3.0 * sqrt(6.0 * biMolData.Dtot * params.timeStep) + 1.0 }; // TODO: here use bindRadius = 1.0
    //                             double R1 { 0.0 };
    //                             double dx = moleculeList[pro1Index].interfaceList[relIface1].coord.x - moleculeList[pro2Index].interfaceList[relIface2].coord.x;
    //                             double dy = moleculeList[pro1Index].interfaceList[relIface1].coord.y - moleculeList[pro2Index].interfaceList[relIface2].coord.y;
    //                             double dz { (std::abs(complexList[moleculeList[pro1Index].myComIndex].D.z - 0) < 1E-10
    //                                             && std::abs(complexList[moleculeList[pro2Index].myComIndex].D.z - 0) < 1E-10)
    //                                     ? 0
    //                                     : moleculeList[pro1Index].interfaceList[relIface1].coord.z - moleculeList[pro2Index].interfaceList[relIface2].coord.z };
    //                             R1 = sqrt((dx * dx) + (dy * dy) + (dz * dz));
    //                             if (R1 < RMax) {
    //                                 moleculeList[pro1Index].interfaceList[relIface1].excludeVolume = true;
    //                                 moleculeList[pro1Index].interfaceList[relIface1].excludeMolList.push_back(pro2Index);
    //                                 moleculeList[pro1Index].interfaceList[relIface1].excludeInfList.push_back(relIface2);
    //                             }
    //                         } //end else 3D
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
}
