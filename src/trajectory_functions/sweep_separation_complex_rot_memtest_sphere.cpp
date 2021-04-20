/*! \file sweep_separation_complex_rot_memtest_sphere.cpp
 * ### Created on 02/28/2020 by Yiben Fu
 * ### Purpose: to check overlap, the sphere boundary
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */

#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

void sweep_separation_complex_rot_memtest_sphere(int simItr, int pro1Index, Parameters& params, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns, const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject)
{
    /*
      In this version, complex comIndex1 is on the membrane.
      If both proteins are on the membrane (Dz==0), evaluate only xy displacement, not z.
      In this _complex_ version, it tests overlap not just for each protein, but for each complex, so all the proteins
     in a complex, before performing position updates.
      NEW 2018: IN THIS VERSION, IT DOES NOT ATTEMPT TO SOLVE OVERLAP FOR PROTEINS WITHIN THE SAME COMPLEX, SINCE THEY
     CANNOT DIFFUSE RELATIVE TO ONE ANOTHER!
      This routine checks whether protein p1 is overlapping any partners in its reaction
      zone at its new position that is given by its current position +traj. If it does
     overlap, the displacement by traj is rejected and a new position for itself and any overlapping partners are
     selected. Once it no longer overlaps anyone, this protein and its complex are moved and the partners retain their
     stored new displacements. If a protein has already updated its position (done in sequential order) then it cannot
     resample a new position, the current protein must still continue to avoid overlapping, however.
    */
    // TRACE();
    int comIndex1 { moleculeList[pro1Index].myComIndex };
    int startProIndex = pro1Index;
    int com1Size = complexList[comIndex1].memberList.size();

    int maxRows { 1 };
    for (auto memMol : complexList[comIndex1].memberList) {
        if (moleculeList[memMol].crossbase.size() > maxRows)
            maxRows = moleculeList[memMol].crossbase.size();
    }

    int ifaceList[maxRows * com1Size];
    int overlapList[maxRows * com1Size];
    int memCheckList[maxRows * com1Size];

    // int reflectList[complexList.size()]; // if this is 0, we need call reflect_traj; if this is 1, we do not need call reflect_traj
    // for (auto i : reflectList) { // initialize
    //     i = 0;
    // }

    /*The sampled displacement for p1 is stored in traj. the position from the
     previous step is still stored in bases[p1].xcom, etc, and will be updated
     at the end of this routine*/

    /*figure out i2*/
    for (int c { 0 }; c < com1Size; ++c) {
        pro1Index = complexList[comIndex1].memberList[c];
        for (int i { 0 }; i < moleculeList[pro1Index].crossbase.size(); ++i) {
            int p2 { moleculeList[pro1Index].crossbase[i] };
            int k2 { moleculeList[p2].myComIndex };
            if (complexList[k2].D.z < 1E-15) {
                memCheckList[maxRows * c + i] = 1;
            } else
                memCheckList[maxRows * c + i] = 0;
            int i1 { moleculeList[pro1Index].mycrossint[i] };
            std::array<int, 3> rxnItr = moleculeList[pro1Index].crossrxn[i];

            // get the partner interface
            ifaceList[maxRows * c + i] = (forwardRxns[rxnItr[0]].reactantListNew[0].relIfaceIndex == i1)
                ? forwardRxns[rxnItr[0]].reactantListNew[1].relIfaceIndex
                : forwardRxns[rxnItr[0]].reactantListNew[0].relIfaceIndex;
        }
    }

    //determine RS3Dinput
    double RS3Dinput { 0.0 };
    Complex targCom { complexList[comIndex1] };
    for (auto& molIndex : targCom.memberList) {
        for (int RS3Dindex = 0; RS3Dindex < 100; RS3Dindex++) {
            if (std::abs(membraneObject.RS3Dvect[RS3Dindex + 400] - moleculeList[molIndex].molTypeIndex) < 1E-2) {
                RS3Dinput = membraneObject.RS3Dvect[RS3Dindex + 300];
                break;
            }
        }
    }

    int tsave = 0;
    // if (reflectList[comIndex1] == 0) {
    //     reflect_traj_complex_rad_rot(params, moleculeList, complexList[comIndex1], membraneObject, RS3Dinput);
    //     reflectList[comIndex1] = 1;
    // }
    int itr { 0 };
    int maxItr { 10 };
    int saveit { 0 };
    while (itr < maxItr) {
        int numOverlap { 0 };
        bool hasOverlap { false };
        int comIndex2 {};
        double dr2 {};
        for (unsigned memMolItr { 0 }; memMolItr < complexList[comIndex1].memberList.size(); ++memMolItr) {
            pro1Index = complexList[comIndex1].memberList[memMolItr];
            for (int crossMemItr { 0 }; crossMemItr < moleculeList[pro1Index].crossbase.size(); ++crossMemItr) {
                int p2 { moleculeList[pro1Index].crossbase[crossMemItr] };
                if (moleculeList[p2].isImplicitLipid)
                    continue;

                comIndex2 = moleculeList[p2].myComIndex;
                /*Do not sweep for overlap if proteins are in the same complex, they cannot move relative to one
                 * another!
                 */
                if (comIndex1 != comIndex2) {
                    int relIface1 { moleculeList[pro1Index].mycrossint[crossMemItr] };
                    int rxnItr { moleculeList[pro1Index].crossrxn[crossMemItr][0] };
                    int relIface2 { ifaceList[maxRows * memMolItr + crossMemItr] };

                    double dx1, dy1, dz1;
                    if (complexList[comIndex1].D.z < 1E-15) { // complex on sphere surface
                        Coord ifacecrds = moleculeList[pro1Index].interfaceList[relIface1].coord;
                        Coord iface_final = calculate_update_position_interface(complexList[comIndex1], ifacecrds);
                        dx1 = iface_final.x;
                        dy1 = iface_final.y;
                        dz1 = iface_final.z;
                    } else { // complex inside sphere
                        Vector iface1Vec { moleculeList[pro1Index].interfaceList[relIface1].coord - complexList[comIndex1].comCoord };
                        std::array<double, 9> M = create_euler_rotation_matrix(complexList[comIndex1].trajRot);
                        iface1Vec = matrix_rotate(iface1Vec, M);
                        dx1 = complexList[comIndex1].comCoord.x + iface1Vec.x + complexList[comIndex1].trajTrans.x;
                        dy1 = complexList[comIndex1].comCoord.y + iface1Vec.y + complexList[comIndex1].trajTrans.y;
                        dz1 = complexList[comIndex1].comCoord.z + iface1Vec.z + complexList[comIndex1].trajTrans.z;
                    }

                    /*Now complex 2*/
                    // if (reflectList[comIndex2] == 0) {
                    //     reflect_traj_complex_rad_rot(params, moleculeList, complexList[comIndex2], membraneObject, RS3Dinput);
                    //     reflectList[comIndex2] = 1;
                    // }
                    double dx2, dy2, dz2;
                    if (complexList[comIndex2].D.z < 1E-15) { // complex on sphere surface
                        Coord ifacecrds = moleculeList[p2].interfaceList[relIface2].coord;
                        Coord iface_final = calculate_update_position_interface(complexList[comIndex2], ifacecrds);
                        dx2 = iface_final.x;
                        dy2 = iface_final.y;
                        dz2 = iface_final.z;
                    } else { // complex inside sphere
                        Vector iface2Vec { moleculeList[p2].interfaceList[relIface2].coord - complexList[comIndex2].comCoord };
                        std::array<double, 9> M2 = create_euler_rotation_matrix(complexList[comIndex2].trajRot);
                        iface2Vec = matrix_rotate(iface2Vec, M2);
                        dx2 = complexList[comIndex2].comCoord.x + iface2Vec.x + complexList[comIndex2].trajTrans.x;
                        dy2 = complexList[comIndex2].comCoord.y + iface2Vec.y + complexList[comIndex2].trajTrans.y;
                        dz2 = complexList[comIndex2].comCoord.z + iface2Vec.z + complexList[comIndex2].trajTrans.z;
                    }

                    /*separation*/
                    double df1 { dx1 - dx2 };
                    double df2 { dy1 - dy2 };
                    double df3 { dz1 - dz2 };

                    dr2 = (df1 * df1) + (df2 * df2) + (df3 * df3);

                    if (dr2 < forwardRxns[rxnItr].bindRadius * forwardRxns[rxnItr].bindRadius) {
                        // we neglect that the bindRadius on sphere surface is different from planar surface. since they are just a little different
                        /*reselect positions for protein p2*/
                        overlapList[numOverlap] = p2;
                        numOverlap++;
                        hasOverlap = true;
                    }
                } // ignore proteins within the same complex.
            }
        }
        /*Now resample positions of p1 and overlapList, if t>0, otherwise no overlap, so
         break from loop*/
        if (hasOverlap) {
            ++itr;
            if (complexList[comIndex1].D.z < 1E-15) { // complex on sphere surface
                Coord targTrans = create_complex_propagation_vectors_on_sphere(params, complexList[comIndex1]);
                complexList[comIndex1].trajTrans.x = targTrans.x;
                complexList[comIndex1].trajTrans.y = targTrans.y;
                complexList[comIndex1].trajTrans.z = targTrans.z;
                complexList[comIndex1].trajRot.x = sqrt(2.0 * params.timeStep * complexList[comIndex1].Dr.x) * GaussV();
                complexList[comIndex1].trajRot.y = sqrt(2.0 * params.timeStep * complexList[comIndex1].Dr.y) * GaussV();
                complexList[comIndex1].trajRot.z = sqrt(2.0 * params.timeStep * complexList[comIndex1].Dr.z) * GaussV();
            } else {
                complexList[comIndex1].trajTrans.x = sqrt(2.0 * params.timeStep * complexList[comIndex1].D.x) * GaussV();
                complexList[comIndex1].trajTrans.y = sqrt(2.0 * params.timeStep * complexList[comIndex1].D.y) * GaussV();
                complexList[comIndex1].trajTrans.z = sqrt(2.0 * params.timeStep * complexList[comIndex1].D.z) * GaussV();
                complexList[comIndex1].trajRot.x = sqrt(2.0 * params.timeStep * complexList[comIndex1].Dr.x) * GaussV();
                complexList[comIndex1].trajRot.y = sqrt(2.0 * params.timeStep * complexList[comIndex1].Dr.y) * GaussV();
                complexList[comIndex1].trajRot.z = sqrt(2.0 * params.timeStep * complexList[comIndex1].Dr.z) * GaussV();
            }

            // reflectList[comIndex1] = 0;
            reflect_traj_complex_rad_rot(params, moleculeList, complexList[comIndex1], membraneObject, RS3Dinput);
            // reflectList[comIndex1] = 1;

            int resampleList[complexList.size()]; // if this is 0, we need resample
            for (auto& i : resampleList) { // initialize
                i = 0;
            }

            for (int checkMolItr { 0 }; checkMolItr < numOverlap; checkMolItr++) {
                int p2 { overlapList[checkMolItr] };
                comIndex2 = moleculeList[p2].myComIndex;
                if (resampleList[comIndex2] == 0) {
                    if (p2 > startProIndex && (moleculeList[p2].trajStatus == TrajStatus::none || moleculeList[p2].trajStatus == TrajStatus::canBeResampled)) {
                        /*
                     We loop over proteins sequentially, so earlier proteins have already moved and avoided
                     their neighbors and should not be moved again.
                     These new positions selected for proteins not yet moved will be stored and
                     then used when they test for overlap themselves.
                     */

                        /*If p2 just dissociated, also don't try to move again*/
                        if (complexList[comIndex2].D.z < 1E-15) { // complex on sphere surface
                            Coord targTrans = create_complex_propagation_vectors_on_sphere(params, complexList[comIndex2]);
                            complexList[comIndex2].trajTrans.x = targTrans.x;
                            complexList[comIndex2].trajTrans.y = targTrans.y;
                            complexList[comIndex2].trajTrans.z = targTrans.z;
                            complexList[comIndex2].trajRot.x = sqrt(2.0 * params.timeStep * complexList[comIndex2].Dr.x) * GaussV();
                            complexList[comIndex2].trajRot.y = sqrt(2.0 * params.timeStep * complexList[comIndex2].Dr.y) * GaussV();
                            complexList[comIndex2].trajRot.z = sqrt(2.0 * params.timeStep * complexList[comIndex2].Dr.z) * GaussV();
                        } else {
                            complexList[comIndex2].trajTrans.x = sqrt(2.0 * params.timeStep * complexList[comIndex2].D.x) * GaussV();
                            complexList[comIndex2].trajTrans.y = sqrt(2.0 * params.timeStep * complexList[comIndex2].D.y) * GaussV();
                            complexList[comIndex2].trajTrans.z = sqrt(2.0 * params.timeStep * complexList[comIndex2].D.z) * GaussV();
                            complexList[comIndex2].trajRot.x = sqrt(2.0 * params.timeStep * complexList[comIndex2].Dr.x) * GaussV();
                            complexList[comIndex2].trajRot.y = sqrt(2.0 * params.timeStep * complexList[comIndex2].Dr.y) * GaussV();
                            complexList[comIndex2].trajRot.z = sqrt(2.0 * params.timeStep * complexList[comIndex2].Dr.z) * GaussV();
                        }

                        // reflectList[comIndex2] = 0;
                        reflect_traj_complex_rad_rot(params, moleculeList, complexList[comIndex2], membraneObject, RS3Dinput);
                        // reflectList[comIndex2] = 1;
                        resampleList[comIndex2] = 1;
                    }
                }
            }
            tsave = numOverlap;
        } else {
            saveit = itr;
            itr = maxItr; // break from loop
        }

        if (itr == maxItr - 1) {
            if (comIndex1 != comIndex2) {
                std::cout << " WARNING ***************************************************** " << '\n';
                std::cout << "At simItr " << simItr << ", pair check can't solve overlap: " << pro1Index
                          << " max cross: " << moleculeList[pro1Index].crossbase.size() << " n overlap: " << tsave
                          << " pro1: " << overlapList[0] << " D: " << complexList[comIndex1].D.x << " "
                          << complexList[comIndex2].D.x << " Last separation: " << sqrt(dr2) << '\n';

                // write_crds_complex(comIndex1, ind_com, bases);
                // write_crds_complex(comIndex2, ind_com, bases);
                for (int c { 0 }; c < com1Size; ++c) {
                    pro1Index = complexList[comIndex1].memberList[c];
                    // std::cout << "p1: " << pro1Index << ' ' << " nfree and com; "
                    //           << moleculeList[pro1Index].freelist.size() << ' ' << moleculeList[pro1Index].comCoord
                    //           << '\t';
                    // std::cout << "traj 1: " << ' ' << complexList[comIndex2].trajTrans << '\n';
                    for (int crossMolItr { 0 }; crossMolItr < moleculeList[pro1Index].crossbase.size(); ++crossMolItr) {
                        int pro2Index { moleculeList[pro1Index].crossbase[crossMolItr] };
                        comIndex2 = moleculeList[pro2Index].myComIndex;
                        int relIface1 { moleculeList[pro1Index].mycrossint[crossMolItr] };
                        int relIface2 { ifaceList[maxRows * c + crossMolItr] };
                        // std::cout << "cross num: " << crossMolItr << " i1: " << relIface1 << " i2: " << relIface2
                        //           << " pro: " << pro2Index << ' ' << " nfree; "
                        //           << moleculeList[pro2Index].freelist.size() << ' ' << moleculeList[pro2Index].comCoord
                        //           << '\t';
                        // std::cout << "traj: " << ' ' << complexList[comIndex2].trajTrans << '\n';
                    }
                }
            } else {
                std::cout << "Reached max iterations in sweep (memCheckList) without solving overlap, but proteins are "
                             "within the same complex. Last separation: "
                          << sqrt(dr2) << '\n';
            }
            // exit(1);
        }

    } // end maximum iterations

    complexList[comIndex1].propagate(moleculeList, membraneObject, molTemplateList);
    // Reset displacements to zero so distance is measured to your current
    // updated position that won't change again this turn
    complexList[comIndex1].trajTrans.zero_crds();
    complexList[comIndex1].trajRot.zero_crds();
}