#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

void sweep_separation_complex_rot(int simItr, int pro1Index, Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject)
{
    if (membraneObject.isSphere)
        sweep_separation_complex_rot_sphere(simItr, pro1Index, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);
    else
        sweep_separation_complex_rot_box(simItr, pro1Index, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);

    // /*
    //   In this version, complex com1Index is on the membrane.
    //   If both proteins are on the membrane (Dz==0), evaluate only xy displacement, not z.

    //   In this _complex_ version, itr tests overlap not just for each protein, but for each complex, so all the proteins
    //  in a complex, before performing position updates.

    //   NEW 2018: IN THIS VERSION, IT DOES NOT ATTEMPT TO SOLVE OVERLAP FOR PROTEINS WITHIN THE SAME COMPLEX, SINCE THEY
    //  CANNOT DIFFUSE RELATIVE TO ONE ANOTHER!

    //   This routine checks whether protein p1 is overlapping any partners in its reaction
    //   zone at its new position that is given by its current position +traj. If itr does
    //  overlap, the displacement by traj is rejected and a new position for itself and any overlapping partners are
    //  selected. Once itr no longer overlaps anyone, this protein and its complex are moved and the partners retain their
    //  stored new displacements. If a protein has already updated its position (done in sequential order) then itr cannot
    //  resample a new position, the current protein must still continue to avoid overlapping, however.
    // */

    // int com1Index { moleculeList[pro1Index].myComIndex };
    // size_t com1Size { complexList[com1Index].memberList.size() };

    // int maxRows { 1 };
    // for (auto memMol : complexList[com1Index].memberList) {
    //     if (moleculeList[memMol].crossbase.size() > maxRows)
    //         maxRows = moleculeList[memMol].crossbase.size();
    // }

    // int ifaceList[maxRows * com1Size];
    // int overlapList[maxRows * com1Size];

    // /*The sampled displacement for p1 is stored in traj. the position from the
    //  previous step is still stored in bases[p1].xcom, etc, and will be updated
    //  at the end of this routine*/

    // /*figure out i2*/
    // for (int memMolItr { 0 }; memMolItr < com1Size; ++memMolItr) {
    //     pro1Index = complexList[com1Index].memberList[memMolItr];
    //     for (int i { 0 }; i < moleculeList[pro1Index].crossbase.size(); ++i) {
    //         int i1 { moleculeList[pro1Index].mycrossint[i] };
    //         std::array<int, 3> rxnItr = moleculeList[pro1Index].crossrxn[i];

    //         ifaceList[maxRows * memMolItr + i] = (forwardRxns[rxnItr[0]].reactantListNew[0].relIfaceIndex == i1)
    //             ? forwardRxns[rxnItr[0]].reactantListNew[1].relIfaceIndex
    //             : forwardRxns[rxnItr[0]].reactantListNew[0].relIfaceIndex;
    //     }
    // }

    // //determine RS3Dinput
    // double RS3Dinput { 0.0 };
    // Complex targCom { complexList[com1Index] };
    // for (auto& molIndex : targCom.memberList) {
    //     for (int RS3Dindex = 0; RS3Dindex < 100; RS3Dindex++) {
    //         if (std::abs(membraneObject.RS3Dvect[RS3Dindex + 400] - moleculeList[molIndex].molTypeIndex) < 1E-2) {
    //             RS3Dinput = membraneObject.RS3Dvect[RS3Dindex + 300];
    //             break;
    //         }
    //     }
    // }

    // int tsave = 0;
    // std::array<double, 9> M = create_euler_rotation_matrix(complexList[com1Index].trajRot);
    // reflect_traj_complex_rad_rot(params, moleculeList, complexList[com1Index], M, membraneObject, RS3Dinput);
    // int itr { 0 };
    // int maxItr { 50 };
    // int saveit { 0 };
    // while (itr < maxItr) {
    //     int numOverlap { 0 };
    //     bool hasOverlap { false };
    //     int com2Index {};
    //     double dr2 {};
    //     for (unsigned memMolItr { 0 }; memMolItr < com1Size; ++memMolItr) {
    //         pro1Index = complexList[com1Index].memberList[memMolItr];
    //         for (int crossMolItr { 0 }; crossMolItr < moleculeList[pro1Index].crossbase.size(); ++crossMolItr) {
    //             int pro2Index { moleculeList[pro1Index].crossbase[crossMolItr] };
    //             com2Index = moleculeList[pro2Index].myComIndex;
    //             /*Do not sweep for overlap if proteins are in the same complex, they cannot move relative to one
    //              * another!
    //              */
    //             if (com1Index != com2Index) {
    //                 int i1 { moleculeList[pro1Index].mycrossint[crossMolItr] };
    //                 std::array<int, 3> rxnItr = moleculeList[pro1Index].crossrxn[crossMolItr];
    //                 int i2 { ifaceList[maxRows * memMolItr + crossMolItr] };

    //                 Vector iface1Vec { moleculeList[pro1Index].interfaceList[i1].coord
    //                     - complexList[com1Index].comCoord };
    //                 iface1Vec = matrix_rotate(iface1Vec, M);

    //                 double dx1 { complexList[com1Index].comCoord.x + iface1Vec.x + complexList[com1Index].trajTrans.x };
    //                 double dy1 { complexList[com1Index].comCoord.y + iface1Vec.y + complexList[com1Index].trajTrans.y };
    //                 double dz1 { complexList[com1Index].comCoord.z + iface1Vec.z + complexList[com1Index].trajTrans.z };

    //                 /*Now complex 2*/
    //                 std::array<double, 9> M2 = create_euler_rotation_matrix(complexList[com2Index].trajRot);
    //                 reflect_traj_complex_rad_rot(params, moleculeList, complexList[com2Index], M2, membraneObject, RS3Dinput);

    //                 Vector iface2Vec { moleculeList[pro2Index].interfaceList[i2].coord
    //                     - complexList[com2Index].comCoord };
    //                 iface2Vec = matrix_rotate(iface2Vec, M2);
    //                 double dx2 { complexList[com2Index].comCoord.x + iface2Vec.x + complexList[com2Index].trajTrans.x };
    //                 double dy2 { complexList[com2Index].comCoord.y + iface2Vec.y + complexList[com2Index].trajTrans.y };
    //                 double dz2 { complexList[com2Index].comCoord.z + iface2Vec.z + complexList[com2Index].trajTrans.z };

    //                 if (dz2 < -membraneObject.waterBox.z / 2.0) {
    //                     std::cout << "In sweep, current pro2Index protein interface is below the membrane: index: "
    //                               << pro2Index << " ,iface: " << i2 << ", dz: " << dz2 << " iteration: " << simItr
    //                               << '\n';
    //                 }

    //                 /*separation*/
    //                 double df1 { dx1 - dx2 };
    //                 double df2 { dy1 - dy2 };
    //                 double df3 { dz1 - dz2 };

    //                 dr2 = (df1 * df1) + (df2 * df2) + (df3 * df3);

    //                 if (dr2 < forwardRxns[rxnItr[0]].bindRadius * forwardRxns[rxnItr[0]].bindRadius) {
    //                     /*reselect positions for protein pro2Index*/
    //                     overlapList[numOverlap] = pro2Index;
    //                     ++numOverlap;
    //                     hasOverlap = true;
    //                 }
    //             } // ignore proteins within the same complex.
    //         }
    //     }
    //     /*Now resample positions of p1 and overlapList, if numOverlap>0, otherwise no overlap, so
    //      break from loop*/
    //     if (hasOverlap) {
    //         ++itr;
    //         complexList[com1Index].trajTrans.x = sqrt(2.0 * params.timeStep * complexList[com1Index].D.x) * GaussV();
    //         complexList[com1Index].trajTrans.y = sqrt(2.0 * params.timeStep * complexList[com1Index].D.y) * GaussV();
    //         complexList[com1Index].trajTrans.z = sqrt(2.0 * params.timeStep * complexList[com1Index].D.z) * GaussV();
    //         complexList[com1Index].trajRot.x = sqrt(2.0 * params.timeStep * complexList[com1Index].Dr.x) * GaussV();
    //         complexList[com1Index].trajRot.y = sqrt(2.0 * params.timeStep * complexList[com1Index].Dr.y) * GaussV();
    //         complexList[com1Index].trajRot.z = sqrt(2.0 * params.timeStep * complexList[com1Index].Dr.z) * GaussV();

    //         M = create_euler_rotation_matrix(complexList[com1Index].trajRot);
    //         reflect_traj_complex_rad_rot(params, moleculeList, complexList[com1Index], M, membraneObject, RS3Dinput);

    //         for (int checkMolItr { 0 }; checkMolItr < numOverlap; checkMolItr++) {
    //             int p2 { overlapList[checkMolItr] };
    //             com2Index = moleculeList[p2].myComIndex;
    //             if (moleculeList[p2].trajStatus == TrajStatus::none
    //                 || moleculeList[p2].trajStatus == TrajStatus::canBeResampled) {
    //                 /*
    //                  We loop over proteins sequentially, so earlier proteins have already moved and avoided
    //                  their neighbors and should not be moved again.
    //                  These new positions selected for proteins not yet moved will be stored and
    //                  then used when they test for overlap themselves.
    //                  */

    //                 /*If p2 just dissociated, also don'numOverlap try to move again*/
    //                 complexList[com2Index].trajTrans.x
    //                     = sqrt(2.0 * params.timeStep * complexList[com2Index].D.x) * GaussV();
    //                 complexList[com2Index].trajTrans.y
    //                     = sqrt(2.0 * params.timeStep * complexList[com2Index].D.y) * GaussV();
    //                 complexList[com2Index].trajTrans.z
    //                     = sqrt(2.0 * params.timeStep * complexList[com2Index].D.z) * GaussV();
    //                 complexList[com2Index].trajRot.x
    //                     = sqrt(2.0 * params.timeStep * complexList[com2Index].Dr.x) * GaussV();
    //                 complexList[com2Index].trajRot.y
    //                     = sqrt(2.0 * params.timeStep * complexList[com2Index].Dr.y) * GaussV();
    //                 complexList[com2Index].trajRot.z
    //                     = sqrt(2.0 * params.timeStep * complexList[com2Index].Dr.z) * GaussV();
    //                 std::array<double, 9> M2 = create_euler_rotation_matrix(complexList[com2Index].trajRot);
    //                 reflect_traj_complex_rad_rot(params, moleculeList, complexList[com2Index], M2, membraneObject, RS3Dinput);
    //             }
    //         }
    //         tsave = numOverlap;
    //     } else {
    //         saveit = itr;
    //         itr = maxItr; // break from loop
    //     }

    //     if (itr == maxItr - 1) {
    //         if (com1Index != com2Index) {
    //             std::cout << " WARNING ***************************************************** " << '\n';
    //             std::cout << "can't solve overlap: " << pro1Index
    //                       << " max cross: " << moleculeList[pro1Index].crossbase.size() << " n overlap: " << tsave
    //                       << " pro1: " << overlapList[0] << " D: " << complexList[com1Index].D.x << " "
    //                       << complexList[com2Index].D.x << " Last separation: " << sqrt(dr2) << '\n';

    //             // write_crds_complex(com1Index, ind_com, bases);
    //             // write_crds_complex(com2Index, ind_com, bases);
    //             for (int memMolItr { 0 }; memMolItr < com1Size; ++memMolItr) {
    //                 pro1Index = complexList[com1Index].memberList[memMolItr];
    //                 std::cout << "p1: " << pro1Index << ' ' << " nfree and com; "
    //                           << moleculeList[pro1Index].freelist.size() << ' ' << moleculeList[pro1Index].comCoord
    //                           << '\t';
    //                 std::cout << "traj 1: " << ' ' << complexList[com2Index].trajTrans << '\n';
    //                 for (int crossMolItr { 0 }; crossMolItr < moleculeList[pro1Index].crossbase.size(); ++crossMolItr) {
    //                     int pro2Index { moleculeList[pro1Index].crossbase[crossMolItr] };
    //                     com2Index = moleculeList[pro2Index].myComIndex;
    //                     int relIface1 { moleculeList[pro1Index].mycrossint[crossMolItr] };
    //                     int relIface2 { ifaceList[maxRows * memMolItr + crossMolItr] };
    //                     std::cout << "cross num: " << crossMolItr << " relIface1: " << relIface1
    //                               << " relIface2: " << relIface2 << " pro: " << pro2Index << ' ' << " nfree; "
    //                               << moleculeList[pro2Index].freelist.size() << ' ' << moleculeList[pro2Index].comCoord
    //                               << '\t';
    //                     std::cout << "traj: " << ' ' << complexList[com2Index].trajTrans << '\n';
    //                 }
    //             }
    //         } else {
    //             std::cout << "Reached max iterations in sweep (memtest) without solving overlap, but proteins are "
    //                          "within the same complex. Last separation: "
    //                       << sqrt(dr2) << '\n';
    //         }
    //         // exit(1);
    //     }

    // } // end maximum iterations

    // complexList[com1Index].propagate(moleculeList);

    // // Reset displacements to zero so distance is measured to your current
    // // updated position that won't change again this turn
    // complexList[com1Index].trajTrans.zero_crds();
    // complexList[com1Index].trajRot.zero_crds();
}
