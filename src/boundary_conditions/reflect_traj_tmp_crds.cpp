#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"

/*evaluates reflection out of box based on tmpCoords of all proteins. targCOM tmpCOM coords also must be consistent-updated.
 *Does not update molecule traj vectors, just updates the passed in vector traj, to collect for both complexes. 
 *
*/
void reflect_traj_tmp_crds(const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, std::array<double, 3>& traj, const Membrane& membraneObject, double RS3Dinput)
{
    /*This routine updated February 2020 to test if a large complex that spans the box or sphere could extend out in both directions
    if so, it attempts to correct for this by resampling the complex's translational and rotational updates.
    */
    if (membraneObject.isSphere == true)
        reflect_traj_tmp_crds_sphere(params, moleculeList, targCom, traj, membraneObject, RS3Dinput);
    else
        reflect_traj_tmp_crds_box(params, moleculeList, targCom, traj, membraneObject, RS3Dinput);
}

// #include "boundary_conditions/reflect_functions.hpp"
// #include "math/matrix.hpp"

// /*evaluates reflection out of box based on tmpCoords of all proteins. targCOM tmpCOM coords also must be consistent-updated.
//  *Does not update molecule traj vectors, just updates the passed in vector traj, to collect for both complexes.
//  *
// */
// void reflect_traj_tmp_crds(
//     const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, std::array<double, 3>& traj, const Membrane& membraneObject, double RS3Dinput)
// {
//     /*This routine updated October 2019 to test if a large complex that spans the box could extend out in both directions
//     if so, it attempts to correct for this by resampling the complex's translational and rotational updates.
//     */
//     double RS3D;
//     if (targCom.OnSurface) {
//         RS3D = 0;
//     } else {
//         RS3D = RS3Dinput;
//     }

//     // declare the six boundary sides of the system box;
//     double posX = membraneObject.waterBox.x / 2.0;
//     double negX = -membraneObject.waterBox.x / 2.0;
//     double posY = membraneObject.waterBox.y / 2.0;
//     double negY = -membraneObject.waterBox.y / 2.0;
//     double posZ = membraneObject.waterBox.z / 2.0;
//     double negZ = -membraneObject.waterBox.z / 2.0 + RS3D;

//     double currx { targCom.tmpComCoord.x + traj[0] };
//     double curry { targCom.tmpComCoord.y + traj[1] };
//     double currz { targCom.tmpComCoord.z + traj[2] };

//     std::array<double, 9> M {};
//     for (int mm = 0; mm < 9; mm++)
//         M[mm] = 0;
//     M[0] = 1;
//     M[4] = 1;
//     M[8] = 1; //set M to identity, perform no rotations here.

//     /*This is to test based on general size if it is close to boundaries, before doing detailed evaluation below.*/

//     bool canBeOutsideX { false };
//     if ((currx + targCom.radius) > posX || (currx - targCom.radius) < negX)
//         canBeOutsideX = true;

//     bool canBeOutsideY { false };
//     if ((curry + targCom.radius) > posY || (curry - targCom.radius) < negY)
//         canBeOutsideY = true;

//     bool canBeOutsideZ { false };
//     if ((currz + targCom.radius) > posZ || (currz - targCom.radius) < negZ)
//         canBeOutsideZ = true;

//     /*Now evaluate all interfaces distance from boundaries.*/
//     bool recheck { false };
//     if (canBeOutsideX) {
//         bool outside { false };
//         bool outsideNeg { false };
//         bool outsidePos { false };

//         double posWall = negX;
//         double negWall = posX;

//         std::array<double, 3> row {};
//         row[0] = M[0];
//         row[1] = M[1];
//         row[2] = M[2];

//         /*these need to be what current positions
//         due to translation and rotation are*/
//         for (auto& memMol : targCom.memberList) {
//             // measure each protein COM to z plane
//             Vector comVec { moleculeList[memMol].tmpComCoord - targCom.tmpComCoord };
//             double dxrot { row[0] * comVec.x + row[1] * comVec.y + row[2] * comVec.z };
//             currx = targCom.tmpComCoord.x + traj[0] + dxrot;

//             if (currx > posWall)
//                 posWall = currx; // farthest point in +X
//             if (currx < negWall)
//                 negWall = currx; // farthest point in -X

//             // measure each interface to x plane
//             for (int ii = 0; ii < moleculeList[memMol].interfaceList.size(); ii++) {
//                 Vector ifaceVec { moleculeList[memMol].tmpICoords[ii] - targCom.tmpComCoord };
//                 dxrot = row[0] * ifaceVec.x + row[1] * ifaceVec.y + row[2] * ifaceVec.z;
//                 currx = targCom.tmpComCoord.x + traj[0] + dxrot;

//                 if (currx > posWall)
//                     posWall = currx;
//                 if (currx < negWall)
//                     negWall = currx;
//             }
//         }
//         // check whether this complex is out of the box
//         if (posWall > posX) {
//             outside = true;
//             outsidePos = true;
//         }
//         if (negWall < negX) {
//             outside = true;
//             outsideNeg = true;
//         }

//         if (outside) {
//             // Put back inside the box, extended out
//             if (outsideNeg && outsidePos) {
//                 // For a large complex, test if it could be pushed back out the other side
//                 std::cout << " extends in both directions of X. negwall " << negWall
//                           << " poswall: " << posWall << '\n';
//                 recheck = true;
//             }
//             if (outsideNeg && !outsidePos) {
//                 // Also need to check that update will not push you out the other side
//                 traj[0] -= 2.0 * (negWall - negX);
//                 if (posWall - 2.0 * (negWall - negX) > posX) {
//                     recheck = true;
//                     std::cout << " Will push out the other side. negwall " << negWall << " poswall: " << posWall
//                               << "traj[0]: " << traj[0] << '\n';
//                 }
//             }
//             if (outsidePos && !outsideNeg) {
//                 traj[0] -= 2.0 * (posWall - posX);
//                 if (negWall - 2.0 * (posWall - posX) < negX) {
//                     recheck = true;
//                     std::cout << " Will push out the other side. negwall " << negWall << " poswall: " << posWall
//                               << "traj[0]: " << traj[0] << '\n';
//                 }
//             }
//         }
//     }

//     if (canBeOutsideY > 0) {
//         bool outside { false };
//         bool outsideNeg { false };
//         bool outsidePos { false };

//         double posWall { negY };
//         double negWall { posY };

//         std::array<double, 3> row {};
//         row[0] = M[3];
//         row[1] = M[4];
//         row[2] = M[5];

//         /*these need to be what current positions
//         due to translation and rotation are*/
//         for (auto& mp : targCom.memberList) {
//             // measure each protein COM
//             Vector comVec { moleculeList[mp].tmpComCoord - targCom.tmpComCoord };

//             // only need y component
//             double dyrot = row[0] * comVec.x + row[1] * comVec.y + row[2] * comVec.z;
//             curry = targCom.tmpComCoord.y + traj[1] + dyrot;

//             if (curry > posWall)
//                 posWall = curry;
//             if (curry < negWall)
//                 negWall = curry;

//             // measure each interface
//             for (int ii = 0; ii < moleculeList[mp].interfaceList.size(); ii++) {
//                 Vector ifaceVec { moleculeList[mp].tmpICoords[ii] - targCom.tmpComCoord };

//                 // only need y component
//                 dyrot = row[0] * ifaceVec.x + row[1] * ifaceVec.y + row[2] * ifaceVec.z;
//                 curry = targCom.tmpComCoord.y + traj[1] + dyrot;

//                 if (curry > posWall)
//                     posWall = curry;
//                 if (curry < negWall)
//                     negWall = curry;
//             }
//         }
//         // check whether this complex is out of the box
//         if (posWall > posY) {
//             outside = true;
//             outsidePos = true;
//         }
//         if (negWall < negY) {
//             outside = true;
//             outsideNeg = true;
//         }

//         if (outside) {
//             // Put back inside the box
//             if (outsideNeg && outsidePos) {
//                 /*For a large complex, test if it could be pushed back out the other side*/
//                 std::cout << " extends in both directions of Y. negWall " << negWall
//                           << " posWall: " << posWall << '\n';
//                 recheck = true;
//             }
//             if (outsideNeg && !outsidePos) {
//                 traj[1] -= 2.0 * (negWall - negY);
//                 // Also need to check that update will not push you out the other side
//                 if (posWall - 2.0 * (negWall - negY) > posY) {
//                     recheck = true;
//                     std::cout << " Will push out the other side. negwall " << negWall << " poswall: " << posWall
//                               << " traj[1]: " << traj[1] << '\n';
//                 }
//             }
//             if (outsidePos && !outsideNeg) {
//                 traj[1] -= 2.0 * (posWall - posY);
//                 if (negWall - 2.0 * (posWall - posY) < negY) {
//                     recheck = true;
//                     std::cout << " Will push out the other side. negwall " << negWall << " poswall: " << posWall
//                               << " traj[1]: " << traj[1] << '\n';
//                 }
//             }
//         }
//     }

//     if (canBeOutsideZ > 0) {
//         bool outside { false };
//         bool outsideNeg { false };
//         bool outsidePos { false };

//         double posWall { negZ };
//         double negWall { posZ };

//         std::array<double, 3> row {};
//         row[0] = M[6];
//         row[1] = M[7];
//         row[2] = M[8];

//         /*these need to be what current positions
//           due to translation and rotation are*/
//         for (auto& memMol : targCom.memberList) {
//             // measure each protein COM
//             Vector comVec { moleculeList[memMol].tmpComCoord - targCom.tmpComCoord };
//             double dzrot = row[0] * comVec.x + row[1] * comVec.y + row[2] * comVec.z;
//             currz = targCom.tmpComCoord.z + traj[2] + dzrot;

//             if (currz > posWall)
//                 posWall = currz;
//             if (currz < negWall)
//                 negWall = currz;

//             // measure each interface
//             for (int ii = 0; ii < moleculeList[memMol].interfaceList.size(); ii++) {
//                 Vector ifaceVec { moleculeList[memMol].tmpICoords[ii] - targCom.tmpComCoord };
//                 dzrot = row[0] * ifaceVec.x + row[1] * ifaceVec.y + row[2] * ifaceVec.z;
//                 currz = targCom.tmpComCoord.z + traj[2] + dzrot;

//                 if (currz > posWall)
//                     posWall = currz;
//                 if (currz < negWall)
//                     negWall = currz;
//             }
//         }
//         // check whether this complex is out of the box
//         if (posWall > posZ) {
//             outside = true;
//             outsidePos = true;
//         }
//         if (negWall < negZ) {
//             outside = true;
//             outsideNeg = true;
//         }

//         if (outside) {
//             // Put back inside the box
//             /*
//             if (targCom.D.z == 0)
//                 std::cout << "Putting complex " << targCom.index << " with ztot " << ztot
//                           << ", zchg " << (targCom.tmpComCoord.z + traj[2]) - membraneObject.waterBox.z / 2.0
//                           << ", currz " << targCom.tmpComCoord.z + traj[2] << " and traj[2] "
//                           << traj[2] << " back in box." << '\n';
//             */
//             if (outsideNeg && outsidePos) {
//                 /*For a large complex, test if it could be pushed back out the other side*/
//                 std::cout << "extends in both directions of Z. negwall " << negWall
//                           << " poswall: " << posWall << '\n';
//                 recheck = true;
//             }
//             if (outsideNeg && !outsidePos) {
//                 // Also need to check that update will not push zou out the other side
//                 traj[2] -= 2.0 * (negWall - negZ);
//                 if (posWall - 2.0 * (negWall - negZ) > posZ) {
//                     recheck = true;
//                     std::cout << " Will push out the other side. negwall " << negWall << " poswall: " << posWall
//                               << " traj[2]: " << traj[2] << '\n';
//                 }
//             }
//             if (outsidePos && !outsideNeg) {
//                 traj[2] -= 2.0 * (posWall - posZ);
//                 if (negWall - 2.0 * (posWall - posZ) < negZ) {
//                     recheck = true;
//                     std::cout << " Will push out the other side. negwall " << negWall << " poswall: " << posWall
//                               << " traj[2]: " << traj[2] << '\n';
//                 }
//             }
//         } // updating traj to reflect
//     }
// }
