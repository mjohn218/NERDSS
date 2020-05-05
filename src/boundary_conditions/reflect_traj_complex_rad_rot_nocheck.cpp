#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"

void reflect_traj_complex_rad_rot_nocheck(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, const Membrane& membraneObject, double RS3Dinput)
{
    if (membraneObject.isSphere == true)
        reflect_traj_complex_rad_rot_nocheck_sphere(params, targCom, moleculeList, membraneObject, RS3Dinput);
    else
        reflect_traj_complex_rad_rot_nocheck_box(params, targCom, moleculeList, membraneObject, RS3Dinput);
}

// #include "boundary_conditions/reflect_functions.hpp"

// void reflect_traj_complex_rad_rot_nocheck(
//     const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, std::array<double, 9>& M, const Membrane& membraneObject, double RS3Dinput)
// {
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

//     double currx { targCom.comCoord.x + targCom.trajTrans.x };
//     double curry { targCom.comCoord.y + targCom.trajTrans.y };
//     double currz { targCom.comCoord.z + targCom.trajTrans.z };

//     bool canBeOutsideX { false };
//     if ((currx + targCom.radius) > posX || (currx - targCom.radius) < negX)
//         canBeOutsideX = true;

//     bool canBeOutsideY { false };
//     if ((curry + targCom.radius) > posY || (curry - targCom.radius) < negY)
//         canBeOutsideY = true;

//     bool canBeOutsideZ { false };
//     if ((currz + targCom.radius) > posZ || (currz - targCom.radius) < negZ)
//         canBeOutsideZ = true;

//     double posWallX { negX };
//     double negWallX { posX };
//     double posWallY { negY };
//     double negWallY { posY };
//     double posWallZ { negZ };
//     double negWallZ { posZ };
//     /*Z is separate to allow the interfaces to approach to the membrane
//      but don't need to test if the entire complex is far enough
//      away from the boundary.
//      */
//     if (canBeOutsideX) {
//         bool outside { false };

//         std::array<double, 3> row {};
//         row[0] = M[0];
//         row[1] = M[1];
//         row[2] = M[2];

//         // these need to be what current positions due to translation and rotation are
//         for (auto& memMol : targCom.memberList) {
//             // Measure each protein COM
//             Vector comVec { moleculeList[memMol].comCoord - targCom.comCoord };

//             // only need X component
//             double dxrot = row[0] * comVec.x + row[1] * comVec.y + row[2] * comVec.z;
//             currx = targCom.comCoord.x + targCom.trajTrans.x + dxrot;

//             if (currx > posWallX)
//                 posWallX = currx; // farthest point in +X
//             if (currx < negWallX)
//                 negWallX = currx; // farthest point in -X

//             // measure each interface
//             for (const auto& iface : moleculeList[memMol].interfaceList) {
//                 Vector ifaceVec { iface.coord - targCom.comCoord };

//                 // only need X component
//                 dxrot = row[0] * ifaceVec.x + row[1] * ifaceVec.y + row[2] * ifaceVec.z;
//                 currx = targCom.comCoord.x + targCom.trajTrans.x + dxrot;

//                 if (currx > posWallX)
//                     posWallX = currx;
//                 if (currx < negWallX)
//                     negWallX = currx;
//             }
//         }
//         if (posWallX > posX) {
//             outside = true;
//             targCom.trajTrans.x -= 2.0 * (posWallX - posX);
//         }
//         if (negWallX < negX) {
//             outside = true;
//             targCom.trajTrans.x -= 2.0 * (negWallX - negX);
//         }
//     }

//     if (canBeOutsideY > 0) {
//         bool outside { false };

//         std::array<double, 3> row {};
//         row[0] = M[3];
//         row[1] = M[4];
//         row[2] = M[5];

//         for (auto& memMol : targCom.memberList) {
//             // measure each protein COM
//             Vector comVec { moleculeList[memMol].comCoord - targCom.comCoord };

//             // only need y component
//             double dyrot = row[0] * comVec.x + row[1] * comVec.y + row[2] * comVec.z;
//             curry = targCom.comCoord.y + targCom.trajTrans.y + dyrot;

//             if (curry > posWallY)
//                 posWallY = curry; // farthest point in +Y
//             if (curry < negWallY)
//                 negWallY = curry; // farthest point in -Y

//             // measure each interface
//             for (const auto& iface : moleculeList[memMol].interfaceList) {
//                 Vector ifaceVec { iface.coord - targCom.comCoord };

//                 // only need y component
//                 dyrot = (row[0] * ifaceVec.x) + (row[1] * ifaceVec.y) + (row[2] * ifaceVec.z);
//                 curry = targCom.comCoord.y + targCom.trajTrans.y + dyrot;

//                 if (curry > posWallY)
//                     posWallY = curry;
//                 if (curry < negWallY)
//                     negWallY = curry;
//             }
//         }
//         if (posWallY > posY) {
//             outside = true;
//             targCom.trajTrans.y -= 2.0 * (posWallY - posY);
//         }
//         if (negWallY < negY) {
//             outside = true;
//             targCom.trajTrans.y -= 2.0 * (negWallY - negY);
//         }
//     }

//     if (canBeOutsideZ) {
//         bool outside { false };

//         std::array<double, 3> row {};
//         row[0] = M[6];
//         row[1] = M[7];
//         row[2] = M[8];

//         for (auto& memMol : targCom.memberList) {
//             // measure each protein com
//             Vector comVec { moleculeList[memMol].comCoord - targCom.comCoord };

//             // only need z component
//             double dzrot = row[0] * comVec.x + row[1] * comVec.y + row[2] * comVec.z;
//             currz = targCom.comCoord.z + targCom.trajTrans.z + dzrot;

//             if (currz > posWallZ)
//                 posWallZ = currz; // farthest point in +Z
//             if (currz < negWallZ)
//                 negWallZ = currz; // farthest point in -Z

//             // measure each interface to z plane
//             for (const auto& iface : moleculeList[memMol].interfaceList) {
//                 Vector ifaceVec { iface.coord - targCom.comCoord };

//                 // only need z component
//                 dzrot = row[0] * ifaceVec.x + row[1] * ifaceVec.y + row[2] * ifaceVec.z;
//                 currz = targCom.comCoord.z + targCom.trajTrans.z + dzrot;

//                 if (currz > posWallZ)
//                     posWallZ = currz;
//                 if (currz < negWallZ)
//                     negWallZ = currz;
//             }
//         }
//         if (posWallZ > posZ) {
//             outside = true;
//             targCom.trajTrans.z -= 2.0 * (posWallZ - posZ);
//         }
//         if (negWallZ < negZ) {
//             outside = true;
//             targCom.trajTrans.z -= 2.0 * (negWallZ - negZ);
//         }
//     }
// }
