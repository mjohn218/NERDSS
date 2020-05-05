/*! \file reflect_complex_targCom.radius_rot.cpp
 * ### Created on 02/25/2020 by Yiben Fu
 * ### Purpose
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */
#include "boundary_conditions/reflect_functions.hpp"

void reflect_complex_rad_rot(const Membrane& membraneObject, Complex& targCom, std::vector<Molecule>& moleculeList, double RS3Dinput)
{
    if (membraneObject.isSphere == true)
        reflect_complex_rad_rot_sphere(membraneObject, targCom, moleculeList, RS3Dinput);
    else
        reflect_complex_rad_rot_box(membraneObject, targCom, moleculeList, RS3Dinput);
}

// /*! \file reflect_complex_targCom.radius_rot.cpp

//  * ### Created on 11/8/18 by Matthew Varga
//  * ### Purpose
//  * ***
//  *
//  * ### Notes
//  * ***
//  *
//  * ### TODO List
//  * ***
//  */
// #include "boundary_conditions/reflect_functions.hpp"

// void reflect_complex_rad_rot(const Membrane& membraneObject, Complex& targCom, std::vector<Molecule>& moleculeList, double RS3Dinput)
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

//     double currx = targCom.comCoord.x;
//     double curry = targCom.comCoord.y;
//     double currz = targCom.comCoord.z;

//     bool canBeOutsideX { false };
//     if ((currx + targCom.radius) > posX || (currx - targCom.radius) < negX)
//         canBeOutsideX = true;

//     bool canBeOutsideY { false };
//     if ((curry + targCom.radius) > posY || (curry - targCom.radius) < negY)
//         canBeOutsideY = true;

//     bool canBeOutsideZ { false };
//     if ((currz + targCom.radius) > posZ || (currz - targCom.radius) < negZ)
//         canBeOutsideZ = true;

//     if (canBeOutsideX) {
//         double posWall = negX;
//         double negWall = posX;
//         bool outsideNeg = false;
//         bool outsidePos = false;
//         bool outside = false;
//         double posdx = 0; // to store the largest distance outside positive X side. mark +
//         double negdx = 0; // to store the largest distance outside negative X side. make -
//         // find the farthest point in +X (posWall) and -X (negWall)
//         for (auto& memMol : targCom.memberList) {
//             //measure each protein COM
//             currx = moleculeList[memMol].comCoord.x;
//             if (currx > posWall)
//                 posWall = currx;
//             if (currx < negWall)
//                 negWall = currx;
//             // measure each interface
//             for (auto& iface : moleculeList[memMol].interfaceList) {
//                 currx = iface.coord.x;
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
//             posdx = posWall - posX;
//         }
//         if (negWall < negX) {
//             outside = true;
//             outsideNeg = true;
//             negdx = negWall - negX;
//         }
//         if (outside) {
//             if (outsideNeg && outsidePos) {
//                 // extends out both the front and back.
//                 std::cout << "IN REFLECT COMPLEX RAD ROT, EXTEND in BOTH directions of X . ALREADY UPDATED POSITIONS. EXITING..." << '\n';
//                 exit(1);
//             }
//             if (outsideNeg && !outsidePos) {
//                 if (posWall - 2.0 * negdx > posX) {
//                     std::cout << "PROBLEM: IN REFLECT COMPLEX RAD ROT, EXTEND in NEGATIVE side of X: try to put back in the box, x: "
//                               << -negdx << "BUT will EXTEND again in POSITIVE side of X" << '\n';
//                     exit(1);
//                 }
//                 std::cout << "IN REFLECT COMPLEX RAD ROT, EXTEND in NEGATIVE direction of X: put back in the box, x: " << -negdx << '\n';
//                 // just update positions.Put back inside the box
//                 targCom.comCoord.x -= 2.0 * negdx;
//                 for (auto memMol : targCom.memberList) {
//                     moleculeList[memMol].comCoord.x -= 2.0 * negdx;
//                     for (auto& iface : moleculeList[memMol].interfaceList)
//                         iface.coord.x -= 2.0 * negdx;
//                 }
//             }
//             if (!outsideNeg && outsidePos) {
//                 if (negWall - 2.0 * posdx < negX) {
//                     std::cout << "PROBLEM: IN REFLECT COMPLEX RAD ROT, EXTEND in POSITIVE side of X: try to put back in the box, x: "
//                               << -posdx << "BUT will EXTEND again in NEGATIVE side of X" << '\n';
//                     exit(1);
//                 }
//                 std::cout << "IN REFLECT COMPLEX RAD ROT, EXTEND in POSITIVE direction of X: put back in the box, x: " << -posdx << '\n';
//                 // just update positions.Put back inside the box
//                 targCom.comCoord.x -= 2.0 * posdx;
//                 for (auto memMol : targCom.memberList) {
//                     moleculeList[memMol].comCoord.x -= 2.0 * posdx;
//                     for (auto& iface : moleculeList[memMol].interfaceList)
//                         iface.coord.x -= 2.0 * posdx;
//                 }
//             }
//         } // update traj for reflection
//     }

//     if (canBeOutsideY) {
//         double posWall = negY;
//         double negWall = posY;
//         bool outsideNeg { false };
//         bool outsidePos { false };
//         bool outside { false };
//         double posdy = 0; // to store the largest distance outside positive Y side. mark +
//         double negdy = 0; // to store the largest distance outside negative Y side. make -
//         // find the farthest point in +Y (posWall) and -Y (negWall)
//         for (auto& memMol : targCom.memberList) {
//             //measure each protein COM
//             curry = moleculeList[memMol].comCoord.y;
//             if (curry > posWall)
//                 posWall = curry;
//             if (curry < negWall)
//                 negWall = curry;
//             // measure each interface
//             for (auto& iface : moleculeList[memMol].interfaceList) {
//                 curry = iface.coord.y;
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
//             posdy = posWall - posY;
//         }
//         if (negWall < negY) {
//             outside = true;
//             outsideNeg = true;
//             negdy = negWall - negY;
//         }
//         if (outside) {
//             if (outsideNeg && outsidePos) {
//                 // extends out both the front and back.
//                 std::cout << "IN REFLECT COMPLEX RAD ROT, EXTEND in BOTH directions of Y . ALREADY UPDATED POSITIONS. EXITING..." << '\n';
//                 exit(1);
//             }
//             if (outsideNeg && !outsidePos) {
//                 if (posWall - 2.0 * negdy > posY) {
//                     std::cout << "PROBLEM: IN REFLECT COMPLEX RAD ROT, EXTEND in NEGATIVE side of Y: try to put back in the box, y: "
//                               << -negdy << "BUT will EXTEND again in POSITIVE side of Y" << '\n';
//                     exit(1);
//                 }
//                 std::cout << "IN REFLECT COMPLEX RAD ROT, EXTEND in NEGATIVE direction of Y: put back in the box, y: " << -negdy << '\n';
//                 // just update positions.Put back inside the box
//                 targCom.comCoord.y -= 2.0 * negdy;
//                 for (auto memMol : targCom.memberList) {
//                     moleculeList[memMol].comCoord.y -= 2.0 * negdy;
//                     for (auto& iface : moleculeList[memMol].interfaceList)
//                         iface.coord.y -= 2.0 * negdy;
//                 }
//             }
//             if (!outsideNeg && outsidePos) {
//                 if (negWall - 2.0 * posdy < negY) {
//                     std::cout << "PROBLEM: IN REFLECT COMPLEX RAD ROT, EXTEND in POSITIVE side of Y: try to put back in the box, y: "
//                               << -posdy << "BUT will EXTEND again in NEGATIVE side of Y" << '\n';
//                     exit(1);
//                 }
//                 std::cout << "IN REFLECT COMPLEX RAD ROT, EXTEND in POSITIVE direction of Y: put back in the box, y: " << -posdy << '\n';
//                 // just update positions.Put back inside the box
//                 targCom.comCoord.y -= 2.0 * posdy;
//                 for (auto memMol : targCom.memberList) {
//                     moleculeList[memMol].comCoord.y -= 2.0 * posdy;
//                     for (auto& iface : moleculeList[memMol].interfaceList)
//                         iface.coord.y -= 2.0 * posdy;
//                 }
//             }
//         } // update traj for reflection
//     }

//     if (canBeOutsideZ) {
//         double posWall = negZ;
//         double negWall = posZ;
//         bool outsideNeg { false };
//         bool outsidePos { false };
//         bool outside { false };
//         double posdz = 0; // to store the largest distance outside positive Z side. mark +
//         double negdz = 0; // to store the largest distance outside negative Z side. make -
//         // find the farthest point in +Z (posWall) and -Z (negWall)
//         for (auto& memMol : targCom.memberList) {
//             //measure each protein COM
//             currz = moleculeList[memMol].comCoord.z;
//             if (currz > posWall)
//                 posWall = currz;
//             if (currz < negWall)
//                 negWall = currz;
//             // measure each interface
//             for (auto& iface : moleculeList[memMol].interfaceList) {
//                 currz = iface.coord.z;
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
//             posdz = posWall - posZ;
//         }
//         if (negWall < negZ) {
//             outside = true;
//             outsideNeg = true;
//             negdz = negWall - negZ;
//         }
//         if (outside) {
//             if (outsideNeg && outsidePos) {
//                 // extends out both the front and back.
//                 std::cout << "IN REFLECT COMPLEX RAD ROT, EXTEND in BOTH directions of Z . ALREADY UPDATED POSITIONS. EXITING..." << '\n';
//                 exit(1);
//             }
//             if (outsideNeg && !outsidePos) {
//                 if (posWall - 2.0 * negdz > posZ) {
//                     std::cout << "PROBLEM: IN REFLECT COMPLEX RAD ROT, EXTEND in NEGATIVE side of Z: try to put back in the box, z: "
//                               << -negdz << "BUT will EXTEND again in POSITIVE side of Z" << '\n';
//                     exit(1);
//                 }
//                 std::cout << "IN REFLECT COMPLEX RAD ROT, EXTEND in NEGATIVE direction of Z: put back in the box, z: " << -negdz << '\n';
//                 // just update positions.Put back inside the box
//                 targCom.comCoord.z -= 2.0 * negdz;
//                 for (auto memMol : targCom.memberList) {
//                     moleculeList[memMol].comCoord.z -= 2.0 * negdz;
//                     for (auto& iface : moleculeList[memMol].interfaceList)
//                         iface.coord.z -= 2.0 * negdz;
//                 }
//             }
//             if (!outsideNeg && outsidePos) {
//                 if (negWall - 2.0 * posdz < negZ) {
//                     std::cout << "PROBLEM: IN REFLECT COMPLEX RAD ROT, EXTEND in POSITIVE side of Z: try to put back in the box, z: "
//                               << -posdz << "BUT will EXTEND again in NEGATIVE side of Z" << '\n';
//                     exit(1);
//                 }
//                 std::cout << "IN REFLECT COMPLEX RAD ROT, EXTEND in POSITIVE direction of Z: put back in the box, z: " << -posdz << '\n';
//                 // just update positions.Put back inside the box
//                 targCom.comCoord.z -= 2.0 * posdz;
//                 for (auto memMol : targCom.memberList) {
//                     moleculeList[memMol].comCoord.z -= 2.0 * posdz;
//                     for (auto& iface : moleculeList[memMol].interfaceList)
//                         iface.coord.z -= 2.0 * posdz;
//                 }
//             }
//         } // update traj for reflection
//     }
// }
