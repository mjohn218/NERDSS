#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"

void reflect_traj_complex_rad_rot(
    const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, const Membrane& membraneObject, double RS3Dinput)
{
    if (membraneObject.isSphere == true)
        reflect_traj_complex_rad_rot_sphere(params, moleculeList, targCom, membraneObject, RS3Dinput);
    else
        reflect_traj_complex_rad_rot_box(params, moleculeList, targCom, membraneObject, RS3Dinput);

    // // NOTE: it only works for a box system with the membrane surface located on the Z-bottom.

    // // for implicit-lipid model, the boundary surface must consider the reflecting-surface RS3D
    // // for explicit-lipid model, RS3D = 0.
    // double RS3D = RS3Dinput;
    // //    std::cout <<" RS3D: "<<membraneObject.RS3D<<std::endl;
    // /*This routine updated March 2017 to test if a large complex that spans the box could extend out in both directions
    // if so, it attempts to correct for this by resampling the complex's translational and rotational updates.
    // */
    // double xtot { 0 };
    // double ytot { 0 };
    // double ztot { 0 };
    // double currx { targCom.comCoord.x + targCom.trajTrans.x };
    // double curry { targCom.comCoord.y + targCom.trajTrans.y };
    // double currz { targCom.comCoord.z + targCom.trajTrans.z };

    // /*This is to test based on general size if it is close to boundaries, before doing detailed evaluation below.*/

    // bool canBeOutsideX { false };
    // double xchg = currx - membraneObject.waterBox.x / 2.0;
    // if ((xchg + targCom.radius) > 0 || (xchg - targCom.radius) < -membraneObject.waterBox.x)
    //     canBeOutsideX = true;

    // bool canBeOutsideY { false };
    // double ychg = curry - membraneObject.waterBox.y / 2.0;
    // if ((ychg + targCom.radius) > 0 || (ychg - targCom.radius) < -membraneObject.waterBox.y)
    //     canBeOutsideY = true;

    // bool canBeOutsideZ { false };
    // double zchg = currz - membraneObject.waterBox.z / 2.0;
    // if ((zchg + targCom.radius) > 0 || (zchg - targCom.radius) < -(membraneObject.waterBox.z - RS3D))
    //     canBeOutsideZ = true;

    // /*Now evaluate all interfaces distance from boundaries.*/
    // bool recheck { false };
    // if (canBeOutsideX) {
    //     bool outside { false };
    //     bool outsideNeg { false };
    //     bool outsidePos { false };

    //     double posWall = membraneObject.waterBox.x;
    //     double negWall = membraneObject.waterBox.x;

    //     std::array<double, 3> row {};
    //     row[0] = M[0];
    //     row[1] = M[1];
    //     row[2] = M[2];

    //     /*these need to be what current positions
    //     due to translation and rotation are*/
    //     for (auto& memMol : targCom.memberList) {
    //         /*measure each protein COM to z plane*/
    //         Vector comVec { moleculeList[memMol].comCoord - targCom.comCoord };

    //         double dxrot { row[0] * comVec.x + row[1] * comVec.y + row[2] * comVec.z };
    //         currx = targCom.comCoord.x + targCom.trajTrans.x + dxrot;

    //         if (membraneObject.waterBox.x / 2.0 - currx < posWall)
    //             posWall = membraneObject.waterBox.x / 2.0 - currx; // shortest distance in x from the +wall
    //         if (currx + membraneObject.waterBox.x / 2.0 < negWall)
    //             negWall = currx + membraneObject.waterBox.x / 2.0; // shortest distance in x from the -wall

    //         xchg = currx - membraneObject.waterBox.x / 2.0;
    //         if (xchg > 0) {
    //             outsidePos = true;
    //             outside = true;
    //             if (-xchg < xtot)
    //                 xtot = -xchg;
    //         } else if (xchg < -membraneObject.waterBox.x) {
    //             outsideNeg = true;
    //             outside = true;
    //             if (-(xchg + membraneObject.waterBox.x) > xtot)
    //                 xtot = -(xchg + membraneObject.waterBox.x);
    //         }
    //         /*measure each interface to z plane*/
    //         for (auto& iface : moleculeList[memMol].interfaceList) {
    //             Vector ifaceVec { iface.coord - targCom.comCoord };

    //             dxrot = row[0] * ifaceVec.x + row[1] * ifaceVec.y + row[2] * ifaceVec.z;
    //             currx = targCom.comCoord.x + targCom.trajTrans.x + dxrot;
    //             if (membraneObject.waterBox.x / 2.0 - currx < posWall)
    //                 posWall = membraneObject.waterBox.x / 2.0 - currx; // shortest distance in x from the +wall
    //             if (currx + membraneObject.waterBox.x / 2.0 < negWall)
    //                 negWall = currx + membraneObject.waterBox.x / 2.0; // shortest distance in x from the -wall

    //             xchg = currx - membraneObject.waterBox.x / 2.0;
    //             if (xchg > 0) {
    //                 outsidePos = true;
    //                 outside = true;
    //                 if (-xchg < xtot)
    //                     xtot = -xchg;
    //             } else if (xchg < -membraneObject.waterBox.x) {
    //                 outsideNeg = true;
    //                 outside = true;
    //                 if (-(xchg + membraneObject.waterBox.x) > xtot)
    //                     xtot = -(xchg + membraneObject.waterBox.x);
    //             }
    //         }
    //     }

    //     if (outside) {
    //         /*Put back inside the box, extended out*/
    //         targCom.trajTrans.x += 2.0 * xtot;
    //         if (outsideNeg > 0 && outsidePos > 0) {
    //             /*For a large complex, test if it could be pushed back out the other side*/
    //             std::cout << "xtot: " << xtot << " but also extends in other direction. negwall " << negWall
    //                       << " poswall: " << posWall << '\n';
    //             recheck = true;
    //         } else if (outsideNeg && posWall < (2.0 * xtot)) {
    //             /*Also need to check that update will not push you out the other side*/
    //             // xtot is positive.
    //             recheck = true;
    //             std::cout << " Will push out the other side. negwall " << negWall << " poswall: " << posWall
    //                       << " 2*xtot: " << 2.0 * xtot << '\n';
    //         } else if (outsidePos && negWall < (-2.0 * xtot)) {
    //             // fpos>0, which means xtot is negative.
    //             recheck = true;
    //             std::cout << " Will push out the other side. negwall " << negWall << " poswall: " << posWall
    //                       << " xtot: " << xtot << '\n';
    //         }
    //     }
    // }
    // if (canBeOutsideY > 0) {
    //     bool outside { false };
    //     bool outsideNeg { false };
    //     bool outsidePos { false };

    //     double poswall { membraneObject.waterBox.y };
    //     double negwall { membraneObject.waterBox.y };

    //     std::array<double, 3> row {};
    //     row[0] = M[3];
    //     row[1] = M[4];
    //     row[2] = M[5];

    //     /*these need to be what current positions
    //     due to translation and rotation are*/
    //     // for (i = 0; i < s1; i++) {
    //     for (auto& mp : targCom.memberList) {
    //         /*measure each protein COM to y plane*/
    //         Vector comVec { moleculeList[mp].comCoord - targCom.comCoord };

    //         /*only need y component*/
    //         double dyrot = row[0] * comVec.x + row[1] * comVec.y + row[2] * comVec.z;
    //         curry = targCom.comCoord.y + targCom.trajTrans.y + dyrot;
    //         if (membraneObject.waterBox.y / 2.0 - curry < poswall)
    //             poswall = membraneObject.waterBox.y / 2.0 - curry; // shortest distance in y from the +wall
    //         if (curry + membraneObject.waterBox.y / 2.0 < negwall)
    //             negwall = curry + membraneObject.waterBox.y / 2.0; // shortest distance in y from the -wall

    //         ychg = curry - membraneObject.waterBox.y / 2.0;
    //         if (ychg > 0) {
    //             outsidePos = true;
    //             outside = true;
    //             if (-ychg < ytot)
    //                 ytot = -ychg;
    //         } else if (ychg < -membraneObject.waterBox.y) {
    //             outside = true;
    //             outsideNeg = true;
    //             if (-(ychg + membraneObject.waterBox.y) > ytot)
    //                 ytot = -(ychg + membraneObject.waterBox.y);
    //         }
    //         /*measure each interface to y plane*/
    //         for (auto& iface : moleculeList[mp].interfaceList) {
    //             Vector ifaceVec { iface.coord - targCom.comCoord };

    //             /*only need y component*/
    //             dyrot = row[0] * ifaceVec.x + row[1] * ifaceVec.y + row[2] * ifaceVec.z;
    //             curry = targCom.comCoord.y + targCom.trajTrans.y + dyrot;
    //             if (membraneObject.waterBox.y / 2.0 - curry < poswall)
    //                 poswall = membraneObject.waterBox.y / 2.0 - curry; // shortest distance in y from the +wall
    //             if (curry + membraneObject.waterBox.y / 2.0 < negwall)
    //                 negwall = curry + membraneObject.waterBox.y / 2.0; // shortest distance in y from the -wall

    //             ychg = curry - membraneObject.waterBox.y / 2.0;
    //             if (ychg > 0) {
    //                 outsidePos = true;
    //                 outside = true;
    //                 if (-ychg < ytot)
    //                     ytot = -ychg;
    //             } else if (ychg < -membraneObject.waterBox.y) {
    //                 outside = true;
    //                 outsideNeg = true;
    //                 if (-(ychg + membraneObject.waterBox.y) > ytot)
    //                     ytot = -(ychg + membraneObject.waterBox.y);
    //             }
    //         }
    //     }
    //     if (outside) {
    //         /*Put back inside the box*/
    //         targCom.trajTrans.y += 2.0 * ytot;
    //         if (outsideNeg && outsidePos) {
    //             /*For a large complex, test if it could be pushed back out the other side*/
    //             std::cout << "ytot: " << ytot << " but also extends in other direction. negwall " << negwall
    //                       << " poswall: " << poswall << '\n';
    //             recheck = true;
    //         } else if (outsideNeg && poswall < (2.0 * ytot)) {
    //             /*Also need to check that update will not push you out the other side*/
    //             // ytot is positive.
    //             recheck = true;
    //             std::cout << " Will push out the other side. negwall " << negwall << " poswall: " << poswall
    //                       << " ytot: " << ytot << '\n';
    //         } else if (outsidePos && negwall < (-2.0 * ytot)) {
    //             // outsidePos>0, which means ytot is negative.
    //             recheck = true;
    //             std::cout << " Will push out the other side. negwall " << negwall << " poswall: " << poswall
    //                       << " ytot: " << ytot << '\n';
    //         }
    //     }
    // }

    // if (canBeOutsideZ > 0) {
    //     bool outside { false };
    //     bool outsideNeg { false };
    //     bool outsidePos { false };

    //     double poswall { membraneObject.waterBox.z };
    //     double negwall { membraneObject.waterBox.z };

    //     std::array<double, 3> row {};
    //     row[0] = M[6];
    //     row[1] = M[7];
    //     row[2] = M[8];

    //     /*these need to be what current positions
    //       due to translation and rotation are*/
    //     // for (i = 0; i < s1; i++) {
    //     for (auto& memMol : targCom.memberList) {
    //         /*measure each protein COM to z plane*/
    //         Vector comVec { moleculeList[memMol].comCoord - targCom.comCoord };

    //         double dzrot = row[0] * comVec.x + row[1] * comVec.y + row[2] * comVec.z;
    //         currz = targCom.comCoord.z + targCom.trajTrans.z + dzrot;
    //         if (membraneObject.waterBox.z / 2.0 - currz < poswall)
    //             poswall = membraneObject.waterBox.z / 2.0 - currz; // shortest distance in z from the +wall
    //         if (currz + membraneObject.waterBox.z / 2.0 - RS3D < negwall)
    //             negwall = currz + membraneObject.waterBox.z / 2.0 - RS3D; // shortest distance in z from the -wall

    //         zchg = currz - membraneObject.waterBox.z / 2.0;
    //         if (zchg > 0) {
    //             outside = true;
    //             outsidePos = true;
    //             if (-zchg < ztot)
    //                 ztot = -zchg;
    //         } else if (zchg < -(membraneObject.waterBox.z - RS3D)) {
    //             outside = true;
    //             outsideNeg = true;
    //             if (-(zchg + membraneObject.waterBox.z - RS3D) > ztot)
    //                 ztot = -(zchg + membraneObject.waterBox.z - RS3D);
    //         }
    //         /*measure each interface to z plane*/
    //         for (auto& iface : moleculeList[memMol].interfaceList) {
    //             Vector ifaceVec { iface.coord - targCom.comCoord };

    //             dzrot = row[0] * ifaceVec.x + row[1] * ifaceVec.y + row[2] * ifaceVec.z;
    //             currz = targCom.comCoord.z + targCom.trajTrans.z + dzrot;

    //             if (membraneObject.waterBox.z / 2.0 - currz < poswall)
    //                 poswall = membraneObject.waterBox.z / 2.0 - currz; // shortest distance in z from the +wall
    //             if (currz + membraneObject.waterBox.z / 2.0 - RS3D < negwall)
    //                 negwall = currz + membraneObject.waterBox.z / 2.0 - RS3D; // shortest distance in z from the -wall

    //             zchg = currz - membraneObject.waterBox.z / 2.0;
    //             if (zchg > 0) {
    //                 outside = true;
    //                 outsidePos = true;
    //                 if (-zchg < ztot)
    //                     ztot = -zchg;
    //             } else if (zchg < -(membraneObject.waterBox.z - RS3D)) {
    //                 outside = true;
    //                 outsideNeg = true;
    //                 if (-(zchg + membraneObject.waterBox.z - RS3D) > ztot)
    //                     ztot = -(zchg + membraneObject.waterBox.z - RS3D);
    //             }
    //         }
    //     }
    //     if (outside) {
    //         /*Put back inside the box*/
    //         if (targCom.D.z == 0)
    //             std::cout << "Putting complex " << targCom.index << " with ztot " << ztot
    //                       << ", zchg " << (targCom.comCoord.z + targCom.trajTrans.z) - membraneObject.waterBox.z / 2.0
    //                       << ", currz " << targCom.comCoord.z + targCom.trajTrans.z << " and trajTrans.z "
    //                       << targCom.trajTrans.z << " back in box." << '\n';
    //         targCom.trajTrans.z += 2.0 * ztot;
    //         if (outsideNeg && outsidePos) {
    //             /*For a large complex, test if it could be pushed back out the other side*/
    //             std::cout << "ztot: " << ztot << " but also extends in other direction. negwall " << negwall
    //                       << " poswall: " << poswall << '\n';
    //             recheck = true;
    //         } else if (outsideNeg && poswall < (2.0 * ztot)) {
    //             /*Also need to check that update will not push zou out the other side*/
    //             // ztot is positive.
    //             recheck = true;
    //             std::cout << " Will push out the other side. negwall " << negwall << " poswall: " << poswall
    //                       << " ztot: " << ztot << '\n';
    //         } else if (outsidePos && negwall < (-2.0 * ztot)) {
    //             // outsidePos>0, which means ztot is negative.
    //             recheck = true;
    //             std::cout << " Will push out the other side. negwall " << negwall << " poswall: " << poswall
    //                       << " ztot: " << ztot << '\n';
    //         }
    //     } // updating traj to reflect
    // }
    // if (recheck) {
    //     std::cout << "rechecking span " << '\n';
    //     /*Test that new coordinates have not pushed you out of the box for a very large complex, if so, resample
    //      * rotation matrix.*/
    //     std::cout << "RECHECK THAT COMPLEX DOES NOT SPAN BOX IN SUBROUTINE REFLECT TRAJ COMPLEX RAD ROT. Complex: "
    //               << targCom.index << " size:" << targCom.memberList.size() << '\n';
    //     reflect_traj_check_span(xtot, ytot, ztot, params, targCom, moleculeList, M, membraneObject, RS3Dinput);
    // }
}
