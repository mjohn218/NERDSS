/*! \file reflect_complex_rad_rot_sphere.cpp
 * ### Created on 02/25/2020 by Yiben Fu
 * ### Purpose: only works for complex already bound on the spherical surface
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */
#include "boundary_conditions/reflect_functions.hpp"
#include "tracing.hpp"

void reflect_complex_rad_rot_sphere(const Membrane& membraneObject, Complex& targCom, std::vector<Molecule>& moleculeList, double RS3Dinput)
{
    // TRACE();
    // only works for the complex after association or diffusion, spherical system
    // For the sphere system, many times of reflections may need to move the complex back inside the sphere!!

    // declare the boundary
    double sphereR;
    if (targCom.D.z < 1E-8) {
        sphereR = membraneObject.sphereR;
    } else {
        sphereR = membraneObject.sphereR - RS3Dinput;
    }

    Coord curr = targCom.comCoord;
    double currR = curr.get_magnitude();

    bool canBeOutsideX { false };
    bool outside { false };

    if ((currR + targCom.radius) > sphereR)
        canBeOutsideX = true;

    if (canBeOutsideX == true) {
        // find the farthest point
        Coord targcrds { 0, 0, sphereR }; // to store the furthest point' coords.
        double rtmp = sphereR; // to store the furthest point' distance to the sphere center.
        for (auto& memMol : targCom.memberList) {
            //measure each protein COM
            curr = moleculeList[memMol].comCoord;
            currR = curr.get_magnitude();
            if (currR > sphereR && currR > rtmp) {
                targcrds = curr;
                rtmp = currR;
                outside = true;
            }

            // measure each interface
            for (auto& iface : moleculeList[memMol].interfaceList) {
                curr = iface.coord;
                currR = curr.get_magnitude();
                if (currR > sphereR && currR > rtmp) {
                    targcrds = curr;
                    rtmp = currR;
                    outside = true;
                }
            }
        }

        int times = 0; // to count the loop-times of 'while'
        while (outside == true) {
            times++;
            rtmp = targcrds.get_magnitude();
            double lamda = -2.0 * (rtmp - sphereR) / rtmp;
            Coord dtrans = lamda * targcrds;
            targCom.comCoord += dtrans;
            for (auto memMol : targCom.memberList) {
                moleculeList[memMol].comCoord += dtrans;
                for (auto& iface : moleculeList[memMol].interfaceList)
                    iface.coord += dtrans;
            }
            // reflecting may make the complex outside the sphere in other direction,
            // thus we need to recheck whether outside
            // initialize outside, targcrds;
            outside = false;
            targcrds = Coord(0, 0, sphereR); // to store the furthest point' coords.
            rtmp = sphereR; // to store the furthest point' distance to the sphere center.
            for (auto& memMol : targCom.memberList) {
                //measure each protein COM
                curr = moleculeList[memMol].comCoord;
                currR = curr.get_magnitude();
                if (currR > sphereR && currR > rtmp) {
                    targcrds = curr;
                    rtmp = currR;
                    outside = true;
                }

                // measure each interface
                for (auto& iface : moleculeList[memMol].interfaceList) {
                    curr = iface.coord;
                    currR = curr.get_magnitude();
                    if (currR > sphereR && currR > rtmp) {
                        targcrds = curr;
                        rtmp = currR;
                        outside = true;
                    }
                }
            }

            if (times > 100) {
                // so many times reflection still cannot make the complex back inside the sphere, thus we may need report 'WRONG!!'
                std::cout << "ALREADY UPDATED POSITIONS.BUT, IN REFLECT COMPLEX_RAD_ROT_SPHERE, 100 TIMES REFLECTIONS CAN'T MOVE THE COMPLEX BACK INSIDE SPHERE." << '\n';
                std::cout << "COMPLEX " << targCom.index << ", COMPLEX COM " << targCom.comCoord.x << ", " << targCom.comCoord.y << ", " << targCom.comCoord.z << '\n';
                std::cout << "COMPLEX RADIUS " << targCom.radius << ", COMPLEX SIZE " << targCom.memberList.size() << '\n';
                std::cout << "EXITING..." << '\n';
                exit(1);
            }
        } // end of while-loop
    } // update reflection

    // ALSO, the position update may cause lipids off sphere, so adjust lipids back onto surface,
    int molnumber = -1;
    double dr = 0.0;
    Coord dtrans;
    for (auto& memMol : targCom.memberList) {
        if (moleculeList[memMol].isLipid == true) {
            Coord targ = moleculeList[memMol].comCoord;
            double rtmp = targ.get_magnitude();
            double drtmp = std::abs(targ.get_magnitude() - membraneObject.sphereR);
            if (drtmp > 1E-4 && drtmp > dr) {
                dr = drtmp;
                molnumber = memMol;
                dtrans = (membraneObject.sphereR / targ.get_magnitude()) * targ - targ;
            }
        }
    }
    if (molnumber != -1) {
        for (auto& mol : targCom.memberList) {
            moleculeList[mol].comCoord += dtrans;
            for (auto& iface : moleculeList[mol].interfaceList) {
                iface.coord += dtrans;
            }
        }
    }
}