/*! \file reflect_traj_check_span_sphere.cpp
 * ### Created on 02/27/2020 by Yiben Fu
 * ### Purpose: to check whether the complex.trajTrans and complex.trajRot will make the complex outside the sphere boundary
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

void reflect_traj_check_span_sphere(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, const Membrane& membraneObject, double RS3Dinput)
{
    // TRACE();
    bool needsRecheck { true };
    int maxItr { 50 };
    int checkItr { 0 };

    double RS3D;
    if (targCom.OnSurface) {
        RS3D = 0;
    } else {
        RS3D = RS3Dinput;
    }
    double sphereR = membraneObject.sphereR - RS3D;

    std::array<double, 9> M;
    M = create_euler_rotation_matrix(targCom.trajRot);

    if (targCom.D.z < 1E-14 || targCom.OnSurface) { // for the complex on the sphere surface
        // in this case, the movement only involves theta and phi, and R doesn't change,
        // so it won't make the complex outside the sphere.

        needsRecheck = false;

    } else { // for the complex inside the sphere

        while (checkItr < maxItr && needsRecheck) {
            needsRecheck = false;

            M = create_euler_rotation_matrix(targCom.trajRot);
            // find the farthest point
            Coord targcrds { 0, 0, sphereR }; // to store the furthest point' coords.
            double rtmp = sphereR; // to store the furthest point' distance to the sphere center.
            for (auto& memMol : targCom.memberList) {
                //measure each protein COM
                Vector comVec { moleculeList[memMol].comCoord - targCom.comCoord };
                Vector rotComVec { matrix_rotate(comVec, M) };

                Coord curr { targCom.comCoord + targCom.trajTrans + rotComVec };
                double currR = curr.get_magnitude();
                if (currR > rtmp) {
                    targcrds = curr;
                    rtmp = currR;
                }

                // measure each interface
                for (auto& iface : moleculeList[memMol].interfaceList) {
                    Vector ifaceVec { iface.coord - targCom.comCoord };
                    Vector rotIfaceVec { matrix_rotate(ifaceVec, M) };
                    curr = targCom.comCoord + targCom.trajTrans + rotIfaceVec;
                    currR = curr.get_magnitude();
                    if (currR > rtmp) {
                        targcrds = curr;
                        rtmp = currR;
                    }
                }
            }
            // check whether this complex is out of the box, if so, change trajTrans by considering the reflection
            if (rtmp > sphereR + 1E-15) {
                double lamda = -2.0 * (rtmp - sphereR) / rtmp;
                Coord dtrans = lamda * targcrds;
                targCom.trajTrans += dtrans;
                // check whether the reflection make the complex inside the sphere
                targcrds = Coord(0, 0, sphereR); // to store the furthest point' coords.
                rtmp = sphereR; // to store the furthest point' distance to the sphere center.
                for (auto& memMol : targCom.memberList) {
                    //measure each protein COM
                    Vector comVec { moleculeList[memMol].comCoord - targCom.comCoord };
                    Vector rotComVec { matrix_rotate(comVec, M) };

                    Coord curr { targCom.comCoord + targCom.trajTrans + rotComVec };
                    double currR = curr.get_magnitude();
                    if (currR > rtmp) {
                        targcrds = curr;
                        rtmp = currR;
                    }

                    // measure each interface
                    for (auto& iface : moleculeList[memMol].interfaceList) {
                        Vector ifaceVec { iface.coord - targCom.comCoord };
                        Vector rotIfaceVec { matrix_rotate(ifaceVec, M) };
                        curr = targCom.comCoord + targCom.trajTrans + rotIfaceVec;
                        currR = curr.get_magnitude();
                        if (currR > rtmp) {
                            targcrds = curr;
                            rtmp = currR;
                        }
                    }
                }
            }
            // recheck whether this complex is still out sphere, if so, regenerate trajTrans
            if (rtmp > sphereR + 1E-15) {
                targCom.trajTrans.x = sqrt(2.0 * params.timeStep * targCom.D.x) * GaussV();
                targCom.trajTrans.y = sqrt(2.0 * params.timeStep * targCom.D.y) * GaussV();
                targCom.trajTrans.z = sqrt(2.0 * params.timeStep * targCom.D.z) * GaussV();
                targCom.trajRot.x = sqrt(2.0 * params.timeStep * targCom.Dr.x) * GaussV();
                targCom.trajRot.y = sqrt(2.0 * params.timeStep * targCom.Dr.y) * GaussV();
                targCom.trajRot.z = sqrt(2.0 * params.timeStep * targCom.Dr.z) * GaussV();

                reflect_traj_complex_rad_rot_nocheck(params, targCom, moleculeList, membraneObject, RS3Dinput);
                ++checkItr;
                needsRecheck = true; // will need to recheck after resampling traj and trajR
            }
        } // end of while-loop
    }

    //    std::cout << "ITERATIONS TO CONVERGE POSITION WITHIN SPHERE: " << checkItr
    //        << " flag at end: true=success: " << !needsRecheck << '\n';
    if (needsRecheck) {
        // std::cout << "WARNING: DID NOT CONVERGE POSITION, NEW POS: " << '\n';
        for (auto memMol : targCom.memberList) {
            Vector comVec { moleculeList[memMol].comCoord - targCom.comCoord };
            Vector rotComVec { matrix_rotate(comVec, M) };

            // first would make xcom=targCom.comCoord.x+vr, then would also add dx
            // std::cout << "i: " << checkItr << " P: " << memMol
            //           << " com:" << targCom.comCoord.x + targCom.trajTrans.x + rotComVec.x << ' '
            //           << targCom.comCoord.y + targCom.trajTrans.y + rotComVec.y << ' '
            //           << targCom.comCoord.z + targCom.trajTrans.z + rotComVec.z << '\n';

            // update interface coords
            for (const auto& iface : moleculeList[memMol].interfaceList) {
                Vector ifaceVec { iface.coord - targCom.comCoord };
                Vector rotIfaceVec { matrix_rotate(ifaceVec, M) };

                /*first would make xcom=targCom.comCoord.x+vr, then would also add dx */
                // std::cout << targCom.comCoord.x + targCom.trajTrans.x + rotIfaceVec.x << ' '
                //           << targCom.comCoord.y + targCom.trajTrans.y + rotIfaceVec.y << ' '
                //           << targCom.comCoord.z + targCom.trajTrans.z + rotIfaceVec.z << '\n';
            }
        }
    } // DID NOT CONVERGE
}
