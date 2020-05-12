#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "tracing.hpp"

void reflect_traj_complex_rad_rot_nocheck_sphere(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, const Membrane& membraneObject, double RS3Dinput)
{
    // TRACE();
    double RS3D;
    if (targCom.OnSurface) {
        RS3D = 0;
    } else {
        RS3D = RS3Dinput;
    }

    double sphereR = membraneObject.sphereR - RS3D;

    std::array<double, 9> M = create_euler_rotation_matrix(targCom.trajRot);

    if (targCom.D.z < 1E-14 || targCom.OnSurface) { // for the complex on the sphere surface
        // in this case, the movement only involves theta and phi, and R doesn't change,
        // so it won't make the complex outside the sphere.

    } else { // for the complex inside the sphere

        Coord curr { targCom.comCoord + targCom.trajTrans };

        if (curr.get_magnitude() + targCom.radius <= sphereR)
            return;

        // for the outside sphere situation:
        Coord targcrds { 0, 0, sphereR }; // to store the furthest point' coords.
        double rtmp = sphereR; // to store the furthest point' distance to the sphere center.
        for (auto& memMol : targCom.memberList) {
            //measure each protein COM
            Vector comVec { moleculeList[memMol].comCoord - targCom.comCoord };
            Vector rotComVec { matrix_rotate(comVec, M) };

            curr = targCom.comCoord + targCom.trajTrans + rotComVec;
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
        double lamda = -2.0 * (rtmp - sphereR) / rtmp;
        Coord dtrans = lamda * targcrds;
        targCom.trajTrans += dtrans;
    }
}