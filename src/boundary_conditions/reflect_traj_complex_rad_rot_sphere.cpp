/*! \file reflect_traj_complex_rad_rot_sphere.cpp
 * ### Created on 02/25/2020 by Yiben Fu
 * ### Purpose:  works for complex inside a sphere 
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
#include "tracing.hpp"

void reflect_traj_complex_rad_rot_sphere(const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, const Membrane& membraneObject, double RS3Dinput)
{
    // TRACE();

    // for implicit-lipid model, the boundary surface must consider the reflecting-surface RS3D
    // for explicit-lipid model, membraneObject.RS3D = 0. And, for those proteins that are bound on surface, they are not allowed to reflect along Z-axis
    // for the complex on the surface, its propagation is different from those inside the sphere
    double RS3D;
    if (targCom.OnSurface) {
        RS3D = 0;
    } else {
        RS3D = RS3Dinput;
    }
    double sphereR = membraneObject.sphereR - RS3D;

    std::array<double, 9> M;
    M = create_euler_rotation_matrix(targCom.trajRot);

    bool recheck = false;
    /*calculate distance to the center of the sphere. */
    /*first just test Complex COM+ radius, if it fits inside membraneObject.sphereR-RS3D.
      if yes, then test if all interfaces inside sphereR.
      then, if they are outside, what is max displacement beyond the sphereR.
      move it along the radial direction inside by the displacement.
      Also, check if it fits inside the sphere.
     */

    if (targCom.D.z < 1E-14 || targCom.OnSurface) { // for the complex on the sphere surface

    } else { // for the complex inside the sphere
        Coord curr;
        curr.x = targCom.comCoord.x + targCom.trajTrans.x;
        curr.y = targCom.comCoord.y + targCom.trajTrans.y;
        curr.z = targCom.comCoord.z + targCom.trajTrans.z;
        /*assume the origin of the sphere is at zero. */
        double rtmp = curr.get_magnitude();
        if (rtmp + targCom.radius > sphereR) {
            /*Now evaluate all molecules and interfaces distance from boundaries.*/
            bool outside = false;
            Coord targcrds;
            double targR = sphereR;
            for (auto& memMol : targCom.memberList) {
                /*measure each protein COM to origin*/
                Vector comVec { moleculeList[memMol].comCoord - targCom.comCoord };
                double dxrot { M[0] * comVec.x + M[1] * comVec.y + M[2] * comVec.z };
                double dyrot { M[3] * comVec.x + M[4] * comVec.y + M[5] * comVec.z };
                double dzrot { M[6] * comVec.x + M[7] * comVec.y + M[8] * comVec.z };
                Vector rot { dxrot, dyrot, dzrot };
                curr = Coord(targCom.comCoord + targCom.trajTrans + rot);
                rtmp = curr.get_magnitude();
                if (rtmp > targR) {
                    outside = true;
                    targR = rtmp;
                    targcrds = curr;
                }
                /*measure each interface to z plane*/
                for (auto& iface : moleculeList[memMol].interfaceList) {
                    Vector ifaceVec { iface.coord - targCom.comCoord };
                    dxrot = M[0] * ifaceVec.x + M[1] * ifaceVec.y + M[2] * ifaceVec.z;
                    dyrot = M[3] * ifaceVec.x + M[4] * ifaceVec.y + M[5] * ifaceVec.z;
                    dzrot = M[6] * ifaceVec.x + M[7] * ifaceVec.y + M[8] * ifaceVec.z;
                    rot = Vector(dxrot, dyrot, dzrot);
                    curr = Coord(targCom.comCoord + targCom.trajTrans + rot);
                    rtmp = curr.get_magnitude();
                    if (rtmp > targR) {
                        outside = true;
                        targR = rtmp;
                        targcrds = curr;
                    }
                }
            }
            if (outside == true) {
                recheck = true;
                double lamda = -2.0 * (targR - sphereR) / targR;
                targCom.trajTrans = Vector(targCom.trajTrans + lamda * targcrds);
            }
        }
    }

    if (recheck) {
        //std::cout << "rechecking span sphere " << '\n';
        //Test that new coordinates have not pushed you out of the box for a very large complex, if so, resample  rotation matrix.
        //std::cout << "RECHECK THAT COMPLEX DOES NOT SPAN SPHERE IN SUBROUTINE REFLECT_TRAJ_COMPLEX_RAD_ROT. Complex: "
        //          << targCom.index << " size:" << targCom.memberList.size() << '\n';
        reflect_traj_check_span(params, targCom, moleculeList, membraneObject, RS3Dinput);
    }
}