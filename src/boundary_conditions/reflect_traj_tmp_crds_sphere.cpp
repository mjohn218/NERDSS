#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"
#include "tracing.hpp"

/*evaluates reflection out of box based on tmpCoords of all proteins. targCOM tmpCOM coords also must be consistent-updated.
 *Does not update molecule traj vectors, just updates the passed in vector traj, to collect for both complexes. 
 *
*/
void reflect_traj_tmp_crds_sphere(
    const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, std::array<double, 3>& traj, const Membrane& membraneObject, double RS3Dinput)
{
    // TRACE();
    /*This routine updated April 2020 to test if a large complex that spans the sphere could extend out in both directions
    if so, it attempts to correct for this by resampling the complex's translational updates.
    */
    // if (targCom.D.z < 1E-15) // for the complex on the sphere, no need to check
    //    return;

    double RS3D;
    if (targCom.OnSurface || targCom.tmpOnSurface) {
        RS3D = 0;
    } else {
        RS3D = RS3Dinput;
    }
    double sphereR = membraneObject.sphereR - RS3D;
    Coord curr;
    curr.x = targCom.tmpComCoord.x + traj[0];
    curr.y = targCom.tmpComCoord.y + traj[1];
    curr.z = targCom.tmpComCoord.z + traj[2];

    std::array<double, 9> M {};
    for (int mm = 0; mm < 9; mm++)
        M[mm] = 0;
    M[0] = 1;
    M[4] = 1;
    M[8] = 1; //set M to identity, perform no rotations here.

    /*This is to test based on general size if it is close to boundaries, before doing detailed evaluation below.*/

    bool canBeOutsphere { false };
    if ((curr.get_magnitude() + targCom.radius) > sphereR)
        canBeOutsphere = true;

    /*Now evaluate all interfaces distance from boundaries.*/
    bool recheck { false };
    double dr = 0.0;
    Coord targcrds;
    if (canBeOutsphere == true) {

        bool outside { false };

        /*these need to be what current positions
        due to translation and rotation are*/
        for (auto& memMol : targCom.memberList) {
            // measure each protein COM to z plane
            Vector comVec { moleculeList[memMol].tmpComCoord - targCom.tmpComCoord };
            double dxrot { M[0] * comVec.x + M[1] * comVec.y + M[2] * comVec.z };
            double dyrot { M[3] * comVec.x + M[4] * comVec.y + M[5] * comVec.z };
            double dzrot { M[6] * comVec.x + M[7] * comVec.y + M[8] * comVec.z };
            curr.x = targCom.tmpComCoord.x + traj[0] + dxrot;
            curr.y = targCom.tmpComCoord.y + traj[1] + dyrot;
            curr.z = targCom.tmpComCoord.z + traj[2] + dzrot;
            double drtmp = curr.get_magnitude() - sphereR;
            if (drtmp > dr) {
                outside = true;
                dr = drtmp;
                targcrds = curr;
            }

            // measure each interface to x plane
            for (int ii = 0; ii < moleculeList[memMol].interfaceList.size(); ii++) {
                Vector ifaceVec { moleculeList[memMol].tmpICoords[ii] - targCom.tmpComCoord };
                dxrot = M[0] * ifaceVec.x + M[1] * ifaceVec.y + M[2] * ifaceVec.z;
                dyrot = M[3] * ifaceVec.x + M[4] * ifaceVec.y + M[5] * ifaceVec.z;
                dzrot = M[6] * ifaceVec.x + M[7] * ifaceVec.y + M[8] * ifaceVec.z;
                curr.x = targCom.tmpComCoord.x + traj[0] + dxrot;
                curr.y = targCom.tmpComCoord.y + traj[1] + dyrot;
                curr.z = targCom.tmpComCoord.z + traj[2] + dzrot;
                drtmp = curr.get_magnitude() - sphereR;
                if (drtmp > dr) {
                    outside = true;
                    dr = drtmp;
                    targcrds = curr;
                }
            }
        }

        if (outside) {
            // Put back inside the box
            double lamda = -2.0 * (targcrds.get_magnitude() - sphereR) / targcrds.get_magnitude();
            traj[0] = lamda * targcrds.x;
            traj[1] = lamda * targcrds.y;
            traj[2] = lamda * targcrds.z;
        }
    }
}