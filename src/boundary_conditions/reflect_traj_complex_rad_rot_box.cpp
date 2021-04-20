#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "tracing.hpp"

void reflect_traj_complex_rad_rot_box(const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, const Membrane& membraneObject, double RS3Dinput)
{
    // TRACE();
    // NOTE: it only works for a box system with the membrane surface located on the Z-bottom.

    // for implicit-lipid model, the boundary surface must consider the reflecting-surface RS3D
    // for explicit-lipid model, membraneObject.RS3D = 0. And, for those proteins that are bound on surface, they are not allowed to reflect along Z-axis
    double RS3D;
    if (targCom.OnSurface || targCom.D.z < 1E-8) {
        RS3D = 0;
    } else {
        RS3D = RS3Dinput;
    }

    std::array<double, 9> M;
    M = create_euler_rotation_matrix(targCom.trajRot);

    // declare the six boundary sides of the system box;
    double PosSideX = membraneObject.waterBox.x / 2.0;
    double NegSideX = -membraneObject.waterBox.x / 2.0;
    double PosSideY = membraneObject.waterBox.y / 2.0;
    double NegSideY = -membraneObject.waterBox.y / 2.0;
    double PosSideZ = membraneObject.waterBox.z / 2.0;
    double NegSideZ = -membraneObject.waterBox.z / 2.0 + RS3D;

    // his routine updated March 2017 to test if a large complex that spans the box could extend out in both directions
    // if so, it attempts to correct for this by resampling the complex's translational and rotational updates.
    double xtot { 0 };
    double ytot { 0 };
    double ztot { 0 };
    double currx { targCom.comCoord.x + targCom.trajTrans.x };
    double curry { targCom.comCoord.y + targCom.trajTrans.y };
    double currz { targCom.comCoord.z + targCom.trajTrans.z };

    // This is to test based on general size if it is close to boundaries, before doing detailed evaluation below.

    bool canBeOutsideX { false };
    //double xchg = currx - membraneObject.waterBox.x / 2.0;
    if ((currx + targCom.radius) > PosSideX || (currx - targCom.radius) < NegSideX)
        canBeOutsideX = true;

    bool canBeOutsideY { false };
    //double ychg = curry - membraneObject.waterBox.y / 2.0;
    if ((curry + targCom.radius) > PosSideY || (curry - targCom.radius) < NegSideY)
        canBeOutsideY = true;

    bool canBeOutsideZ { false };
    //double zchg = currz - membraneObject.waterBox.z / 2.0;
    if ((currz + targCom.radius) > PosSideZ || (currz - targCom.radius) < NegSideZ)
        canBeOutsideZ = true;

    // Now evaluate all interfaces distance from boundaries.
    bool recheck { false };
    /////////////////////////////////////////////////////////////////////////////////////////////////
    if (canBeOutsideX) {
        bool outside { false };
        bool outsideNeg { false };
        bool outsidePos { false };

        double Posdx = 0; // to store the largest distance outside positive X side. mark +
        double Negdx = 0; // to store the largest distance outside negative X side. make -
        double PosWallX = NegSideX; // to store the farthest position in positive X direction.
        double NegWallX = PosSideX; // to store the farthest position in negative X direction.

        std::array<double, 3> row {};
        row[0] = M[0];
        row[1] = M[1];
        row[2] = M[2];

        //these need to be what current positions due to translation and rotation are
        for (auto& memMol : targCom.memberList) {
            //measure each protein COM to x sides
            Vector comVec { moleculeList[memMol].comCoord - targCom.comCoord };
            double dxrot { row[0] * comVec.x + row[1] * comVec.y + row[2] * comVec.z };
            currx = targCom.comCoord.x + targCom.trajTrans.x + dxrot;

            if (currx > PosSideX) {
                outsidePos = true;
                outside = true;
                if (currx - PosSideX > Posdx)
                    Posdx = currx - PosSideX; // outside x, longest distance to +wall
            }
            if (currx < NegSideX) {
                outsideNeg = true;
                outside = true;
                if (currx - NegSideX < Negdx)
                    Negdx = currx - NegSideX; // outside x, longest distance to -wall
            }
            if (currx > PosWallX)
                PosWallX = currx;
            if (currx < NegWallX)
                NegWallX = currx;

            // measure each interface to x sides
            for (auto& iface : moleculeList[memMol].interfaceList) {
                Vector ifaceVec { iface.coord - targCom.comCoord };
                dxrot = row[0] * ifaceVec.x + row[1] * ifaceVec.y + row[2] * ifaceVec.z;
                currx = targCom.comCoord.x + targCom.trajTrans.x + dxrot;

                if (currx > PosSideX) {
                    outsidePos = true;
                    outside = true;
                    if (currx - PosSideX > Posdx)
                        Posdx = currx - PosSideX; // outside x, longest distance to +wall
                }
                if (currx < NegSideX) {
                    outsideNeg = true;
                    outside = true;
                    if (currx - NegSideX < Negdx)
                        Negdx = currx - NegSideX; // outside x, longest distance to -wall
                }
                if (currx > PosWallX)
                    PosWallX = currx;
                if (currx < NegWallX)
                    NegWallX = currx;
            }
        }
        if (outside) {
            // Put back inside the box, extended out
            if (outsidePos) {
                targCom.trajTrans.x -= 2.0 * Posdx;
                xtot = -Posdx;
            } else if (outsideNeg) {
                targCom.trajTrans.x -= 2.0 * Negdx;
                xtot = -Negdx;
            }

            if (outsideNeg && outsidePos) {
                // For a large complex, test if it could be pushed back out the other side
                // std::cout << "Complex extends both X sides of the box. Out posX side " << Posdx
                //           << ", and out negX side " << Negdx << '\n';
                recheck = true;
            } else if (outsideNeg && PosWallX - 2.0 * Negdx > PosSideX) {
                // Also need to check that update will not push you out the other side
                recheck = true;
                // std::cout << " Will push out the other side. negwallx " << NegWallX << ", poswallx: " << PosWallX
                //           << ". 2*xtot: " << -2.0 * Negdx << '\n';
            } else if (outsidePos && NegWallX - 2.0 * Posdx < NegSideX) {
                recheck = true;
                // std::cout << " Will push out the other side. negwallx " << NegWallX << ", poswallx: " << PosWallX
                //           << ". 2*xtot: " << -2.0 * Posdx << '\n';
            }
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////
    if (canBeOutsideY) {
        bool outside { false };
        bool outsideNeg { false };
        bool outsidePos { false };

        double Posdy = 0; // to store the largest distance outside positive Y side. mark +
        double Negdy = 0; // to store the largest distance outside negative Y side. make -
        double PosWallY = NegSideY; // to store the farthest position in positive Y direction.
        double NegWallY = PosSideY; // to store the farthest position in negative Y direction.

        std::array<double, 3> row {};
        row[0] = M[3];
        row[1] = M[4];
        row[2] = M[5];

        //these need to be what current positions due to translation and rotation are
        for (auto& memMol : targCom.memberList) {
            //measure each protein COM to y sides
            Vector comVec { moleculeList[memMol].comCoord - targCom.comCoord };
            double dyrot { row[0] * comVec.x + row[1] * comVec.y + row[2] * comVec.z };
            curry = targCom.comCoord.y + targCom.trajTrans.y + dyrot;

            if (curry > PosSideY) {
                outsidePos = true;
                outside = true;
                if (curry - PosSideY > Posdy)
                    Posdy = curry - PosSideY; // outside Y, longest distance to +wall
            }
            if (curry < NegSideY) {
                outsideNeg = true;
                outside = true;
                if (curry - NegSideY < Negdy)
                    Negdy = curry - NegSideY; // outside Y, longest distance to -wall
            }
            if (curry > PosWallY)
                PosWallY = curry;
            if (curry < NegWallY)
                NegWallY = curry;

            // measure each interface to y sides
            for (auto& iface : moleculeList[memMol].interfaceList) {
                Vector ifaceVec { iface.coord - targCom.comCoord };
                dyrot = row[0] * ifaceVec.x + row[1] * ifaceVec.y + row[2] * ifaceVec.z;
                curry = targCom.comCoord.y + targCom.trajTrans.y + dyrot;

                if (curry > PosSideY) {
                    outsidePos = true;
                    outside = true;
                    if (curry - PosSideY > Posdy)
                        Posdy = curry - PosSideY; // outside Y, longest distance to +wall
                }
                if (curry < NegSideY) {
                    outsideNeg = true;
                    outside = true;
                    if (curry - NegSideY < Negdy)
                        Negdy = curry - NegSideY; // outside Y, longest distance to -wall
                }
                if (curry > PosWallY)
                    PosWallY = curry;
                if (curry < NegWallY)
                    NegWallY = curry;
            }
        }
        if (outside) {
            // Put back inside the box, extended out
            if (outsidePos) {
                targCom.trajTrans.y -= 2.0 * Posdy;
                ytot = -Posdy;
            } else if (outsideNeg) {
                targCom.trajTrans.y -= 2.0 * Negdy;
                ytot = -Negdy;
            }

            if (outsideNeg && outsidePos) {
                // For a large complex, test if it could be pushed back out the other side
                // std::cout << "Complex extends both Y sides of the box. Out posY side " << Posdy
                //           << ", and out negY side " << Negdy << '\n';
                recheck = true;
            } else if (outsideNeg && PosWallY - 2.0 * Negdy > PosSideY) {
                // Also need to check that update will not push you out the other side
                recheck = true;
                // std::cout << " Will push out the other side. negwally " << NegWallY << ", poswally: " << PosWallY
                //           << ". 2*ytot: " << -2.0 * Negdy << '\n';
            } else if (outsidePos && NegWallY - 2.0 * Posdy < NegSideY) {
                recheck = true;
                // std::cout << " Will push out the other side. negwally " << NegWallY << ", poswally: " << PosWallY
                //           << ". 2*ytot: " << -2.0 * Posdy << '\n';
            }
        }
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////
    if (canBeOutsideZ) {
        bool outside { false };
        bool outsideNeg { false };
        bool outsidePos { false };

        double Posdz = 0; // to store the largest distance outside positive Z side. mark +
        double Negdz = 0; // to store the largest distance outside negative Z side. make -
        double PosWallZ = NegSideZ; // to store the farthest position in positive Z direction.
        double NegWallZ = PosSideZ; // to store the farthest position in negative Z direction.

        std::array<double, 3> row {};
        row[0] = M[6];
        row[1] = M[7];
        row[2] = M[8];

        //these need to be what current positions due to translation and rotation are
        for (auto& memMol : targCom.memberList) {
            //measure each protein COM to z sides
            Vector comVec { moleculeList[memMol].comCoord - targCom.comCoord };
            double dzrot { row[0] * comVec.x + row[1] * comVec.y + row[2] * comVec.z };
            currz = targCom.comCoord.z + targCom.trajTrans.z + dzrot;

            if (currz > PosSideZ) {
                outsidePos = true;
                outside = true;
                if (currz - PosSideZ > Posdz)
                    Posdz = currz - PosSideZ; // outside Z, longest distance to +wall
            }
            if (currz < NegSideZ) {
                outsideNeg = true;
                outside = true;
                if (currz - NegSideZ < Negdz)
                    Negdz = currz - NegSideZ; // outside Z, longest distance to -wall
            }
            if (currz > PosWallZ)
                PosWallZ = currz;
            if (currz < NegWallZ)
                NegWallZ = currz;

            // measure each interface to z sides
            for (auto& iface : moleculeList[memMol].interfaceList) {
                Vector ifaceVec { iface.coord - targCom.comCoord };
                dzrot = row[0] * ifaceVec.x + row[1] * ifaceVec.y + row[2] * ifaceVec.z;
                currz = targCom.comCoord.z + targCom.trajTrans.z + dzrot;

                if (currz > PosSideZ) {
                    outsidePos = true;
                    outside = true;
                    if (currz - PosSideZ > Posdz)
                        Posdz = currz - PosSideZ; // outside Z, longest distance to +wall
                }
                if (currz < NegSideZ) {
                    outsideNeg = true;
                    outside = true;
                    if (currz - NegSideZ < Negdz)
                        Negdz = currz - NegSideZ; // outside Z, longest distance to -wall
                }
                if (currz > PosWallZ)
                    PosWallZ = currz;
                if (currz < NegWallZ)
                    NegWallZ = currz;
            }
        }
        if (outside) {
            // Put back inside the box, extended out
            if (outsidePos) {
                targCom.trajTrans.z -= 2.0 * Posdz;
                ztot = -Posdz;
            } else if (outsideNeg) {
                targCom.trajTrans.z -= 2.0 * Negdz;
                ztot = -Negdz;
            }

            if (outsideNeg && outsidePos) {
                // For a large complex, test if it could be pushed back out the other side
                // std::cout << "Complex extends both Z sides of the box. Out posZ side " << Posdz
                //           << ", and out negZ side " << Negdz << '\n';
                recheck = true;
            } else if (outsideNeg && PosWallZ - 2.0 * Negdz > PosSideZ) {
                // Also need to check that update will not push you out the other side
                recheck = true;
                // std::cout << " Will push out the other side. negwallz " << NegWallZ << ", poswallz: " << PosWallZ
                //           << ". 2*ztot: " << -2.0 * Negdz << '\n';
            } else if (outsidePos && NegWallZ - 2.0 * Posdz < NegSideZ) {
                recheck = true;
                // std::cout << " Will push out the other side. negwallz " << NegWallZ << ", poswallz: " << PosWallZ
                //           << ". 2*ztot: " << -2.0 * Posdz << '\n';
            }
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    if (recheck) {
        //std::cout << "rechecking span " << '\n';
        //Test that new coordinates have not pushed you out of the box for a very large complex, if so, resample  rotation matrix.
        //std::cout << "RECHECK THAT COMPLEX DOES NOT SPAN BOX IN SUBROUTINE REFLECT_TRAJ_COMPLEX_RAD_ROT. Complex: "
        //          << targCom.index << " size:" << targCom.memberList.size() << '\n';
        reflect_traj_check_span_box(params, targCom, moleculeList, membraneObject, RS3Dinput);
    }
}