#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "tracing.hpp"

void reflect_traj_check_span_box(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, const Membrane& membraneObject, double RS3Dinput)
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

    std::array<double, 9> M;
    M = create_euler_rotation_matrix(targCom.trajRot);

    // declare the six boundary sides of the system box;
    double posX = membraneObject.waterBox.x / 2.0;
    double negX = -membraneObject.waterBox.x / 2.0;
    double posY = membraneObject.waterBox.y / 2.0;
    double negY = -membraneObject.waterBox.y / 2.0;
    double posZ = membraneObject.waterBox.z / 2.0;
    double negZ = -membraneObject.waterBox.z / 2.0 + RS3D;

    // double tol = 1E-14;
    while (checkItr < maxItr && needsRecheck) {
        needsRecheck = false; // without double span, this will stay 0

        bool moveFailed { false };
        bool outsideBox { false };
        bool outsidePosX { false };
        bool outsideNegX { false };
        bool outsidePosY { false };
        bool outsideNegY { false };
        bool outsidePosZ { false };
        bool outsideNegZ { false };

        double posWallX { negX };
        double negWallX { posX };
        double posWallY { negY };
        double negWallY { posY };
        double posWallZ { negZ };
        double negWallZ { posZ };

        // these need to be what current positions
        // due to translation and rotation are
        for (auto& memMol : targCom.memberList) {
            // measure protein COM to plane
            Vector comVec { moleculeList[memMol].comCoord - targCom.comCoord };
            Vector rotComVec { matrix_rotate(comVec, M) };

            double currX { targCom.comCoord.x + targCom.trajTrans.x + rotComVec.x };
            double currY { targCom.comCoord.y + targCom.trajTrans.y + rotComVec.y };
            double currZ { targCom.comCoord.z + targCom.trajTrans.z + rotComVec.z };

            if (currX > posWallX)
                posWallX = currX; // farthest point in +X
            if (currX < negWallX)
                negWallX = currX; // farthest point in -X

            if (currY > posWallY)
                posWallY = currY; // farthest point in +Y
            if (currY < negWallY)
                negWallY = currY; // farthest point in -Y

            if (currZ > posWallZ)
                posWallZ = currZ; // farthest point in +Z
            if (currZ < negWallZ)
                negWallZ = currZ; // farthest point in -Z

            // measure each interface
            for (auto& iface : moleculeList[memMol].interfaceList) {
                Vector ifaceVec { iface.coord - targCom.comCoord };
                Vector rotIfaceVec { matrix_rotate(ifaceVec, M) };

                currX = targCom.comCoord.x + targCom.trajTrans.x + rotIfaceVec.x;
                currY = targCom.comCoord.y + targCom.trajTrans.y + rotIfaceVec.y;
                currZ = targCom.comCoord.z + targCom.trajTrans.z + rotIfaceVec.z;

                if (currX > posWallX)
                    posWallX = currX;
                if (currX < negWallX)
                    negWallX = currX;

                if (currY > posWallY)
                    posWallY = currY;
                if (currY < negWallY)
                    negWallY = currY;

                if (currZ > posWallZ)
                    posWallZ = currZ;
                if (currZ < negWallZ)
                    negWallZ = currZ;
            } // loop over interfaces
        } // loop over proteins in complex.

        if (posWallX > posX) {
            outsidePosX = true;
            outsideBox = true;
        }
        if (negWallX < negX) {
            outsideNegX = true;
            outsideBox = true;
        }
        if (std::abs(posWallX - negWallX) > 1.0 / 2.0 * std::abs(posX - negX))
            maxItr = 20;
        if (std::abs(posWallX - negWallX) > 2.0 / 3.0 * std::abs(posX - negX))
            maxItr = 10;
        if (std::abs(posWallX - negWallX) > 4.0 / 5.0 * std::abs(posX - negX))
            maxItr = 5;

        if (posWallY > posY) {
            outsidePosY = true;
            outsideBox = true;
        }
        if (negWallY < negY) {
            outsideNegY = true;
            outsideBox = true;
        }
        if (std::abs(posWallY - negWallY) > 1.0 / 2.0 * std::abs(posY - negY))
            maxItr = 20;
        if (std::abs(posWallY - negWallY) > 2.0 / 3.0 * std::abs(posY - negY))
            maxItr = 10;
        if (std::abs(posWallY - negWallY) > 4.0 / 5.0 * std::abs(posY - negY))
            maxItr = 5;

        if (posWallZ > posZ) {
            outsidePosZ = true;
            outsideBox = true;
        }
        if (negWallZ < negZ) {
            outsideNegZ = true;
            outsideBox = true;
        }
        if (std::abs(posWallZ - negWallZ) > 1.0 / 2.0 * std::abs(posZ - negZ))
            maxItr = 20;
        if (std::abs(posWallZ - negWallZ) > 2.0 / 3.0 * std::abs(posZ - negZ))
            maxItr = 10;
        if (std::abs(posWallZ - negWallZ) > 4.0 / 5.0 * std::abs(posZ - negZ))
            maxItr = 5;

        if (outsideBox) {
            // Check the x-dimension
            if (outsideNegX && outsidePosX) {
                // For a large complex, test if checkItr could be pushed back out the other side
                // std::cout << " extends in both directions of X. negwall " << negWallX
                //           << " poswall: " << posWallX << '\n';
                // std::cout << "Pushed back out of box X, resample M and traj " << '\n';
                moveFailed = true;
            }
            if (outsideNegX && !outsidePosX) {
                targCom.trajTrans.x -= 2.0 * (negWallX - negX);
                // Also need to check that update will not push you out the other side
                // std::cout << "Pushed back out of box X, resample M and traj " << '\n';
                if (posWallX - 2.0 * (negWallX - negX) > posX) {
                    moveFailed = true;
                    // std::cout << " Will push out the other side. negwall " << negWallX << " poswall: " << posWallX
                    //           << " Trans.x: " << targCom.trajTrans.x << '\n';
                }
            }
            if (!outsideNegX && outsidePosX) {
                // std::cout << "Pushed back out of box X, resample M and traj " << '\n';
                targCom.trajTrans.x -= 2.0 * (posWallX - posX);
                if (negWallX - 2.0 * (posWallX - posX) < negX) {
                    moveFailed = true;
                    // std::cout << " Will push out the other side. negwall " << negWallX << " poswall: " << posWallX
                    //           << " Trans.x: " << targCom.trajTrans.x << '\n';
                }
            }

            // Check the y-dimension
            if (outsideNegY && outsidePosY) {
                // For a large complex, test if checkItr could be pushed back out the other side
                // std::cout << " extends in both directions of Y. negwall " << negWallY
                //           << " poswall: " << posWallY << '\n';
                // std::cout << "Pushed back out of box Y, resample M and traj " << '\n';
                moveFailed = true;
            }
            if (outsideNegY && !outsidePosY) {
                targCom.trajTrans.y -= 2.0 * (negWallY - negY);
                // Also need to check that update will not push you out the other side
                // std::cout << "Pushed back out of box Y, resample M and traj " << '\n';
                if (posWallY - 2.0 * (negWallY - negY) > posY) {
                    moveFailed = true;
                    // std::cout << " Will push out the other side. negwall " << negWallY << " poswall: " << posWallY
                    //           << " Trans.y: " << targCom.trajTrans.y << '\n';
                }
            }
            if (!outsideNegY && outsidePosY) {
                // std::cout << "Pushed back out of box Y, resample M and traj " << '\n';
                targCom.trajTrans.y -= 2.0 * (posWallY - posY);
                if (negWallY - 2.0 * (posWallY - posY) < negY) {
                    moveFailed = true;
                    // std::cout << " Will push out the other side. negwall " << negWallY << " poswall: " << posWallY
                    //           << " Trans.y: " << targCom.trajTrans.y << '\n';
                }
            }

            // Check the z-dimension
            if (outsideNegZ && outsidePosZ) {
                // For a large complex, test if checkItr could be pushed back out the other side
                // std::cout << " extends in both directions of Z. negwall " << negWallZ
                //           << " poswall: " << posWallZ << '\n';
                // std::cout << "Pushed back out of box Z, resample M and traj " << '\n';
                moveFailed = true;
            }
            if (outsideNegZ && !outsidePosZ) {
                targCom.trajTrans.z -= 2.0 * (negWallZ - negZ);
                // Also need to check that update will not push you out the other side
                // std::cout << "Pushed back out of box Z, resample M and traj " << '\n';
                if (posWallZ - 2.0 * (negWallZ - negZ) > posZ) {
                    moveFailed = true;
                    // std::cout << " Will push out the other side. negwall " << negWallZ << " poswall: " << posWallZ
                    //           << " Trans.z: " << targCom.trajTrans.z << '\n';
                }
            }
            if (!outsideNegZ && outsidePosZ) {
                // std::cout << "Pushed back out of box Z, resample M and traj " << '\n';
                targCom.trajTrans.z -= 2.0 * (posWallZ - posZ);
                if (negWallZ - 2.0 * (posWallZ - posZ) < negZ) {
                    moveFailed = true;
                    // std::cout << " Will push out the other side. negwall " << negWallZ << " poswall: " << posWallZ
                    //           << " Trans.z: " << targCom.trajTrans.z << '\n';
                }
            }
        } // recheck span

        if (moveFailed == true) {
            // Resample, extends in x, y, and/or z
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
    } // loop over iterations and flag condition

    // std::cout << "ITERATIONS TO CONVERGE POSITION WITHIN BOX: " << checkItr
    //           << " flag at end: 0 = success: " << needsRecheck << '\n';
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