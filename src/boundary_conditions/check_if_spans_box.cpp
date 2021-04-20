/*! \file check_box_span.cpp
 * ### Created on 2018-11-30 by Matthew Varga
 * ### TODO List
 * ***
 * can make this more efficient by just kicking out whenever cancelAssoc = true
 */
#include "boundary_conditions/reflect_functions.hpp"
#include "tracing.hpp"

#include <iostream>

void check_if_spans_box(bool& cancelAssoc, const Parameters& params, Complex& reactCom1, Complex& reactCom2,
    std::vector<Molecule>& moleculeList, const Membrane& membraneObject)
{
    // TRACE();
    // Associating proteins have been moved to contact. Before assigning them to the complexsame complex,
    // test to see if the complex is too big to fit in the box.
    // declare the six boundary sides of the system box;
    double posX = membraneObject.waterBox.x / 2.0;
    double negX = -membraneObject.waterBox.x / 2.0;
    double posY = membraneObject.waterBox.y / 2.0;
    double negY = -membraneObject.waterBox.y / 2.0;
    double posZ = membraneObject.waterBox.z / 2.0;
    double negZ = -membraneObject.waterBox.z / 2.0; // no considering the reflecting-surface, because here we are checking whether to span the box

    bool canBeOutsideX { ((reactCom1.radius + reactCom2.radius) > (membraneObject.waterBox.x / 2.0)) };
    bool canBeOutsideY { ((reactCom1.radius + reactCom2.radius) > (membraneObject.waterBox.y / 2.0)) };
    bool canBeOutsideZ { ((reactCom1.radius + reactCom2.radius) > (membraneObject.waterBox.z / 2.0)) };

    if (canBeOutsideZ) {
        // std::cout << "CHECK BOX SPAN in Z. Complex 1 Radius, Complex 2 radius " << reactCom1.radius << ' ' << reactCom2.radius << '\n';
        // The approximate size of the complex in z (max size) puts it as outside, now test interface positions.
        bool outsidePos { false };
        bool outsideNeg { false };
        double posWall = -membraneObject.waterBox.z / 2.0; //initialize to negative value.
        double negWall = -membraneObject.waterBox.z / 2.0;
        double zplus = posWall;
        double zneg = negWall;

        double currz;
        // find the farthest position of complex1 both in +Z and -Z
        for (int memMol : reactCom1.memberList) {
            currz = moleculeList[memMol].tmpComCoord.z;
            zplus = currz - posZ;
            zneg = negZ - currz;
            if (zplus > posWall)
                posWall = zplus; //distance to +z
            if (zneg > negWall)
                negWall = zneg; //distance to -z

            // measure each interface
            for (const auto& iface : moleculeList[memMol].tmpICoords) {
                currz = iface.z;
                zplus = currz - posZ;
                zneg = negZ - currz;
                if (zplus > posWall)
                    posWall = zplus; //distance to +z
                if (zneg > negWall)
                    negWall = zneg; //distance to -z
            }
        }
        // find the farthest position of complex2 both in +Z and -Z
        for (int memMol : reactCom2.memberList) {
            currz = moleculeList[memMol].tmpComCoord.z;
            zplus = currz - posZ;
            zneg = negZ - currz;
            if (zplus > posWall)
                posWall = zplus; //distance to +z
            if (zneg > negWall)
                negWall = zneg; //distance to -z

            // measure each interface
            for (const auto& iface : moleculeList[memMol].tmpICoords) {
                currz = iface.z;
                zplus = currz - posZ;
                zneg = negZ - currz;
                if (zplus > posWall)
                    posWall = zplus; //distance to +z
                if (zneg > negWall)
                    negWall = zneg; //distance to -z
            }
        }

        // check whether to span the box
        if (posWall > 0)
            outsidePos = true;
        if (negWall > 0)
            outsideNeg = true;
        // translation or cancel
        if (outsideNeg && outsidePos) {
            // std::cout << "STICKS OUT BOTH SIZES in Z, CANCEL ASSOCIATION " << '\n';
            cancelAssoc = true;
        }
        /*Also check if it sticks out far enough in one direction, that pushing back in will cause 
	       it to stick out the other side.
	     */
        if (posWall + negWall > 0) {
            // std::cout << "WILL STICK OUT BOTH SIZES in Z, CANCEL ASSOCIATION " << '\n';
            cancelAssoc = true;
        }
        if (outsideNeg && !outsidePos) {
            // put back in the box. put at edge, rather than bouncing off.
            double transZ { -negWall };
            reactCom1.tmpComCoord.z -= transZ;
            for (int memMol : reactCom1.memberList) {
                moleculeList[memMol].tmpComCoord.z -= transZ;
                // update interface coords
                for (auto& iface : moleculeList[memMol].tmpICoords)
                    iface.z -= transZ;
            }
            reactCom2.tmpComCoord.z -= transZ;
            for (int memMol : reactCom2.memberList) {
                moleculeList[memMol].tmpComCoord.z -= transZ;
                // update interface coords
                for (auto& iface : moleculeList[memMol].tmpICoords)
                    iface.z -= transZ;
            }
        }
        if (!outsideNeg && outsidePos) {
            // put back in the box. put at edge, rather than bouncing off.
            double transZ { posWall };
            reactCom1.tmpComCoord.z -= transZ;
            for (int memMol : reactCom1.memberList) {
                moleculeList[memMol].tmpComCoord.z -= transZ;
                // update interface coords
                for (auto& iface : moleculeList[memMol].tmpICoords)
                    iface.z -= transZ;
            }
            reactCom2.tmpComCoord.z -= transZ;
            for (int memMol : reactCom2.memberList) {
                moleculeList[memMol].tmpComCoord.z -= transZ;
                // update interface coords
                for (auto& iface : moleculeList[memMol].tmpICoords)
                    iface.z -= transZ;
            }
        }
    }

    if (canBeOutsideY) {
        // std::cout << "CHECK BOX SPAN in Y. Complex 1 Radius, Complex 2 radius " << reactCom1.radius << ' ' << reactCom2.radius << '\n';
        // The approximate size of the complex in y (max size) puts it as outside, now test interface positions.
        bool outsidePos { false };
        bool outsideNeg { false };
        //double posWall{ negY };
        // double negWall{ posY };
        double curry;
        double posWall = -membraneObject.waterBox.y / 2.0; //initialize to negative value.
        double negWall = -membraneObject.waterBox.y / 2.0;
        double yplus = posWall;
        double yneg = negWall;
        // find the farthest position of complex1 both in +Y and -Y
        for (int memMol : reactCom1.memberList) {
            curry = moleculeList[memMol].tmpComCoord.y;
            /*These will be positive if there is extension beyond the box*/
            yplus = curry - posY;
            yneg = negY - curry;
            if (yplus > posWall)
                posWall = yplus;
            if (yneg > negWall)
                negWall = yneg;
            // measure each interface
            for (const auto& iface : moleculeList[memMol].tmpICoords) {
                curry = iface.y;
                yplus = curry - posY;
                yneg = negY - curry;
                if (yplus > posWall)
                    posWall = yplus; //distance to +y
                if (yneg > negWall)
                    negWall = yneg; //distance to -y
            }
        }
        // find the farthest position of complex2 both in +y and -y
        for (int memMol : reactCom2.memberList) {
            curry = moleculeList[memMol].tmpComCoord.y;
            yplus = curry - posY;
            yneg = negY - curry;
            if (yplus > posWall)
                posWall = yplus; //distance to +y
            if (yneg > negWall)
                negWall = yneg; //distance to -y

            // measure each interface
            for (const auto& iface : moleculeList[memMol].tmpICoords) {
                curry = iface.y;
                yplus = curry - posY;
                yneg = negY - curry;
                if (yplus > posWall)
                    posWall = yplus; //distance to +y
                if (yneg > negWall)
                    negWall = yneg; //distance to -y
            }
        }
        // check whether to span the box
        if (posWall > 0)
            outsidePos = true;
        if (negWall > 0)
            outsideNeg = true;
        // translation or cancel
        if (outsideNeg && outsidePos) {
            // std::cout << "STICKS OUT BOTH SIZES in Y, CANCEL ASSOCIATION " << '\n';
            cancelAssoc = true;
        }
        /*Also check if it sticks out far enough in one direction, that pushing back in will cause 
	  it to stick out the other side.
	 */
        if (posWall + negWall > 0) {
            // std::cout << "WILL STICK OUT BOTH SIZES in Y, CANCEL ASSOCIATION " << posWall << ' ' << negWall << '\n';
            cancelAssoc = true;
        }

        if (outsideNeg && !outsidePos) {
            // put back in the box. put at edge, rather than bouncing off.
            double transy { -negWall };
            reactCom1.tmpComCoord.y -= transy;
            for (int memMol : reactCom1.memberList) {
                moleculeList[memMol].tmpComCoord.y -= transy;
                // update interface coords
                for (auto& iface : moleculeList[memMol].tmpICoords)
                    iface.y -= transy;
            }
            reactCom2.tmpComCoord.y -= transy;
            for (int memMol : reactCom2.memberList) {
                moleculeList[memMol].tmpComCoord.y -= transy;
                // update interface coords
                for (auto& iface : moleculeList[memMol].tmpICoords)
                    iface.y -= transy;
            }
        }
        if (!outsideNeg && outsidePos) {
            // put back in the box. put at edge, rather than bouncing off.
            double transy { posWall };
            reactCom1.tmpComCoord.y -= transy;
            for (int memMol : reactCom1.memberList) {
                moleculeList[memMol].tmpComCoord.y -= transy;
                // update interface coords
                for (auto& iface : moleculeList[memMol].tmpICoords)
                    iface.y -= transy;
            }
            reactCom2.tmpComCoord.y -= transy;
            for (int memMol : reactCom2.memberList) {
                moleculeList[memMol].tmpComCoord.y -= transy;
                // update interface coords
                for (auto& iface : moleculeList[memMol].tmpICoords)
                    iface.y -= transy;
            }
        }
    }

    if (canBeOutsideX) {
        // std::cout << "CHECK BOX SPAN in X. Complex 1 Radius, Complex 2 radius " << reactCom1.radius << ' ' << reactCom2.radius << '\n';
        // The approximate size of the complex in x (max size) puts it as outside, now test interface positions.
        bool outsidePos { false };
        bool outsideNeg { false };
        double posWall = -membraneObject.waterBox.x / 2.0; //initialize to negative value.
        double negWall = -membraneObject.waterBox.x / 2.0;
        double xplus = posWall;
        double xneg = negWall;

        double currx;
        // find the farthest position of complex1 both in +X and -X
        for (int memMol : reactCom1.memberList) {
            currx = moleculeList[memMol].tmpComCoord.x;
            xplus = currx - posX;
            xneg = negX - currx;
            if (xplus > posWall)
                posWall = xplus; //distance to +x
            if (xneg > negWall)
                negWall = xneg; //distance to -x
                    // measure each interface
            for (const auto& iface : moleculeList[memMol].tmpICoords) {
                currx = iface.x;
                xplus = currx - posX;
                xneg = negX - currx;
                if (xplus > posWall)
                    posWall = xplus; //distance to +x
                if (xneg > negWall)
                    negWall = xneg; //distance to -x
            }
        }
        // find the farthest position of complex2 both in +x and -x
        for (int memMol : reactCom2.memberList) {
            currx = moleculeList[memMol].tmpComCoord.x;
            xplus = currx - posX;
            xneg = negX - currx;
            if (xplus > posWall)
                posWall = xplus; //distance to +x
            if (xneg > negWall)
                negWall = xneg; //distance to -x

            // measure each interface
            for (const auto& iface : moleculeList[memMol].tmpICoords) {
                currx = iface.x;
                xplus = currx - posX;
                xneg = negX - currx;
                if (xplus > posWall)
                    posWall = xplus; //distance to +x
                if (xneg > negWall)
                    negWall = xneg; //distance to -x
            }
        }
        // check whether to span the box
        if (posWall > 0)
            outsidePos = true;
        if (negWall > 0)
            outsideNeg = true;
        // translation or cancel
        if (outsideNeg && outsidePos) {
            // std::cout << "STICKS OUT BOTH SIZES in X, CANCEL ASSOCIATION " << '\n';
            cancelAssoc = true;
        }
        /*Also check if it sticks out far enough in one direction, that pushing back in will cause 
	  it to stick out the other side.
	 */
        if (posWall + negWall > 0) {
            // std::cout << "WILL STICK OUT BOTH SIZES in X, CANCEL ASSOCIATION " << '\n';
            cancelAssoc = true;
        }

        if (outsideNeg && !outsidePos) {
            // put back in the box. put at edge, rather than bouncing off.
            double transx { -negWall };
            reactCom1.tmpComCoord.x -= transx;
            for (int memMol : reactCom1.memberList) {
                moleculeList[memMol].tmpComCoord.x -= transx;
                // update interface coords
                for (auto& iface : moleculeList[memMol].tmpICoords)
                    iface.x -= transx;
            }
            reactCom2.tmpComCoord.x -= transx;
            for (int memMol : reactCom2.memberList) {
                moleculeList[memMol].tmpComCoord.x -= transx;
                // update interface coords
                for (auto& iface : moleculeList[memMol].tmpICoords)
                    iface.x -= transx;
            }
        }
        if (!outsideNeg && outsidePos) {
            // put back in the box. put at edge, rather than bouncing off.
            double transx { posWall };
            reactCom1.tmpComCoord.x -= transx;
            for (int memMol : reactCom1.memberList) {
                moleculeList[memMol].tmpComCoord.x -= transx;
                // update interface coords
                for (auto& iface : moleculeList[memMol].tmpICoords)
                    iface.x -= transx;
            }
            reactCom2.tmpComCoord.x -= transx;
            for (int memMol : reactCom2.memberList) {
                moleculeList[memMol].tmpComCoord.x -= transx;
                // update interface coords
                for (auto& iface : moleculeList[memMol].tmpICoords)
                    iface.x -= transx;
            }
        }
    }
}