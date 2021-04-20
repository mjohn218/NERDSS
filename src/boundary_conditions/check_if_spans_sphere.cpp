/*! \file check_if_spans_sphere.cpp
 * ### Created on 2020-02-23 by Yiben Fu
 */
#include "boundary_conditions/reflect_functions.hpp"
#include "reactions/association/association.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"
#include "tracing.hpp"

#include <iostream>

void check_if_spans_sphere(bool& cancelAssoc, const Parameters& params, Complex& reactCom1, Complex& reactCom2,
    std::vector<Molecule>& moleculeList, const Membrane& membraneObject)
{
    // TRACE();
    // Associating proteins have been moved to contact. Before assigning them to the complexsame complex,
    // test to see if the complex is too big to fit in the box.

    // declare the boundary side of the system;
    double sphereR = membraneObject.sphereR; // no considering the reflecting-surface, because here we are checking whether to span the box

    // find the new COM and new radius, to check whether the new radius is larger than the sphere radius
    Coord newCom;
    double newRadius = 0.0;
    com_of_two_tmp_complexes(reactCom1, reactCom2, newCom, moleculeList);
    for (auto& memMol : reactCom1.memberList) {
        if (moleculeList[memMol].isImplicitLipid)
            continue;

        Vector disVec { moleculeList[memMol].tmpComCoord - newCom };
        disVec.calc_magnitude();
        if (disVec.magnitude > newRadius)
            newRadius = disVec.magnitude;
        for (const auto& iface : moleculeList[memMol].tmpICoords) {
            disVec = Vector(iface - newCom);
            disVec.calc_magnitude();
            if (disVec.magnitude > newRadius)
                newRadius = disVec.magnitude;
        }
    }
    for (auto& memMol : reactCom2.memberList) {
        if (moleculeList[memMol].isImplicitLipid)
            continue;

        Vector disVec { moleculeList[memMol].tmpComCoord - newCom };
        disVec.calc_magnitude();
        if (disVec.magnitude > newRadius)
            newRadius = disVec.magnitude;
        for (const auto& iface : moleculeList[memMol].tmpICoords) {
            disVec = Vector(iface - newCom);
            disVec.calc_magnitude();
            if (disVec.magnitude > newRadius)
                newRadius = disVec.magnitude;
        }
    }

    if (newRadius > sphereR) {
        // std::cout << "STICKS OUT THE SPHERE, CANCEL ASSOCIATION " << '\n';
        cancelAssoc = true;
        return;
    }

    // std::cout << "CHECK SPHERE SPAN. Complex 1 Radius, Complex 2 radius " << reactCom1.radius << ' ' << reactCom2.radius << '\n';
    // The approximate size of the complex (max size) puts it as outside, now test interface positions.
    bool outside { false };

    Coord curr;
    double dr = 0.0;
    Coord targcrds;
    // find the farthest position of complex1
    for (int memMol : reactCom1.memberList) {
        if (moleculeList[memMol].isImplicitLipid)
            continue;

        curr = moleculeList[memMol].tmpComCoord;
        double drtmp = curr.get_magnitude() - sphereR;
        if (drtmp > dr) {
            outside = true;
            dr = drtmp;
            targcrds = curr;
        }
        // measure each interface
        for (const auto& iface : moleculeList[memMol].tmpICoords) {
            curr = iface;
            drtmp = curr.get_magnitude() - sphereR;
            if (drtmp > dr) {
                outside = true;
                dr = drtmp;
                targcrds = curr;
            }
        }
    }
    // find the farthest position of complex2
    for (int memMol : reactCom2.memberList) {
        if (moleculeList[memMol].isImplicitLipid)
            continue;

        curr = moleculeList[memMol].tmpComCoord;
        double drtmp = curr.get_magnitude() - sphereR;
        if (drtmp > dr) {
            outside = true;
            dr = drtmp;
            targcrds = curr;
        }
        // measure each interface
        for (const auto& iface : moleculeList[memMol].tmpICoords) {
            curr = iface;
            drtmp = curr.get_magnitude() - sphereR;
            if (drtmp > dr) {
                outside = true;
                dr = drtmp;
                targcrds = curr;
            }
        }
    }

    // put back in the box. put at edge, rather than bouncing off.
    if (outside == true) {
        double lamda = -dr / targcrds.get_magnitude();
        Coord trans = lamda * targcrds;
        reactCom1.tmpComCoord += trans;
        for (int memMol : reactCom1.memberList) {
            moleculeList[memMol].tmpComCoord += trans;
            // update interface coords
            for (auto& iface : moleculeList[memMol].tmpICoords)
                iface += trans;
        }
        reactCom2.tmpComCoord += trans;
        for (int memMol : reactCom2.memberList) {
            moleculeList[memMol].tmpComCoord += trans;
            // update interface coords
            for (auto& iface : moleculeList[memMol].tmpICoords)
                iface += trans;
        }
    }
}