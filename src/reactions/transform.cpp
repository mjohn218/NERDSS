#include "reactions/association/association.hpp"

void transform(Coord& reactIface, Molecule& reactMol1, Molecule& reactMol2, const Vector& axis)
{
    Vector alignAxis {};
//     if (isOnMembrane)
//         alignAxis = Vector { 0, 1, 0 };
//     else
    alignAxis = Vector { 0, 0, 1 };//z-axis

    alignAxis.magnitude = 1.0;

    // return without transformation if axis is already on the z-axis
    double ang { alignAxis.dot_theta(axis) };
    if (ang == M_PI || ang == 0)
        return;

    Vector rotAxis { alignAxis.cross(axis) };
    rotAxis.normalize();
    double theta { axis.dot_theta(alignAxis) };

    // I cant remember why the angles are negative
    // TODO: might need to have a signs check
    Quat rotQuat(
        cos(-theta / 2), sin(-theta / 2) * rotAxis.x, sin(-theta / 2) * rotAxis.y, sin(-theta / 2) * rotAxis.z);

    { // base1
        // rotate COM
        Vector comVec { reactMol1.tmpComCoord - reactIface };
        rotQuat.rotate(comVec);
        reactMol1.tmpComCoord = Coord(comVec.x, comVec.y, comVec.z) + reactIface;

        // rotate all the other interfaces
        for (auto& coord : reactMol1.tmpICoords) {
            Vector ifaceVec { coord - reactIface };
            // rotate the vector
            rotQuat.rotate(ifaceVec);
            coord = Coord(ifaceVec.x, ifaceVec.y, ifaceVec.z) + reactIface;
        }
    }
    { // base2
        Vector comVec { reactMol2.tmpComCoord - reactIface };
        rotQuat.rotate(comVec);
        reactMol2.tmpComCoord = Coord(comVec.x, comVec.y, comVec.z) + reactIface;

        // rotate the ifaces
        for (auto& coord : reactMol2.tmpICoords) {
            Vector ifaceVec { coord - reactIface };
            rotQuat.rotate(ifaceVec);
            coord = Coord(ifaceVec.x, ifaceVec.y, ifaceVec.z) + reactIface;
        }
    }
}
