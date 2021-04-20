#include "reactions/association/association.hpp"

void rotate(Coord& rotOrigin, Quat& rotQuat, Complex& targCom,
    std::vector<Molecule>& moleculeList)
{
    // First rotate all points in complex 1 around the reacting interface of
    // protein p1, and translate them by vector
    for (auto& mol : targCom.memberList) {
        Vector comVec { moleculeList[mol].tmpComCoord - rotOrigin };
        rotQuat.rotate(comVec);
        moleculeList[mol].tmpComCoord = Coord(comVec.x, comVec.y, comVec.z) + Coord(rotOrigin.x, rotOrigin.y, rotOrigin.z);

        // now rotate each member molecule of the complex
        for (auto& iface : moleculeList[mol].tmpICoords) {
            // get the vector from the interface to the target interface
            Vector ifaceVec { iface - rotOrigin };
            // rotate
            rotQuat.rotate(ifaceVec);
            iface = Coord(ifaceVec.x, ifaceVec.y, ifaceVec.z) + Coord(rotOrigin.x, rotOrigin.y, rotOrigin.z);
        }
    }
}