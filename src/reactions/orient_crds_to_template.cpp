#include "reactions/association/association.hpp"
#include "tracing.hpp"

#include <cmath>

Quat orient_crds_to_template(const MolTemplate& oneTemplate, Molecule& targMol)
{
    // TRACE();
    Quat firstRot {};
    Quat secondRot {};
    { // first rotation
        // determine the interface-center of mass vector (v0, v1), the rotation vector (u), and the angle to rotate
        // (angle)
        // arbitrarily the center of mass to first interface of the MolTemplate
        Vector vec1 { oneTemplate.interfaceList[0].iCoord - oneTemplate.comCoord };
        // the center of mass to first interface of the target Molecule
        Vector vec2 { targMol.tmpICoords[0] - targMol.tmpComCoord };
        // the rotation vector formed by v0 cross v1
        vec1.calc_magnitude();
        vec2.calc_magnitude();
        Vector rotAxis { vec1.cross(vec2) };
        double angle { vec1.dot_theta(vec2) };
        // TODO: the above vectors are correct
        // determine if we need to flip the sign
        //        if (requiresSignFlip(rotAxis, vec1, vec2))
        //            angle = -angle;
        double sa = sin(angle / 2);
        firstRot = Quat { cos(angle / 2), sa * rotAxis.x, sa * rotAxis.y,
            sa * rotAxis.z };
        firstRot = firstRot.unit();
        {
            Vector tmpVec { targMol.tmpICoords[0] - targMol.tmpComCoord };
            firstRot.rotate(tmpVec);
            if (std::abs(tmpVec.x - oneTemplate.interfaceList[0].iCoord.x) > 1E-8 || std::abs(tmpVec.y - oneTemplate.interfaceList[0].iCoord.y) > 1E-8 || std::abs(tmpVec.z - oneTemplate.interfaceList[0].iCoord.z) > 1E-8) {
                angle = -angle;
                sa = sin(angle / 2);
                firstRot = Quat { cos(angle / 2), sa * rotAxis.x, sa * rotAxis.y,
                    sa * rotAxis.z };
                firstRot = firstRot.unit();
            }
        }

        for (auto& iface : targMol.tmpICoords) {
            Vector tmpVec { iface - targMol.tmpComCoord };
            firstRot.rotate(tmpVec);
            iface = Coord { tmpVec.x, tmpVec.y, tmpVec.z };
        }
    }

    if (targMol.interfaceList.size() > 1) {
        // if the protein has more than one interface, use a second one to make sure all the interfaces line up
        // First check to make sure they're not in a line. If so, use the last interface
        size_t ifaceIndex { 1 };
        {
            Vector ifaceVec1 { oneTemplate.interfaceList[0].iCoord - oneTemplate.comCoord };
            Vector ifaceVec2 { oneTemplate.interfaceList[1].iCoord - oneTemplate.comCoord };
            ifaceVec1.calc_magnitude();
            ifaceVec2.calc_magnitude();

            double ang1 { ifaceVec1.dot_theta(ifaceVec2) };
            if ((ang1 == 0 || ang1 == M_PI) && !oneTemplate.isRod && oneTemplate.interfaceList.size() > 2) {
                size_t tmpIndex { oneTemplate.interfaceList.size() - 1 };
                Vector ifaceVec3 { oneTemplate.interfaceList[tmpIndex].iCoord - oneTemplate.comCoord };
                ifaceVec3.calc_magnitude();
                double ang2 { ifaceVec1.dot_theta(ifaceVec3) };
                if (ang2 == 0 || ang2 == M_PI) {
                    // continue on or quit?
                    ifaceIndex = tmpIndex; // TODO: ONLY FOR NOW
                } else {
                    ifaceIndex = tmpIndex;
                }
                ifaceIndex = tmpIndex;
            }
        }
        Vector v0 { oneTemplate.interfaceList[ifaceIndex].iCoord - oneTemplate.comCoord };
        Vector v1 { targMol.tmpICoords[ifaceIndex] - targMol.tmpComCoord };
        Vector rotAxis { targMol.tmpICoords[0] - targMol.tmpComCoord };
        v0.calc_magnitude();
        v1.calc_magnitude();
        rotAxis.normalize();

        // project the current and desired iface-com vectors
        // onto a plane of which the first iface-com vector is normal to
        Vector projVec0(v0.vector_projection(rotAxis));
        Vector projVec1(v1.vector_projection(rotAxis));
        projVec0.calc_magnitude();
        projVec1.calc_magnitude();
        double angle { projVec0.dot_theta(projVec1) };

        if (requiresSignFlip(rotAxis, projVec0, projVec1))
            angle = -angle;

        double sa { std::sin(angle / 2) };
        secondRot = Quat(cos(angle / 2), sin(angle / 2) * rotAxis.x, sin(angle / 2) * rotAxis.y, sin(angle / 2) * rotAxis.z);
        secondRot = secondRot.unit();

        //        {
        //            Vector tmpVec { targMol.tmpICoords[ifaceIndex] - targMol.tmpComCoord };
        //            firstRot.rotate(tmpVec);
        //            if(std::abs(tmpVec.x - oneTemplate.interfaceList[ifaceIndex].iCoord.x) > 1E-8 || std::abs(tmpVec.y - oneTemplate.interfaceList[ifaceIndex].iCoord.y) > 1E-8 || std::abs(tmpVec.z - oneTemplate.interfaceList[ifaceIndex].iCoord.z) > 1E-8) {
        //                angle = -angle;
        //                sa = sin(angle/2);
        //                secondRot = Quat { cos(angle / 2), sa * rotAxis.x, sa * rotAxis.y,
        //                                  sa * rotAxis.z };
        //                secondRot = secondRot.unit();
        //            }
        //        }

        // might not need this, but good for a check after
        // to make sure it worked
        for (auto& iface : targMol.tmpICoords) {
            Vector tmpVec { iface - targMol.tmpComCoord };
            secondRot.rotate(tmpVec);
            iface = Coord { tmpVec.x, tmpVec.y, tmpVec.z };
        }
    }

    // now use the inverse quat product to rotate the norm to
    // its actual position
    return (targMol.tmpICoords.size() > 1) ? (secondRot * firstRot) : firstRot;
}
