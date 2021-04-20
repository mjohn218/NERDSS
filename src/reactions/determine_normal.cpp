#include "reactions/association/association.hpp"
#include "tracing.hpp"

Vector determine_normal(Vector normal, const MolTemplate& molTemplate, Molecule oneMol)
{
    // TRACE();
    if (oneMol.interfaceList.empty()) {
        // std::cout << "Molecule is a point, has no normal." << std::endl;
        return { 0, 0, 0 };
    }

    for (auto& coord : oneMol.tmpICoords)
        coord -= oneMol.tmpComCoord;
    oneMol.tmpComCoord -= oneMol.tmpComCoord;
    Molecule tmpMol { oneMol }; // tmp molecule for rot check later

    // see if it's already oriented, this is not necessary.
    /*  int numUnmatched { 0 };
    for (unsigned ifaceItr { 0 }; ifaceItr < oneMol.interfaceList.size(); ++ifaceItr) {
        Vector diffVec {oneMol.tmpICoords[ifaceItr] - molTemplate.interfaceList[ifaceItr].iCoord};
        diffVec.calc_magnitude();
        if (std::abs(diffVec.magnitude) > 1E-8)
            ++numUnmatched;
	    }*/
    //    if (numUnmatched != 0) {
    Quat totalRotQuat { orient_crds_to_template(molTemplate, oneMol) };

    totalRotQuat = totalRotQuat.unit();
    totalRotQuat = totalRotQuat.inverse();
    totalRotQuat.rotate(normal);

    { // check to make sure the rotations were successful
        for (unsigned ifaceIndex { 0 }; ifaceIndex < oneMol.tmpICoords.size(); ifaceIndex++) {
            Vector tmpVec { oneMol.tmpICoords[ifaceIndex] - oneMol.tmpComCoord };
            totalRotQuat.rotate(tmpVec);
            oneMol.tmpICoords[ifaceIndex] = Coord(tmpVec.x, tmpVec.y, tmpVec.z);
            if (oneMol.tmpICoords[ifaceIndex] != tmpMol.tmpICoords[ifaceIndex]) {
                // std::cout << "Backwards rotation unsuccessful on interface " << ifaceIndex << std::endl;
                return { 0, 0, 0 };
            }
        }
    }
    //    }

    normal.normalize(); // just in case
    return normal;
}
