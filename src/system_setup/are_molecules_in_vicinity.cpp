#include "system_setup/system_setup.hpp"

bool are_molecules_in_vicinity(const Molecule& mol1, const Molecule& mol2, const std::vector<MolTemplate>& molTemplates)
{
    Vector tmpVec { mol2.comCoord - mol1.comCoord };
    tmpVec.calc_magnitude();
    double distance = tmpVec.magnitude;

    // Use the largest binding radius of each molecule interface and add a safety margin of 2.
    float radiusSum = molTemplates[mol1.molTypeIndex].radius + molTemplates[mol2.molTypeIndex].radius + 2.0;

    return distance < radiusSum;
}
