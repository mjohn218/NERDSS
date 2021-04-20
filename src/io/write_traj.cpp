#include "io/io.hpp"
#include "tracing.hpp"
#include <chrono>
#include <ctime>
#include <iomanip>

void write_traj(long long int iter, std::ofstream& trajFile, const Parameters& params, const std::vector<Molecule>& moleculeList,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject)
{
    // TRACE();
    // TODO: I need to change this so that it fills up the trajectory with empties so that molecules are always in the
    // same line in the trajectory
    std::string chain { "ABCDEFGHIJKLMNOPQRSTUVWXYZ" };
    std::vector<std::string> molTypeNames;
    for (auto& molTemp : molTemplateList)
        molTypeNames.push_back(molTemp.molName.substr(0, 2));

    trajFile << params.numTotalUnits << '\n';
    trajFile << "iteration: " << iter << std::endl;
    int numWritten { 0 };
    for (auto& mol : moleculeList) {
        if (mol.isEmpty || mol.isImplicitLipid)
            continue;
        {
            trajFile << std::setw(4) << molTypeNames[mol.molTypeIndex] << ' ' << std::fixed << mol.comCoord << '\n';
            ++numWritten;
            for (auto& iface : mol.interfaceList) {
                trajFile << std::setw(4) << molTypeNames[mol.molTypeIndex] << ' ' << std::fixed << iface.coord << '\n';
                ++numWritten;
            }
            trajFile << std::flush;
        }
    }

    while (numWritten < params.numTotalUnits) {
        trajFile << std::setw(4) << "EMTY" << ' ' << std::fixed << membraneObject.waterBox.x / 2.0 << std::fixed
                 << membraneObject.waterBox.y / 2.0 << std::fixed << membraneObject.waterBox.z / 2.0 << '\n';
        ++numWritten;
        trajFile << std::flush;
    }
}
