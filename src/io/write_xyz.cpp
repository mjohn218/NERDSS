#include "io/io.hpp"
#include "tracing.hpp"
#include <chrono>
#include <ctime>
#include <iomanip>

void write_xyz(std::string filename, const Parameters& params, const std::vector<Molecule>& moleculeList,
    const std::vector<MolTemplate>& molTemplateList)
{
    // TRACE();
    std::ofstream out(filename);
    std::string lett { "ABCDEFGHIJKLMNOPQRSTUVWXYZ" };
    std::vector<std::string> names;
    for (auto& molTemp : molTemplateList)
        names.push_back(molTemp.molName);
    out << params.numTotalUnits << std::endl;
    out << "mol output final" << std::endl;
    int numWritten { 0 };
    for (auto& mol : moleculeList) {
        if (molTemplateList[mol.molTypeIndex].isImplicitLipid) {
            continue;
        }
        out << std::setw(4) << names[mol.molTypeIndex] << ' ' << std::fixed << mol.comCoord << std::endl;
        ++numWritten;
        for (auto& iface : mol.interfaceList) {
            out << std::setw(4) << names[mol.molTypeIndex] << ' ' << std::fixed << iface.coord << std::endl;
            ++numWritten;
        }
    }
    while (numWritten < params.numTotalUnits) {
        out << std::setw(4) << "EMTY" << ' ' << std::fixed << Coord { 150.0, 150.0, 150.0 } << std::endl;
        ++numWritten;
    }
}
