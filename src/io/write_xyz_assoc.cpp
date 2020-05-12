#include "io/io.hpp"
#include "tracing.hpp"
#include <chrono>
#include <ctime>
#include <iomanip>

void write_xyz_assoc(
    std::string filename, const Complex& reactCom1, const Complex& reactCom2, const std::vector<Molecule>& moleculeList)
{
    // TRACE();
    std::ofstream out(filename);

    int totAtoms { int(reactCom1.memberList.size() + reactCom2.memberList.size()) };
    for (auto memMol : reactCom1.memberList)
        totAtoms += (int)moleculeList[memMol].interfaceList.size();
    for (auto memMol : reactCom2.memberList)
        totAtoms += (int)moleculeList[memMol].interfaceList.size();

    out << totAtoms << std::endl;
    out << "mol output final" << std::endl;
    for (auto& memMol : reactCom1.memberList) {
        out << 'A' << ' ' << std::fixed << moleculeList[memMol].tmpComCoord << std::endl;
        for (auto& iface : moleculeList[memMol].tmpICoords) {
            out << 'A' << ' ' << std::fixed << iface << std::endl;
        }
    }

    for (auto& memMol : reactCom2.memberList) {
        out << 'B' << ' ' << std::fixed << moleculeList[memMol].tmpComCoord << std::endl;
        for (auto& iface : moleculeList[memMol].tmpICoords) {
            out << 'B' << ' ' << std::fixed << iface << std::endl;
        }
    }
}
