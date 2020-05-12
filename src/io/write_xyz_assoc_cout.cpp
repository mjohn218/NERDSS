#include "io/io.hpp"
#include "tracing.hpp"
#include <chrono>
#include <ctime>
#include <iomanip>

void write_xyz_assoc_cout(
    const Complex& reactCom1, const Complex& reactCom2, const std::vector<Molecule>& moleculeList)
{
    // TRACE();
    //std::ofstream out(filename);

    int totAtoms { int(reactCom1.memberList.size() + reactCom2.memberList.size()) };
    for (auto memMol : reactCom1.memberList)
        totAtoms += (int)moleculeList[memMol].interfaceList.size();
    for (auto memMol : reactCom2.memberList)
        totAtoms += (int)moleculeList[memMol].interfaceList.size();

    //std::cout << totAtoms << std::endl;
    std::cout << "mol std::coutput final" << std::endl;
    for (auto& memMol : reactCom1.memberList) {
        //std::cout << 'A' << ' ' << std::fixed << moleculeList[memMol].tmpComCoord << std::endl;
        std::cout << std::fixed << moleculeList[memMol].tmpComCoord << std::endl;
        for (auto& iface : moleculeList[memMol].tmpICoords) {
            //std::cout << 'A' << ' ' << std::fixed << iface << std::endl;
            std::cout << std::fixed << iface << std::endl;
        }
    }
    std::cout << " Complex : " << reactCom2.index << std::endl;
    for (auto& memMol : reactCom2.memberList) {
        if (moleculeList[memMol].isImplicitLipid == true) {
            std::cout << "Implicit Lipid" << std::endl;
            continue;
        }
        //std::cout << 'B' << ' ' << std::fixed << moleculeList[memMol].tmpComCoord << std::endl;
        std::cout << std::fixed << moleculeList[memMol].tmpComCoord << std::endl;
        for (auto& iface : moleculeList[memMol].tmpICoords) {
            //std::cout << 'B' << ' ' << std::fixed << iface << std::endl;
            std::cout << std::fixed << iface << std::endl;
        }
    }
}
