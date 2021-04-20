#include "io/io.hpp"
#include "tracing.hpp"
#include <chrono>
#include <ctime>
#include <iomanip>

void write_xyz_assoc_cout(
    const Complex& reactCom1, const Complex& reactCom2, const std::vector<Molecule>& moleculeList)
{
    std::cout << "Complex : " << reactCom1.index << std::endl;
    for (auto& memMol : reactCom1.memberList) {
        std::cout << "Mol type : " << moleculeList[memMol].molTypeIndex << std::endl;
        std::cout << std::fixed << moleculeList[memMol].tmpComCoord << std::endl;
        for (auto& iface : moleculeList[memMol].tmpICoords) {
            std::cout << std::fixed << iface << std::endl;
        }
    }
    std::cout << "Complex : " << reactCom2.index << std::endl;
    for (auto& memMol : reactCom2.memberList) {
        std::cout << "Mol type : " << moleculeList[memMol].molTypeIndex << std::endl;
        std::cout << std::fixed << moleculeList[memMol].tmpComCoord << std::endl;
        for (auto& iface : moleculeList[memMol].tmpICoords) {
            std::cout << std::fixed << iface << std::endl;
        }
    }
}
