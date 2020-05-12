#include "io/io.hpp"
#include "tracing.hpp"
#include <chrono>
#include <ctime>

void write_complex_crds(std::string name, const Complex& complex1, const Complex& complex2, std::vector<Molecule>& moleculeList)
{
    // TRACE();
    for (auto& mp : complex1.memberList) {
        std::ofstream out("out/c" + std::to_string(complex1.index) + "_p" + std::to_string(mp) + "_" + name + ".dat");
        out << moleculeList[mp].molTypeIndex << ' ' << moleculeList[mp].myComIndex << std::endl;
        moleculeList[mp].write_crd_file(out);
    }
    for (auto& mp : complex2.memberList) {
        std::ofstream out("out/c" + std::to_string(complex2.index) + "_p" + std::to_string(mp) + "_" + name + ".dat");
        out << moleculeList[mp].molTypeIndex << ' ' << moleculeList[mp].myComIndex << std::endl;
        moleculeList[mp].write_crd_file(out);
    }
}
