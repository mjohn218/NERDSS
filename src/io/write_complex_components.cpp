#include "io/io.hpp"
#include "tracing.hpp"
#include <chrono>
#include <ctime>

void write_complex_components(long long int simItr, std::ofstream& complexFile, const Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList)
{
    // TRACE();
    complexFile << "iter: " << simItr << " Ncomplexes: " << Complex::numberOfComplexes << '\n';
    int totMols { 0 };
    for (auto& complex : complexList) {
        if (!complex.isEmpty) {
            complexFile << complex.index << ' ';
            for (int i { 0 }; i < molTemplateList.size(); ++i) {
                complexFile << molTemplateList[i].molName << ' ' << complex.numEachMol[i] << ' ';
                complexFile << '\n';
                totMols += complex.numEachMol[i];
            }
        } else
            complexFile << complex.index << " EMPTY" << std::endl;
    }
}
