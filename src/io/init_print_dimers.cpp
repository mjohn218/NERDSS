#include "classes/class_Molecule_Complex.hpp"
#include "io/io.hpp"
#include "tracing.hpp"
#include <iostream>

using namespace std;

/*first line description of the MONODIMER output file.*/
void init_print_dimers(std::ofstream& outfile, Parameters params, std::vector<MolTemplate>& molTemplateList)
{
    // TRACE();
    int i, j;
    int size;
    int p1;

    outfile << "TIME (s) " << '\t';
    for (j = 0; j < molTemplateList.size(); j++) {
        if (molTemplateList[j].isImplicitLipid)
            continue;
        outfile << "MONO:" << molTemplateList[j].molName << '\t' << "DIMERS W:" << molTemplateList[j].molName << '\t';
    }
    outfile << endl;
}
