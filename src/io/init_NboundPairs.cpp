#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_copyCounters.hpp"

using namespace std;

/*N bounds pairs are the count of all directly bound pairs of protein A and its partner protein B. The fact that A or B
may also be bound to other proteins does not matter, so it counts all A-B bonds that exist in the system*/

void init_NboundPairs(
    copyCounters& counterArrays, ofstream& outfile, Parameters params, std::vector<MolTemplate>& molTemplateList)
{
    int i, j;
    bool bindSurface = false;
    int index;
    int p1, p2;
    outfile << "TIME(us)" << '\t';
    for (p1 = 0; p1 < molTemplateList.size(); p1++) {
        for (j = 0; j < molTemplateList[p1].rxnPartners.size(); j++) {

            p2 = molTemplateList[p1].rxnPartners[j];
            if (p2 >= p1) {

                index = p1 * params.numMolTypes
                    + p2; // only store pair of proteins once, so first index (p1) must be <= second index.
                counterArrays.proPairlist.push_back(index);

                outfile << "'" << molTemplateList[p1].molName << "," << molTemplateList[p2].molName << "'" << '\t';
                cout << "Pro pair: "
                     << " p1: " << p1 << " p2: " << p2 << " index: " << index << ' ' << molTemplateList[p1].molName
                     << "," << molTemplateList[p2].molName << "'" << '\t';
            }
        }
        if (molTemplateList[p1].bindToSurface == true) {
            bindSurface = true; // YES
            std::cout << " SURFACE TREATED AS AN IMPLICIT-LIPID MODEL " << std::endl;
            index = params.numMolTypes * params.numMolTypes + p1;
            counterArrays.proPairlist.push_back(index);
            outfile << "'" << molTemplateList[p1].molName << ","
                    << "MEM"
                    << "'" << '\t';
            cout << "Pro pair: "
                 << " p1: " << p1 << " to surface: "
                 << " index: " << index << ' ' << molTemplateList[p1].molName << '\n';
        }
    }
    outfile << "Nloops" << endl;

    /*Initialize copy numbers of bound pairs to zero*/

    int pairSize = params.numMolTypes * params.numMolTypes;
    cout << "size of proPairlist: " << counterArrays.proPairlist.size() << " Number of pairs: " << pairSize
         << " numMolTypes: " << params.numMolTypes << endl;

    if (bindSurface == true)
        pairSize
            = params.numMolTypes * params.numMolTypes + params.numMolTypes; // for each protein binding to the surface
    // initiaize bound copies to zero. For restart, this will need to be updated.
    counterArrays.nBoundPairs.reserve(pairSize);
    for (int i = 0; i < pairSize; i++)
        counterArrays.nBoundPairs.push_back(0);
}
