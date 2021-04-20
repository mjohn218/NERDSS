#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_copyCounters.hpp"
#include "tracing.hpp"

using namespace std;

/*N bounds pairs are the count of all directly bound pairs of protein A and its partner protein B. The fact that A or B
may also be bound to other proteins does not matter, so it counts all A-B bonds that exist in the system*/

void init_NboundPairs(
    copyCounters& counterArrays, ofstream& outfile, Parameters params, std::vector<MolTemplate>& molTemplateList, std::vector<Molecule>& moleculeList)
{
    // TRACE();
    int i, j;
    bool bindSurface = false;
    int index;
    int p1, p2;
    outfile << "TIME(s)" << '\t';
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
        // if (molTemplateList[p1].bindToSurface == true) {
        //     bindSurface = true; // YES
        //     std::cout << " SURFACE TREATED AS AN IMPLICIT-LIPID MODEL " << std::endl;
        //     index = params.numMolTypes * params.numMolTypes + p1;
        //     counterArrays.proPairlist.push_back(index);
        //     outfile << "'" << molTemplateList[p1].molName << ","
        //             << "MEM"
        //             << "'" << '\t';
        //     cout << "Pro pair: "
        //          << " p1: " << p1 << " to surface: "
        //          << " index: " << index << ' ' << molTemplateList[p1].molName << '\n';
        // }
    }
    outfile << "Nloops" << '\t' << "nOverlapPartner" << '\t' << "nOverlapSystem" << '\t' << "nOverlapSpanBox" << '\t' << "nDisplace2D" << '\t' << "nDisplace3D" << '\t' << "nDisplace3Dto2D" << '\t' << "nAssocSuccess" << endl;

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

    // For restart, initialize bound copies by loop the moleculeList
    if (params.fromRestart == true) {
        for (auto& molTemp : moleculeList) {
            int bounder1 { molTemp.index };
            for (auto& bounder2 : molTemp.bndpartner) {
                // make sure bound1 <= bound2
                if (bounder1 >= bounder2) {
                    // find out the index in proPairList
                    int p1 { moleculeList[bounder1].molTypeIndex };
                    int p2 { moleculeList[bounder2].molTypeIndex };
                    // here, make sure p2 >= p1
                    if (p2 < p1) {
                        int temp { p2 };
                        p2 = p1;
                        p1 = temp;
                    }
                    int index { p1 * params.numMolTypes + p2 };
                    counterArrays.nBoundPairs[index]++;
                }
            }
            // if (molTemp.linksToSurface > 0) {
            //     int p1 { moleculeList[bounder1].molTypeIndex };
            //     int index { params.numMolTypes * params.numMolTypes + p1 };
            //     counterArrays.nBoundPairs[index]++;
            // }
        }
    }
}
