#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_copyCounters.hpp"
#include "tracing.hpp"

using namespace std;

void write_NboundPairs(copyCounters& counterArrays, ofstream& outfile, int it, const Parameters& params, std::vector<Molecule>& moleculeList)
{
    // TRACE();
    // if (params.implicitLipid == true) {
    //     //update NboundPairs with Implicit Lipid
    //     for (int typeIndex = 0; typeIndex < params.numMolTypes; typeIndex++) {
    //         counterArrays.nBoundPairs[params.numMolTypes * params.numMolTypes + typeIndex] = 0;
    //     }
    //     for (auto& molTemp : moleculeList) {
    //         if (molTemp.isImplicitLipid == true)
    //             continue;
    //         int bounder1 { molTemp.index };

    //         if (molTemp.linksToSurface > 0) {
    //             int p1 { moleculeList[bounder1].molTypeIndex };
    //             int index { params.numMolTypes * params.numMolTypes + p1 };
    //             counterArrays.nBoundPairs[index]++;
    //         }
    //     }
    // }

    int i, j;

    int index;
    outfile << (it - params.itrRestartFrom) * params.timeStep * 1E-6 + params.timeRestartFrom << '\t';
    for (i = 0; i < counterArrays.proPairlist.size(); i++) {
        index = counterArrays.proPairlist[i];
        // outfile<<" i in loop: "<<i<<" index: "<<index<<'\t';
        outfile << counterArrays.nBoundPairs[index] << '\t';
    }
    outfile << counterArrays.nLoops << '\t' << counterArrays.nCancelOverlapPartner << '\t' << counterArrays.nCancelOverlapSystem << '\t' << counterArrays.nCancelSpanBox << '\t' << counterArrays.nCancelDisplace2D << '\t' << counterArrays.nCancelDisplace3D << '\t' << counterArrays.nCancelDisplace3Dto2D << '\t' << counterArrays.nAssocSuccess << endl;
}
