#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_copyCounters.hpp"

using namespace std;

void write_NboundPairs(copyCounters& counterArrays, ofstream& outfile, int it, const Parameters& params)
{
    int i, j;

    int index;
    outfile << it * params.timeStep << '\t';
    for (i = 0; i < counterArrays.proPairlist.size(); i++) {
        index = counterArrays.proPairlist[i];
        // outfile<<" i in loop: "<<i<<" index: "<<index<<'\t';
        outfile << counterArrays.nBoundPairs[index] << '\t';
    }
    outfile << counterArrays.nLoops << endl;
}
