#include "classes/class_Parameters.hpp"
#include "classes/class_copyCounters.hpp"
#include "tracing.hpp"

void update_Nboundpairs(int ptype1, int ptype2, int chg, const Parameters& params, copyCounters& counterArrays)
{
    // TRACE();
    /*After a reaction occurs, change the count of protein-protein pair bonds */
    /*To avoid storing the same pair in two bins, ptype1 must be <=ptype2*/

    int Npro = params.numMolTypes;
    int ind;
    if (ptype2 == -1) { // this is a binding to a surface reaction
        ind = Npro * Npro + ptype1;
    } else {
        ind = ptype1 * Npro + ptype2;
        if (ptype1 > ptype2)
            ind = ptype2 * Npro + ptype1;
        // cout <<"Change Nbound pairs from : "<<nBoundPairs[ind]<<" by: "<<chg<<"  index: "<<ind<<endl;
    }

    counterArrays.nBoundPairs[ind] += chg; // update by adding or subtracting 'chg' number of copies.
}
