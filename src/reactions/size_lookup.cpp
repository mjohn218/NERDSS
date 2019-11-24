#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

size_t size_lookup(double bindRadius, double Dtot, const Parameters& params, double Rmax)
{
    int ctr{ 0 };
    double stepSize {sqrt(Dtot * params.timeStep) / 50};

    double RIndex{ bindRadius };
    while (RIndex <= Rmax + stepSize) {
        ctr += 1;
        RIndex += stepSize;
    }

    return static_cast<size_t>(ctr);
}
