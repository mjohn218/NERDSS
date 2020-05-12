#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "tracing.hpp"

double calc_pirr(gsl_matrix* pirMatrix, gsl_matrix* survMatrix, double RStepSize, double r, double r0, double a)
{
    // TRACE();
    // This is a 2D lookup table, decide the value of an arbitrary point (r,r0) using 2D linear interpolation
    int indexr0 { static_cast<int>(floor((r0 - a) / RStepSize)) };
    int indexr { static_cast<int>(floor((r - a) / RStepSize)) };
    double prob;

    if (indexr < 0)
        indexr = 0;
    if (indexr0 < 0)
        indexr0 = 0;

    if (indexr0 == indexr) {
        // otherwise you follow the algorithm below and get divide by zero for prob
        prob = gsl_matrix_get(pirMatrix, indexr0, indexr);
    } else {
        // find the prob of coord (r0,r) that falls inside the square (Ar0,Ar)-(Br0,Br)-(Cr0,Cr)-(Dr0,Dr)
        double Ar0 = gsl_matrix_get(survMatrix, 0, indexr0);
        double Ar = gsl_matrix_get(survMatrix, 0, indexr);
        double valA = gsl_matrix_get(pirMatrix, indexr0, indexr);
        double distA = sqrt(pow(Ar0 - r0, 2) + pow(Ar - r, 2));

        double Br0 = gsl_matrix_get(survMatrix, 0, indexr0);
        double Br = gsl_matrix_get(survMatrix, 0, indexr + 1);
        double valB = gsl_matrix_get(pirMatrix, indexr0, indexr + 1);
        double distB = sqrt(pow(Br0 - r0, 2) + pow(Br - r, 2));

        double Cr0 = gsl_matrix_get(survMatrix, 0, indexr0 + 1);
        double Cr = gsl_matrix_get(survMatrix, 0, indexr + 1);
        double valC = gsl_matrix_get(pirMatrix, indexr0 + 1, indexr + 1);
        double distC = sqrt(pow(Cr0 - r0, 2) + pow(Cr - r, 2));

        double Dr0 = gsl_matrix_get(survMatrix, 0, indexr0 + 1);
        double Dr = gsl_matrix_get(survMatrix, 0, indexr);
        double valD = gsl_matrix_get(pirMatrix, indexr0 + 1, indexr);
        double distD = sqrt(pow(Dr0 - r0, 2) + pow(Dr - r, 2));

        prob = ((valA / distA) + (valB / distB) + (valC / distC) + (valD / distD))
            / ((1 / distA) + (1 / distB) + (1 / distC) + (1 / distD));
    }
    return prob;
}
