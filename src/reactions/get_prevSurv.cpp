#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

double get_prevSurv(const gsl_matrix* survMatrix, double Dtot, double deltaT, double r0, double bindRadius)
{
    /*RstepSize is coupled to the DDmatrixcreate RstepSize, they must be the SAME definition,
     *so RstepSize should not be defined any other way!
     */
    const double RstepSize = sqrt(Dtot * deltaT) / 50;
    int index{ static_cast<int>(floor((r0 - bindRadius) / RstepSize)) };
    if (index < 0)
        index = 0;

    double val01 = gsl_matrix_get(survMatrix, 1, index);
    double val02 = gsl_matrix_get(survMatrix, 1, index + 1);
    double r01 = gsl_matrix_get(survMatrix, 0, index);
    double r02 = gsl_matrix_get(survMatrix, 0, index + 1);

    return (val01 * (r02 - r0) + val02 * (r0 - r01)) / RstepSize;
}
