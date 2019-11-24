#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

double get_prevNorm(gsl_matrix *normMatrix, double RStepSize, double r0, double bindRadius)
{
    int index{ static_cast<int>(floor((r0 - bindRadius) / RStepSize)) };
    if (index < 0)
        index = 0;

    double val01 {gsl_matrix_get(normMatrix, 1, index)};
    double val02 {gsl_matrix_get(normMatrix, 1, index + 1)};
    double r01 {gsl_matrix_get(normMatrix, 0, index)};
    double r02 {gsl_matrix_get(normMatrix, 0, index + 1)};

    return (val01 * (r02 - r0) + val02 * (r0 - r01)) / RStepSize;
}
