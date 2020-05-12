#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "tracing.hpp"
#include <gsl/gsl_sf_bessel.h>

double DDpirr_pfree_ratio_ps(gsl_matrix* pirMatrix, gsl_matrix* survMatrix, gsl_matrix* normMatrix, double r, double Dtot, double deltaT, double r0, double ps_prev, double rTol, double bindRadius)
{
    // TRACE();
    double pFree {};
    double pNormVal {};
    double pirrVal {};
    double temp { r * r0 / (2.0 * Dtot * deltaT) };
    double temp2 { 1 / (4.0 * M_PI * deltaT * Dtot) };
    double temp3 { (exp(temp - (r0 * r0 + r * r) / (4.0 * deltaT * Dtot))) };
    /*RstepSize is coupled to the DDmatrixcreate RstepSize, they must be the SAME definition,
     *so RstepSize should not be defined any other way!
     */
    const double RstepSize { sqrt(Dtot * deltaT) / 50 };

    pFree = temp2 * temp3 * gsl_sf_bessel_I0_scaled(temp);
    pNormVal = get_prevNorm(normMatrix, RstepSize, r0, bindRadius);
    pirrVal = calc_pirr(pirMatrix, survMatrix, RstepSize, r, r0, bindRadius);

    double pfreeN = pFree / pNormVal; // NORMALIZES TO ONE

    double ratio;
    if (std::abs(pirrVal - pfreeN * ps_prev) < rTol)
        ratio = 1.0;
    else {
        ratio = pirrVal / (pfreeN * ps_prev);
    }
    return ratio;
}
