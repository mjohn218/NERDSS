#include "math/Faddeeva.hpp"
#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "tracing.hpp"

double pirr_pfree_ratio_psF_1D(
    double rCurr, double r0, double tCurr, double Dtot, double bindrad, double ka, double ps_prev)
{
    // combined variables
    const double fDt = 4.0 * Dtot * tCurr;
    const double sq_fDt = sqrt(fDt);
    const double sq_pifDt = sqrt(4.0 * M_PI * Dtot * tCurr);
    const double dx_m = rCurr - r0;
    const double dx_p = rCurr + r0;
    const double dx_sigma = rCurr + r0 - 2.0 * bindrad;
    const double x_sigma_m = r0 - bindrad;
    const double x_sigma_p = r0 + bindrad;
    const double sq_scale_t{ka * sqrt(tCurr / Dtot)};

    // normalized p_free
    double ep1 = exp(-dx_m * dx_m / fDt);
    double ep2 = exp(-dx_p * dx_p / fDt);
    double I_free = sqrt(M_PI * Dtot * tCurr) *
               (erf(x_sigma_m / sq_fDt) + erf(x_sigma_p / sq_fDt));
    double p_free_norm = (ep1 - ep2) / I_free;

    // p_irr
    double ep3 = exp(-dx_sigma * dx_sigma / fDt);
    double sep = dx_sigma / sq_fDt + sq_scale_t;
    double p_irr = (ep1 + ep3) / sq_pifDt - ka / Dtot * ep3 * Faddeeva::erfcx(sep);

    // reweighting ratio
    double ratio = p_irr / ps_prev / p_free_norm;

    return ratio;
}
