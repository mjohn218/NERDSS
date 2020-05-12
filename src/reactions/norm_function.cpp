#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "tracing.hpp"

#include <gsl/gsl_sf_bessel.h>

double norm_function(double x, void* p)
{
    // TRACE();
    IntegrandParams& params = *reinterpret_cast<IntegrandParams*>(p);
    double temp { x * params.r0 / (2.0 * params.D * params.t) };
    double temp2 { x / (2.0 * params.D * params.t) };
    double temp3 = (exp(temp - (params.r0 * params.r0 + (x * x)) / (4.0 * params.t * params.D)));
    return temp2 * temp3 * gsl_sf_bessel_I0_scaled(temp);
}
