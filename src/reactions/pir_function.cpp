#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "tracing.hpp"

double pir_function(double x, void* p)
{
    // TRACE();
    double alp, bet, tet, P, T, h, f;
    IntegrandParams& params = *reinterpret_cast<IntegrandParams*>(p);

    h = (2.0 * M_PI * params.a * params.D);

    if (params.k < 1.0 / 0.0) {
        alp = h * x * j1(x * params.a) + params.k * j0(x * params.a);
        bet = h * x * y1(x * params.a) + params.k * y0(x * params.a);
        tet = sqrt(alp * alp + bet * bet);

        P = (j0(x * params.r) * bet - y0(x * params.r) * alp) / tet;
        T = (j0(x * params.r0) * bet - y0(x * params.r0) * alp) / tet;

        f = x * exp(-params.D * params.t * x * x) * P * T / (2.0 * M_PI);
    } else { // absorbing boundary conditions
        alp = j0(x * params.a);
        bet = y0(x * params.a);
        tet = sqrt(alp * alp + bet * bet);

        P = (j0(x * params.r) * bet - y0(x * params.r) * alp) / tet;
        T = (j0(x * params.r0) * bet - y0(x * params.r0) * alp) / tet;

        f = x * exp(-params.D * params.t * x * x) * P * T / (2.0 * M_PI);
    }

    return f;
}
