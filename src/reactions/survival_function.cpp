#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

double survival_function(double x, void* p)
{
    double f{};
    IntegrandParams& params = *reinterpret_cast<IntegrandParams*>(p);

    if (params.k < 1.0 / 0.0) { // This gives the association probability
        double h = (2.0 * M_PI * params.a * params.D);

        double alp = h * x * j1(x * params.a) + params.k * j0(x * params.a);
        double bet = h * x * y1(x * params.a) + params.k * y0(x * params.a);
        double tet = alp * alp + bet * bet;
        double P = j0(x * params.a) * y1(x * params.a) - j1(x * params.a) * y0(x * params.a);
        double T = (j0(x * params.r0) * bet - y0(x * params.r0) * alp) / tet;

        f = T * P * (1 - exp(-params.D * params.t * x * x));
        f = f * params.a * params.k;
    } else { // absorbing boundary conditions... this gives the survival probability
        double alp = j0(x * params.a);
        double bet = y0(x * params.a);
        double tet = alp * alp + bet * bet;
        double P = j0(x * params.a) * y0(x * params.r0) - j0(x * params.r0) * y0(x * params.a);
        double T = 1.0 / tet;

        f = (2.0 / M_PI) * T * P * exp(-params.D * params.t * x * x) / x;
    }

    return f;
}
