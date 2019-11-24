#include "math/Faddeeva.hpp"
#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

double passocF(double r0, double tCurr, double Dtot, double bindRadius, double alpha, double cof)
{
    const double fDt { 4.0 * Dtot * tCurr };
    const double sqrtfDt { sqrt(fDt) };

    const double f1 { cof * bindRadius / r0 };

    const double sqrttCurr { sqrt(tCurr) };
    const double a2 { alpha * alpha };
    const double sep { (r0 - bindRadius) / sqrtfDt }; // a

    const double e1 { 2.0 * sep * sqrttCurr * alpha + a2 * tCurr };
    const double ef1 { sep + alpha * sqrttCurr }; // a+b
    const double ep1 = exp(e1);

    double term1 { erfc(sep) };
    double term2 {};
    if (std::isinf(ep1)) {
        std::complex<double> z;
        z.real(0.0); // = 0.0;
        z.imag(ef1); // = ef1;

        double relerr { 0.0 };
        std::complex<double> value { Faddeeva::w(z, relerr) };
        double ea2 { exp(-sep * sep) };
        term2 = ea2 * real(value);

    } else {
        term2 = exp(e1) * erfc(ef1);
    }

    return (term1 - term2) * f1;
}
