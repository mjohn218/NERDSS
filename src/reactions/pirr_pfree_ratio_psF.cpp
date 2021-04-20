#include "math/Faddeeva.hpp"
#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "tracing.hpp"

double pirr_pfree_ratio_psF(
    double rCurr, double r0, double tCurr, double Dtot, double bindrad, double alpha, double ps_prev, double rtol)
{
    // TRACE();

    double fDt = 4.0 * Dtot * tCurr;
    double sq_fDt = sqrt(fDt);

    double f1 = 1.0 / (sqrt(4.0 * M_PI * tCurr));
    double f2 = 1.0 / (4.0 * M_PI * r0 * sqrt(Dtot));

    double sep, dist;
    double sqrt_t = sqrt(tCurr);
    double a2 = alpha * alpha;
    double r1, term1, term2, e1, ef1, sum;

    r1 = rCurr;
    sep = r1 + r0 - 2.0 * bindrad;
    dist = r1 - r0;
    term1 = f1 * (exp(-dist * dist / fDt) + exp(-sep * sep / fDt));
    double a = sep / sq_fDt;
    double b = sqrt_t * alpha;
    e1 = 2.0 * a * b + a2 * tCurr;
    ef1 = a + b;
    double ep1 = exp(e1);
    if (std::isinf(ep1)) {
        std::complex<double> z;
        z.real(0.0); // = 0.0;
        z.imag(ef1); // = ef1;
        // cout <<"Complex number: "<<z<<endl;
        std::complex<double> value;
        double relerr = 0;
        value = Faddeeva::w(z, relerr);
        double ea2 = exp(-a * a);
        term2 = ea2 * real(value);

    } else {
        term2 = exp(e1) * erfc(ef1);
    }

    // term2=alpha*exp(e1)*erfc(ef1);

    sum = term1 - alpha * term2;
    sum *= f2 / r1;
    double pirr = sum;

    /*Normalization for no diffusion inside binding radius!*/
    double cof = f1 * f2; // equal to 1/( 8pir_0 sqrt(piDt))
    double c1 = 4.0 * M_PI * cof;
    double ndist = bindrad - r0;
    double nadist = bindrad + r0;
    double sq_P = sqrt(M_PI);
    term1 = -0.5 * fDt * exp(-ndist * ndist / fDt) - 0.5 * sq_fDt * sq_P * r0 * erf(-ndist / sq_fDt);
    term2 = 0.5 * fDt * exp(-nadist * nadist / fDt) + 0.5 * sq_fDt * sq_P * r0 * erf(nadist / sq_fDt);
    double pnorm = 1.0 - c1 * (term1 + term2); // this is then the normlization from sigma->inf
    //  cout <<"Normlization: integratl from sigma to infinity: "<<pnorm<<endl;

    double adist = r1 + r0;
    term1 = exp(-dist * dist / fDt) - exp(-adist * adist / fDt);
    double pfree = cof / r1 * term1 / pnorm; // NORMALIZES TO ONE

    double ratio;
    if (std::abs(pirr - pfree * ps_prev) < rtol)
        ratio = 1.0;
    else {
        ratio = pirr / (pfree * ps_prev);
        //    cout <<"pirr: "<<pirr<<" pfree*ps: "<<pfree*ps_prev<<" ratio: "<<ratio<<endl;
    }

    if (std::isnan(ratio) || std::isinf(ratio)) {
        // std::cout << ratio << '\n';
    }

    return ratio;
}
