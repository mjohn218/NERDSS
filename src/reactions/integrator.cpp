#include "io/io.hpp"
#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "tracing.hpp"

#include <iomanip>

double integrator(gsl_function F, IntegrandParams params, gsl_integration_workspace* w, double r0, double bindrad,
    double Dtot, double kr, double deltat, char* funcID, double (*f)(double, void*))
{
    // TRACE();
    double result {};
    double error {};
    const double xlow { 0 };
    const double epsabs { 1e-7 };
    const double epsrel { 1e-7 };
    double publicationcrit = 1e-6;

    /*For successful qagiu, takes ~0.01:0.3s, depending on R0.
      using qags even just twice takes ~4s-10s. qags seems necessary for large ka values. .*/
    std::chrono::steady_clock::time_point start;
    int status { gsl_integration_qagiu(&F, xlow, epsabs, epsrel, 1e6, w, &result, &error) };

    int it { 1 };
    /*Try again with a weaker criterion*/
    if (status != GSL_SUCCESS)
        status = gsl_integration_qagiu(&F, xlow, publicationcrit, publicationcrit, 1e6, w, &result, &error);

    /*Now try with a different upper integration bound. Issues can arise due to oscillations of Bessel functions,
      convergence can vary with choice of upper integration.
     */
    if (status != GSL_SUCCESS) {
        start = std::chrono::steady_clock::now();
        /*This value of umax is not optimized*/
        double umax { 10000.0 }; // sqrt(-log(DBL_MIN) / (Dtot * deltat));
        while (std::abs((*f)(umax, &params)) > 1e-10)
            umax = umax * 1.2;

        // std::cout << linebreak;
        // std::cout << std::setprecision(12);
        // std::cout << funcID << " qagiu failed with status: " << status << '\n';
        // std::cout << "No solution found with rel/abs error smaller than: " << publicationcrit << '\n';
        // std::cout << "   umax: " << umax << '\n';
        // std::cout << "   @time: " << deltat << '\n';
        // std::cout << "   ka: " << kr << '\n';
        // std::cout << "   D: " << Dtot << '\n';
        // std::cout << "   r0: " << r0 << '\n';
        // std::cout << "Truncation will be performed on the semi-infinite domain..." << '\n';
        // std::cout << std::flush;

        it = 0;
        while (status != GSL_SUCCESS) {
            /*Can we use qag, instead of qags?*/
            int key { 2 };
            status = gsl_integration_qag(&F, xlow, umax, epsabs, publicationcrit, 1e6, key, w, &result, &error);
            umax = umax * 0.9;
            it += 1;
        }

        // std::cout << "New integration upper bound " << umax / 0.9 << " found after " << it << " iterations." << '\n';
        // std::cout << "   err: " << error << '\n';
        // std::cout << linebreak;
        // std::cout << "Time of 2D table generation: "
        //           << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << '\n';
    }
    return result;
}
