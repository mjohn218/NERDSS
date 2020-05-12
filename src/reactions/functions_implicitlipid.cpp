#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "tracing.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_bessel.h>
#include <math.h>

// unbinding probability
// h is the time-step; sigma is the bind_radius, Na is the number of proteins in solution,
// Nlipid is the number of lipids on the membrane surface, A is the area of membrane surface
double dissociate2D(paramsIL& parameters2D)
{
    double h = parameters2D.dt;
    double D = parameters2D.Dtot;
    double sigma = parameters2D.sigma;
    double ka = parameters2D.ka;
    double kb = parameters2D.kb / 1.0e6;
    int Na = parameters2D.Na;
    int Nlipid = parameters2D.Nlipid;
    double A = parameters2D.area;
    if (kb < 1E-15) {
        return 0.0;
    }

    double KD = kb / ka;
    double maxNaNlipid = 0;
    if (Na > Nlipid) {
        maxNaNlipid = Na;
    } else {
        maxNaNlipid = Nlipid;
    }

    double b = 2.0 * sqrt(A / M_PI / maxNaNlipid + sigma * sigma);
    double kon = 1.0 / (1.0 / ka + 1.0 / (8.0 * M_PI * D) * (4.0 * log(b / sigma) / pow(1.0 - pow(sigma / b, 2.0), 2.0) - 2.0 / (1.0 - pow(sigma / b, 2.0)) - 1.0));
    double koff = kon * KD;
    double out = 1.0 - exp(-koff * h);
    return out;
}

// a function that is necessary for other caculation
double function2D(double u, void* parameter)
{
    struct paramsIL* params = (struct paramsIL*)parameter;
    double sigma = (params->sigma);
    double D = (params->Dtot);
    double r = (params->R2D);
    double ka = (params->ka);
    double h = (params->dt);

    double Rmax = sigma + 3.0 * sqrt(4.0 * D * h);

    double H = 2.0 * M_PI * sigma * D;
    double rmax = 5 * Rmax;
    double a, b, alpha, peta;
    alpha = H * u * gsl_sf_bessel_Y1(sigma * u) + ka * gsl_sf_bessel_Y0(sigma * u);
    peta = H * u * gsl_sf_bessel_J1(sigma * u) + ka * gsl_sf_bessel_J0(sigma * u);
    a = u * rmax * gsl_sf_bessel_J1(rmax * u) - u * r * gsl_sf_bessel_J1(r * u);
    b = u * rmax * gsl_sf_bessel_Y1(rmax * u) - u * r * gsl_sf_bessel_Y1(r * u);
    double out = 1.0 / pow(u, 3.0) * (exp(-D * u * u * h) - 1.0) / (alpha * alpha + peta * peta) * (alpha * a - peta * b);
    return out;
}

// the block-distance
double integral_for_blockdistance2D(paramsIL& parameters2D)
{
    paramsIL params = parameters2D;

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1e6);
    double result, error;
    double eps1 = 1.0e-5;
    double eps2 = eps1;
    gsl_function F;
    F.function = &function2D;
    F.params = &params;
    gsl_set_error_handler_off();
    int status = gsl_integration_qagiu(&F, 0, eps1, eps2, 1000000, w, &result, &error);
    if (status != GSL_SUCCESS) {
        double u1 = 0;
        double u2 = 1.0e4;
        while (std::abs(function2D(u2, F.params)) > 1.0e-5) {
            u2 = u2 * 1.5;
        }
        while (status != GSL_SUCCESS) {
            status = gsl_integration_qags(&F, u1, u2, eps1, eps1, 1000000, w, &result, &error);
            u2 = u2 * 0.9;
        }
    }
    gsl_integration_workspace_free(w);
    gsl_set_error_handler(NULL);
    return result;
}

void block_distance(paramsIL& parameters2D)
{
    double kb = parameters2D.kb / 1.0e6;
    double sigma = parameters2D.sigma;
    double D = parameters2D.Dtot;
    double h = parameters2D.dt;
    double Rmax = sigma + 3.0 * sqrt(4.0 * D * h);
    double left = dissociate2D(parameters2D);
    double criterion = 1e-5;
    double rmin = sigma;
    double rmax = Rmax;
    double rmean, right;
    while (std::abs(rmax - rmin) > criterion) {
        rmean = 0.5 * (rmax + rmin);
        parameters2D.R2D = rmean;
        right = 4 * kb * integral_for_blockdistance2D(parameters2D);
        if (right > left) {
            rmin = rmean;
        } else {
            rmax = rmean;
        }
    }
    parameters2D.R2D = rmean;
    //std::cout<<left<<std::endl;
    //std::cout<<rmean<<", "<<Rmax<<std::cin.get();
}

// binding probability, but must time the lipid density
double pimplicitlipid_2D(paramsIL& parameters2D)
{
    double ka = parameters2D.ka;
    if (ka < 1E-15) {
        return 0.0;
    }
    block_distance(parameters2D);
    // std::cout<<parameters2D.R2D<<std::endl;
    paramsIL params = parameters2D;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1e6);
    double result, error;
    double eps1 = 1.0e-5;
    double eps2 = eps1;
    gsl_function F;
    F.function = &function2D;
    F.params = &params;
    gsl_set_error_handler_off();
    int status = gsl_integration_qagiu(&F, 0, eps1, eps2, 1000000, w, &result, &error);
    if (status != GSL_SUCCESS) {
        double u1 = 0;
        double u2 = 1.0e4;
        while (std::abs(function2D(u2, F.params)) > 1.0e-5) {
            u2 = u2 * 1.5;
        }
        while (status != GSL_SUCCESS) {
            status = gsl_integration_qags(&F, u1, u2, eps1, eps1, 1000000, w, &result, &error);
            u2 = u2 * 0.9;
        }
    }
    gsl_integration_workspace_free(w);
    gsl_set_error_handler(NULL);

    //double ka = parameters2D.ka;
    return result * 4 * ka;
}

///////////////////////////////////////////////////////////////
// 3D
double dissociate3D(double h, double D, double sigma, double ka, double kbsecond)
{
    //double h = parameters3D.dt;
    //double D = parameters3D.Dtot;
    //double sigma = parameters3D.sigma;
    //double ka = parameters3D.ka;
    //double kb = parameters3D.kb;
    double kb = kbsecond / 1.0e6; // change the unit S into us.
    if (kb < 1E-15) {
        return 0.0;
    }
    double KD = 2.0 * kb / ka;
    double kon = 0.5 / (1.0 / ka + 1.0 / (4.0 * M_PI * D * sigma));
    double koff = kon * KD;
    double out = 1.0 - exp(-koff * h);
    return out;
}

// binding probability, but must time the lipid density
double pimplicitlipid_3D(double z, paramsIL& parameters3D)
{
    double h = parameters3D.dt;
    double D = parameters3D.Dtot;
    double sigma = parameters3D.sigma;
    double ka = parameters3D.ka;
    if (ka < 1E-15) {
        return 0.0;
    }

    double out;
    if (z > sigma) {
        double alpha = sqrt(D) / sigma * (1.0 + ka / (4.0 * M_PI * sigma * D));
        double conf = 2.0 * M_PI * sigma * sigma * ka * (4.0 * M_PI * sigma * D) / (ka + 4.0 * M_PI * sigma * D) / (ka + 4.0 * M_PI * sigma * D);
        double a = (z - sigma) / sqrt(4.0 * D * h);
        double b = alpha * sqrt(h);
        if (std::isinf(exp(2.0 * a * b + b * b))) {
            out = conf * (exp(-a * a) / sqrt(M_PI) / (a + b) - (2.0 * a * b + 1.0) * erfc(a) + 2.0 * alpha * sqrt(h / M_PI) * exp(-a * a));
        } else {
            out = conf * (exp(2.0 * a * b + b * b) * erfc(a + b) - (2.0 * a * b + 1.0) * erfc(a) + 2.0 * alpha * sqrt(h / M_PI) * exp(-a * a));
        }
    } else {
        z = sigma;
        double alpha = sqrt(D) / sigma * (1.0 + ka / (4.0 * M_PI * sigma * D));
        double conf = 2.0 * M_PI * sigma * ka * sqrt(D) / alpha / (ka + 4.0 * M_PI * sigma * D);
        double a = 0;
        double b = alpha * sqrt(h);
        if (std::isinf(exp(b * b))) {
            out = conf * (1.0 / sqrt(M_PI) / b - 1.0 + 2.0 * alpha * sqrt(h / M_PI));
        } else {
            out = conf * (exp(b * b) * erfc(b) - 1.0 + 2.0 * alpha * sqrt(h / M_PI));
        }
    }
    return out;
}
