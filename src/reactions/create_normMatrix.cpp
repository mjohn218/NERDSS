#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "tracing.hpp"

void create_normMatrix(gsl_matrix*& normMatrix, double bindRadius, double Dtot, double kr, double comRMax,
    double RStepSize, const Parameters& params)
{
    // TRACE();
    //    double result, error;
    const double xLowB { bindRadius };
    int itr { 0 };
    const double epsAbs { 1e-6 };
    const double epsRel { 1e-6 };

    gsl_function F;
    F.function = &norm_function;
    IntegrandParams intParams;
    intParams.a = bindRadius;
    intParams.D = Dtot;
    intParams.k = kr;
    intParams.t = params.timeStep;

    double RIndex { intParams.a };
    //    for (double Rindex = intParams.a; Rindex <= Rmax + RStepSize; Rindex += RStepSize) {
    while (RIndex <= comRMax + RStepSize) {
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(10000000);
        intParams.r0 = RIndex;
        gsl_matrix_set(normMatrix, 0, itr, RIndex);
        F.params = reinterpret_cast<void*>(&intParams);

        double result {};
        double error {};
        gsl_integration_qagiu(&F, xLowB, epsAbs, epsRel, 10000000, w, &result, &error);
        gsl_matrix_set(normMatrix, 1, itr, result);
        gsl_integration_workspace_free(w);
        ++itr;
        RIndex += RStepSize;
    }
    //    std::cout << "ELEMENT 318: " << gsl_matrix_get(normMatrix, 0, 318) << ' ' << gsl_matrix_get(normMatrix, 1, 318)
    //              << '\n';
    //    std::cout << "ELEMENT 319: " << gsl_matrix_get(normMatrix, 0, 319) << ' ' << gsl_matrix_get(normMatrix, 1, 319)
    //              << '\n';
}
