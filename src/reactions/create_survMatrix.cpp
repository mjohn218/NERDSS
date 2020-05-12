#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "tracing.hpp"

void create_survMatrix(gsl_matrix*& survMatrix, double bindRadius, double Dtot, double kr, double comRMax,
    double RStepSize, const Parameters& params)
{
    // TRACE();
    size_t ctr { 0 };
    //    double result; //, RstepSize = sqrt(Dtot * deltat) / 50;
    //	double Rmax = 3.0 * sqrt(4.0 * Dtot * deltat) + bindrad;
    char funcID[] = "tblsur";

    gsl_function F;
    F.function = &survival_function;
    IntegrandParams intParams;
    intParams.a = bindRadius;
    intParams.D = Dtot;
    intParams.k = kr;
    intParams.t = params.timeStep;

    gsl_set_error_handler_off();
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1e6);

    double RIndex { intParams.a };
    while (RIndex <= comRMax + RStepSize) {
        //    for (double RIndex = intParams.a; RIndex <= com + RStepSize; RIndex += RStepSize) {
        intParams.r0 = RIndex;
        gsl_matrix_set(survMatrix, 0, ctr, RIndex);
        F.params = reinterpret_cast<void*>(&intParams);

        if (kr < 1.0 / 0.0) {
            double result { integrator(
                F, intParams, w, RIndex, bindRadius, Dtot, kr, params.timeStep, funcID, survival_function) };
            gsl_matrix_set(survMatrix, 1, ctr, result);
        } else {
            double result { integrator(
                F, intParams, w, RIndex, bindRadius, Dtot, kr, params.timeStep, funcID, survival_function) };
            gsl_matrix_set(survMatrix, 1, ctr, 1.0 - result);
        }
        ctr = ctr + 1;
        RIndex += RStepSize;
    }
    gsl_integration_workspace_free(w);
    gsl_set_error_handler(nullptr);
}
