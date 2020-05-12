#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "tracing.hpp"

void create_pirMatrix(gsl_matrix*& pirMatrix, double bindRadius, double Dtot, double kr, double comRMax,
    double RStepSize, const Parameters& params)
{
    // TRACE();
    int itr1 { 0 };
    int itr2 { 0 };
    double result;
    // double RStepSize = sqrt(Dtot * deltat) / 50;
    // double Rmax = 3.0 * sqrt(4.0 * Dtot * deltat) + bindRadius;
    char funcID[] = "tblpirr";

    gsl_function F;
    F.function = &pir_function;
    IntegrandParams intParams;
    intParams.a = bindRadius;
    intParams.D = Dtot;
    intParams.k = kr;
    intParams.t = params.timeStep;

    gsl_set_error_handler_off();
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1e6);

    //    for (double Rindex = intParams.a; Rindex <= Rmax + RStepSize; Rindex += RStepSize) {
    double RIndex { intParams.a };
    while (RIndex <= comRMax + RStepSize) {
        intParams.r = RIndex;
        double R0Index { intParams.a };
        //        for (double R0index = intParams.a; R0index <= Rmax + RStepSize; R0index += RStepSize) {
        while (R0Index <= comRMax + RStepSize) {
            intParams.r0 = R0Index;
            F.params = reinterpret_cast<void*>(&intParams);
            result = integrator(F, intParams, w, R0Index, bindRadius, Dtot, kr, params.timeStep, funcID, pir_function);
            gsl_matrix_set(pirMatrix, itr1, itr2, result);
            itr1 += 1;
            R0Index += RStepSize;
        }
        itr1 = 0;
        itr2 += 1;
        RIndex += RStepSize;
    }
    gsl_integration_workspace_free(w);
    gsl_set_error_handler(nullptr);
}
