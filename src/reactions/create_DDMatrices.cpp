#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

void create_DDMatrices(gsl_matrix*& survMatrix, gsl_matrix*& normMatrix, gsl_matrix*& pirMatrix, double bindRadius,
    double Dtot, double comRMax, double ktemp, const Parameters& params)
{
    //    char fname[100];
    double RStepSize { sqrt(Dtot * params.timeStep) / 50 };

    create_normMatrix(normMatrix, bindRadius, Dtot, ktemp, comRMax, RStepSize, params);
    // std::cout << "Norm matrix successfully created for Dtot: " << Dtot << " timeStep "
    //           << ", bindRadius: " << bindRadius << '\n';

    create_survMatrix(survMatrix, bindRadius, Dtot, ktemp, comRMax, RStepSize, params);
    // std::cout << "Survival matrix successfully created for Dtot: " << Dtot << " timeStep "
    //           << ", bindRadius: " << bindRadius << '\n';

    create_pirMatrix(pirMatrix, bindRadius, Dtot, ktemp, comRMax, RStepSize, params);
    // std::cout << "Pir matrix successfully created for Dtot: " << Dtot << " timeStep "
    //           << ", bindRadius: " << bindRadius << '\n';
}
