#include "gsl/gsl_rng.h"
#include "math/rand_gsl.hpp"
#include <cmath>
#include <iostream>

//static gsl_rng* the_generator = nullptr;

double rand_gsl()
{
    //++randNum;
    //if (!the_generator)
    //    the_generator = gsl_rng_alloc(gsl_rng_taus);
    return gsl_rng_uniform(r) + 1.0 / (gsl_rng_max(r) + 1.0) * gsl_rng_uniform(r);
}

void srand_gsl(int num)
{
    //if (!the_generator)
    //    the_generator = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(r, num);
}

void write_rng_state()
{
    FILE* stateOut = fopen("rng_state", "w");
    if (ferror(stateOut)) {
        std::cerr << "ERROR: Could not open RNG state file for writing. Exiting.\n";
        exit(1);
    }

    int stateWriteState = gsl_rng_fwrite(stateOut, r);
    if (stateWriteState == GSL_EFAILED) {
        std::cerr << "ERROR: Could not write RNG state file. Exiting.\n";
        exit(1);
    }
    fclose(stateOut);
}

void write_rng_state_simItr(int simItr)
{
    char fnameProXYZ[100];
    sprintf(fnameProXYZ, "rng_state%d", simItr);
    FILE* stateOut = fopen(fnameProXYZ, "w");
    if (ferror(stateOut)) {
        std::cerr << "ERROR: Could not open RNG state file for writing. Exiting.\n";
        exit(1);
    }

    int stateWriteState = gsl_rng_fwrite(stateOut, r);
    if (stateWriteState == GSL_EFAILED) {
        std::cerr << "ERROR: Could not write RNG state file. Exiting.\n";
        exit(1);
    }
    fclose(stateOut);
}

void read_rng_state()
{
    //std::cout << "Reading RNG state file.\n";
    FILE* stateIn = fopen("rng_state", "r");
    if (stateIn == nullptr || ferror(stateIn)) {
        std::cerr << "Could not find RNG state file, initializing new RNG..\n";
        fclose(stateIn);
        return;
    }

    int stateReadStatus = gsl_rng_fread(stateIn, r);
    if (stateReadStatus == GSL_EFAILED) {
        std::cerr << "Could not read RNG state file, initializing new RNG..\n";
        fclose(stateIn);
        return;
    }
    fclose(stateIn);
}

double GaussV()
{
    double R { 2.0 };
    double V1 {};

    while (R >= 1.0) {
        V1 = 2.0 * gsl_rng_uniform(r) - 1.0;
        double V2 = 2.0 * gsl_rng_uniform(r) - 1.0;
        R = (V1 * V1) + (V2 * V2);
    }
    return (V1 * sqrt(-2.0 * log(R) / R));
}
