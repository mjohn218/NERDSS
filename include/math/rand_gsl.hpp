#pragma once

#include "gsl/gsl_rng.h"

extern gsl_rng* r;
extern long long randNum;

/*!
 * \brief Uses the previously initialized random number generator to return a random number
 * \param[out] double Uniformly distributed random double.
 */
double rand_gsl();

/*!
 * \brief Initializes the GSL random number generator.
 */
void srand_gsl(int);

/*!
 * \brief Wrapper for the internal GSL RNG state read function.
 *
 * Reads a previously written binary file with the current status of the RNG, so restarting will give the same random
 * numbers as a continuous run would.
 */
void read_rng_state();

/*!
 * \brief Wrapper for the internal GSL RNG state write function.
 *
 * Writes the current state of the RNG to a binary file.
 */
void write_rng_state();

/*!
 * \brief Wrapper for the internal GSL RNG state write function.
 *
 * Writes the current state of the RNG to a binary file for check point.
 */
void write_rng_state_simItr(int simItr);

/*!
 * \brief Uses Box-Mueller method to greate Gaussian-distributed random numbers from a uniform random number generator.
 * \param[out] double Gaussian-distributed random double.
 */
double GaussV();
