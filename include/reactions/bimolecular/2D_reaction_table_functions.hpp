/*! \file 2D_reaction_table_functions.hpp

 * ### Created on 11/12/18 by Matthew Varga
 * ### Purpose
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */
#pragma once

#include "classes/class_Molecule_Complex.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>

#include <chrono>

/*! \defgroup 2DReactions
 * \brief
 */

/*! \ingroup 2DReactions
 * \brief Parameters for integration during 2D matrix creation
 */
struct IntegrandParams { // Structure for integrands
    double a{}; //!< sigma
    double D{}; //!< Diff. Coeff
    double k{}; //!< rate of association
    double r0{}; //!<  initial loc.
    double r{}; //!<  initial radius
    double t{}; //!<  init. t.
    double rho{}; //!<  current lipid density

    IntegrandParams() = default;
    IntegrandParams(double _a, double _D, double _k, double _t)
        : a(_a)
        , D(_D)
        , k(_k)
        , t(_t)
    {
    }
};

size_t size_lookup(double bindRadius, double Dtot, const Parameters& params, double Rmax);

/*! \ingroup 2DReactions
 * \brief Main function for 2D reaction tables
 */
void create_DDMatrices(gsl_matrix *& survMatrix, gsl_matrix *& normMatrix, gsl_matrix *& pirMatrix, double bindRadius,
                       double Dtot, double comRMax, double ktemp, const Parameters& params);

/* MATH FUNCTIONS */
/*! \ingroupd 2DReactions
 * \brief Survival function for integration.
 */
double survival_function(double x, void* p);

/*! \ingroup 2DReactions
 * \brief Norm function for integration
 */
double norm_function(double x, void* p);

/*! \ingroup 2DReactions
 * \brief Pir function for integration
 */
double pir_function(double x, void* p);

/*! \ingroup 2DReactions
 * \brief 2D lookup table for pnorm
 */
double get_prevNorm(gsl_matrix *normMatrix, double RStepSize, double r0, double bindRadius);

/*! \ingroup 2DReactions
 * \brief 2D lookup table for psurv
 */
double get_prevSurv(const gsl_matrix* survMatrix, double Dtot, double deltaT, double r0, double bindRadius);

/*! \ingroupd 2DReactions
 * \brief Calculates the ratio of pirr and pfree
 */
double DDpirr_pfree_ratio_ps(gsl_matrix *pirMatrix, gsl_matrix *survMatrix, gsl_matrix *normMatrix,
                             double r, double Dtot, double deltaT, double r0, double ps_prev, double rTol,
                             double bindRadius);

double pirr_pfree_ratio_psF(
    double rCurr, double r0, double tCurr, double Dtot, double bindrad, double alpha, double ps_prev, double rtol);

/*!
 * \brief 2D lookup table for pirr
 */
double calc_pirr(
    gsl_matrix *pirMatrix, gsl_matrix *survMatrix, double RStepSize, double r, double r0, double a);

/* INTEGRATION FUNCTIONS */
/*! \ingroupd 2DReactions
 * \brief Integrator for the functions involved in the creation of the 2D reaxtion matrix creation
 */
double integrator(gsl_function F, IntegrandParams params, gsl_integration_workspace* w, double r0, double bindrad,
    double Dtot, double kr, double deltat, char* funcID, double (*f)(double, void*));

/* 2D MATRIX CREATION FUNCTIONS */
/*! \inground 2DReactions
 * \brief Function to create the survival matrix.
 */
void create_survMatrix(gsl_matrix *& survMatrix, double bindRadius, double Dtot, double kr, double comRMax,
                       double RStepSize, const Parameters& params);

/*! \ingroup 2DReactions
 * \brief Function to create the norm matrix.
 */
void create_normMatrix(gsl_matrix *& normMatrix, double bindRadius, double Dtot, double kr, double comRMax,
                       double RStepSize, const Parameters& params);

/*! \ingroup 2DReactions
 * \brief  Function to create the pir matrix.
 */
void create_pirMatrix(gsl_matrix *& pirMatrix, double bindRadius, double Dtot, double kr, double comRMax,
                      double RStepSize, const Parameters& params);
