/*! \file debug.hpp

 * ### Created on 10/18/18 by Matthew Varga
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

#include <cstring>

#include "classes/class_Molecule_Complex.hpp"
#include "split.cpp"

using namespace std;

/*! \defgroup error
 * \brief Functions for detecting errors.
 */

/*! \ingroup error
 * \brief print errors and exit
 */
void error(string errorString);

void error(MpiContext &mpiContext, string errorString);

void error(MpiContext &mpiContext, Molecule &mol, string errorString);
