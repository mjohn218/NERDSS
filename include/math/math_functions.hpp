/*! \file math_functions.hpp
 * ### Created on 2019-02-06 by Matthew Varga
 */

#pragma once

namespace MathFuncs {
/*!
     * \brief Computes the factorial
     */
long double factorial(unsigned n);

// Following from Numerical Recipes, ch. 6
double gammln(double n);

double gammFactorial(int n);
};