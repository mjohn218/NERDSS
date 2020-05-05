/*! \file check_if_spans.cpp
 * ### Created on 2020-02-23 by Yiben
 */
#include "boundary_conditions/reflect_functions.hpp"

#include <iostream>

void check_if_spans(bool& cancelAssoc, const Parameters& params, Complex& reactCom1, Complex& reactCom2,
    std::vector<Molecule>& moleculeList, const Membrane& membraneObject)
{
    // Associating proteins have been moved to contact. Before assigning them to the complexsame complex,
    // test to see if the complex is too big to fit in the box.
    if (membraneObject.isSphere == true)
        check_if_spans_sphere(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);
    else
        check_if_spans_box(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);
}