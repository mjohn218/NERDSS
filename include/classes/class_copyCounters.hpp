/*! \file class_Membrane.hpp
 *
 */

#pragma once

#include "classes/class_copyCounters.hpp"
#include <string>
#include <vector>

class copyCounters {
public:
    std::vector<int> nBoundPairs; // How many of each type of pair of reactive Molecules (sums over all interfaces) have
                                  // links to each other.
  std::vector<int> proPairlist; // index of pair in the array.
  std::vector<int> copyNumSpecies; // array keeping track of counts of all species.
  int nLoops { 0 }; //!< Counter of closed loops formed (e.g. hexagons).
  std::vector<int> singleDouble;//is the species a single molecule, or a bound complex?
};
