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
  int nCancelOverlapPartner {0}; //!< Counter of number of association events rejected by overlap, of binders with partners in complex
  int nCancelOverlapSystem {0}; //!< Counter of number of association events rejected by overlap, of binders with non-bindiners
  int nCancelDisplace2D {0}; //!< Counter of number of association events rejected by displacement.
  int nCancelDisplace3D {0}; //!< Counter of number of association events rejected by displacement.
  int nCancelDisplace3Dto2D {0}; //!< Counter of number of association events rejected by displacement.
  int nCancelSpanBox {0}; //!< Counter of number of association events rejected by spanning the box.
  int nAssocSuccess {0}; //!< Counter of number of association events successful
    std::vector<int> singleDouble; //is the species a single molecule, or a bound complex?
    std::vector<bool> implicitDouble; // is the bound complex an implicit?
    std::vector<bool> canDissociate; // can the specie dissociate?
    std::vector<std::vector<int>> bindPairList; // one mol index of the specie (the first forward reactant), used for the dissociation pool
    std::vector<std::vector<int>> bindPairListIL2D; // one mol index of the specie, used for the dissociation pool, 2D implicit lipid bound, complex.linksToSurface > 1
    std::vector<std::vector<int>> bindPairListIL3D; // one mol index of the specie, used for the dissociation pool, 3D implicit lipid bound, complex.linksToSurface == 1
  int eventArraySize {20};
  std::vector<int> events3D;
  std::vector<int> events2D;
  std::vector<int> events3Dto2D;
};
