/*! \file class_Membrane.hpp
 *
 */

#pragma once

#include <string>
#include <vector>

#include "classes/mpi_functions.hpp"

#include "classes/class_copyCounters.hpp"

class copyCounters {
 public:
  std::vector<int>
      nBoundPairs;  // How many of each type of pair of reactive Molecules (sums
                    // over all interfaces) have links to each other.
  std::vector<int> proPairlist;  // index of pair in the array.
  std::vector<int>
      copyNumSpecies;  // array keeping track of counts of all species.
  int nLoops{0};       //!< Counter of closed loops formed (e.g. hexagons).
  int nCancelOverlapPartner{
      0};  //!< Counter of number of association events rejected by overlap, of
           //!< binders with partners in complex
  int nCancelOverlapSystem{
      0};  //!< Counter of number of association events rejected by overlap, of
           //!< binders with non-bindiners
  int nCancelDisplace2D{0};      //!< Counter of number of association events
                                 //!< rejected by displacement.
  int nCancelDisplace3D{0};      //!< Counter of number of association events
                                 //!< rejected by displacement.
  int nCancelDisplace3Dto2D{0};  //!< Counter of number of association events
                                 //!< rejected by displacement.
  int nCancelSpanBox{0};  //!< Counter of number of association events rejected
                          //!< by spanning the box.
  int nAssocSuccess{0};  //!< Counter of number of association events successful
  std::vector<int>
      singleDouble;  // is the species a single molecule, or a bound complex?
  std::vector<bool> implicitDouble;  // is the bound complex an implicit?
  std::vector<bool> canDissociate;   // can the specie dissociate?

  int eventArraySize{20};
  std::vector<int> events3D;
  std::vector<int> events2D;
  std::vector<int> events3Dto2D;

  /*
  Function serialize serializes copyCounters
  into array of bytes.
  */

  void temp_serialize_primitive_vector(std::vector<bool> to_serialize,
                                       unsigned char *arrayRank,
                                       int &nArrayRank) {
    // serialize to_serialize vector of primitive type bool
    // Store vector size in first bytes needed for an int
    // and set current number of bytes read to the size of int:
    *((int *)(arrayRank + nArrayRank)) = to_serialize.size();
    nArrayRank += sizeof(int);
    for (auto it : to_serialize) {  // for all vector elements
      *((bool *)&(arrayRank[nArrayRank])) =
          it;                      // read one primitive type bool, e.g. int
      nArrayRank += sizeof(bool);  // increase number of bytes read
    }
  }
  void serialize(unsigned char *arrayRank, int &nArrayRank) {
    serialize_primitive_vector<int>(nBoundPairs, arrayRank, nArrayRank);
    serialize_primitive_vector<int>(proPairlist, arrayRank, nArrayRank);
    serialize_primitive_vector<int>(copyNumSpecies, arrayRank, nArrayRank);
    PUSH(nLoops);
    PUSH(nCancelOverlapPartner);
    PUSH(nCancelOverlapSystem);
    PUSH(nCancelDisplace2D);
    PUSH(nCancelDisplace3D);
    PUSH(nCancelDisplace3Dto2D);
    PUSH(nCancelSpanBox);
    PUSH(nAssocSuccess);
    serialize_primitive_vector<int>(singleDouble, arrayRank, nArrayRank);
    temp_serialize_primitive_vector(implicitDouble, arrayRank, nArrayRank);
    temp_serialize_primitive_vector(canDissociate, arrayRank, nArrayRank);
    //        serialize_primitive_vector<bool>(implicitDouble, arrayRank+start);
    //        serialize_primitive_vector<bool>(canDissociate, arrayRank+start);
    PUSH(eventArraySize);
    serialize_primitive_vector<int>(events3D, arrayRank, nArrayRank);
    serialize_primitive_vector<int>(events2D, arrayRank, nArrayRank);
    serialize_primitive_vector<int>(events3Dto2D, arrayRank, nArrayRank);
  }
  /*
  Function deserialize deserializes copyCounters
  from array of bytes.
  */
  void deserialize(unsigned char *arrayRank, int &nArrayRank) {
    deserialize_primitive_vector<int>(nBoundPairs, arrayRank, nArrayRank);
    deserialize_primitive_vector<int>(proPairlist, arrayRank, nArrayRank);
    deserialize_primitive_vector<int>(copyNumSpecies, arrayRank, nArrayRank);
    POP(nLoops);
    POP(nCancelOverlapPartner);
    POP(nCancelOverlapSystem);
    POP(nCancelDisplace2D);
    POP(nCancelDisplace3D);
    POP(nCancelDisplace3Dto2D);
    POP(nCancelSpanBox);
    POP(nAssocSuccess);
    deserialize_primitive_vector<int>(singleDouble, arrayRank, nArrayRank);
    deserialize_primitive_vector<bool>(implicitDouble, arrayRank, nArrayRank);
    deserialize_primitive_vector<bool>(canDissociate, arrayRank, nArrayRank);
    POP(eventArraySize);
    deserialize_primitive_vector<int>(events3D, arrayRank, nArrayRank);
    deserialize_primitive_vector<int>(events2D, arrayRank, nArrayRank);
    deserialize_primitive_vector<int>(events3Dto2D, arrayRank, nArrayRank);
  }
};
