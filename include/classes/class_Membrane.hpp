/*! \file class_Membrane.hpp
 *
 */

#pragma once

#include <string>
#include <vector>

#include "classes/mpi_functions.hpp"

/*! \enum BoundaryKeywords
 * \ingroup Parser
 * \brief Boundary parameters read in from the command line
 */
enum class BoundaryKeyword : int {
  implicitLipid = 0,  //!< use implicit lipid model
  waterBox = 1,       //!< use a rectangular box, specify x y and z lengths
  xBCtype = 2,        //!< reflecting or periodic?
  yBCtype = 3,        //!< reflecting or periodic?
  zBCtype = 4,        //!< reflecting or periodic?
  isSphere = 5,  //!< use a sphere boundary, if provide sphereR, this will be
                 //!< set to true.
  sphereR = 6,   //!< sphere radius
};

struct Membrane {
  // public:
  struct WaterBox {
    /*!
     * \brief Just a container for the water box dimensions
     * Only cubic at the moment. Not a Coord because then it'll be a circular
     * include (since Coord needs Parameters)
     */

    double x{0};
    double y{0};
    double z{0};
    double xLeft{0.0};   // left bound of the rank
    double xRight{0.0};  // right bound of the rank
    double volume{0};
    WaterBox() = default;
    explicit WaterBox(std::vector<double> vals)
        : x(vals[0]), y(vals[1]), z(vals[2]) {
      volume = x * y * z;
      xLeft = -x / 2.0;
      xRight = x / 2.0;
    }
    /*
    Function serialize serializes the WaterBox
    into array of bytes.
    */
    void serialize(unsigned char *arrayRank, int &nArrayRank) {
      PUSH(x);
      PUSH(y);
      PUSH(z);
      PUSH(xLeft);
      PUSH(xRight);
      PUSH(volume);
    }
    /*
    Function deserialize deserializes the WaterBox
    from array of bytes.
    */
    void deserialize(unsigned char *arrayRank, int &nArrayRank) {
      POP(x);
      POP(y);
      POP(z);
      POP(xLeft);
      POP(xRight);
      POP(volume);
    }
  };
  WaterBox waterBox;   //!< water box x, y, z. used to be xboxl, yboxl, zboxl
  double sphereR = 0;  //!< for sphere, value of radius in nm.
  double sphereVol = 0;
  int nSites;
  int nStates{0};  // number of the states of implicit lipid
  int No_free_lipids;
  std::vector<int> numberOfFreeLipidsEachState{};  // record the free lipids of
                                                   // each state for IL, updated
                                                   // each step in main function
  int No_protein;  // use for implicit-lipid model;
  std::vector<int>
      numberOfProteinEachState{};  // record the number of proteins that can
                                   // bound to each state for IL
  int implicitlipidIndex{-1};
  std::vector<double> RS3Dvect;  // this is the look-up table for RS3D, which is
                                 // the reflecting-surface for 3D-->2D reaction
                                 // of implicit-lipid case

  //    double RD2D = 0; // block-distance for 2D->2D reaction of implicit-lipid
  //    case
  double totalSA;
  double Dx;
  double Dy;
  double Dz;
  double Drx;
  double Dry;
  double Drz;
  double offset;
  double lipidLength{0.0};
  bool implicitLipid = false;
  bool TwoD = false;
  bool isBox = false;
  bool isSphere = false;
  std::string xBCtype;  // allow reflect, or pbc
  std::string yBCtype;
  std::string zBCtype;

  /*set_value_BC is defined in src/parser/parse_input.cpp
    And the map to BoundaryKeyword keywords is also defined in that file.
    BoundaryKeyword keywords are defined above.
    ParameterKeywords are in include/classes/class_Parameters.hpp
   */

  void set_value_BC(std::string value, BoundaryKeyword keywords);
  /*In here, we could also store coordinate vector
    for a single representative lipid
  */

  void display();  // display the information for the boundary, define in the
                   // src/parse/parse_input.cpp

  void create_water_box();  // create box for sphere boundary, define in the
                            // src/parse/parse_input.cpp

  /*
  Function serialize serializes the Molecule
  into array of bytes.
  */
  void serialize(unsigned char *arrayRank, int &nArrayRank) {
    // std::cout << "+Membrane serialization starts here..." << std::endl;
    waterBox.serialize(arrayRank, nArrayRank);  // serialize starting
    // from beginning of arrayRank
    // increased by the number of bytes already serialized
    PUSH(sphereR);
    PUSH(sphereVol);
    PUSH(nSites);
    PUSH(nStates);
    PUSH(No_free_lipids);
    serialize_primitive_vector<int>(numberOfFreeLipidsEachState, arrayRank,
                                    nArrayRank);
    PUSH(No_protein);
    serialize_primitive_vector<int>(numberOfProteinEachState, arrayRank,
                                    nArrayRank);
    PUSH(implicitlipidIndex);
    serialize_primitive_vector<double>(RS3Dvect, arrayRank, nArrayRank);
    //    double RD2D = 0; // block-distance for 2D->2D reaction of
    //    implicit-lipid case
    PUSH(totalSA);
    PUSH(Dx);
    PUSH(Dy);
    PUSH(Dz);
    PUSH(Drx);
    PUSH(Dry);
    PUSH(Drz);
    PUSH(offset);
    PUSH(lipidLength);
    PUSH(implicitLipid);
    PUSH(TwoD);
    PUSH(isBox);
    PUSH(isSphere);
    serialize_string(xBCtype, arrayRank, nArrayRank);
    serialize_string(yBCtype, arrayRank, nArrayRank);
    serialize_string(zBCtype, arrayRank, nArrayRank);
    // std::cout << "+Total Membrane size in bytes: " << nArrayRank <<
    // std::endl;
  }
  void deserialize(unsigned char *arrayRank, int &nArrayRank) {
    waterBox.deserialize(arrayRank, nArrayRank);
    POP(sphereR);
    POP(sphereVol);
    POP(nSites);
    POP(nStates);
    POP(No_free_lipids);
    deserialize_primitive_vector<int>(numberOfFreeLipidsEachState, arrayRank,
                                      nArrayRank);
    POP(No_protein);
    deserialize_primitive_vector<int>(numberOfProteinEachState, arrayRank,
                                      nArrayRank);
    POP(implicitlipidIndex);
    deserialize_primitive_vector<double>(RS3Dvect, arrayRank, nArrayRank);
    //    double RD2D = 0; // block-distance for 2D->2D reaction of
    //    implicit-lipid case
    POP(totalSA);
    POP(Dx);
    POP(Dy);
    POP(Dz);
    POP(Drx);
    POP(Dry);
    POP(Drz);
    POP(offset);
    POP(lipidLength);
    POP(implicitLipid);
    POP(TwoD);
    POP(isBox);
    POP(isSphere);
    deserialize_string(xBCtype, arrayRank, nArrayRank);
    deserialize_string(yBCtype, arrayRank, nArrayRank);
    deserialize_string(zBCtype, arrayRank, nArrayRank);
  }
};
