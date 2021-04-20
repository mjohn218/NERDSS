/*! \file class_Membrane.hpp
 * 
 */

#pragma once

#include "classes/class_Membrane.hpp"

#include <string>
#include <vector>

/*! \enum BoundaryKeywords
 * \ingroup Parser
 * \brief Boundary parameters read in from the command line
 */
enum class BoundaryKeyword : int {
    implicitLipid = 0, //!< use implicit lipid model
    waterBox = 1, //!< use a rectangular box, specify x y and z lengths
    xBCtype = 2, //!< reflecting or periodic?
    yBCtype = 3, //!< reflecting or periodic?
    zBCtype = 4, //!< reflecting or periodic?
    isSphere = 5, //!< use a sphere boundary, if provide sphereR, this will be set to true.
    sphereR = 6, //!< sphere radius
};

struct Membrane {
    //public:
    struct WaterBox {
        /*!
         * \brief Just a container for the water box dimensions
         * Only cubic at the moment. Not a Coord because then it'll be a circular include (since Coord needs Parameters)
         */

        double x { 0 };
        double y { 0 };
        double z { 0 };
        double volume { 0 };
        WaterBox() = default;
        explicit WaterBox(std::vector<double> vals)
            : x(vals[0])
            , y(vals[1])
            , z(vals[2])
        {
            volume = x * y * z;
        }
    };
    WaterBox waterBox; //!< water box x, y, z. used to be xboxl, yboxl, zboxl
    double sphereR = 0; //!< for sphere, value of radius in nm.
    double sphereVol = 0;
    int nSites;
    int nStates { 0 }; // number of the states of implict lipid
    int No_free_lipids;
    std::vector<int> numberOfFreeLipidsEachState {}; // record the free lipids of each state for IL, updated each step in main function
    int No_protein; // use for implicit-lipid model;
    std::vector<int> numberOfProteinEachState {}; // record the number of proteins that can bound to each state for IL
    int implicitlipidIndex { -1 };
    std::vector<double> RS3Dvect; //this is the look-up table for RS3D, which is the reflecting-surface for 3D-->2D reaction of implicit-lipid case

    //    double RD2D = 0; // block-distance for 2D->2D reaction of implicit-lipid case
    double totalSA;
    double Dx;
    double Dy;
    double Dz;
    double Drx;
    double Dry;
    double Drz;
    double offset;
    double lipidLength { 0.0 };
    bool implicitLipid = false;
    bool TwoD = false;
    bool isBox = false;
    bool isSphere = false;
    std::string xBCtype; //allow reflect, or pbc
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

    void display(); // display the information for the boundary, define in the src/parse/parse_input.cpp

    void create_water_box(); // create box for sphere boundary, define in the src/parse/parse_input.cpp
};
