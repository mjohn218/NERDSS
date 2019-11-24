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
    sphere  = 5, //!< use a sphere, specify R.
};

struct Membrane{
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
    int nSites;
    int No_free_lipids;
    int No_protein; // use for implicit-lipid model;
  int implicitlipidIndex {-1};
    double RS3D = 0; // reflecting-surface for 3D->2D reaction of implicit-lipid case
    double Punbinding3D = 0; // unbinding probability for 2D particle unbinding from implicit-lipid surface
    double RD2D = 0; // block-distance for 2D->2D reaction of implicit-lipid case 
    double Pbinding2D = 0; // binding probability for 2D implicit-lipid surface.
    double Punbinding2D = 0; // unbinding probability for 2D implicit-lipid surface.
    double totalSA;
    double Dx;
    double Dy;
    double Dz;
    double Drx;
    double Dry;
    double Drz;
    double offset;
    bool implicitLipid = false;
    bool TwoD = false;
    bool isBox=false;
    bool isSphere=false;
    double sphereRadius=0;
    std::string xBCtype;//allow reflect, or pbc
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
};
