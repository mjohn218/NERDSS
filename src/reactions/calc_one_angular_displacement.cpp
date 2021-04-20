#include "reactions/association/association.hpp"
#include "tracing.hpp"
/*What is the angular displacement between the site-COM of an associating molecule
  calculate the dot-product between the two vectors, the tmpICoords and the Coords.
*/

double calc_one_angular_displacement(int ifaceIndex1, Molecule& reactMol1, Complex& reactCom1)
{

    double zero = 0.0;
    Vector v1 { reactMol1.tmpICoords[ifaceIndex1] - reactMol1.tmpComCoord }; //temporary vector due to association moves.

    Vector v2 { reactMol1.interfaceList[ifaceIndex1].coord - reactMol1.comCoord }; //original vector
    v1.calc_magnitude();
    if (v1.magnitude < 1E-12) {
        // std::cout << "No rotation for a POINT particle \n";
        return zero;
    }
    v2.calc_magnitude();

    double currTheta = v2.dot_theta(v1);
    Vector test { v2.cross(v1) }; //original.cross.final
    if (test.z > 0 && std::abs(currTheta) > 1E-12 && (M_PI - std::abs(currTheta)) > 1E-12) //positive z, flip currTheta
        currTheta = -currTheta;

    // std::cout << "Angle between Molecule1 interface to COM, temporary vs original: " << currTheta << std::endl;
    return currTheta;
}
