#include "reactions/association/association.hpp"
#include "tracing.hpp"
/*What is the angular displacement between the site-COM of each associating molecule
  calculate the dot-product between the two vectors, the tmpICoords and the Coords.
*/
/*What is the angular displacement between the site-complexCOM of each associating molecule*/

void calc_angular_displacement(int ifaceIndex1, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2,
    Complex& reactCom1, Complex& reactCom2, std::vector<Molecule>& moleculeList)
{

    Vector v1 { reactMol1.tmpICoords[ifaceIndex1] - reactMol1.tmpComCoord }; //temporary vector due to association moves.
    Vector v2 { reactMol1.interfaceList[ifaceIndex1].coord - reactMol1.comCoord }; //original vector
    v1.calc_magnitude();
    if (v1.magnitude < 1E-12) {
        // std::cout << "No rotation for a POINT particle \n";
        return;
    }
    v2.calc_magnitude();

    double currTheta = v2.dot_theta(v1);
    // std::cout << "Angle between Molecule1 interface to COM, temporary vs original: " << currTheta << std::endl;

    /*molecule 2*/
    Vector v3 { reactMol2.tmpICoords[ifaceIndex2] - reactMol2.tmpComCoord }; //temporary vector due to association moves.
    Vector v4 { reactMol2.interfaceList[ifaceIndex2].coord - reactMol2.comCoord }; //original vector
    v3.calc_magnitude();
    if (v3.magnitude < 1E-12) {
        // std::cout << "No rotation for a POINT particle \n";
        return;
    }
    v4.calc_magnitude();

    currTheta = v4.dot_theta(v3);
    // std::cout << "Angle between Molecule2 interface to COM, temporary vs original: " << currTheta << std::endl;

    /*calculate angle between interface and the whole complex COM.*/
    Vector v5 { reactMol1.tmpICoords[ifaceIndex1] - reactCom1.tmpComCoord }; //temporary vector due to association moves.
    Vector v6 { reactMol1.interfaceList[ifaceIndex1].coord - reactCom1.comCoord }; //original vector
    v5.calc_magnitude();
    if (v5.magnitude < 1E-12) {
        // std::cout << "No rotation for a POINT particle \n";
        return;
    }
    v6.calc_magnitude();

    currTheta = v6.dot_theta(v5);
    // std::cout << "Angle between Molecule1 interface to ComplexCOM, temporary vs original: " << currTheta << std::endl;

    Vector v7 { reactMol2.tmpICoords[ifaceIndex2] - reactCom2.tmpComCoord }; //temporary vector due to association moves.
    Vector v8 { reactMol2.interfaceList[ifaceIndex2].coord - reactCom2.comCoord }; //original vector
    v7.calc_magnitude();
    if (v7.magnitude < 1E-12) {
        // std::cout << "No rotation for a POINT particle \n";
        return;
    }
    v8.calc_magnitude();

    currTheta = v8.dot_theta(v7);
    // std::cout << "Angle between Molecule2 interface to ComplexCOM, temporary vs original: " << currTheta << std::endl;
}
