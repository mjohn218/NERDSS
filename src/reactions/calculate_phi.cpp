#include "reactions/association/association.hpp"
#include "tracing.hpp"

double calculate_phi(Coord reactIface1, int ifaceIndex2, Molecule reactMol1, Molecule reactMol2, const Vector& normal,
    Vector axis, const ForwardRxn& currRxn, const std::vector<MolTemplate>& molTemplateList)
{
    // TRACE();
    // coordinate transform along com-iface vector
    transform(reactIface1, reactMol1, reactMol2, axis);

    // orthographic projection onto xy-plane
    Vector vec1 { reactIface1 - reactMol2.tmpICoords[ifaceIndex2] }; //iface1-iface2= sigma
    Vector vec2 { determine_normal(normal, molTemplateList[reactMol1.molTypeIndex], reactMol1) };

    // remove z coordinates
    Vector projVec1 { vec1.x, vec1.y, 0 }; //sigma vector
    Vector projVec2 { vec2.x, vec2.y, 0 }; //normal vector of reactMol1
    projVec1.calc_magnitude();
    projVec2.calc_magnitude();
    double phi = projVec1.dot_theta(projVec2);

    /*angle between sigma and normal is phi. 
      as a dihedral, ranges from -pi to pi, 
      since phi here is between 0 and pi, use defined convention to decide sign flip.
    */
    //  if (angleSignIsCorrect(projVec1, projVec2))
    //         return phi;
    //     else
    //         return -phi;
    /*if their normal points in -z, keep theta, otherwise flip sign                                                                                                                                             
    ADDED: Do not flip the sign if it is PI or Zero.                                                                                                                                                          
    */
    Vector test { projVec1.cross(projVec2) };
    //double *test=new double[3];
    //crossproduct(sigma, normal, test);
    if (test.z > 0 && std::abs(phi) > 1E-12 && (M_PI - std::abs(phi)) > 1E-12) //positive z, flip phi
        phi = -phi;

    return phi;
}
