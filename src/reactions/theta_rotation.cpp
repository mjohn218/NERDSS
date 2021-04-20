#include "reactions/association/association.hpp"
#include "tracing.hpp"

void theta_rotation(Coord& reactIface1, Coord& reactIface2, Molecule& reactMol1, Molecule& reactMol2, double targAngle,
    Complex& reactCom1, Complex& reactCom2, std::vector<Molecule>& moleculeList)
{
    // TRACE();
    Vector v1 { reactIface1 - reactMol1.tmpComCoord };
    Vector sigma { reactIface1 - reactIface2 };
    v1.calc_magnitude();
    if (v1.magnitude < 1E-12) {
        // std::cout << "No theta rotation for a POINT particle \n";
        return;
    }
    sigma.calc_magnitude();

    double currTheta { sigma.dot_theta(v1) };
    // std::cout << "Desired theta: " << targAngle << " Current theta: " << currTheta << std::endl;

    // Determine if we even need to rotate (i.e. if the theta is already aligned )
    if (std::abs(targAngle - currTheta) < 1E-8) {
        // std::cout << "No theta rotation needed" << std::endl;
    } else {

        // if sigma and v1 are parallel, rotation axis can't be determined by the cross product
        // so we need to find an arbitrary vector orthogonal to either. Choose either x or y axis,
        //whichever one is not parallel.
        Vector rotAxis = (areParallel(currTheta)) ? create_arbitrary_vector(v1) : sigma.cross(v1);
        rotAxis.normalize();

        // Determine how far to rotate each Molecule, based on their respective rotation diffusion constants
        double rotAngPos {};
        double rotAngNeg {};
        determine_rotation_angles(targAngle, currTheta, rotAngPos, rotAngNeg, reactCom1, reactCom2);

        // std::cout << "Positive half angle: " << rotAngPos << "\nNegative half angle: " << rotAngNeg << std::endl;

        // Create rotation quaternions
        Quat rotQuatPos { cos(rotAngPos / 2), sin(rotAngPos / 2) * rotAxis.x, sin(rotAngPos / 2) * rotAxis.y,
            sin(rotAngPos / 2) * rotAxis.z };
        Quat rotQuatNeg { cos(rotAngNeg / 2), sin(rotAngNeg / 2) * rotAxis.x, sin(rotAngNeg / 2) * rotAxis.y,
            sin(rotAngNeg / 2) * rotAxis.z };
        // make them unit quaternions
        rotQuatPos = rotQuatPos.unit();
        rotQuatNeg = rotQuatNeg.unit();

        // rotate the molecules
        rotate(reactIface1, rotQuatPos, reactCom1, moleculeList);
        rotate(reactIface1, rotQuatNeg, reactCom2, moleculeList);

        v1 = Vector { reactIface1 - reactMol1.tmpComCoord };
        sigma = Vector { reactIface1 - reactIface2 };
        v1.calc_magnitude();
        sigma.calc_magnitude();
        currTheta = sigma.dot_theta(v1);

        // if the angle between sigma and the rotated Molecule are correct, return.
        if ((areSameAngle(targAngle, M_PI) || areSameAngle(targAngle, 0)) && areSameAngle(currTheta, targAngle)) {
            // std::cout << "Theta After: " << currTheta << std::endl;
            return;
        }
        // std::cout << "Theta After: " << currTheta << std::endl;
    }
}
