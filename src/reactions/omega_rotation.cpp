#include "reactions/association/association.hpp"
#include "tracing.hpp"

void omega_rotation(Coord& reactIface1, Coord& reactIface2, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2,
    Complex& reactCom1, Complex& reactCom2, double targOmega, const ForwardRxn& currRxn,
    std::vector<Molecule>& moleculeList, const std::vector<MolTemplate>& molTemplateList)
{
    // TRACE();
    /*! \ingroup Associate
     * \brief Performs a rotation of two molecules to some target omega angle, the com-iface1 to sigma to com-iface2
     * dihedral angle
     *
     * \param[in] reactIface1 first reactant interface (order is arbitrary)
     * \param[in] reactIface2 second reactant interface
     * \param[in] reactMol1 first reactant molecule
     * \param[in] reactMol2 second reactant molecule
     * \param[in] reactCom1 parent complex of first reactant molecule
     * \param[in] reactCom2 parent complex of second reactant molecule
     * \param[in] currRxn current association reaction
     * \param[in] moleculeList list of all Molecules in the system
     * \param[in] molTemplateList list of all MolTemplates
     */
    Vector rotAxis { reactIface1 - reactIface2 };
    rotAxis.normalize();

    double currOmega { calculate_omega(reactIface1, ifaceIndex2, rotAxis, currRxn, reactMol1, reactMol2, molTemplateList) };

    // double omega {dotproduct(norm1, norm2)};
    // std::cout << "Desired omega: " << targOmega << " current omega: " << currOmega << std::endl;

    if (std::abs(targOmega - currOmega) < 1E-8) {
        // std::cout << "No omega rotation needed" << std::endl;
        return;
    } else if ((areSameAngle(targOmega, M_PI) || areSameAngle(targOmega, 0)) && areSameAngle(targOmega, currOmega)) {
        //std::cout << "No omega rotation needed" << std::endl;
    } else {

        // both complexes rotate around the axis of rotation
        // proportional to their respective rotational diffusion constants
        double posOmegaRotAngle { 0 };
        double negOmegaRotAngle { 0 };
        determine_rotation_angles(targOmega, currOmega, posOmegaRotAngle, negOmegaRotAngle, reactCom1, reactCom2);
        // std::cout << "Positive half angle: " << posOmegaRotAngle << "\nNegative half angle: " << negOmegaRotAngle
        //           << std::endl;

        // Create the quaternions
        Quat omegaPos(cos(posOmegaRotAngle / 2), sin(posOmegaRotAngle / 2) * rotAxis.x,
            sin(posOmegaRotAngle / 2) * rotAxis.y, sin(posOmegaRotAngle / 2) * rotAxis.z);
        Quat omegaNeg(cos(negOmegaRotAngle / 2), sin(negOmegaRotAngle / 2) * rotAxis.x,
            sin(negOmegaRotAngle / 2) * rotAxis.y, sin(negOmegaRotAngle / 2) * rotAxis.z);
        // make them unit quaternions
        omegaPos = omegaPos.unit();
        omegaNeg = omegaNeg.unit();

        // rotate the two molecules
        rotate(reactIface2, omegaNeg, reactCom2, moleculeList);
        rotate(reactIface2, omegaPos, reactCom1, moleculeList);

        rotAxis = Vector { reactIface1 - reactIface2 };
        rotAxis.normalize();

        currOmega = calculate_omega(reactIface1, ifaceIndex2, rotAxis, currRxn, reactMol1, reactMol2, molTemplateList);

        if ((areSameAngle(targOmega, M_PI) || areSameAngle(targOmega, 0)) && areSameAngle(targOmega, currOmega)) {
            // std::cout << "Omega After: " << currOmega << std::endl;
            return;
        }
        if (areSameAngle(targOmega, M_PI) && areSameAngle(-M_PI, currOmega)) {
            // std::cout << "Omega After (-pi==pi): " << currOmega << std::endl;
            return;
        }

        double tol = 1e-8;
        if (std::abs(currOmega - targOmega) > tol) {
            // std::cout << "Rotated to omega in the wrong direction. " << currOmega << " Reversing and rotating the other way." << '\n';
            reverse_rotation(reactIface1, reactMol1, reactMol2, reactCom1, reactCom2, omegaPos, omegaNeg, moleculeList);

            currOmega = calculate_omega(reactIface1, ifaceIndex2, rotAxis, currRxn, reactMol1, reactMol2, molTemplateList);
            // std::cout << "Omega after reset: " << currOmega << std::endl;

            determine_rotation_angles(targOmega, currOmega, posOmegaRotAngle, negOmegaRotAngle, reactCom2, reactCom1);

            // std::cout << "Reversed, new rotation with\n";
            // std::cout << "Positive half angle: " << posOmegaRotAngle << "\nNegative half angle: " << negOmegaRotAngle
            //           << std::endl;
            omegaPos = Quat(cos(posOmegaRotAngle / 2), sin(posOmegaRotAngle / 2) * rotAxis.x,
                sin(posOmegaRotAngle / 2) * rotAxis.y, sin(posOmegaRotAngle / 2) * rotAxis.z);
            omegaNeg = Quat(cos(negOmegaRotAngle / 2), sin(negOmegaRotAngle / 2) * rotAxis.x,
                sin(negOmegaRotAngle / 2) * rotAxis.y, sin(negOmegaRotAngle / 2) * rotAxis.z);
            omegaPos = omegaPos.unit();
            omegaNeg = omegaNeg.unit();

            rotate(reactIface2, omegaPos, reactCom2, moleculeList);
            rotate(reactIface2, omegaNeg, reactCom1, moleculeList);

            // check angle again, if not correct, cancel association

            rotAxis = Vector { reactIface1 - reactIface2 };
            rotAxis.normalize();

            currOmega = calculate_omega(reactIface1, ifaceIndex2, rotAxis, currRxn, reactMol1, reactMol2, molTemplateList);
            if ((areSameAngle(targOmega, M_PI) || areSameAngle(targOmega, 0)) && areSameAngle(targOmega, currOmega)) {
                // std::cout << "Omega After: " << currOmega << std::endl;
                return;
            } else {
                // std::cout << "WARNING: OMEGA NOT CONVERGED TO CORRECT VALUE, CURR VALUE: " << currOmega << std::endl;
            }
        }
    }
}
