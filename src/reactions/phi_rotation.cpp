#include "io/io.hpp"
#include "reactions/association/association.hpp"
#include "tracing.hpp"

void phi_rotation(Coord& reactIface1, Coord& reactIface2, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2,
    Complex& reactCom1, Complex& reactCom2, const Vector& normal, const double& targPhi, const ForwardRxn& currRxn,
    std::vector<Molecule>& moleculeList, const std::vector<MolTemplate>& molTemplateList)
{
    // TRACE();
    /*! \ingroup Associate
     * \brief Function to rotate a Molecule to some a target dihedral phi, relative to sigma (the
     (sigma)-(interface)-(normal) dihedral)
     *
     * @param[in] iface1 associating Molecule Interface of mol1
     * @param[in] iface2 associating Molecule Interface of mol2
     * @param[in] ifaceIndex2 the index of iface2 in mol2.assocICoords (see below for why)
     * @param[in] mol1 associating Molecule which we are rotating to the target theta angle
     * @param[in] mol2 the other associating Molecule
     * @param[in] targAngle target theta angle
     *
     * ### MAP
     *   1. Calculate current phi using calculate_phi(). See for details.
     *   2. If current phi is not the target angle, within tolerance, rotation is needed.
     *   3. Determine how far to rotate each protein according to their respective diffusion constants
     *   4. Determine rotation quaternions using positive rotation for ldrinant protein and negative for fferior
     *   5. Rotate.
     *   6. If rotated phi is incorrect, reverse rotation and rotate in the opposite direction (negative for ldrinant
     *      positive for fferior)
     *
     * Why the index of the fferior interface?
     *
     * The index of iface2 is needed because iface1 and iface2 are references to references to (yes, references to
     * references) their respective coordinates in their respective Molecule's assocICoords vector. Since iface1,
     * iface2, mol1, and mol2 are passed by value to calculate_phi(), and through it also to transform(). The
     * transform() function does a coordinate transformation which changes the mol1 and mol2 coordinates but does
     * not change iface1 or iface2, since they are references to different objects. Because of this, the creation of
     * the projected vectors in calculate_phi() will fail unless the index of iface2 is passed to it.
     */

    // axis of rotation is com-iface vector
    //    Vector rotAxis = Vector { reactMol1.tmpComCoord-reactIface1 };
    Vector rotAxis = Vector { reactIface1 - reactMol1.tmpComCoord };
    rotAxis.calc_magnitude();
    //    std::cout <<" in phi, rotAxis: "<<rotAxis.x<<' '<<rotAxis.y<<' '<<rotAxis.z<<std::endl;
    double currPhi { calculate_phi(
        reactIface1, ifaceIndex2, reactMol1, reactMol2, normal, rotAxis, currRxn, molTemplateList) };
    // std::cout << "Desired phi: " << targPhi << " Current phi: " << currPhi << std::endl;

    // quit out if the angles are damn near the same, or if the target is 0/M_PI, if the current angle is -0 or -M_PI
    if ((std::abs(targPhi - currPhi) < 1E-8)) {
        // std::cout << "No phi rotation needed" << std::endl;
    } else if ((areSameAngle(targPhi, M_PI) || areSameAngle(targPhi, 0)) && areSameAngle(targPhi, currPhi)) {
        // std::cout << "No phi rotation needed" << std::endl;
    } else {
        rotAxis.normalize(); // rotation axis must be normalized

        // Determine how much to rotate each complex by based on their respective contributions to the overall
        // rotational diffusion constant of the two molecules
        double posPhiRotAng {};
        double negPhiRotAng {};
        determine_rotation_angles(targPhi, currPhi, posPhiRotAng, negPhiRotAng, reactCom1, reactCom2);

        // create the rotation quaternions
        // std::cout << "rot_phi_pos: " << posPhiRotAng << " rot_phi_neg: " << negPhiRotAng << std::endl;
        Quat rotQuatPos(cos(posPhiRotAng / 2), sin(posPhiRotAng / 2) * rotAxis.x, sin(posPhiRotAng / 2) * rotAxis.y,
            sin(posPhiRotAng / 2) * rotAxis.z);
        Quat rotQuatNeg(cos(negPhiRotAng / 2), sin(negPhiRotAng / 2) * rotAxis.x, sin(negPhiRotAng / 2) * rotAxis.y,
            sin(negPhiRotAng / 2) * rotAxis.z);
        rotQuatPos = rotQuatPos.unit();
        rotQuatNeg = rotQuatNeg.unit();

        // rotate the two complexes
        rotate(reactIface1, rotQuatPos, reactCom1, moleculeList);
        rotate(reactIface1, rotQuatNeg, reactCom2, moleculeList);

        rotAxis = Vector(reactIface1 - reactMol1.tmpComCoord);
        rotAxis.normalize();
        currPhi
            = calculate_phi(reactIface1, ifaceIndex2, reactMol1, reactMol2, normal, rotAxis, currRxn, molTemplateList);

        //write_xyz_assoc("phi_forward2.xyz", reactCom1, reactCom2, moleculeList);
        if ((areSameAngle(targPhi, M_PI) || areSameAngle(targPhi, 0)) && areSameAngle(targPhi, currPhi)) {
            // std::cout << "Phi After: " << currPhi << std::endl;
            return;
        }

        // if it rotated the wrong way, reverse the rotation and swap quaternions
        if (std::abs(currPhi - targPhi) > 1E-11) {
            // reverse the rotation
            // std::cout << "Reversing rotation, current phi: " << currPhi << ", target phi: " << targPhi << '\n';
            reverse_rotation(
                reactIface1, reactMol1, reactMol2, reactCom1, reactCom2, rotQuatPos, rotQuatNeg, moleculeList);

            // make sure the reversal was successful
            currPhi = calculate_phi(
                reactIface1, ifaceIndex2, reactMol1, reactMol2, normal, rotAxis, currRxn, molTemplateList);
            // std::cout << "Phi after reversal: " << currPhi << std::endl;

            determine_rotation_angles(targPhi, currPhi, posPhiRotAng, negPhiRotAng, reactCom2, reactCom1);

            // get the new rotation quaternions
            // std::cout << "rot_phi_pos: " << posPhiRotAng << " rot_phi_neg: " << negPhiRotAng << std::endl;
            rotQuatPos = Quat(cos(posPhiRotAng / 2), sin(posPhiRotAng / 2) * rotAxis.x,
                sin(posPhiRotAng / 2) * rotAxis.y, sin(posPhiRotAng / 2) * rotAxis.z);
            rotQuatNeg = Quat(cos(negPhiRotAng / 2), sin(negPhiRotAng / 2) * rotAxis.x,
                sin(negPhiRotAng / 2) * rotAxis.y, sin(negPhiRotAng / 2) * rotAxis.z);
            rotQuatPos = rotQuatPos.unit();
            rotQuatNeg = rotQuatNeg.unit();

            // rotate the two complexes
            rotate(reactIface1, rotQuatNeg, reactCom1, moleculeList);
            rotate(reactIface1, rotQuatPos, reactCom2, moleculeList);

            rotAxis = Vector { reactIface1 - reactMol1.tmpComCoord };
            rotAxis.normalize();

            currPhi = calculate_phi(
                reactIface1, ifaceIndex2, reactMol1, reactMol2, normal, rotAxis, currRxn, molTemplateList);

            if ((areSameAngle(targPhi, M_PI) || areSameAngle(targPhi, 0)) && areSameAngle(targPhi, currPhi)) {
                // std::cout << "Phi After: " << currPhi << std::endl;
                return;
            }
            //write_xyz_assoc("phi_reversal2.xyz", reactCom1, reactCom2, moleculeList);

            // If the rotation still didn't give the desired angle
            if (std::abs(currPhi - targPhi) > 1E-4) {
                std::cerr << "Cannot resolve phi angle for protein " << reactMol1.index << " in complex "
                          << reactCom1.index << " relative to protein " << reactMol2.index << " in complex "
                          << reactCom2.index << ". Final angle: " << currPhi << ". Exiting..." << std::endl;
                exit(1);
            }
        }
    }
}
