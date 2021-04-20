#include "reactions/association/association.hpp"
#include "tracing.hpp"

#include <iomanip>

void check_bases(bool& cancelAssoc, const Coord& reactIface1, const Coord& reactIface2, int ifaceIndex1,
    int ifaceIndex2, const Molecule& reactMol1, const Molecule& reactMol2, const Complex& reactCom1,
    const Complex& reactCom2, const ForwardRxn& currRxn, const std::vector<Molecule>& moleculeList,
    const std::vector<MolTemplate>& molTemplateList)
{
    // TRACE();
    Vector sigma { reactIface1 - reactIface2 };
    Vector v1 { reactIface1 - reactMol1.tmpComCoord };
    Vector v2 { reactIface2 - reactMol2.tmpComCoord };
    sigma.calc_magnitude();
    v1.calc_magnitude();
    v2.calc_magnitude();

    if (std::abs(currRxn.assocAngles.theta1 - sigma.dot_theta(v1)) < 1E-6) {
        //        && roundv(currRxn.assocAngles.theta1) != roundv((2 * M_PI) - sigma.dot_theta(v1))) {
        // std::cerr << "Error: Angle theta1 is " << sigma.dot_theta(v1) << " but should be "
        //           << roundv(currRxn.assocAngles.theta1) << std::endl;
        cancelAssoc = true;
        return;
    }

    Vector tmpSigma { -sigma }; // if you change sign of original sigma, screws up the calculation of omega below
    tmpSigma.calc_magnitude();
    if (roundv(currRxn.assocAngles.theta2) != roundv(tmpSigma.dot_theta(v2))
        && roundv(currRxn.assocAngles.theta2) != roundv((2 * M_PI) - tmpSigma.dot_theta(v2))) {
        // std::cerr << "Error: Angle theta_2 is " << tmpSigma.dot_theta(v2) << " but should be "
        //           << currRxn.assocAngles.theta2 << std::endl;
        cancelAssoc = true;
        return;
    }

    if (!std::isnan(currRxn.assocAngles.phi1)) {
        double phi { roundv(calculate_phi(reactIface1, ifaceIndex2, reactMol1, reactMol2, currRxn.norm1,
            Vector { reactIface1 - reactMol1.tmpComCoord }, currRxn, molTemplateList)) };

        if (areSameAngle(currRxn.assocAngles.phi1, M_PI) && ((std::abs(phi) > 3.1414) && std::abs(phi) < 3.1417))
            phi = std::abs(phi);

        if (roundv(currRxn.assocAngles.phi1) != roundv(phi)) {
            // std::cerr << "Error: Reactant 1 dihedral is " << phi << " but should be " << currRxn.assocAngles.phi1
            //           << std::endl;
            cancelAssoc = true;
            return;
        }
    }

    if (!std::isnan(currRxn.assocAngles.phi2)) {
        double phi { roundv(calculate_phi(reactIface2, ifaceIndex1, reactMol2, reactMol1, currRxn.norm2,
            Vector { reactIface2 - reactMol2.tmpComCoord }, currRxn, molTemplateList)) };

        if (areSameAngle(currRxn.assocAngles.phi2, M_PI) && ((std::abs(phi) > 3.1414) && std::abs(phi) < 3.1417))
            phi = std::abs(phi);

        if (roundv(currRxn.assocAngles.phi2) != roundv(phi)) {
            // std::cerr << "Error: Reactant 2 dihedral is " << phi << " but should be " << currRxn.assocAngles.phi2
            //           << std::endl;
            cancelAssoc = true;
            return;
        }
    }

    // if either molecule is a rod and has theta 0 or M_PI to sigma, there is no omega angle
    //    bool hasOmegaMol1{ molTemplateList[reactMol1.molTypeIndex].isRod
    //        && (areEqual(currRxn.assocAngles.theta1, 0) || areSameAngle(currRxn.assocAngles.theta1, M_PI)) };
    //    bool hasOmegaMol2{ molTemplateList[reactMol2.molTypeIndex].isRod
    //        && (areEqual(currRxn.assocAngles.theta2, 0) || areSameAngle(currRxn.assocAngles.theta2, M_PI)) };
    //    if (!(hasOmegaMol1 || hasOmegaMol2)) {
    if (!std::isnan(currRxn.assocAngles.omega)) {
        double omega { roundv(
            calculate_omega(reactIface1, ifaceIndex2, sigma, currRxn, reactMol1, reactMol2, molTemplateList)) };

        // check this weirdness
        if (areSameAngle(currRxn.assocAngles.omega, M_PI) && ((std::abs(omega) > 3.1414) && std::abs(omega) < 3.1417))
            omega = std::abs(omega);

        if (roundv(currRxn.assocAngles.omega) != roundv(omega)) {
            // std::cerr << "Error: Complex dihedral is " << omega << " but should be " << currRxn.assocAngles.omega
            //           << std::endl;
            cancelAssoc = true;
            return;
        }
    }

    // check to make sure vector magnitudes were conserved (TODO: remove for final version)
    if (!conservedMags(reactCom1, moleculeList)) {
        cancelAssoc = true;
        return;
    }
    if (!conservedMags(reactCom2, moleculeList)) {
        cancelAssoc = true;
        return;
    }

    // check to make sure the rigidity of the molecule didn't get fucked up
    if (!conservedRigid(reactCom1, moleculeList)) {
        cancelAssoc = true;
        return;
    }
    if (!conservedRigid(reactCom2, moleculeList)) {
        cancelAssoc = true;
        return;
    }

    // std::cout << "Angles and vector magnitudes are correct, association successful." << std::endl;
    // std::cout << std::setw(30) << std::setfill('-') << ' ' << std::setfill(' ') << std::endl;
}
