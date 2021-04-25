/*! \file association_functions.hpp
 * ### Created on 5/20/18 by Matthew Varga
 * ### Purpose
 * ***
 *
 * ### Notes
 * ***
 */

#pragma once

#include "classes/class_Quat.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"

extern unsigned numAssoc;

/* MAIN FUNCTION */
/*! \ingroup Associate
 * \brief Main association function, which puts the two complexes at sigma and then performs the rotations.
 *
 * \param[in] ifaceIndex1 Index of first reactant interface in reactMol1.interfaceList.
 * \param[in] ifaceIndex2 Index of second reactant interface in reactMol2.interfaceList.
 * \param[in] reactMol1 Molecule object for the parent Molecule of the first reacting interface
 * \param[in] reactMol2 Molecule object for the parent Molecule of the second reacting interface
 * \param[in] reactCom1 Parent Complex object for reactMol1
 * \param[in] reactCom2 Parent Complex object for reactMol2
 *
 * Delegates to theta_rotation(), phi_rotation(), and omega_rotation() to do the actual rotations themselves
 *
 * ### MAP
 *   1. Set up temporary association coordinates
 *   2. Translate the reacting interfaces such that they overlap
 *   3. Translate the reacting interfaces to sigma along the center of mass to center of mass vector
 *   4. Rotate theta, starting with theta_2
 *   5. Rotate phi, starting with phi_2
 *   6. Rotate omega
 *   7. Check angles, vector magnitudes to make sure everything worked properly
 *
 */
void associate(long long int iter, int ifaceIndex1, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2, Complex& reactCom1,
    Complex& reactCom2, const Parameters& params, ForwardRxn& currRxn, std::vector<Molecule>& moleculeList,
    std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays, std::vector<Complex>& complexList, Membrane& membraneObject, const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns);
void associate_box(long long int iter, int ifaceIndex1, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2, Complex& reactCom1,
    Complex& reactCom2, const Parameters& params, ForwardRxn& currRxn, std::vector<Molecule>& moleculeList,
    std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays, std::vector<Complex>& complexList, Membrane& membraneObject, const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns);
void associate_sphere(long long int iter, int ifaceIndex1, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2, Complex& reactCom1,
    Complex& reactCom2, const Parameters& params, ForwardRxn& currRxn, std::vector<Molecule>& moleculeList,
    std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays, std::vector<Complex>& complexList, Membrane& membraneObject, const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns);

void associate_implicitlipid(int ifaceIndex1, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2, Complex& reactCom1,
    Complex& reactCom2, const Parameters& params, ForwardRxn& currRxn, std::vector<Molecule>& moleculeList,
    std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays, std::vector<Complex>& complexList, Membrane& membraneObject, const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns);
void associate_implicitlipid_box(int ifaceIndex1, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2, Complex& reactCom1,
    Complex& reactCom2, const Parameters& params, ForwardRxn& currRxn, std::vector<Molecule>& moleculeList,
    std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays, std::vector<Complex>& complexList, Membrane& membraneObject, const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns);
void associate_implicitlipid_sphere(int ifaceIndex1, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2, Complex& reactCom1,
    Complex& reactCom2, const Parameters& params, ForwardRxn& currRxn, std::vector<Molecule>& moleculeList,
    std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays, std::vector<Complex>& complexList, Membrane& membraneObject, const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns);

/* BOOLEANS */
/*! \ingroup Associate
 * \brief This checks if two angles are equal, accounting for minus signs
 */
bool areSameAngle(double ang1, double ang2);

/*! \ingroup Associate
 * \brief Checks if two vectors are parallel
 */
bool areParallel(const double& angle);

/*! \ingroup Associate
 * \brief Checks if the sign of the angle is the correct sign demanded by the desired angle.
 */
bool angleSignIsCorrect(const Vector& vec1, const Vector& vec2);

/*!
 * \brief Checks the center of mass to interface vectors of the parent complex of an associating molecule to make sure
 * all magnitudes were conserved during the association.
 *
 * Just a safety check. Can be removed, if need be.
 */
bool conservedMags(const Complex& targCom, const std::vector<Molecule>& moleculeList);

/*!
 * \brief Checks the center of mass to interface vectors of the parent Complex of an associating Molecule to ensure the
 * angles were conserved during association.
 *
 * Just a safety check. Can be removed, if need be.
 */
bool conservedRigid(const Complex& targCom, const std::vector<Molecule>& moleculeList);

/*! \ingroup Associate
 * \brief A generic version of dihedral_projection and phi_projection, to project any arbitray axes onto a plane
 * TODO: Rename this
 */
bool requiresSignFlip(Vector axis, Vector v1, Vector v2);

/* MATH FUNCTIONS */
/*! \ingroup Associate
 * \brief Create an arbitrary vector orthogonal another vector, base on the cross product of it and either the x or
 * y axis.
 */
Vector create_arbitrary_vector(Vector& vec);

/*! \ingroup Associate
 * \brief Transforms a molecule's interface and center of coordinates, used only by projection functions (project(),
 * project_omega(), project_phi()).
 *
 * \param[in] reactIface reacting interface of the current molecule being rotated and transformed
 * \param[in] reactMol1 current molecule being rotated and transformed
 * \param[in] reactMol2 other molecule in the association event
 * \param[in] axis axis around which the molecules are rotating
 *
 * NOTE: Only use this function if the calling function is passed reactMol1 and reactMol2 by value, NOT by
 * reference. This function alters the coordinates in a way that is not desired by association, but required in the
 * course of it.
 */
void transform(Coord& reactIface, Molecule& reactMol1, Molecule& reactMol2, const Vector& axis);

/*! \ingroup Parser
 * \brief This subroutine determines the normal of the protein based on
 *  some previously determined vector w.r.t. the internal coordinates
 *  of the protein in question.
 *
 *  This is probably super inefficient.
 * Procedure:
 *   1. center molecule with center of mass at (0,0,0) (for this to work, all proteins must have their center of
 *   mass at (0,0,0))
 *   2. First rotation: determine quat and rotate the protein so that interface 1 is in the correct position
 *   3. Second rotation:
 *       1. only if the protein has more than one interface (obviously)
 *       2. use the first interface (already rotated) as the axis of rotation
 *       3. project the current second (interface)-(center of mass) vector and the desired
 *       (interface)-(center of mass) vector onto a plane normal to the axis of rotation
 *       (otherwise the angle of rotation will be wrong)
 *       4. determine quaternion and rotate
 *   4. combine quaternions, if interfaces > 1, and rotate the internal coordinate normal in the reverse
 *   direction to determine the normal for the real coordinates of the protein
 *   5. normalize and pass back
 */
Vector determine_normal(Vector normal, const MolTemplate& molTemplate, Molecule oneMol);

/*! \ingroup Associate
 * \brief Determine the rotation angles for each Complex, according to their diffusion constants.
 *
 * If both complexes have no rotational diffusion (i.e. both complexes are on the membrane), it uses translational
 * diffusion constants, in either the x- or z-component, depending on if either are larger than zero.
 */
void determine_rotation_angles(double targAngle, double currAngle, double& rotAngPos, double& rotAngNeg,
    const Complex& reactCom1, const Complex& reactCom2);

/* ROTATION FUNCTIONS */
/*! \ingroup Associate
 * \brief This function is meant to determine a quaternion for the rotation of the target Molecule to its
 * MolTemplate coordinates, based on the interface to center of mass vectors of the MolTemplate.
 *
 * If the Molecule is not a rod, two interfaces to center of mass vectors are used to determine the quaternion.
 * TODO: Write how/why
 */
Quat orient_crds_to_template(const MolTemplate& oneTemplate, Molecule& targMol);

/*! \ingroup Associate
 * \brief This function is meant to rotate a target Complex with a rotation quaternion determined from the
 * associating Molecule, and then translate it to a specific position.
 *
 * \param[in] rotOrigin The origin for all vectors to rotate around.
 * \param[in] rotQuat Rotation Quaternion determined in the angle-relevant rotation function
 * \param[in] targCom The complex being rotated.
 * \param[in] moleculeList List of all Molecules in the system.
 *
 * TODO: Write how
 *
 */
void rotate(Coord& rotOrigin, Quat& rotQuat, Complex& targCom, std::vector<Molecule>& moleculeList);

/*! \ingroup Associate
 * \brief This function is meant to rotate a Molecule in the reverse direction.
 *
 * - It does so by determining the inverse of the two quaterions (from the half angles), and rotating with these new
 * quaternions.
 * - \f$ Q^{-1} (Q(vec)Q^{-1})Q \f$
 *
 */
void reverse_rotation(Coord& reactIface1, Molecule& reactMol1, Molecule& reactMol2, Complex& reactCom1,
    Complex& reactCom2, Quat rotQuatPos, Quat rotQuatNeg, std::vector<Molecule>& moleculeList);

/*! \ingroup Associate
ForwardRxn& rxn, Complex& domCom, Complex& infCom, std::vector<Molecule>& molList)
 * \brief Function to rotate a Molecule to some a target angle theta, relative to sigma (the angle between the two
 * associating interfaces).
 *
 * \param[in] reactIface1 The associating Molecule Interface of domMol
 * \param[in] reactIface2 The associating Molecule Interface of infMol
 * \param[in] reactMol1 The associating Molecule which we are rotating to the target theta angle
 * \param[in] reactMol2 The other associating Molecule
 * \param[in] targAngle Target theta angle between the two Molecules.
 * \param[in] currRxn Foward reaction (ForwardRxn) currently being performed.
 * \param[in] reactCom1 The parent Complex of the associating dominant Molecule (the Molecule currently being rotated)
 * \param[in] reactCom2 The parent Complex of the associating inferior Molecule (the other Molecule in the reaction)
 * \param[in] moleculeList List of all molecules in the system.
 *
 * ### MAP
 *   1. Calculate current theta by taking the dot product of the interface to center of mass vector and sigma:
 *      \f$ \theta = acos((com-iface) \cdot \sigma) \f$
 *   2. If the current theta is not equal to the target theta, within a tolerance, rotation is needed.
 *   3. Determine the axis of rotation
 *     - If the target angle is either \f$\pi\f$ or 0, we can't use the cross product of the (com)-(iface) and
 *     \f$\sigma\f$ as the axis of rotation, so we use the one of the cartesian coordinate axes.
 *   4. Get quaternions for the rotation
 *     - Positive rotation for the dominant protein and negative rotation for the inferior protein.
 *     - Based on each complex's diffusion constants
 *   5. Rotate
 *   6. Check rotation, if incorrect, reverse rotation and rotate the other way (i.e. negative rotation for dominant
 *      and positive for inferior).
 *     - This is a safeguard, and it's not often needed. Probably due to missing a sign check somewhere.
 *
 * ### Uses quaternion rotation:
 *   - \f$Q (vec) Q^{-1}\f$ --> rotation of vector using quaternion \f$Q\f$. \f$Q^{-1}\f$ is the inverse of \f$Q\f$,
 *  \f$Q = (C, w)\f$ where C is the rotation angle and w is the axis of rotation \f$(ai + bj +
 *  ck)\f$.
 *   - Q must be normalized for rotation. Normalizes entire thing, including C
 *   - For rotation we use half angles (ang/2), \f$Q = (sin(ang/2), cos(ang/2) * w.x, cos(ang/2) * w.y,
 * cos(ang/2) * w.z)\f$
 */
void theta_rotation(Coord& reactIface1, Coord& reactIface2, Molecule& reactMol1, Molecule& reactMol2, double targAngle,
    Complex& reactCom1, Complex& reactCom2, std::vector<Molecule>& moleculeList);

/*! \ingroup Associate
 * \brief Function to rotate a Molecule to some a target dihedral phi, relative to sigma (the
 (sigma)-(interface)-(normal) dihedral)
 *
 * \param[in] reactIface1 The associating Molecule Interface of domMol
 * \param[in] reactIface2 The associating Molecule Interface of infMol
 * \param[in] reactMol1 The associating Molecule which we are rotating to the target theta angle
 * \param[in] reactMol2 The other associating Molecule
 * \param[in] targAngle Target theta angle between the two Molecules.
 * \param[in] currRxn Foward reaction (ForwardRxn) currently being performed.
 * \param[in] reactCom1 The parent Complex of the associating dominant Molecule (the Molecule currently being rotated)
 * \param[in] reactCom2 The parent Complex of the associating inferior Molecule (the other Molecule in the reaction)
 * \param[in] moleculeList List of all molecules in the system.
 *
 * ### MAP
 *   1. Calculate current phi using calculate_phi(). See for details.
 *   2. If current phi is not the target angle, within tolerance, rotation is needed.
 *   3. Determine how far to rotate each protein according to their respective diffusion constants
 *   4. Determine rotation quaternions using positive rotation for dominant protein and negative for inferior
 *   5. Rotate.
 *   6. If rotated phi is incorrect, reverse rotation and rotate in the opposite direction (negative for dominant
 *      positive for inferior)
 *
 * Why the index of the inferior interface?
 *
 * The index of infIface is needed because domIface and infIface are references to references to (yes, references to
 * references) their respective coordinates in their respective Molecule's assocICoords vector. Since domIface,
 * infIface, domMol, and infMol are passed by value to calculate_phi(), and through it also to transform(). The
 * transform() function does a coordinate transformation which changes the domMol and infMol coordinates but does
 * not change domIface or infIface, since they are references to different objects. Because of this, the creation of
 * the projected vectors in calculate_phi() will fail unless the index of infIface is passed to it.
 */
void phi_rotation(Coord& reactIface1, Coord& reactIface2, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2,
    Complex& reactCom1, Complex& reactCom2, const Vector& normal, const double& targPhi, const ForwardRxn& currRxn,
    std::vector<Molecule>& moleculeList, const std::vector<MolTemplate>& molTemplateList);

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
void omega_rotation(Coord& reactIface1, Coord& reactIface2, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2,
    Complex& reactCom1, Complex& reactCom2, double targOmega, const ForwardRxn& currRxn,
    std::vector<Molecule>& moleculeList, const std::vector<MolTemplate>& molTemplateList);

/* PROJECTIONS */
/*! \ingroup Associate
 * \brief this should get omega via an orthographic projected onto the xy-plane.
 *  why? because if the com-iface/norm vectors are not parallel,
 *  i.e. if theta1, theta2 != M_PI the code gets confused
 *
 *   1. rotate the two proteins such that sigma is on the z axis
 *   2. do orthographic projection onto the xy-plane.
 *   3. determine omega between the com-iface, or the norm vectors, if one of the theta angles is M_PI
 *   4. return omega.
 */
double calculate_omega(Coord reactIface1, int reactIface2, Vector& sigma,
    const ForwardRxn& currRxn, Molecule reactMol1, Molecule reactMol2, const std::vector<MolTemplate>& molTemplateList);

/*! \ingroup Associate
 * \brief this should get phi via an orthographic projected onto the xy-plane.
 *
 *   1. first rotation the two proteins such the com-iface vector is on the z axis
 *   2. do orthographic projection onto the xy-plane.
 *   3. determine phi between sigma and norm.
 *   4. return phi.
 */
double calculate_phi(Coord reactIface1, int ifaceIndex2, Molecule reactMol1, Molecule reactMol2, const Vector& normal,
    Vector axis, const ForwardRxn& currRxn, const std::vector<MolTemplate>& molTemplateList);

/* OTHER */
/*! \ingroup Associate
 * \brief Checks to make sure association was successful.
 *
 * Checks if
 *   1. all angles are at their desired value
 *   2. magnitude of the interface-center of mass vectors was conserved
 *   3. the rigitidy of the molecules was conserved
 */
void check_bases(bool& cancelAssoc, const Coord& reactIface1, const Coord& reactIface2, int ifaceIndex1,
    int ifaceIndex2, const Molecule& reactMol1, const Molecule& reactMol2, const Complex& reactCom1,
    const Complex& reactCom2, const ForwardRxn& currRxn, const std::vector<Molecule>& moleculeList,
    const std::vector<MolTemplate>& molTemplateList);

/*! \ingroup Associate
 * \brief Checks to see if the centers of masses of any of the molecules that are undergoing physical association 
 * to prevent steric overlap. Only checks molecules that are flagged with checkOverlap=1
 *
 * If so, cancels association.
 */
void check_for_structure_overlap(bool& cancelAssoc, const Complex& reactCom1, const Complex& reactCom2,
    const std::vector<Molecule>& moleculeList, const Parameters& params,
    const std::vector<MolTemplate>& molTemplateList);

/*! \ingroup Associate
 * \brief Checks to see if the centers of masses of any of the molecules that are undergoing physical association
 * overlap with any of the other molecules in the system! Only checks molecules that are flagged with checkOverlap=1
 *
 * If so, cancels association.
 */

void check_for_structure_overlap_system(bool& flag, const Complex& reactCom1, const Complex& reactCom2,
    std::vector<Molecule>& moleculeList, const Parameters& params,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns);

/*! \ingroup Associate
 * \brief Checks to see if the centers of masses of any of the molecules that are undergoing physical association
 * have displaced by a large distance, returns currently a warning if the displacement is larger than LARGE_DISP.
 *
 * Can be set to cancel association if the (flag=true) if the displacement of either complex is too large,
 *as it is unphysical within a single time-step. Can occur for large complexes re-orienting to bind to one another, causing
 *some proteins to move long distances as the complex rotates.
 */
void measure_complex_displacement(bool& flag, Complex& reactCom1, Complex& reactCom2,
    std::vector<Molecule>& moleculeList, const Parameters& params,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<Complex>& complexList);

/*! \ingroup Associate
 * \brief  Routine to calculate how far a vector has rotated (in radians) from its position store in original coordinates to its position store in tmpCoords. Used during associate to evaluate the extent to which rotation into the proper orientation has caused re-alignment of the interface-COM vectors of associating interfaces.
 *  */
void calc_angular_displacement(int ifaceIndex1, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2, Complex& reactCom1, Complex& reactCom2, std::vector<Molecule>& moleculeList);

/*! \ingroup Associate
 * \brief  Routine to calculate how far a vector has rotated (in radians) from its position store in original coordinates to its position store in tmpCoords. Returns angle of rotation in radians.
 *  */

double calc_one_angular_displacement(int ifaceIndex1, Molecule& reactMol1, Complex& reactCom1);


/*! \ingroup Associate
 * \brief measures separation between interfaces of two proteins, one that is associating (use tmpCoords) one that is in the system 
 *(use full coords).  does not test whether they are capable of binding with each other--all interface pairs evaluated.
 *is called within check_for_structure_overlap_system for two proteins that checkOverlap. 
 *
 * If interfaces overlap, cancels association.
 */

void measure_overlap_protein_interfaces(Molecule base1, Molecule baseTmp, bool& flagCancel);

void measure_overlap_free_protein_interfaces(Molecule base1, Molecule baseTmp, bool& flagCancel,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns);

/*! \ingroup Associate
 * \brief Store Calculate rotation matrix for orienting one molecule to itself (at another timepoint, e.g.)
 * First molecule is tru coordinate, second molecule--uses temporarary coordinates
 */

Quat save_mem_orientation(Molecule baseTarget, Molecule base1, MolTemplate onePro);

/*! \ingroup Associate
 * \brief calculate the COM of two complexes prior to their association
 *
 */

void com_of_two_tmp_complexes(Complex& reactCom1, Complex& reactCom2, Coord& vectorCOM, std::vector<Molecule>& moleculeList);

/*! \ingroup Associate
 * \brief update the COM of a complex, based on its tmpComCoord, and the members tmpComCoords as well (during association only!)
 *
 */

void update_complex_tmp_com_crds(Complex& reactCom, std::vector<Molecule>& moleculeList);

double get_geodesic_distance(Coord intFace1, Coord intFace2);

/*! \ingroup Associate
 * \brief keep track of any association event, based on the size of the smaller complex being added to the larger complex. If it is 1+1, specifically denotes dimerizaiton.
 *
 */

void track_association_events(Complex& reactCom1, Complex& reactCom2, bool transitionToSurface, bool isOnMembrane, copyCounters& counterArrays);

void print_association_events( copyCounters& counterArrays, std::ofstream & outfile, int it, Parameters params);
void init_association_events( copyCounters& counterArrays);

