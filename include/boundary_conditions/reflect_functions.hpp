/*! \file old_funcs.hpp
 * \brief Functions for enforcing reflecting boundary conditions
 * ### Created on 11/1/18 by Matthew Varga
 * ### TODO List
 * ***
 */
#pragma once

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Membrane.hpp"

/*! \defgroup BoundaryConditions
 * \brief Functions associated with enforcing various boundary conditions
 */

/* DISSOCIATION REFLECTION */

/*! \ingroup BoundaryConditions
 * \brief Enforces reflecting boundary conditions during dissociation reactions.
 */
void reflect_complex_rad_rot(const Membrane& membraneObject, Complex& targCom, std::vector<Molecule>& moleculeList);
void reflect_complex_rad_rot_inx(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList);

/* TRAJ PROPAGATION FUNCTIONS */

/*! \ingroup BoundaryConditions
 * \brief Enforces reflecting boundary conditions during complex propagation.
 */
void reflect_traj_complex_rad_rot(
				  const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, std::array<double, 9>& M, const Membrane &membraneObject);
void reflect_traj_complex_rad_rot_new(
				  const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, std::array<double, 9>& M, const Membrane &membraneObject);

/*! \ingroup BoundaryConditions
 * \brief Checks to make sure the complex doesn't span the box (go out of the box in both sides).
 *
 * Child function of reflect_traj_complex_rad_rot
 */
void reflect_traj_check_span(double xtot, double ytot, double ztot, const Parameters& params, Complex& targCom,
			     std::vector<Molecule>& moleculeList, std::array<double, 9>& M, const Membrane &membraneObject);

/*!
 * \brief Enforces reflecting boundary conditions by placing Complex and component Molecules back into the box without
 * checking bounds after doing so.
 *
 * Child function of reflect_traj_check_span
 */
void reflect_traj_complex_rad_rot_nocheck(
					  const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, std::array<double, 9>& M, const Membrane  &membraneObject);

/*!
 * \brief Checks if the newly associated complex spans the simulation volume.
 *
 * If it does, cancels association. If complex is out of the box in either the positive or negative direction (of any x,
 * y, z dimensions), it puts the complex back into the box at the edge, not bounced off, as in
 * reflect_traj_complex_rad_rot, and its children.
 */

void check_if_spans_box(bool &cancelAssoc, const Parameters& params, Complex& reactCom1, Complex& reactCom2,
                        std::vector<Molecule>& moleculeList, const Membrane &membraneObject);

/*! \ingroup BoundaryConditions
 * \brief evaluates size of reflection off of walls during association, stores in temporary vector traj. 
 *based on tmpCoords, does not update complex.traj vectors..
 */
void reflect_traj_tmp_crds(const Parameters& params, std::vector<Molecule>& moleculeList, Complex& targCom, std::array<double, 3> &traj, const Membrane &membraneObject);

// if the complex dissociate from the implicit-lipid on membrane surface, then it must be put at RS3D, according to the rules of implicit-lipid model 
void reflect_complex_RS3D(const Membrane& membraneObject, Complex& targCom, std::vector<Molecule>& moleculeList);
