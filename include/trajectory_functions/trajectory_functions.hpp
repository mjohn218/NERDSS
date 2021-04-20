/*! \file trajectory_functions.hpp

 * ### Created on 10/25/18 by Matthew Varga
 * ### Purpose
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */
#pragma once

//#include "classes/class_mol_containers.hpp"
#include "classes/class_Rxns.hpp"

/*!
 * \brief Simply translates a Complex and its member Molecules along some translation Vector.
 * \param[in] transVec Vector along which the Complex will be translated.
 * \param[in] targCom Complex to be translated.
 * \param[in] moleculeList List of all Molecules in the system.
 */
void translate_complex_tmpCoords(const Vector& transVec, Complex& targCom, std::vector<Molecule>& moleculeList);

/*!
 * \brief Propagate the molecule with random translational and rotational motion
 * \param[in] params Simulation parameters as provided by user.
 * \param[in] targCom The complex to be propagated.
 * \param[in] moleculeList List of all Molecules in the system.
 *
 * Only occurs if the molecule has not dissociated or associated during the timestep
 */
void create_complex_propagation_vectors(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject);
Coord create_complex_propagation_vectors_on_sphere(const Parameters& params, Complex& targCom);

bool complexSpansBox(Vector& transVec, const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList);

/*!
 * \brief Check if propagation will cause the complex to span the entire box.
 * \param[in] params Simulation parameters as provided by user.
 * \param[in] targCom The complex to be propagated.
 * \param[in] moleculeList List of all Molecules in the system.
 *
 * If so, resample the translational and rotational propagation. Replaces reflect_traj_complex_rat_rot.cpp
 */
void reflect_complex(Vector& transVec, const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList,
    const std::vector<MolTemplate>& molTemplateList);

/*!
 * \brief Simply sets the list of previous reweighting vectors with the current reweighting vectors and clears the
 * current reweighting vectors
 */
void clear_reweight_vecs(Molecule& oneMol);

/*!
 * \brief Checks for overlap of proteins on the membrane.
 *
 * Membrane version of sweep_separation_complex_rot
 */
void sweep_separation_complex_rot_memtest(int simItr, int pro1Index, Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject);
void sweep_separation_complex_rot_memtest_box(int simItr, int pro1Index, Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject);
void sweep_separation_complex_rot_memtest_sphere(int simItr, int pro1Index, Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject);
void sweep_separation_complex_rot_memtest_cluster(int simItr, int pro1Index, Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject);
void sweep_separation_complex_rot_memtest_cluster_box(int simItr, int pro1Index, Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject);
void sweep_separation_complex_rot_memtest_cluster_sphere(int simItr, int pro1Index, Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject);

/*!
 * \brief Checks for overlap of proteins in solution.
 *
 * 3D version of sweep_separation_complex_rot_memtest
 */
void sweep_separation_complex_rot(int simItr, int pro1Index, Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject);
void sweep_separation_complex_rot_box(int simItr, int pro1Index, Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject);
void sweep_separation_complex_rot_sphere(int simItr, int pro1Index, Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject);
// void sweep_separation_complex_rot_cluster(int simItr, int pro1Index, const Parameters& params,
//     std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
//     const std::vector<ForwardRxn>& forwardRxns,
//     const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject);
// void sweep_separation_complex_rot_cluster_box(int simItr, int pro1Index, const Parameters& params,
//     std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
//     const std::vector<ForwardRxn>& forwardRxns,
//     const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject);
// void sweep_separation_complex_rot_cluster_sphere(int simItr, int pro1Index, const Parameters& params,
//     std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
//     const std::vector<ForwardRxn>& forwardRxns,
//     const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject);
