/*! \file system_setup.hpp
 * \brief Functions for initial system creation
 *
 * ### Created on 10/5/18 by Matthew Varga
 */
#pragma once

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"

/*! \defgroup SystemSetup SystemSetup
 */

/*!
 * \ingroup SystemSetup
 * \brief Check for overlap of two molecules' bounding spheres.
 * \param[in] mol1 Molecule 1 in the pairwise check.
 * \param[in] mol2 Molecule 2 in the pairwise check.
 * \param[in] molTemplateList List of provided MolTemplates.
 * \param[out] Boolean for if the bounding spheres overlap.
 *
 * Returns boolean for if the mol1 COM - mol2 COM vector is shorter than the combined radii of the two
 * molecules before spending the time determining if their interfaces ovelap.
 */
bool areInVicinity(const Molecule& mol1, const Molecule& mol2, const std::vector<MolTemplate>& molTemplateList);

/*!
 * \ingroup SystemSetup
 * \brief Function to check for overlap between proteins.
 * \param[in] params Simulation Parameters as provided by user.
 * \param[in] moleculeList List of all Molecules in the system.
 * \param[in] complexList List of all Complexes in the system.
 * \param[in] molTemplateList List of all provided MolTemplates.
 * \param[in] forwardRxns List of all provided ForwardRxns.
 *
 * ## Algorithm
 *   1. Extrapolate along the mol1 iface - mol2 iface vector using each molecule's interface, to see if they even
 * can be overlapping
 */
void generate_coordinates(const Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, const std::vector<MolTemplate>& molTemplateList,
			  const std::vector<ForwardRxn>& forwardRxns, const Membrane &membraneObject);

/*!
 * \ingroup SystemSetup
 * \brief Sets the RMax for the simulation, e.g. the maximum distance any protein can diffuse per timestep.
 * \param[in] params Simulation Parameters as provided by user.
 * \param[in] moleculeList List of all Molecules in the system.
 * \param[in] forwardRxns List of all provided ForwardRxns.
 *
 * Used to create the dimensions for the box cells.
 */
void set_rMaxLimit(Parameters& params, const std::vector<MolTemplate>& molTemplateList,
                   const std::vector<ForwardRxn>& forwardRxns);



//void create(const MolTemplate& oneTemp, std::vector<int>& emptyMolList, std::vector<int>& emptyComList,
//	    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList);
Molecule initialize_molecule(int comIndex, const Parameters& params, const MolTemplate& molTemplate, const Membrane &membraneObject);
Complex initialize_complex(const Molecule& mol, const MolTemplate& molTemp);

// set up some important parameters for implicit-lipid model;
void initialize_paramters_for_implicitlipid_model(int& implicitlipidIndex, const Parameters& params, std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList, 
                                             std::vector<MolTemplate>& molTemplateList, std::vector<Complex>& complexList, Membrane& membraneObject);
