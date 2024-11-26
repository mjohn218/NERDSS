/*! \file system_setup.hpp
 * \brief Functions for initial system creation
 *
 * ### Created on 10/5/18 by Matthew Varga
 */
#pragma once

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"

/**
 * Determines whether two molecules are in the vicinity of each other.
 *
 * @param mol1         The first molecule to check.
 * @param mol2         The second molecule to check.
 * @param molTemplates A vector of MolTemplate objects specifying the binding radii of each molecule type.
 *
 * @return             True if the distance between the two molecules is less than the sum of their radii plus a safety margin of 2.
 */
bool are_molecules_in_vicinity(const Molecule& mol1, const Molecule& mol2, const std::vector<MolTemplate>& molTemplates);
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
void generate_coordinates(const Parameters &params,
                          std::vector<Molecule> &moleculeList,
                          std::vector<Complex> &complexList,
                          std::vector<MolTemplate> &molTemplateList,
                          const std::vector<ForwardRxn> &forwardRxns,
                          const Membrane &membraneObject,
                          std::string &coordinateFileName);

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
    const std::vector<ForwardRxn>& forwardRxns, int numDoubleBeforeAdd, int numMolTemplateBeforeAdd);

//void create(const MolTemplate& oneTemp, std::vector<int>& emptyMolList, std::vector<int>& emptyComList,
//	    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList);
Molecule initialize_molecule(int comIndex, const Parameters& params, const MolTemplate& molTemplate, const Membrane& membraneObject);
Complex initialize_complex(const Molecule& mol, const MolTemplate& molTemp);

// set up some important parameters for implicit-lipid model;
void initialize_paramters_for_implicitlipid_and_compartment_model(int& implicitlipidIndex, const Parameters& params, std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns, std::vector<TransmissionRxn>& transmissionRxns, std::vector<Molecule>& moleculeList,
    std::vector<MolTemplate>& molTemplateList, std::vector<Complex>& complexList, Membrane& membraneObject);

void initialize_paramters_for_implicitlipid_model(
    int& implicitlipidIndex, const Parameters& params,
    std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<Molecule>& moleculeList,
    std::vector<MolTemplate>& molTemplateList,
    std::vector<Complex>& complexList, Membrane& membraneObject);

//functions to generate new added molecules and complexes fo a restart simulation
void generate_coordinates_for_restart(Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, std::vector<MolTemplate>& molTemplateList,
    const std::vector<ForwardRxn>& forwardRxns, const Membrane& membraneObject, int numMolTemplateBeforeAdd, int numForwardRxnBdeforeAdd);

void create_molecule_and_complex_for_restart(MolTemplate& createdMolTemp, Parameters& params, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns, const Membrane& membraneObject);

Molecule initialize_molecule_for_restart(
    int index, Parameters& params, MolTemplate& molTemplate, const Membrane& membraneObject);

bool moleculeOverlapsForRestart(const Parameters& params, Molecule& createdMol,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject);

// function to determine molecule's isPoint and isRod
void determine_shape_molecule(std::vector<MolTemplate>& molTemplateList);

// function to initialize starting copy number for each state
void initialize_states(std::vector<Molecule>& moleculeList, std::vector<MolTemplate>& molTemplateList, Membrane& membraneObject);

void set_exclude_volume_bound(std::vector<ForwardRxn>& forwardRxns, std::vector<MolTemplate>& molTemplateList);