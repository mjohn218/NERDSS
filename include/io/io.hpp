/*! \file io.hpp

 * ### Created on 10/18/18 by Matthew Varga
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
#include "classes/class_Observable.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_SimulVolume.hpp"
#include "classes/class_copyCounters.hpp"

/*! \defgroup IO
 * \brief Functions for input/output of coordinates, parameters, trajectories, etc.
 */

#if defined(__APPLE__) || defined(__linux__)
std::ostream& bon(std::ostream& os);
std::ostream& boff(std::ostream& os);
#endif

std::ostream& linebreak(std::ostream& os);
std::ostream& llinebreak(std::ostream& os);
void write_psf(const Parameters& params, const std::vector<Molecule>& moleculeList,
    const std::vector<MolTemplate>& molTemplateList);
void write_crds(std::string name, const std::vector<Complex>& Complexlist, const std::vector<Molecule>& bases);
void write_crds(const std::vector<Complex>& Complexlist,
    const std::vector<Molecule>& bases); // overloaded crd dump for errors

/*!
 * \brief Writes the coordinates of all Molecules in the system to an XYZ coordinate file.
 */
void write_xyz(std::string filename, const Parameters& params, const std::vector<Molecule>& moleculeList,
    const std::vector<MolTemplate>& molTemplateList);

/*! \ingroup IO
 * \brief Writes/appends the coordinates of all Molecules in the system to the trajectory file.
 */
void write_traj(long long int iter, std::ofstream& trajFile, const Parameters& params, const std::vector<Molecule>& moleculeList,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject);

/*! \ingroup IO
 * \brief Debugging IO function to write an xyz file for a associating complexes.
 */
void write_xyz_assoc(std::string filename, const Complex& reactCom1, const Complex& reactCom2,
    const std::vector<Molecule>& moleculeList);

/*! \ingroup IO
 * \brief Debugging IO function to write an xyz file for a associating complexes to std cout.
 */
void write_xyz_assoc_cout(const Complex& reactCom1, const Complex& reactCom2,
    const std::vector<Molecule>& moleculeList);

/*! \ingroup IO
 * \brief Writes the coordinates of a complex to a file
 */
void write_complex_crds(
    std::string name, const Complex& complex1, const Complex& complex2, std::vector<Molecule>& moleculeList);

/*! \ingroup IO
 * \brief Writes out timestep information, such as the different molecular types in the system.
 *
 * TODO: This description needs more information
 */
void write_timestep_information(long long int simItr, std::ofstream& outFile, std::ofstream& molecTypesFile,
    std::ofstream& textTimeStatFile, const Parameters& params, std::vector<std::vector<int>>& molecTypesList,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, std::vector<MolTemplate>& molTemplateList);

/*! \ingroup IO
 * \brief Writes the number of each component Molecule in all complexes at the current timestep
 */
void write_complex_components(long long int simItr, std::ofstream& complexFile, const Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList);

/*! \ingroup IO
 * \brief Writes a plain text restart file at intervals specified in the Parameters file.
 *
 * This is a formatted text file, which is essentially illegible to the user, but it's not like they'd need to look at
 * it anyway.
 */
void write_restart(long long int simItr, std::ofstream& restartFile, const Parameters& params, const SimulVolume& simulVolume,
    const std::vector<Molecule>& moleculeList, const std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, const std::vector<CreateDestructRxn>& createDestructRxns,
    const std::map<std::string, int>& observablesList, const Membrane& membraneObject, const copyCounters& counterArrays);

/*! \ingroup IO
 * \brief Reads a restart file and sets up the simulation
 */
void read_restart(long long int& simItr, std::ifstream& restartFile, Parameters& params, SimulVolume& simulVolume,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<CreateDestructRxn>& createDestructRxns,
    std::map<std::string, int>& observablesList, Membrane& membraneObject, copyCounters& counterArrays);

/*! \ingroup IO
 * \ingroup SpeciesTracker
 * \brief Writes current number of each Observable to a CSV formatted file
 */
void write_observables(
    double simTime, std::ofstream& observablesFile, const std::map<std::string, int>& observablesList);

/*! \ingroup IO
 * \brief Writes a pdb file for the current frame
 */
void write_pdb(long long int simItr, unsigned frameNum, const Parameters& params, const std::vector<Molecule>& moleculeList,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject);

/*! \ingroup IO
 * \brief first line description of the MONODIMER output file.
 */
void init_print_dimers(std::ofstream& outfile, Parameters params, std::vector<MolTemplate>& molTemplateList);
void print_dimers(std::vector<Complex>& complexList, std::ofstream& outfile, int it, Parameters params,
    std::vector<MolTemplate>& molTemplateList);

/*! \ingroup IO
 * \brief Nbound pairs are counting all directly bound pairs of protein A and partner B. Does not matter if A or B are
 * also bound via other interfaces to other proteins.
 */
void init_NboundPairs(
    copyCounters& counterArray, std::ofstream& outfile, Parameters params, std::vector<MolTemplate>& molTemplateList, std::vector<Molecule>& moleculeList);

/*! \ingroup IO
 * \brief write: Nbound pairs are counting all directly bound pairs of protein A and partner B. Does not matter if A or
 * B are also bound via other interfaces to other proteins.
 */
void write_NboundPairs(copyCounters& counterArrays, std::ofstream& outfile, int it, const Parameters& params, std::vector<Molecule>& moleculeList);

/*! \ingroup IO
 * \brief print calculated histogram of distinct types of complexes, based on protein/lipid composition.
 */
double print_complex_hist(std::vector<Complex>& complexList, std::ofstream& outfile, int it, Parameters params,
    std::vector<MolTemplate>& molTemplateList, int nImplicitLipids);

/*! \ingroup IO
 * \brief initialize array of counterArrays.copyNumSpecies, based on the initial molecule Species and their interface
 * States.
 */
void init_counterCopyNums(copyCounters& counterArrays, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject, int totalSpeciesNum, Parameters& params);

/*! \ingroup IO
 * \brief writes out the names of species in the all_species.dat file.
 */

int init_speciesFile(std::ofstream& speciesFile, copyCounters& counterArrays, std::vector<MolTemplate>& molTemplateList, std::vector<ForwardRxn>& forwardRxns, Parameters& params);
/*! \ingroup IO
 * \brief Writes all the species in the system, from a copyCounter object
 */
void write_all_species(double simTime, std::ofstream& speciesFile, const copyCounters& counterArray);

/*! \ingroup IO
 * \brief Prints out information about all species in the system (Molecules and Complexes)
 */
void print_system_information(long long int simItr, std::ofstream& systemFile, const std::vector<Molecule>& moleculeList,
    const std::vector<Complex>& complexList, const std::vector<MolTemplate>& molTemplateList);
