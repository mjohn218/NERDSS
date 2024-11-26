#pragma once

#include <cstring>
#include <set>
#include <vector>

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Observable.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_SimulVolume.hpp"
#include "classes/class_copyCounters.hpp"
#include "mpi.h"
#include "split.cpp"

using namespace std;

/*! \defgroup mpi_function
 * \brief Functions for mpi implementation.
 */

/*! \ingroup mpi_function
 * \brief move complex
 */
void move_complexes_based_on_propagation(MpiContext &mpiContext,
                                         vector<Molecule> &moleculeList,
                                         SimulVolume &simulVolume,
                                         Membrane &membraneObject);

/*! \ingroup mpi_function
 * \brief Search for the molecule with the same ID within moleculeList
 */
int find_molecule(vector<Molecule> &moleculeList, int id);

/*! \ingroup mpi_function
 * \brief Search for the complex with the same ID within complexList
 */
int find_complex(vector<Complex> &complexList, int id);

/*! \ingroup mpi_function
 * \brief update the properties of molecule when it is in the shared zone
 */
inline void in_shared_zone(Molecule &mol) {
  // TODO: introduce shared complex if needed
}

/*! \ingroup mpi_function
 * \brief update the properties of the molecule when it is not in the share zone
 */
inline void not_in_shared_zone(Molecule &mol) {
  // TODO: remove complex sharing if needed
}

/*! \ingroup mpi_function
 * \brief update the copycounters when complex enter ghosted zone
 */
void update_copyCounters_enter_ghosted_zone(
    MpiContext &mpiContext, Molecule &mol, vector<MolTemplate> &molTemplateList,
    copyCounters &counterArrays, vector<Molecule> &moleculeList,
    SimulVolume &simulVolume);

/*! \ingroup mpi_function
 * \brief update the copycounters when complex exit ghosted zone
 */
void update_copyCounters_exit_ghosted_zone(MpiContext &mpiContext,
                                           Molecule &mol,
                                           vector<MolTemplate> &molTemplateList,
                                           copyCounters &counterArrays,
                                           vector<Molecule> &moleculeList,
                                           SimulVolume &simulVolume);

/*! \ingroup mpi_function
 * \brief update the properties of the molecule when it enters the ghosted zone
 */
void enter_ghosted_zone(MpiContext &mpiContext, Molecule &mol,
                        vector<MolTemplate> &molTemplateList,
                        copyCounters &counterArrays,
                        vector<Molecule> &moleculeList,
                        SimulVolume &simulVolume);

/*! \ingroup mpi_function
 * \brief update the properties of the molecule when it exits the ghosted zone
 */
void exit_ghosted_zone(MpiContext &mpiContext, Molecule &mol,
                       vector<MolTemplate> &molTemplateList,
                       copyCounters &counterArrays,
                       vector<Molecule> &moleculeList,
                       SimulVolume &simulVolume);

/*! \ingroup mpi_function
 * \brief update the properties of the molecule when it moves across the shared
 * zones
 */
void check_molecule_coordinates(
    MpiContext &mpiContext, vector<Molecule> &moleculeList,
    vector<int> &molIndexList, vector<Complex> &complexList,
    vector<MolTemplate> &molTemplateList, SimulVolume &simulVolume,
    Membrane &membraneObject, copyCounters &counterArrays);

/*! \ingroup mpi_function
 * \brief map local index to unique global id
 */
void indices_to_IDs(Molecule &mol, vector<Molecule> &moleculeList,
                    vector<Complex> &complexList);

/*! \ingroup mpi_function
 * \brief map mySubVolIndex of one molecule
 */
void update_mySubVolIndex(MpiContext &mpiContext, Molecule &mol,
                          Membrane &membraneObject, SimulVolume &simulVolume);

/*! \ingroup mpi_function
 * \brief check if a recieved mol enter the ghosted zone
 */
bool check_enter_ghosted_zone(MpiContext &mpiContext, int xBinOld, int xBinNew,
                              vector<int> &indicesEnteringMolecules,
                              SimulVolume &simulVolume);

/*! \ingroup mpi_function
 * \brief check if a recieved mol exit the ghosted zone
 */
bool check_exit_ghosted_zone(MpiContext &mpiContext, int xBinOld, int xBinNew,
                             vector<int> &indicesEnteringMolecules,
                             SimulVolume &simulVolume);

/*! \ingroup mpi_function
 * \brief serialize molecules
 */
void serialize_molecules(MpiContext &mpiContext, SimulVolume &simulVolume,
                         vector<Molecule> &moleculeList,
                         vector<Complex> &complexList, unsigned char *arrayRank,
                         int &nArrayRank, set<int> &complexesSet,
                         bool isLeft);

void serialize_ids(vector<int>& ids, char *arrayRank, int &nArrayRank);
void deserialize_ids(vector<int>& ids, char *arrayRank, int &nArrayRank);

void serialize_molecules_after_updating_subboxes(MpiContext &mpiContext, SimulVolume &simulVolume,
                         vector<Molecule> &moleculeList,
                         vector<Complex> &complexList, unsigned char *arrayRank,
                         int &nArrayRank, set<int> &complexesSet,
                         bool isLeft);

/*! \ingroup mpi_function
 * \brief serialize some molecules
 */
void serialize_molecules_partial(MpiContext &mpiContext,
                                 SimulVolume &simulVolume, vector<int> &indices,
                                 vector<Molecule> &moleculeList,
                                 vector<Complex> &complexList,
                                 unsigned char *arrayRank, int &nArrayRank,
                                 set<int> &complexesSet, int cellX);

/*! \ingroup mpi_function
 * \brief deserialize molecules
 */
void deserialize_molecules(MpiContext &mpiContext, SimulVolume &simulVolume,
                           vector<Molecule> &moleculeList,
                           vector<Complex> &complexList,
                           vector<MolTemplate> &molTemplateList,
                           Membrane &membraneObject,
                           copyCounters &counterArrays, vector<int> &indices,
                           vector<int> &indicesEnteringMolecules,
                           vector<int> &indicesExitingMolecules,
                           unsigned char *arrayRank, int &nArrayRank);

void deserialize_molecules_right(
    MpiContext &mpiContext, SimulVolume &simulVolume,
    vector<Molecule> &moleculeList, vector<Complex> &complexList,
    vector<MolTemplate> &molTemplateList, Membrane &membraneObject,
    copyCounters &counterArrays, vector<int> &indices,
    vector<int> &indicesEnteringMolecules, vector<int> &indicesExitingMolecules,
    unsigned char *arrayRank, int &nArrayRank);

/*! \ingroup mpi_function
 * \brief Serialize complexes from complexList with indices within complexesSet
 */
void serialize_complexes(MpiContext &mpiContext, SimulVolume &simulVolume,
                         Membrane &membraneObject, set<int> &complexesSet,
                         vector<Complex> &complexList,
                         vector<Molecule> &moleculeList,
                         unsigned char *arrayRank, int &nArrayRank);

/*! \ingroup mpi_function
 * \brief Deserialize complexes into complexList if they don't exist
 */
void deserialize_complexes(MpiContext &mpiContext,
                           vector<Molecule> &moleculeList,
                           vector<Complex> &complexList,
                           unsigned char *arrayRank, int &nArrayRank);

void deserialize_complexes_right(MpiContext &mpiContext,
                                 vector<Molecule> &moleculeList,
                                 vector<Complex> &complexList,
                                 unsigned char *arrayRank, int &nArrayRank);

/*! \ingroup mpi_function
 * \brief Map global unique id to local index
 */
void IDs_to_indices(MpiContext &mpiContext, vector<Molecule> &moleculeList,
                    vector<int> &indices);

/*! \ingroup mpi_function
 * \brief disconnect molecule partners
 */
void disconnect_molecule_partners(unsigned &targMolIndex, Molecule &mol,
                                  vector<Molecule> &moleculeList);

/*! \ingroup mpi_function
 * \brief Delete all molecules that were in the shared zone before sending to
 * the neighbor rank, but were not received back:
 */
void delete_disappeared_molecules(
    MpiContext &mpiContext, SimulVolume &simulVolume, Membrane &membraneObject,
    vector<Molecule> &moleculeList, vector<Complex> &complexList,
    vector<int> &region, bool isLeft, vector<MolTemplate> &molTemplateList);

/*! \ingroup mpi_function
 * \brief Delete all complexes that were in the shared zone before sending to
 * the neighbor rank, but were not received back:
 */
void delete_disappeared_complexes_partial(MpiContext &mpiContext,
                                          vector<Molecule> &moleculeList,
                                          vector<Complex> &complexList,
                                          bool left);

/*! \ingroup mpi_function
 * \brief Delete all complexes that were in the shared zone before sending to
 * the neighbor rank, but were not received back:
 */
void delete_disappeared_complexes(MpiContext &mpiContext,
                                  vector<Molecule> &moleculeList,
                                  vector<Complex> &complexList);

/*! \ingroup mpi_function
 * \brief Send data to neighbor rank
 */
void send_data_to_neighboring_ranks(MpiContext &mpiContext,
                                    long long int simItr, Parameters &params,
                                    SimulVolume &simulVolume,
                                    vector<Molecule> &moleculeList,
                                    vector<Complex> &complexList,
                                    vector<MolTemplate> &molTemplateList,
                                    Membrane &membraneObject);

/*! \ingroup mpi_function
 * \brief Send data to left neighbor rank
 */
void send_data_to_left_neighboring_ranks(
    MpiContext &mpiContext, long long int simItr, Parameters &params,
    SimulVolume &simulVolume, vector<int> &left, vector<Molecule> &moleculeList,
    vector<Complex> &complexList, vector<MolTemplate> &molTemplateList,
    Membrane &membraneObject);

void send_molids_to_left_neighboring_ranks(MpiContext &mpiContext, SimulVolume &simulVolume, vector<Molecule> &moleculeList);

void send_data_to_left_neighboring_ranks_after_updating_subboxes(
    MpiContext &mpiContext, long long int simItr, Parameters &params,
    SimulVolume &simulVolume, vector<int> &left, vector<Molecule> &moleculeList,
    vector<Complex> &complexList, vector<MolTemplate> &molTemplateList,
    Membrane &membraneObject);

/*! \ingroup mpi_function
 * \brief Send data to right neighbor rank
 */
void send_data_to_right_neighboring_ranks(
    MpiContext &mpiContext, long long int simItr, Parameters &params,
    SimulVolume &simulVolume, vector<int> &right,
    vector<Molecule> &moleculeList, vector<Complex> &complexList,
    vector<MolTemplate> &molTemplateList, Membrane &membraneObject);

/*! \ingroup mpi_function
 * \brief Receive shared zones from both left and right neighboring rank,
 * deserialize received molecules, update interfaces to match molecules on this
 * rank, and delete those molecules that were sent, but not received back
 */
void receive_neighborhood_zones(
    MpiContext &mpiContext, long long int simItr, long long int startSimItr,
    long long int stopSimItr, SimulVolume &simulVolume,
    vector<Molecule> &moleculeList, vector<Complex> &complexList,
    vector<MolTemplate> &molTemplateList, Membrane &membraneObject,
    copyCounters &counterArrays);

/*! \ingroup mpi_function
 * \brief Receive shared zones from right neighboring rank, deserialize received
 * molecules, update interfaces to match molecules on this rank, and delete
 * those molecules that were sent, but not received back
 */
void receive_right_neighborhood_zones(
    MpiContext &mpiContext, long long int simItr, SimulVolume &simulVolume,
    vector<int> &right, vector<Molecule> &moleculeList,
    vector<Complex> &complexList, vector<MolTemplate> &molTemplateList,
    Membrane &membraneObject, copyCounters &counterArrays);

void receive_ids(MpiContext &mpiContext, SimulVolume &simulVolume, vector<Molecule> &moleculeList);

void receive_right_neighborhood_zones_after_updating_subboxes(
    MpiContext &mpiContext, long long int simItr, SimulVolume &simulVolume,
    vector<int> &right, vector<Molecule> &moleculeList,
    vector<Complex> &complexList, vector<MolTemplate> &molTemplateList,
    Membrane &membraneObject, copyCounters &counterArrays);

/*! \ingroup mpi_function
 * \brief Receive shared zones from left neighboring rank, deserialize received
 * molecules, update interfaces to match molecules on this rank, and delete
 * those molecules that were sent, but not received back
 */
void receive_left_neighborhood_zones(
    MpiContext &mpiContext, long long int simItr, SimulVolume &simulVolume,
    vector<int> &left, vector<Molecule> &moleculeList,
    vector<Complex> &complexList, vector<MolTemplate> &molTemplateList,
    Membrane &membraneObject, copyCounters &counterArrays);

/*! \ingroup mpi_function
 * \brief Prepare rank data
 */
int prepare_rank_data(int tempRank, vector<Molecule> &moleculeList,
                      SimulVolume &simulVolume, Membrane &membraneObject,
                      vector<MolTemplate> &molTemplateList, Parameters &params,
                      vector<ForwardRxn> &forwardRxns,
                      vector<BackRxn> &backRxns,
                      vector<CreateDestructRxn> &createDestructRxns,
                      copyCounters &counterArrays, MpiContext &mpiContext,
                      vector<Complex> &complexList, ofstream &pairOutfile,
                      unsigned char *arrayRank);

/*! \ingroup mpi_function
 * \brief Prepare data structures for parallel execution
 */
void prepare_data_structures_for_parallel_execution(
    vector<Molecule> &moleculeList, SimulVolume &simulVolume,
    Membrane &membraneObject, vector<MolTemplate> &molTemplateList,
    Parameters &params, vector<ForwardRxn> &forwardRxns,
    vector<BackRxn> &backRxns, vector<CreateDestructRxn> &createDestructRxns,
    copyCounters &counterArrays, MpiContext &mpiContext,
    vector<Complex> &complexList, ofstream &pairOutfile);

/*! \ingroup mpi_function
 * \brief Used for the writing of output files
 */
bool is_ghosted(Molecule &mol, MpiContext &mpiContext,
                SimulVolume &simulVolume);

bool is_owned_by_processor(Molecule &mol, MpiContext &mpiContext,
                           SimulVolume &simulVolume);