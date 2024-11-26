#include <iostream>

#include "debug/debug.hpp"
#include "macro.hpp"
#include "mpi/mpi_function.hpp"
#include "classes/mpi_functions.hpp"

using namespace std;

void send_data_to_left_neighboring_ranks(
    MpiContext &mpiContext, long long int simItr, Parameters &params,
    SimulVolume &simulVolume, vector<int> &left, vector<Molecule> &moleculeList,
    vector<Complex> &complexList, vector<MolTemplate> &molTemplateList,
    Membrane &membraneObject) {
  // Send all the left ghost and left edge molecules and complexes to Left neighbor processor
  
  if (VERBOSE) cout << "send_data_to_left_neighboring_ranks begins" << endl;
  if (mpiContext.rank) {  // if not first rank
    // Store data for left rank into the array:
    set<int> complexesSet;
    mpiContext.nMPIArrayToLeft = 0;
    serialize_molecules(mpiContext, simulVolume, moleculeList, complexList,
                        mpiContext.MPIArrayToLeft, mpiContext.nMPIArrayToLeft,
                        complexesSet, true);
    mpiContext.check_buffer_size(mpiContext.nMPIArrayToLeft);
    serialize_complexes(mpiContext, simulVolume, membraneObject, complexesSet,
                        complexList, moleculeList, mpiContext.MPIArrayToLeft,
                        mpiContext.nMPIArrayToLeft);

    mpiContext.check_buffer_size(mpiContext.nMPIArrayToLeft);
    MPI_Isend(mpiContext.MPIArrayToLeft, mpiContext.nMPIArrayToLeft, MPI_CHAR,
              mpiContext.rank - 1, 0, MPI_COMM_WORLD,
              &mpiContext.requestSendToLeft);

    // Wait for sending pieces to left ranks to finish
    MPI_Status status;
    MPI_Wait(&mpiContext.requestSendToLeft, &status);
    if (COUNT_COMMUNICATIONS) {
      int count;  // Number of elements (bytes) send
      MPI_Get_count(&status, MPI_CHAR, &count);
      double sizeInMegabytes = count / 1024.0 / 1024.0;
      std::cout << "Rank " << mpiContext.rank << " sent " << sizeInMegabytes
                << " MB to rank " << mpiContext.rank - 1 << endl;
    }
  }

  mpiContext.increase_send_buffers();
  if (VERBOSE) cout << "send_data_to_left_neighboring_ranks ends" << endl;
}

void send_data_to_right_neighboring_ranks(
    MpiContext &mpiContext, long long int simItr, Parameters &params,
    SimulVolume &simulVolume, vector<int> &right,
    vector<Molecule> &moleculeList, vector<Complex> &complexList,
    vector<MolTemplate> &molTemplateList, Membrane &membraneObject) {
  if (VERBOSE) cout << "send_data_to_right_neighboring_ranks begins" << endl;
  if (mpiContext.rank < mpiContext.nprocs - 1) {  // if not last rank
    // Store data for rank into the array:
    set<int> complexesSet;
    mpiContext.nMPIArrayToRight = 0;
    serialize_molecules(mpiContext, simulVolume, moleculeList, complexList,
                        mpiContext.MPIArrayToRight, mpiContext.nMPIArrayToRight,
                        complexesSet, false);
    mpiContext.check_buffer_size(mpiContext.nMPIArrayToRight);

    serialize_complexes(mpiContext, simulVolume, membraneObject, complexesSet,
                        complexList, moleculeList, mpiContext.MPIArrayToRight,
                        mpiContext.nMPIArrayToRight);

    mpiContext.check_buffer_size(mpiContext.nMPIArrayToRight);
    MPI_Isend(mpiContext.MPIArrayToRight, mpiContext.nMPIArrayToRight, MPI_CHAR,
              mpiContext.rank + 1, 0, MPI_COMM_WORLD,
              &mpiContext.requestSendToRight);

    MPI_Status status;
    MPI_Wait(&mpiContext.requestSendToRight,
             &status);  // waiting for sending to finish
    if (COUNT_COMMUNICATIONS) {
      int count;  // Number of elements (bytes) send
      MPI_Get_count(&status, MPI_CHAR, &count);
      double sizeInMegabytes = count / 1024.0 / 1024.0;
      std::cout << "Rank " << mpiContext.rank << " sent " << sizeInMegabytes
                << " MB to rank " << mpiContext.rank + 1 << endl;
    }
  }

  mpiContext.increase_send_buffers();
  if (VERBOSE) cout << "send_data_to_right_neighboring_ranks ends" << endl;
}

void receive_right_neighborhood_zones(
    MpiContext &mpiContext, long long int simItr, SimulVolume &simulVolume,
    vector<int> &right, vector<Molecule> &moleculeList,
    vector<Complex> &complexList, vector<MolTemplate> &molTemplateList,
    Membrane &membraneObject, copyCounters &counterArrays) {
  if (VERBOSE) cout << "receive_right_neighborhood_zones begins" << endl;

  // If not last rank, receive data from the right rank,
  // deserialize received molecules,
  // update interfaces to match molecules on this rank,
  // and delete those molecules that were sent, but not received back:
  if (mpiContext.rank < mpiContext.nprocs - 1) {  // if not last rank
    MPI_Irecv(mpiContext.MPIArrayFromRight, mpiContext.recvBufferSize, MPI_CHAR,
              mpiContext.rank + 1, 0, MPI_COMM_WORLD,
              &mpiContext.requestRecvFromRight);
    MPI_Wait(
        &mpiContext.requestRecvFromRight,
        &mpiContext.statusRecvFromRight);  // waiting for receiving to finish
    MPI_Get_count(&mpiContext.statusRecvFromRight, MPI_CHAR,
                  &mpiContext.nMPIArrayFromRight);
    if (COUNT_COMMUNICATIONS) {
      double sizeInMegabytes = mpiContext.nMPIArrayFromRight / 1024.0 / 1024.0;
      std::cout << "Rank " << mpiContext.rank << " receive " << sizeInMegabytes
                << " MB from rank " << mpiContext.rank + 1 << endl;
    }
    if (VERBOSE) {
      mpiContext.print_spaces();
      cout << "-- simItr= " << simItr << ", Rank " << mpiContext.rank
           << " received from right: " << mpiContext.nMPIArrayFromRight << endl;
    }
    mpiContext.nMPIArrayFromRight = 0;
    vector<int> indices;
    vector<int> indicesEnteringMolecules;
    vector<int> indicesExitingMolecules;
    // Deserialize molecules and complexes from shared zone and update
    // interfaces, indices, myComIndex-es, etc. Updating molecules myComIndex-es
    // requires deserializating new complexes.
    deserialize_molecules(
        mpiContext, simulVolume, moleculeList, complexList, molTemplateList,
        membraneObject, counterArrays, indices, indicesEnteringMolecules,
        indicesExitingMolecules, mpiContext.MPIArrayFromRight,
        mpiContext.nMPIArrayFromRight);
    
    deserialize_complexes(mpiContext, moleculeList, complexList,
                                mpiContext.MPIArrayFromRight,
                                mpiContext.nMPIArrayFromRight);

    IDs_to_indices(mpiContext, moleculeList, indices);
  }

  mpiContext.increase_recv_buffers();
  if (VERBOSE) {
    mpiContext.print_spaces();
    cout << "-- simItr= " << simItr << ", Rank " << mpiContext.rank
         << " STOP receiving..." << endl;
  }
  if (VERBOSE) cout << "receive_right_neighborhood_zones ends" << endl;
}

void receive_left_neighborhood_zones(
    MpiContext &mpiContext, long long int simItr, SimulVolume &simulVolume,
    vector<int> &left, vector<Molecule> &moleculeList,
    vector<Complex> &complexList, vector<MolTemplate> &molTemplateList,
    Membrane &membraneObject, copyCounters &counterArrays) {
  if (VERBOSE) cout << "receive_left_neighborhood_zones begins" << endl;

  // If not first rank, receive data from the left rank,
  // deserialize received molecules,
  // update interfaces to match molecules on this rank,
  // and delete those molecules that were sent, but not received back:
  if (mpiContext.rank) {
    MPI_Irecv(mpiContext.MPIArrayFromLeft, mpiContext.recvBufferSize, MPI_CHAR,
              mpiContext.rank - 1, 0, MPI_COMM_WORLD,
              &mpiContext.requestRecvFromLeft);
    MPI_Wait(
        &mpiContext.requestRecvFromLeft,
        &mpiContext.statusRecvFromLeft);  // waiting for receiving to finish
    MPI_Get_count(&mpiContext.statusRecvFromLeft, MPI_CHAR,
                  &mpiContext.nMPIArrayFromLeft);
    if (COUNT_COMMUNICATIONS) {
      double sizeInMegabytes = mpiContext.nMPIArrayFromLeft / 1024.0 / 1024.0;
      std::cout << "Rank " << mpiContext.rank << " receive " << sizeInMegabytes
                << " MB from rank " << mpiContext.rank - 1 << endl;
    }
    if (VERBOSE) {
      mpiContext.print_spaces();
      cout << "-- simItr= " << simItr << ", Rank " << mpiContext.rank
           << " received from left: " << mpiContext.nMPIArrayFromLeft << endl;
    }
    mpiContext.nMPIArrayFromLeft = 0;
    vector<int> indices;
    vector<int> indicesEnteringMolecules;
    vector<int> indicesExitingMolecules;
    // Deserialize molecules and complexes from shared zone and update
    // interfaces, indices, myComIndex-es, etc. Updating molecules myComIndex-es
    // requires deserializating new complexes.
    deserialize_molecules(mpiContext, simulVolume, moleculeList, complexList,
                          molTemplateList, membraneObject, counterArrays,
                          indices, indicesEnteringMolecules,
                          indicesExitingMolecules, mpiContext.MPIArrayFromLeft,
                          mpiContext.nMPIArrayFromLeft);

    deserialize_complexes(mpiContext, moleculeList, complexList,
                          mpiContext.MPIArrayFromLeft,
                          mpiContext.nMPIArrayFromLeft);

    IDs_to_indices(mpiContext, moleculeList, indices);
  }

  mpiContext.increase_recv_buffers();
  if (VERBOSE) {
    mpiContext.print_spaces();
    cout << "-- simItr= " << simItr << ", Rank " << mpiContext.rank
         << " STOP receiving..." << endl;
  }
  if (VERBOSE) cout << "receive_left_neighborhood_zones ends" << endl;
}
