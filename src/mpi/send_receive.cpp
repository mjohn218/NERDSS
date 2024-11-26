#include <iostream>

#include "debug/debug.hpp"
#include "macro.hpp"
#include "mpi/mpi_function.hpp"

using namespace std;

void send_data_to_left_neighboring_ranks(
    MpiContext &mpiContext, long long int simItr, Parameters &params,
    SimulVolume &simulVolume, vector<int> &left, vector<Molecule> &moleculeList,
    vector<Complex> &complexList, vector<MolTemplate> &molTemplateList,
    Membrane &membraneObject) {
  if (VERBOSE) cout << "send_data_to_left_neighboring_ranks begins" << endl;
  // simulVolume.update_memberMolLists(params, moleculeList, complexList,
  //                                   molTemplateList, membraneObject, simItr,
  //                                   mpiContext);

  // Send pieces to left ranks
  if (mpiContext.rank) {  // if not first rank
    // Store portion of molecules for left rank into the array:
    set<int> complexesSet;
    set<int> moleculesSet;

    int xItr = 0;
    for (int j = 0; j < simulVolume.numSubCells.y; j++) {
      for (int k = 0; k < simulVolume.numSubCells.z; k++) {
        int currBin = xItr + j * simulVolume.numSubCells.x +
                      k * simulVolume.numSubCells.x * simulVolume.numSubCells.y;
        for (auto& molIdx : simulVolume.subCellList[currBin].memberMolList) {
          auto& mol = moleculeList[molIdx];
          if (mol.myComIndex == -1 || mol.isImplicitLipid == true || mol.isEmpty == true) continue;
          mol.isGhosted = true;
          moleculesSet.insert(molIdx);
          complexesSet.insert(mol.myComIndex);
        }
      }
    }

    xItr = 1;
    for (int j = 0; j < simulVolume.numSubCells.y; j++) {
      for (int k = 0; k < simulVolume.numSubCells.z; k++) {
        int currBin = xItr + j * simulVolume.numSubCells.x +
                      k * simulVolume.numSubCells.x * simulVolume.numSubCells.y;
        for (auto& molIdx : simulVolume.subCellList[currBin].memberMolList) {
          auto& mol = moleculeList[molIdx];
          if (mol.myComIndex == -1 || mol.isImplicitLipid == true || mol.isEmpty == true) continue;
          mol.isShared = true;
          moleculesSet.insert(molIdx);
          complexesSet.insert(mol.myComIndex);
          for (auto& memIdx : complexList[mol.myComIndex].memberList) {
            moleculeList[memIdx].isShared = true;
            moleculesSet.insert(memIdx);
          }
        }
      }
    }

    mpiContext.nMPIArrayToLeft = 0;
    serialize_molecules(mpiContext, simulVolume, moleculeList, complexList,
                        mpiContext.MPIArrayToLeft, mpiContext.nMPIArrayToLeft,
                        complexesSet, moleculesSet);
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
  // simulVolume.update_memberMolLists(params, moleculeList, complexList,
  //                                   molTemplateList, membraneObject, simItr,
  //                                   mpiContext);

  // Send pieces to right ranks
  if (mpiContext.rank < mpiContext.nprocs - 1) {  // if not last rank
    // Store portion of molecules for rank into the array:
    set<int> complexesSet;
    set<int> moleculesSet;

    int xItr = simulVolume.numSubCells.x - 1;
    for (int j = 0; j < simulVolume.numSubCells.y; j++) {
      for (int k = 0; k < simulVolume.numSubCells.z; k++) {
        int currBin = xItr + j * simulVolume.numSubCells.x +
                      k * simulVolume.numSubCells.x * simulVolume.numSubCells.y;
        for (auto& molIdx : simulVolume.subCellList[currBin].memberMolList) {
          auto& mol = moleculeList[molIdx];
          if (mol.myComIndex == -1 || mol.isImplicitLipid == true || mol.isEmpty == true) continue;
          mol.isGhosted = true;
          moleculesSet.insert(molIdx);
          complexesSet.insert(mol.myComIndex);
        }
      }
    }

    xItr = simulVolume.numSubCells.x - 2;
    for (int j = 0; j < simulVolume.numSubCells.y; j++) {
      for (int k = 0; k < simulVolume.numSubCells.z; k++) {
        int currBin = xItr + j * simulVolume.numSubCells.x +
                      k * simulVolume.numSubCells.x * simulVolume.numSubCells.y;
        for (auto& molIdx : simulVolume.subCellList[currBin].memberMolList) {
          auto& mol = moleculeList[molIdx];
          if (mol.myComIndex == -1 || mol.isImplicitLipid == true || mol.isEmpty == true) continue;
          mol.isShared = true;
          moleculesSet.insert(molIdx);
          complexesSet.insert(mol.myComIndex);
          for (auto& memIdx : complexList[mol.myComIndex].memberList) {
            moleculeList[memIdx].isShared = true;
            moleculesSet.insert(memIdx);
          }
        }
      }
    }

    mpiContext.nMPIArrayToRight = 0;
    serialize_molecules(mpiContext, simulVolume, moleculeList, complexList,
                        mpiContext.MPIArrayToRight, mpiContext.nMPIArrayToRight,
                        complexesSet, moleculesSet);
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

  if (DEBUG) {
    debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                     "//before deserializing()");
  }
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
    if (DEBUG) {
      DEBUG_FIND_MOL("0_ (before right deserialize_molecules)");
    }
    deserialize_molecules_right(
        mpiContext, simulVolume, moleculeList, complexList, molTemplateList,
        membraneObject, counterArrays, indices, indicesEnteringMolecules,
        indicesExitingMolecules, mpiContext.MPIArrayFromRight,
        mpiContext.nMPIArrayFromRight);  // last parameter is ghostStripe

    //        DEBUG_FIND_MOL("3");
    // Remove molecules that were not updated (received back from the neighbor
    // rank):
    if (DEBUG) {
      DEBUG_FIND_MOL("1_ (after right deserialize_molecules)");
    }
    delete_disappeared_molecules(
        mpiContext, simulVolume, membraneObject, moleculeList, complexList,
        right, molTemplateList);
    if (DEBUG) {
      DEBUG_FIND_MOL("1.4_ (after delete_disappeared_molecules)");
      DEBUG_FIND_COMPLEX("1.4_ (after delete_disappeared_molecules)");
    }
    deserialize_complexes_right(mpiContext, moleculeList, complexList,
                                mpiContext.MPIArrayFromRight,
                                mpiContext.nMPIArrayFromRight);
    if (DEBUG) {
      DEBUG_FIND_MOL("1.5_ (after deserialize_complexes)");
      DEBUG_FIND_COMPLEX("1.5_ (after deserialize_complexes)");
    }
    IDs_to_indices(mpiContext, moleculeList, indices);

    // for (auto i : indicesEnteringMolecules) {
    //   update_copyCounters_enter_ghosted_zone(mpiContext, moleculeList[i],
    //                                          molTemplateList, counterArrays,
    //                                          moleculeList, simulVolume);
    // }
    // for (auto i : indicesExitingMolecules) {
    //   update_copyCounters_exit_ghosted_zone(mpiContext, moleculeList[i],
    //                                         molTemplateList, counterArrays,
    //                                         moleculeList, simulVolume);
    // }

    if (DEBUG) {
      DEBUG_FIND_MOL("2_ (after IDs_to_indices)");
      DEBUG_FIND_COMPLEX("2_ (after IDs_to_indices)");
    }
  }
  // Delete all complexes that were just in the shared zone before sending to
  // the neighbor rank, but were not received back: only process the right bound
  // of the rank
  delete_disappeared_complexes_partial(mpiContext, moleculeList, complexList,
                                       false);

  if (DEBUG) {
      DEBUG_FIND_MOL("2_ (after delete_disappeared_complexes_partial)");
      DEBUG_FIND_COMPLEX("2_ (after delete_disappeared_complexes_partial)");
    }

  if (DEBUG) {
    debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                     "//after delete com partial()");
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
  if (DEBUG) {
    debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                     "//before deserializing()");
  }
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
    if (DEBUG) {
      DEBUG_FIND_MOL("0 (before left deserialize_molecules)");
    }
    deserialize_molecules(mpiContext, simulVolume, moleculeList, complexList,
                          molTemplateList, membraneObject, counterArrays,
                          indices, indicesEnteringMolecules,
                          indicesExitingMolecules, mpiContext.MPIArrayFromLeft,
                          mpiContext.nMPIArrayFromLeft);
    // Remove molecules that were not updated (received back from the neighbor
    // rank): DEBUG_FIND_MOL("1 (after left deserialize_molecules)");
    if (DEBUG) {
      DEBUG_FIND_MOL("1 (after left deserialize_molecules)");
    }
    delete_disappeared_molecules(mpiContext, simulVolume, membraneObject,
                                 moleculeList, complexList, left,
                                 molTemplateList);
    // DEBUG_FIND_MOL("1.4 (after delete_disappeared_molecules)");
    if (DEBUG) {
      DEBUG_FIND_MOL("1.4 (after delete_disappeared_molecules)");
    }
    if (DEBUG) {
      DEBUG_FIND_COMPLEX("1.4 (before deserialize_complexes)");
    }
    deserialize_complexes(mpiContext, moleculeList, complexList,
                          mpiContext.MPIArrayFromLeft,
                          mpiContext.nMPIArrayFromLeft);
    if (DEBUG) {
      DEBUG_FIND_MOL("1.5 (after deserialize_complexes)");
      DEBUG_FIND_COMPLEX("1.5 (after deserialize_complexes)");
    }
    IDs_to_indices(mpiContext, moleculeList, indices);

    // for (auto i : indicesEnteringMolecules) {
    //   update_copyCounters_enter_ghosted_zone(mpiContext, moleculeList[i],
    //                                          molTemplateList, counterArrays,
    //                                          moleculeList, simulVolume);
    // }
    // for (auto i : indicesExitingMolecules) {
    //   update_copyCounters_exit_ghosted_zone(mpiContext, moleculeList[i],
    //                                         molTemplateList, counterArrays,
    //                                         moleculeList, simulVolume);
    // }

    if (DEBUG) {
      DEBUG_FIND_MOL("2 (after IDs_to_indices)");
    }
  }

  // Delete all complexes that were just in the shared zone before sending to
  // the neighbor rank, but were not received back: only process the left bound
  // of the rank
  delete_disappeared_complexes_partial(mpiContext, moleculeList, complexList,
                                       true);
  if (DEBUG) {
    DEBUG_FIND_COMPLEX("1.6 (after delete_disappeared_complexes_partial)");
  }

  mpiContext.increase_recv_buffers();
  if (VERBOSE) {
    mpiContext.print_spaces();
    cout << "-- simItr= " << simItr << ", Rank " << mpiContext.rank
         << " STOP receiving..." << endl;
  }
  if (VERBOSE) cout << "receive_left_neighborhood_zones ends" << endl;
}
