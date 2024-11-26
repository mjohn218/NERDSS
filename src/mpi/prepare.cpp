#include <cmath>
#include <iostream>

#include "debug/debug.hpp"
#include "io/io.hpp"
#include "macro.hpp"
#include "mpi/mpi_function.hpp"

using namespace std;

/**
 * @brief Prepares data structures for parallel execution
 *
 * Following data structures need updating before spreading data to ranks:
 * Molecules have to be divided onto ranks.
 * Molecule indices (index field) have to be updated.
 * Cell indices of molecules have to be updated, as cell indices differ in the
 * serial and the parallel version. Cells have to be divided onto ranks. Cell
 * neighborhood lists have to be updated, as indices at ranks differ from those
 * from the serial version. Cell members have to be updated, as indices of
 * molecules differ from those from the serial version. Formula (in Simulation
 * volume) for x-offset on ranks have to be updated, so that neighborhood cells
 * are addressed properly. List of complexes have to be updated, as portions of
 * complexes will be sent to ranks. Member list of complexes have to be updated,
 * as molecule indices differ. Simulation volume number of cells per x axis for
 * ranks have to be updated, as only portion of x-domain will be received by any
 * rank. Simulation volume number of cells (total) have to be updated. If
 * implicit lipid, then molecule 0 should be transferred to all ranks. Complexes
 * that are shared between ranks have to be divided for each rank. IDs of
 * molecules and complexes have to be introduced and must be unique system-wide.
 *
 * @param moleculeList Reference to the list of Molecule objects
 * @param simulVolume Reference to the SimulVolume object
 * @param membraneObject Reference to the Membrane object
 * @param molTemplateList Reference to the list of MolTemplate objects
 * @param params Reference to the Parameters object
 * @param forwardRxns Reference to the list of ForwardRxn objects
 * @param backRxns Reference to the list of BackRxn objects
 * @param createDestructRxns Reference to the list of CreateDestructRxn objects
 * @param counterArrays Reference to the copyCounters object
 * @param mpiContext Reference to the MpiContext object
 * @param complexList Reference to the list of Complex objects
 * @param pairOutfile Reference to the output file stream for writing pairs
 */
void prepare_data_structures_for_parallel_execution(
    vector<Molecule> &moleculeList, SimulVolume &simulVolume,
    Membrane &membraneObject, vector<MolTemplate> &molTemplateList,
    Parameters &params, vector<ForwardRxn> &forwardRxns,
    vector<BackRxn> &backRxns, vector<CreateDestructRxn> &createDestructRxns,
    copyCounters &counterArrays, MpiContext &mpiContext,
    vector<Complex> &complexList, ofstream &pairOutfile) {
  cout << "*************** PREPARING DATA STRUCTURES FOR PARALLEL EXECUTION "
          "**************** rank:"
       << mpiContext.rank << endl;

  // Reserve memory for an array for serialized data to be sent from rank 0 to
  // others
  unsigned char *arrayRank =
      (unsigned char *)malloc(RANK0_BUFFER_SIZE);  // TODO: determine size

  // Testing serialization and deserialization is done just for debugging
  // purposes, solely on rank 0:
  test_serialization(mpiContext, moleculeList, simulVolume, membraneObject,
                     molTemplateList, params, forwardRxns, backRxns,
                     createDestructRxns, counterArrays, complexList, arrayRank);
  debug_molecule_complex_missmatch(
      mpiContext, moleculeList, complexList,
      "//before prepare_data_structures_for_parallel_execution()");

  int totalxBins = simulVolume.numSubCells.x;

  // Rank 0 prepares data for each rank (tempRank), including itself
  if (!mpiContext.rank) {
    // Create customized simulation volume, cell list, molecule list, complex
    // list, etc. for each rank. Loop from the last to rank 0, so that prepared
    // array for rank 0 doesn't need to be sent for deserialization, nor it gets
    // overwritten
    for (int tempRank = mpiContext.nprocs - 1; tempRank >= 0; tempRank--) {
      int nArrayRank = prepare_rank_data(
          tempRank, moleculeList, simulVolume, membraneObject, molTemplateList,
          params, forwardRxns, backRxns, createDestructRxns, counterArrays,
          mpiContext, complexList, pairOutfile, arrayRank);

      // Send all data prepared for rank tempRank
      if (tempRank) {  // send only to other ranks, not from rank 0 to rank 0
        MPI_Send(&nArrayRank, 1, MPI_INT, tempRank, 0,
                 MPI_COMM_WORLD);  // send the length of the array
        MPI_Send(arrayRank, nArrayRank, MPI_CHAR, tempRank, 0,
                 MPI_COMM_WORLD);  // sending the array
      }
    }
  }

  // Receive the array from rank 0 (if needed) with packed data and deserialize
  // it
  int nBytes = 0;         // length of the array to receive
  if (mpiContext.rank) {  // if rank is 0, no need for receiving - the arrayRank
                          // is already prepared
    MPI_Recv(&nBytes, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);  // receive the length of the array
    MPI_Recv(arrayRank, nBytes, MPI_CHAR, 0, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);  // receive the array
  }

  // Deserialize data from rank 0 on each rank, including 0;
  // Rank 0 holds initialized data for the whole system, i.e. all cells,
  // molecules, complexes, etc. Clear existing arrays, so that deserialized ones
  // are the only ones:

  // The order of deserialization must be exactly the same as the order of
  // serialization. Starting from char at address 0, all vectors and objects
  // should be deserialized from array arrayRank, where they are kept one after
  // another; therefore, each deserialization increases nArrayRank address.
  int nArrayRank = 0;
  complexList.clear();

  // arrayRank is an array, and acts as a pointer to the first byte of the array
  // arrayRank, nArrayRank acts as a pointer to nArrayRank byte of the array,
  deserialize_abstract_vector<Complex>(complexList, arrayRank,
                                       nArrayRank);  // deserialize complexes

  moleculeList.clear();
  deserialize_abstract_vector<Molecule>(moleculeList, arrayRank,
                                        nArrayRank);  // deserialize molecules

  SimulVolume simulVolumeClean;
  simulVolume = simulVolumeClean;
  simulVolume.deserialize(arrayRank, nArrayRank);

  mpiContext.deserialize(arrayRank, nArrayRank);

  Membrane membraneObjectClean;
  membraneObject = membraneObjectClean;
  membraneObject.deserialize(arrayRank, nArrayRank);

  molTemplateList.clear();
  deserialize_abstract_vector<MolTemplate>(molTemplateList, arrayRank,
                                           nArrayRank);

  Parameters paramsClean;
  params = paramsClean;
  params.deserialize(arrayRank, nArrayRank);

  forwardRxns.clear();
  deserialize_abstract_vector<ForwardRxn>(forwardRxns, arrayRank, nArrayRank);

  backRxns.clear();
  deserialize_abstract_vector<BackRxn>(backRxns, arrayRank, nArrayRank);

  createDestructRxns.clear();
  deserialize_abstract_vector<CreateDestructRxn>(createDestructRxns, arrayRank,
                                                 nArrayRank);

  copyCounters counterArraysClean;
  counterArrays = counterArraysClean;
  counterArrays.deserialize(arrayRank, nArrayRank);

  // Recalculate the membraneObject.totalSA and membraneObject.nSites and
  // membraneObject.numberOfFreeLipidsEachState divide
  // membraneObject.numberOfFreeLipidsEachState onto ranks proportional to
  // number of bins per rank
  double ratio =
      1.0 * (mpiContext.endCell - mpiContext.startCell + 1) / totalxBins;

  if (membraneObject.isSphere) {
  } else {
    membraneObject.totalSA =
        membraneObject.waterBox.x * membraneObject.waterBox.y * ratio;
    membraneObject.waterBox.volume = membraneObject.waterBox.volume * ratio;
    membraneObject.waterBox.xLeft =
        mpiContext.xLeft * (membraneObject.waterBox.x) -
        (membraneObject.waterBox.x / 2.0);
    membraneObject.waterBox.xRight =
        mpiContext.xRight * (membraneObject.waterBox.x) -
        (membraneObject.waterBox.x / 2.0);
  }

  if (membraneObject.implicitLipid == true) {
    membraneObject.nSites = round(membraneObject.nSites * ratio);
    if (params.fromRestart == false) {
      for (auto &iface :
           molTemplateList[moleculeList[membraneObject.implicitlipidIndex]
                               .molTypeIndex]
               .interfaceList) {
        for (auto &state : iface.stateList) {
          if (&state - &iface.stateList[0] ==
              0) {  // if is the first state, free lipids is initial copies of
                    // IL
            membraneObject.numberOfFreeLipidsEachState[0] =
                membraneObject.nSites;
          } else {  // others are zero
          }
        }
      }
    }
  }

  free(arrayRank);

  // Store addresses of vectors and objects to be potentially used for debugging
  // purposes:
  mpiContext.membraneObject = &membraneObject;
  mpiContext.simulVolume = &simulVolume;
  mpiContext.moleculeList = &moleculeList;
  mpiContext.complexList = &complexList;

  // Update isGhosted for each molecule:
  for (auto &mol : moleculeList) {
    if (mol.isImplicitLipid == true) {
      mol.isGhosted = false;
      continue;
    }
    int xBin = get_x_bin(mpiContext, mol);
    mol.isGhosted = false;
    if ((mpiContext.rank) &&  // if not first rank
        (xBin == 0))
      mol.isGhosted = true;
    if ((mpiContext.rank < mpiContext.nprocs - 1) &&  // if not first rank
        (xBin == simulVolume.numSubCells.x - 1))
      mol.isGhosted = true;
  }

  // For the left-right division model, in the first step, all ranks (except the
  // first rank) will send their left share zones to the left neighbor rank. So,
  // all the right share zones of each rank (except the last rank) should set
  // the receivedFromNighborRank to be False:
  for (auto &mol : moleculeList) {
    bool receivedFromNeighborRank = true;
    int xBin = get_x_bin(mpiContext, mol);

    if (mpiContext.rank < mpiContext.nprocs - 1) {
      if ((xBin == simulVolume.numSubCells.x - 1) ||
          (xBin == simulVolume.numSubCells.x - 2)) {
        receivedFromNeighborRank = false;
      }
    }

    if (mol.isImplicitLipid) {
      receivedFromNeighborRank = true;
    }

    if (!receivedFromNeighborRank) {
      mol.receivedFromNeighborRank = false;
      complexList[mol.myComIndex].receivedFromNeighborRank = false;
      complexList[mol.myComIndex].deleteIfNotReceivedBack = true;
    }
  }  // end loop moleculeList to update receivedFromNeighborRank

  debug_molecule_complex_missmatch(
      mpiContext, moleculeList, complexList,
      "//prepare_data_structures_for_parallel_execution()");

  //  Count how many interface connections are on each rank by calling
  //  init_NboundPairs on each rank.
  init_NboundPairs(
      counterArrays, pairOutfile, params, molTemplateList,
      moleculeList);  // initializes to zero, re-calculated for a restart!!

  Molecule::maxID = Molecule::maxID + (INT_MAX - Molecule::maxID) /
                                          mpiContext.nprocs * mpiContext.rank;
  Complex::maxID = Complex::maxID + (INT_MAX - Complex::maxID) /
                                        mpiContext.nprocs * mpiContext.rank;
}

int prepare_rank_data(int tempRank, vector<Molecule> &moleculeList,
                      SimulVolume &simulVolume, Membrane &membraneObject,
                      vector<MolTemplate> &molTemplateList, Parameters &params,
                      vector<ForwardRxn> &forwardRxns,
                      vector<BackRxn> &backRxns,
                      vector<CreateDestructRxn> &createDestructRxns,
                      copyCounters &counterArrays, MpiContext &mpiContext,
                      vector<Complex> &complexList, ofstream &pairOutfile,
                      unsigned char *arrayRank) {
  // Initialize domain of cells for each rank, including ghosted ones:
  mpiContext.init_x_domain_and_offset(simulVolume.numSubCells.x, tempRank);

  int nArrayRank = 0;  ///< number of bytes in arrayRank prepared for tempRank
  int nComplexes = 0;  ///< number of complexes to transfer

  // Preparing simulation volume for tempRank
  SimulVolume simulVolumeRank{};  // Custom simulation volume for tempRank
  simulVolumeRank.maxNeighbors = simulVolume.maxNeighbors;
  // Divide number of cells per x axis (simulVolume.numSubCells.x)
  // onto number of ranks (mpiContext.nprocs):
  simulVolumeRank.numSubCells.x =
      mpiContext.endGhosted - mpiContext.startGhosted + 1;
  simulVolumeRank.numSubCells.y = simulVolume.numSubCells.y;
  simulVolumeRank.numSubCells.z = simulVolume.numSubCells.z;
  simulVolumeRank.subCellSize = simulVolume.subCellSize;
  simulVolumeRank.subCellList.clear();

  // Initializing auxiliary data

  // Assign unigue ID to each molecule and complex:
  for (size_t i = 0; i < moleculeList.size(); i++) moleculeList[i].id = i;
  for (size_t i = 0; i < complexList.size(); i++) complexList[i].id = i;

  // Calculate starting IDs that ranks will be assigning based on their rank
  // numbers:
  //..........|----rank0----|----rank1----|----rank2----|----rank3----|
  Molecule::maxID = moleculeList.size();
  Complex::maxID = complexList.size();

  vector<Molecule> moleculeListRank;  ///< list of all molecules for tempRank
  // Init molecule index maping from serial version index to parallel version
  // index; needed for updating numEachMol before transferring each complex
  int mapSerialToParallelMolecule[moleculeList.size()];
  for (auto &it : mapSerialToParallelMolecule)
    it = -1;  // molecule does not belong to temprank

  // Calculate which complexes should be transferred to tempRank:
  bool transferComplex[complexList.size()];  // true <=> complex at least
                                             // partionaly belongs to tempRank
  // Initially, assume no complex should be transferred to tempRank
  // Later, as molecules are found for tempRank, their complexes will be marked
  // for transferring to tempRank
  for (unsigned iComplex = 0; iComplex < complexList.size(); iComplex++)
    transferComplex[iComplex] = false;

  // If implicid lipid exists, it is on index 0;
  // therefore, molecule 0 should be copied to each rank
  if ((moleculeList.size()) && (params.implicitLipid == true)) {
    mapSerialToParallelMolecule[0] = 0;
    moleculeListRank.push_back(moleculeList[0]);
    transferComplex[0] = true;
  }

  // Cells for particular rank have different indices than initial cells
  // Before sending cells to a rank, neighborhood cell indices have to updated
  // Vector mapSerialToParallelCell stores mapping from initial (serial) cells
  // to cells prepared for a tempRank:
  vector<int> mapSerialToParallelCell(simulVolume.subCellList.size());
  // Vector mapParallelToSerialCell stores mapping from cells prepared for a
  // tempRank to initial (serial) cells:
  vector<int> mapParallelToSerialCell(simulVolume.subCellList.size());

  // Preparing cells for tempRank

  // Iterate over all cells, and determine whether each has to be transferred to
  // tempRank, along with molecules and complexes which occupy these cells:
  for (unsigned cellItr{0}; cellItr < simulVolume.subCellList.size();
       ++cellItr) {
    int &xBin = simulVolume.subCellList[cellItr].xIndex;
    // If a cell belongs to tempRank that includes ghosted zones:
    if ((xBin >= mpiContext.startGhosted) && (xBin <= mpiContext.endGhosted)) {
      // The cell at old index cellItr will correspond to
      // new index which is equal to current size of cell list for tempRank:
      mapSerialToParallelCell[cellItr] = simulVolumeRank.subCellList.size();
      // When looping over cells prepared to tempRank,
      // finding the corresponding cell in serial version will be needed
      // for finding molecules that need to be transfered:
      mapParallelToSerialCell[simulVolumeRank.subCellList.size()] = cellItr;
      // Add this cell to cell list for tempRank:
      simulVolumeRank.subCellList.push_back(simulVolume.subCellList[cellItr]);
    } else {
      // Cell not to be serialized and sent to the rank:
      mapSerialToParallelCell[cellItr] = -1;
    }
  }

  // Update cell's neighborhood list by mapping from initial cells to cells
  // prepared for tempRank
  for (auto &it_cell : simulVolumeRank.subCellList) {  // for all tempRank cells
    vector<int> newNeighborList;
    for (auto &it_neigh :
         it_cell.neighborList) {  // for all neighborList elements
      // If it gets transfered to tempRank:
      if (mapSerialToParallelCell[it_neigh] != -1) {
        // Push its mapping instead of serial version index:
        newNeighborList.push_back(mapSerialToParallelCell[it_neigh]);
      }
    }
    it_cell.neighborList = newNeighborList;  // apply the new list
  }

  // Set total number of cells that tempRank will hold (including ghosted) to
  // nCell
  simulVolumeRank.numSubCells.tot = simulVolumeRank.subCellList.size();

  // Preparing molecules for tempRank

  // Determine which molecules have to be transfered to tempRank
  for (unsigned cellItrRank{0};
       cellItrRank < simulVolumeRank.subCellList.size();
       ++cellItrRank) {  // for all cells
    // Loop over all molecules that reside within current cell,
    // prepare them for serialization,
    // and mark their complexes for transferring to tempRank:

    // mapParallelToSerialCell is used for finding corresponding cell from
    // serial version list of cells:
    unsigned cellItr = mapParallelToSerialCell[cellItrRank];

    // MemberMolList for the cell for tempRank has to be reconstructed:
    simulVolumeRank.subCellList[cellItrRank].memberMolList.clear();

    // Loop over all molecules from corresponding cell in serial version:
    for (unsigned memItr{0};
         memItr < simulVolume.subCellList[cellItr].memberMolList.size();
         ++memItr) {
      // Find molecule index:
      unsigned targMolIndex =
          simulVolume.subCellList[cellItr].memberMolList[memItr];

      // Find molecule:
      Molecule &mol = moleculeList[targMolIndex];

      // skip implicit lipid
      if (mol.isImplicitLipid == true) {
        continue;
      }

      // Create a backup of the molecule,
      // since mol.index, mol.mySubVolIndex,... have to be restored after
      // changing, as the same molecule may be shared by two ranks:
      Molecule molBackup = mol;

      // Update mySubVolIndex to match tempRank cell index
      mol.mySubVolIndex = cellItrRank;

      // Each molecule keeps an index equal to its position in moleculeList;
      // update this index:
      mol.index = moleculeListRank.size();

      // Update memberMolList by adding the index of mol
      simulVolumeRank.subCellList[cellItrRank].memberMolList.push_back(
          mol.index);

      // Mark complex that this molecule belongs to for transferring to
      // tempRank:
      transferComplex[mol.myComIndex] = true;

      // Serial myComIndex will have to be changed
      // after the complexList for tempRank has been constructed.

      // Store mapping from serial molecule index to parallel molecule index for
      // later use with complexes, etc.:
      mapSerialToParallelMolecule[targMolIndex] = moleculeListRank.size();

      // Prepare molecule for serialization:
      moleculeListRank.push_back(mol);

      // Restore old values of mol, as one mol can be copied to two ranks, where
      // one copy is a ghost:
      mol = molBackup;
    }
  }

  // Each rank stores the number of molecules it can see (includes ghosted
  // molecules) as the size of the molecule list prepared for the rank:
  Molecule::numberOfMolecules = moleculeListRank.size();

  // Preparing and serializing complexes for tempRank

  unsigned startComplexByte = nArrayRank;
  nArrayRank += sizeof(int);  // reserved for number of complexes

  // Complexes for parallel have different indices than for serial complexList.
  // Before sending molecules to a rank, their myComIndex indices have to be
  // updated: Vector mapSerialToParallelComplex stores mapping from serial
  // complex index to parallel complex index:
  vector<int> mapSerialToParallelComplex(complexList.size());

  // Loop through all complexes and check which have to be transfered:
  for (unsigned iComplex = 0; iComplex < complexList.size(); iComplex++) {
    if (transferComplex[iComplex]) {
      // Reduce the size of complex to match molecules from tempRank only;
      // TODO: connect with neighborhood ranks parts of the same complex
      //       using complex management API.

      // Create a complex for tempRank by coping this serial complex:
      Complex c = complexList[iComplex];

      if ((iComplex == 0) && (params.implicitLipid == true)) {
        c.memberList.clear();
        c.memberList.push_back(0);
      } else {
        // Add only those whose members which belong to tempRank
        // (mapSerialToParallelMolecule[ ] != -1)
        c.memberList.clear();
        for (auto &tempMol : complexList[iComplex].memberList)
          if (mapSerialToParallelMolecule[tempMol] != -1) {
            // Add mapping of initial molecule as member of complex iComplex
            c.memberList.push_back(mapSerialToParallelMolecule[tempMol]);
          }
      }

      // Set index to current size of the complexList for tempRank
      c.index = nComplexes;

      // Serialize this complex into arrayRank:
      c.serialize(arrayRank, nArrayRank);

      // Complex at position iComplex will be at position nComplexes at
      // tempRank:
      mapSerialToParallelComplex[iComplex] = nComplexes;

      // Increase number of complexes to be serialized:
      nComplexes++;
    } else {  // don't transfer complex iComplex to tempRank
      mapSerialToParallelComplex[iComplex] = -1;
    }
  }
  // Store the number of serialized complexes to be sent:
  *((int *)(arrayRank + startComplexByte)) = nComplexes;

  // Update myComIndex values in molecules before serializing molecules:
  for (auto &mol : moleculeListRank)
    mol.myComIndex = mapSerialToParallelComplex[mol.myComIndex];

  // Preparing molTemplateList for tempRank

  vector<MolTemplate> molTemplateListRank;
  for (auto &molTemp : molTemplateList) {
    // This should actualy be done only if destruction reaction is allowed for
    // this molecule template
    MolTemplate molTempRank = molTemp;
    molTemplateListRank.push_back(molTempRank);
  }

  vector<int> nBoundPairs;

  for (auto &it : counterArrays.nBoundPairs) it = 0;

  // arrays tracking bound molecule pairs, and species copy nums for tempRank:
  copyCounters counterArraysRank;

  for (bool it : counterArrays.canDissociate) {
    counterArraysRank.canDissociate.push_back(it);
  }

  for (auto &it : counterArrays.proPairlist) {
    counterArraysRank.proPairlist.push_back(it);
  }

  for (auto &it : counterArrays.copyNumSpecies) {
    counterArraysRank.copyNumSpecies.push_back(it);
  }

  for (auto &it : counterArrays.singleDouble) {
    counterArraysRank.singleDouble.push_back(it);
  }

  for (bool it : counterArrays.implicitDouble) {
    counterArraysRank.implicitDouble.push_back(it);
  }

  // Serializing vectors and objects for tempRank

  serialize_abstract_vector<Molecule>(moleculeListRank, arrayRank, nArrayRank);
  simulVolumeRank.serialize(arrayRank, nArrayRank);
  mpiContext.serialize(arrayRank, nArrayRank);
  membraneObject.serialize(arrayRank, nArrayRank);
  serialize_abstract_vector<MolTemplate>(molTemplateListRank, arrayRank,
                                         nArrayRank);
  params.serialize(arrayRank, nArrayRank);
  serialize_abstract_vector<ForwardRxn>(forwardRxns, arrayRank, nArrayRank);
  serialize_abstract_vector<BackRxn>(backRxns, arrayRank, nArrayRank);
  serialize_abstract_vector<CreateDestructRxn>(createDestructRxns, arrayRank,
                                               nArrayRank);
  counterArraysRank.serialize(arrayRank, nArrayRank);

  return nArrayRank;
}
