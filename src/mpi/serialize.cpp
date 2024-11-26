#include <iostream>

#include "debug/debug.hpp"
#include "error/error.hpp"
#include "macro.hpp"
#include "mpi/mpi_function.hpp"

using namespace std;

void serialize_ids(vector<int>& ids, char *arrayRank, int &nArrayRank){
  if (VERBOSE) cout << "serialize_ids begins" << endl;
  // Store number of ids to be serialized:
  PUSH(ids.size());
  // Loop over set elements:
  for (auto itr = ids.begin(); itr != ids.end(); itr++) {
    // Serialize the id *itr into the arrayRank:
    PUSH(*itr);
  }
  if (VERBOSE) cout << "serialize_ids ends" << endl;
}

void deserialize_ids(vector<int>& ids, char *arrayRank, int &nArrayRank){
  if (VERBOSE) cout << "deserialize_ids begins" << endl;
  // Get number of ids to be deserialized:
  int nIds;
  POP(nIds);
  // Loop over set elements:
  for (int i = 0; i < nIds; i++) {
    // Deserialize the id into the arrayRank:
    int id;
    POP(id);
    ids.push_back(id);
  }
  if (VERBOSE) cout << "deserialize_ids ends" << endl;
}

void serialize_complexes(MpiContext &mpiContext, SimulVolume &simulVolume,
                         Membrane &membraneObject, set<int> &complexesSet,
                         vector<Complex> &complexList,
                         vector<Molecule> &moleculeList,
                         unsigned char *arrayRank, int &nArrayRank) {
  if (VERBOSE) cout << "serialize_complexes begins" << endl;
  // Store number of complexes to be serialized:
  PUSH(complexesSet.size());

  // TODO: AssumING 2 molecule complexes:
  // for (auto &com : complexList) {
  //   com.deleteIfNotReceivedBack = true;
  // }

  // Loop over set elements:
  for (auto itr = complexesSet.begin(); itr != complexesSet.end(); itr++) {
    if (complexList[*itr].isEmpty) continue;
    Complex c = complexList[*itr];
    vector<int> memberList{};
    for (auto &it : c.memberList) {
      memberList.push_back(moleculeList[it].id);
    }
    c.memberList = memberList;

    // Serialize the complex *itr into the arrayRank:
    c.serialize(arrayRank, nArrayRank);

    // complexList[*itr].receivedFromNeighborRank = false;
  }

  if (VERBOSE) cout << "serialize_complexes ends" << endl;
}

void serialize_molecules(MpiContext &mpiContext, SimulVolume &simulVolume,
                         vector<Molecule> &moleculeList,
                         vector<Complex> &complexList, unsigned char *arrayRank,
                         int &nArrayRank, set<int> &complexesSet,
                         bool isLeft) {
  if (VERBOSE) cout << "serialize_molecules begins" << endl;
  // Molecules are serialized after the number of molecules,
  // making it possible to deserialize.
  // Store index in arrayRank where number of molecules will be serialized once
  // known:
  int writeLengthAt = nArrayRank;
  // Start serializing molecules right after the position for the number of
  // molecules:
  nArrayRank += sizeof(int);
  // Count molecules:
  int nMols = 0;

  set<int> moleculesSet{};

  if(isLeft){if (mpiContext.rank > 0) {
    // To left neighbor processor: locate all left ghost and edge molecules and complexes
    for (int yItr{0}; yItr < simulVolume.numSubCells.y; ++yItr) {
      for (int zItr{0}; zItr < simulVolume.numSubCells.z; ++zItr) {
        int leftBoxIndex = yItr * simulVolume.numSubCells.x +
                           zItr * simulVolume.numSubCells.x * simulVolume.numSubCells.y;
        for (auto& molIndex : simulVolume.subCellList[leftBoxIndex].memberMolList) {
          Molecule& mol = moleculeList[molIndex];
          if (mol.isImplicitLipid || mol.isEmpty) continue;
          complexesSet.insert(mol.myComIndex);
          moleculesSet.insert(molIndex);
        }
        leftBoxIndex++;
        for (auto& molIndex : simulVolume.subCellList[leftBoxIndex].memberMolList) {
          Molecule& mol = moleculeList[molIndex];
          if (mol.isImplicitLipid || mol.isEmpty) continue;
          complexesSet.insert(mol.myComIndex);
          moleculesSet.insert(molIndex);
          for (auto& i : complexList[mol.myComIndex].memberList) {
            if (moleculeList[i].isImplicitLipid || moleculeList[i].isEmpty) continue;
            moleculesSet.insert(i);
          }
        }
      }
    }
    // cout << "At simItr " << mpiContext.simItr << endl;
    // cout << "Rank " << mpiContext.rank << " has " << moleculesSet.size()
    //      << " molecules and " << complexesSet.size()
    //      << " complexes to send to the left side." << endl;
  }} else {if (mpiContext.rank < mpiContext.nprocs - 1) {
    for (int yItr{0}; yItr < simulVolume.numSubCells.y; ++yItr) {
      for (int zItr{0}; zItr < simulVolume.numSubCells.z; ++zItr) {
        int rightBoxIndex = yItr * simulVolume.numSubCells.x +
                            zItr * simulVolume.numSubCells.x * simulVolume.numSubCells.y +
                            simulVolume.numSubCells.x - 1;
        for (auto& molIndex : simulVolume.subCellList[rightBoxIndex].memberMolList) {
          Molecule& mol = moleculeList[molIndex];
          if (mol.isImplicitLipid || mol.isEmpty) continue;
          complexesSet.insert(mol.myComIndex);
          moleculesSet.insert(molIndex);
        }
        rightBoxIndex--;
        for (auto& molIndex : simulVolume.subCellList[rightBoxIndex].memberMolList) {
          Molecule& mol = moleculeList[molIndex];
          if (mol.isImplicitLipid || mol.isEmpty) continue;
          complexesSet.insert(mol.myComIndex);
          moleculesSet.insert(molIndex);
          for (auto& i : complexList[mol.myComIndex].memberList) {
            if (moleculeList[i].isImplicitLipid || moleculeList[i].isEmpty) continue;
            moleculesSet.insert(i);
          }
        }
      }
    }
    // cout << "At simItr " << mpiContext.simItr << endl;
    // cout << "Rank " << mpiContext.rank << " has " << moleculesSet.size()
    //      << " molecules and " << complexesSet.size()
    //      << " complexes to send to the right side." << endl;
  }}

  for (auto &molIndex : moleculesSet) {
    Molecule &mol = moleculeList[molIndex];
    int myComIndex = mol.myComIndex;
    indices_to_IDs(mol, moleculeList, complexList); // the index of parent complex is mapped to ID
    mol.serialize(arrayRank, nArrayRank);
    mol.myComIndex = myComIndex; // restore the index of parent complex
    nMols++;
  }

  // Serialize number of molecules before serialized molecules:
  *((int *)&(arrayRank[writeLengthAt])) = nMols;
  if (VERBOSE) cout << "Total serialized number of molecules:" << nMols << endl;
}

void serialize_molecules_after_updating_subboxes(MpiContext &mpiContext, SimulVolume &simulVolume,
                         vector<Molecule> &moleculeList,
                         vector<Complex> &complexList, unsigned char *arrayRank,
                         int &nArrayRank, set<int> &complexesSet,
                         bool isLeft) {
  if (VERBOSE) cout << "serialize_molecules_after_updating_subboxes begins" << endl;
  // Molecules are serialized after the number of molecules,
  // making it possible to deserialize.
  // Store index in arrayRank where number of molecules will be serialized once
  // known:
  int writeLengthAt = nArrayRank;
  // Start serializing molecules right after the position for the number of
  // molecules:
  nArrayRank += sizeof(int);
  // Count molecules:
  int nMols = 0;

  set<int> moleculesSet{};

  if(isLeft){if (mpiContext.rank > 0) {
    for (int yItr{0}; yItr < simulVolume.numSubCells.y; ++yItr) {
      for (int zItr{0}; zItr < simulVolume.numSubCells.z; ++zItr) {
        int leftBoxIndex = yItr * simulVolume.numSubCells.x +
                           zItr * simulVolume.numSubCells.x * simulVolume.numSubCells.y;
        for (auto& molIndex : simulVolume.subCellList[leftBoxIndex].memberMolList) {
          Molecule& mol = moleculeList[molIndex];
          if (mol.isImplicitLipid || mol.isEmpty) continue;
          complexesSet.insert(mol.myComIndex);
          moleculesSet.insert(molIndex);
        }
        leftBoxIndex++;
        for (auto& molIndex : simulVolume.subCellList[leftBoxIndex].memberMolList) {
          Molecule& mol = moleculeList[molIndex];
          if (mol.isImplicitLipid || mol.isEmpty) continue;
          complexesSet.insert(mol.myComIndex);
          moleculesSet.insert(molIndex);
          for (auto& i : complexList[mol.myComIndex].memberList) {
            if (moleculeList[i].isImplicitLipid || moleculeList[i].isEmpty) continue;
            moleculesSet.insert(i);
          }
        }
      }
    }
  }} else {if (mpiContext.rank < mpiContext.nprocs - 1) {
    for (int yItr{0}; yItr < simulVolume.numSubCells.y; ++yItr) {
      for (int zItr{0}; zItr < simulVolume.numSubCells.z; ++zItr) {
        int rightBoxIndex = yItr * simulVolume.numSubCells.x +
                            zItr * simulVolume.numSubCells.x * simulVolume.numSubCells.y +
                            simulVolume.numSubCells.x - 1;
        for (auto& molIndex : simulVolume.subCellList[rightBoxIndex].memberMolList) {
          Molecule& mol = moleculeList[molIndex];
          if (mol.isImplicitLipid || mol.isEmpty) continue;
          complexesSet.insert(mol.myComIndex);
          moleculesSet.insert(molIndex);
        }
        rightBoxIndex--;
        for (auto& molIndex : simulVolume.subCellList[rightBoxIndex].memberMolList) {
          Molecule& mol = moleculeList[molIndex];
          if (mol.isImplicitLipid || mol.isEmpty) continue;
          complexesSet.insert(mol.myComIndex);
          moleculesSet.insert(molIndex);
          for (auto& i : complexList[mol.myComIndex].memberList) {
            if (moleculeList[i].isImplicitLipid || moleculeList[i].isEmpty) continue;
            moleculesSet.insert(i);
          }
        }
      }
    }
  }}

  for (auto &molIndex : moleculesSet) {
    Molecule &mol = moleculeList[molIndex];
    int myComIndex = mol.myComIndex;
    indices_to_IDs(mol, moleculeList, complexList);
    mol.serialize(arrayRank, nArrayRank);
    mol.myComIndex = myComIndex;
    nMols++;
  }

  // Serialize number of molecules before serialized molecules:
  *((int *)&(arrayRank[writeLengthAt])) = nMols;
  if (VERBOSE) cout << "Total serialized number of molecules:" << nMols << endl;
}

void serialize_molecules_partial(MpiContext &mpiContext,
                                 SimulVolume &simulVolume, vector<int> &indices,
                                 vector<Molecule> &moleculeList,
                                 vector<Complex> &complexList,
                                 unsigned char *arrayRank, int &nArrayRank,
                                 set<int> &complexesSet, int cellX) {
  if (VERBOSE) cout << "serialize_molecules begins" << endl;
  // Molecules are serialized after the number of molecules,
  // making it possible to deserialize.
  // Store index in arrayRank where number of molecules will be serialized once
  // known:
  int writeLengthAt = nArrayRank;
  // Start serializing molecules right after the position for the number of
  // molecules:
  nArrayRank += sizeof(int);
  // Count molecules:
  int nMols = 0;
  // Iterate over all cells, and determine whether their molecules have to be
  // transfered
  // TODO: transfer also complexes which occupy these cells:
  for (unsigned cellItr{0}; cellItr < simulVolume.subCellList.size();
       ++cellItr) {  // for all cells
    auto &cell = simulVolume.subCellList[cellItr];
    // cout << "Testing whether cellItr = " << cellItr << " should be
    // serialized" << endl;
    if (cell.xIndex - mpiContext.xOffset == cellX) {
      // cout << "Looping over all molecules from cell: cellItr = " << cellItr
      // << endl;
      for (unsigned memItr{0}; memItr < cell.memberMolList.size(); ++memItr) {
        // Find molecule index:
        unsigned molIndex = cell.memberMolList[memItr];
        if (DEBUG && (molIndex >= moleculeList.size()))
          error("1: molIndex >= moleculeList.size()");
        // Find molecule:
        Molecule &mol = moleculeList[molIndex];
        if (moleculeList[molIndex].myComIndex == -1 ||
            moleculeList[molIndex].isImplicitLipid ==
                true)  // molecule to be deleted or implicit lipid
          continue;
        // cout << "Serializing molecule, index = " << molIndex << endl;
        
        // Mark myComIndex complex to be serialized as well:
        complexesSet.insert(mol.myComIndex);

        // Store current myComIndex,
        // as the one for serialization will be replaced by complexId:
        int myComIndex = mol.myComIndex;

        // Based on interfaces' partnerIndex, assign partnerIDs:
        indices_to_IDs(mol, moleculeList, complexList);

        // Serialize mol into the arrayRank:
        mol.serialize(arrayRank, nArrayRank);

        // Restore original myComIndex, instead of complexId:
        mol.myComIndex = myComIndex;

        // Before sending molecules to the neighbor rank,
        // mark each molecule as not received back from the other rank,
        // so that once molecules from the same zones are received back from
        // other rank, those molecules that were sent, but are not received
        // back, can be deleted:
        mol.receivedFromNeighborRank = false;
        nMols++;
        // TODO: it is a complex => determine whether the complex is broken, and
        // introduce COM rank for the complex.
      }
    }
  }
  // Serialize number of molecules before serialized molecules:
  *((int *)&(arrayRank[writeLengthAt])) = nMols;
  if (VERBOSE) cout << "Total serialized number of molecules:" << nMols << endl;
}
