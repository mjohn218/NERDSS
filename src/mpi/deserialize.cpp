#include <cstdio>
#include <iostream>

#include "debug/debug.hpp"
#include "error/error.hpp"
#include "macro.hpp"
#include "mpi/mpi_function.hpp"

using namespace std;

void deserialize_complexes(MpiContext &mpiContext,
                           vector<Molecule> &moleculeList,
                           vector<Complex> &complexList,
                           unsigned char *arrayRank, int &nArrayRank) {
  if (VERBOSE) cout << "deserialize_complexes begins" << endl;
  // Store number of complexes to be deserialized:
  size_t nComplexes;
  POP(nComplexes);
  int oldCom = 0;
  int newCom = 0;
  if (VERBOSE) cout << "nComplexes = " << nComplexes << endl;
  // Deserialize received complexes:
  for (size_t i = 0; i < nComplexes; i++) {
    if (VERBOSE) cout << "i = " << i << endl;
    // Deserialize i-th received complex:
    Complex c;
    c.deserialize(arrayRank, nArrayRank);

    c.receivedFromNeighborRank = true;
    int complexIndex = find_complex(complexList, c.id);

    // Remap received IDs into moleculeList indices on this rank:
    vector<int> memberList{};
    if (VERBOSE) cout << "members: ";
    for (auto &it : c.memberList) {  // looping over IDs
      if (VERBOSE) cout << it;
      int molIndex = find_molecule(moleculeList, it);
      if (molIndex == -1)
        error(mpiContext, "5: complex member mol not found");
      memberList.push_back(molIndex);
      if (VERBOSE) cout << "(" << molIndex << "); ";
    }
    if (VERBOSE) cout << endl;
    c.memberList = memberList;

    // Extract unique complex identifier in the system,
    // and find index to local complex with the same ID:

    if (VERBOSE)
      cout << "complexIndex = " << complexIndex << "complexId = " << c.id
           << endl;

    if (complexIndex == -1) {  // new complex on this rank
      newCom++;
      if (VERBOSE)
        printf("This is a new complex (id=%d) here\n",
               complexList[complexIndex].id);

      c.index = complexList.size();
      complexIndex = c.index;
      complexList.push_back(c);
    } else {
      oldCom++;
      if (VERBOSE)
        printf("This is a old complex (id=%d) here\n",
               complexList[complexIndex].id);

      auto vecOld = complexList[complexIndex].memberList;
      c.isLeftEdge = complexList[complexIndex].isLeftEdge;
      c.isRightEdge = complexList[complexIndex].isRightEdge;
      c.isLeftGhost = complexList[complexIndex].isLeftGhost;
      c.isRightGhost = complexList[complexIndex].isRightGhost;

      complexList[complexIndex] = c;  // copy all deserialized properties
      complexList[complexIndex].index = complexIndex;
      auto &vecNew = complexList[complexIndex].memberList;
    }

    // Update myComIndex of member molecules to complexIndex:
    complexList[complexIndex].index = complexIndex;
    for (auto &it : complexList[complexIndex].memberList) {
      // cout << "updating member " << it << endl;
      if (DEBUG && ((it >= moleculeList.size()) || (it == -1)))
        error(mpiContext,
              "6: complex member mol not found: it=" + to_string(it));
      Molecule &mol = moleculeList[it];
      mol.myComIndex = complexIndex;
    }
  }
  // cout << "At simItr " << mpiContext.simItr << endl;
  // cout << "Received oldCom = " << oldCom << ", newCom = " << newCom << " total = " << oldCom+newCom << endl;

  if (VERBOSE) cout << "deserialize_complexes ends" << endl;
}

void deserialize_complexes_right(MpiContext &mpiContext,
                                 vector<Molecule> &moleculeList,
                                 vector<Complex> &complexList,
                                 unsigned char *arrayRank, int &nArrayRank) {
  if (VERBOSE) cout << "deserialize_complexes begins..." << endl;
  // Store number of complexes to be deserialized:
  size_t nComplexes;
  POP(nComplexes);
  if (VERBOSE) cout << "nComplexes = " << nComplexes << endl;
  // Deserialize received complexes:
  for (size_t i = 0; i < nComplexes; i++) {
    if (VERBOSE) cout << "i = " << i << endl;
    // Deserialize i-th received complex:
    Complex c;
    c.deserialize(arrayRank, nArrayRank);

    c.receivedFromNeighborRank = true;
    int complexIndex = find_complex(complexList, c.id);

    // Remap received IDs into moleculeList indices on this rank:
    vector<int> memberList{};
    if (VERBOSE) cout << "members: ";
    for (auto &it : c.memberList) {  // looping over IDs
      if (VERBOSE) cout << it;
      int molIndex = find_molecule(moleculeList, it);
      memberList.push_back(molIndex);
      if (VERBOSE) cout << "(" << molIndex << "); ";
    }
    if (VERBOSE) cout << endl;
    c.memberList = memberList;

    if (VERBOSE)
      cout << "complexIndex = " << complexIndex << "complexId = " << c.id
           << endl;

    if (complexIndex == -1) {  // new complex on this rank
      if (VERBOSE)
        printf("This is a new complex (id=%d) here\n",
               complexList[complexIndex].id);

      c.index = complexList.size();
      complexIndex = c.index;
      complexList.push_back(c);
    } else {
      if (VERBOSE)
        printf("This is a old complex (id=%d) here\n",
               complexList[complexIndex].id);

      auto vecOld = complexList[complexIndex].memberList;
      c.isLeftEdge = complexList[complexIndex].isLeftEdge;
      c.isRightEdge = complexList[complexIndex].isRightEdge;
      c.isLeftGhost = complexList[complexIndex].isLeftGhost;
      c.isRightGhost = complexList[complexIndex].isRightGhost;

      complexList[complexIndex] = c;  // copy all deserialized properties
      complexList[complexIndex].index = complexIndex;
      auto &vecNew = complexList[complexIndex].memberList;
    }

    // Update myComIndex of member molecules to complexIndex:
    complexList[complexIndex].index = complexIndex;
    for (auto &it : complexList[complexIndex].memberList) {
      // cout << "updating member " << it << endl;
      if (DEBUG && ((it >= moleculeList.size()) || (it == -1)))
        error(mpiContext,
              "6: complex member mol not found: it=" + to_string(it));
      Molecule &mol = moleculeList[it];
      mol.myComIndex = complexIndex;
    }
  }
  //    debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
  //    "//end deserialize_complexes()");
  if (VERBOSE) cout << "deserialize_complexes ends" << endl;
}

// Deserializes stripe of molecules with cellX x-bin coordinate,
// and updates moleculeList.
// Example stripes x-bin coordinates:
// rank0 |rank1|rank2| rank3  <- ownership
// 0 1 2 | 3                  <- rank0 visibility (includes ghosted 3)
//     2 | 3 4 | 5            <- rank1 visibility (includes ghosted 2 and 5)
//           4 | 5 6 | 7      <- rank2 visibility (includes ghosted 4 and 7)
//                 6 | 7 8 9  <- rank3 visibility (includes ghosted 6)
void deserialize_molecules(MpiContext &mpiContext, SimulVolume &simulVolume,
                           vector<Molecule> &moleculeList,
                           vector<Complex> &complexList,
                           vector<MolTemplate> &molTemplateList,
                           Membrane &membraneObject,
                           copyCounters &counterArrays, vector<int> &indices,
                           vector<int> &indicesEnteringMolecules,
                           vector<int> &indicesExitingMolecules,
                           unsigned char *arrayRank, int &nArrayRank) {
  int nMols;   // number of molecules to deserialize
  POP(nMols);  // deserialize from arrayRank
  int oldMol = 0;
  int newMol = 0;
  if (VERBOSE)
    printf("deserializing (%d) molecules from the right rank...\n", nMols);
  // Deserialize received molecules and update structures:
  for (int iDeserialize = 0; iDeserialize < nMols; iDeserialize++) {
    // Deserialize each received molecule:
    Molecule mol;

    mol.deserialize(arrayRank, nArrayRank);

    if (VERBOSE) printf("deserializing mol (id=%d) ...\n", mol.id);

    int molIndex = find_molecule(moleculeList, mol.id);

    // Check whether the molecule was already in this processor,
    // or it has just arrived (index == -1):
    if (molIndex != -1) {  // if found
      oldMol++;
      if (VERBOSE) printf("mol (id=%d) is already here.\n", mol.id);

      if (DEBUG && (moleculeList[molIndex].id != mol.id))
        error(mpiContext, "7 // deserialize_molecules()");

      // update the index of the molecule to match the current processor
      mol.index = molIndex;

      Molecule oldMol = moleculeList[molIndex];
      // keep the region flags
      mol.isLeftEdge = oldMol.isLeftEdge;
      mol.isRightEdge = oldMol.isRightEdge;
      mol.isLeftGhost = oldMol.isLeftGhost;
      mol.isRightGhost = oldMol.isRightGhost;
      mol.receivedFromNeighborRank = true;
      mol.mySubVolIndex = oldMol.mySubVolIndex;

      moleculeList[molIndex] = mol;
    } else {  // if not found, add this new molecule
      newMol++;
      mol.index = moleculeList.size();
      molIndex = mol.index;
      mol.receivedFromNeighborRank = true;
      moleculeList.push_back(mol);
    }
    indices.push_back(molIndex);
    // if(DEBUG) DEBUG_MOL("1DESERIALIZed A MOLECULE");
  }
  // cout << "At simItr " << mpiContext.simItr << endl;
  // cout << "Received " << oldMol << " old molecules and " << newMol
  //      << " new molecules." << " Total: " << oldMol+newMol << endl;

  if (VERBOSE)
    cout << "Total deserialized number of molecules: " << nMols << endl;
}

void deserialize_molecules_right(
    MpiContext &mpiContext, SimulVolume &simulVolume,
    vector<Molecule> &moleculeList, vector<Complex> &complexList,
    vector<MolTemplate> &molTemplateList, Membrane &membraneObject,
    copyCounters &counterArrays, vector<int> &indices,
    vector<int> &indicesEnteringMolecules, vector<int> &indicesExitingMolecules,
    unsigned char *arrayRank, int &nArrayRank) {
  int nMols;   // number of molecules to deserialize
  POP(nMols);  // deserialize from arrayRank
  if (VERBOSE)
    printf("deserializing (%d) molecules from the right rank...\n", nMols);
  // Deserialize received molecules and update structures:
  for (int iDeserialize = 0; iDeserialize < nMols; iDeserialize++) {
    // Deserialize each received molecule:
    Molecule mol;

    mol.deserialize(arrayRank, nArrayRank);

    if (VERBOSE) printf("deserializing mol (id=%d) ...\n", mol.id);

    int molIndex = find_molecule(moleculeList, mol.id);

    // Check whether the molecule was already in this processor,
    // or it has just arrived (index == -1):
    if (molIndex != -1) {  // if found
      if (VERBOSE) printf("mol (id=%d) is already here.\n", mol.id);

      if (DEBUG && (moleculeList[molIndex].id != mol.id))
        error(mpiContext, "7 // deserialize_molecules()");

      // update the index of the molecule to match the current processor
      mol.index = molIndex;

      Molecule oldMol = moleculeList[molIndex];
      // keep the region flags
      mol.isLeftEdge = oldMol.isLeftEdge;
      mol.isRightEdge = oldMol.isRightEdge;
      mol.isLeftGhost = oldMol.isLeftGhost;
      mol.isRightGhost = oldMol.isRightGhost;
      mol.receivedFromNeighborRank = true;
      mol.mySubVolIndex = oldMol.mySubVolIndex;

      moleculeList[molIndex] = mol;
    } else {  // if not found, add this new molecule
      mol.index = moleculeList.size();
      molIndex = mol.index;
      mol.receivedFromNeighborRank = true;
      moleculeList.push_back(mol);
    }
    indices.push_back(molIndex);
    // if(DEBUG) DEBUG_MOL("1DESERIALIZed A MOLECULE");
  }

  if (VERBOSE)
    cout << "Total deserialized number of molecules: " << nMols << endl;
}

void update_mySubVolIndex(MpiContext &mpiContext, Molecule &mol,
                          Membrane &membraneObject, SimulVolume &simulVolume) {
  if (VERBOSE) cout << "update_mySubVolIndex begins" << endl;
  // The molecule that a rank sends might return with new coordinates/cell
  // index; therefore, searching should be done at least around current cell:
  // ...for(auto &molIndex :
  // simulVolume.subCellList[mol.mySubVolIndex].memberList){... Another possible
  // optimization: receive action requrement to move the molecule to certain
  // cell Another possible optimization: for each molecule, keep left index and
  // right index for finding it in the moleculeList of neighbor rank
  // TODO:ANSW: (based on answer on another question)
  //         int currBin = mol.mySubVolIndex - mpiContext.xOffset;
  int xItr{int((mol.comCoord.x + membraneObject.waterBox.x / 2) /
               simulVolume.subCellSize.x) -
           mpiContext.xOffset};
  int yItr{int((mol.comCoord.y + membraneObject.waterBox.y / 2) /
               simulVolume.subCellSize.y)};
  int zItr;
  if (membraneObject.waterBox.z > 0)
    zItr = int(-(mol.comCoord.z + 1E-6 - membraneObject.waterBox.z / 2.0) /
               simulVolume.subCellSize.z);
  else
    zItr = 0;

  int xItrAdjusted, yItrAdjusted, zItrAdjusted;

  if (xItr < 0) xItrAdjusted = 0;
  else if (xItr > simulVolume.numSubCells.x -1) xItrAdjusted = simulVolume.numSubCells.x - 1;
  else xItrAdjusted = xItr;
  if (yItr < 0) yItrAdjusted = 0;
  else if (yItr > simulVolume.numSubCells.y -1) yItrAdjusted = simulVolume.numSubCells.y - 1;
  else yItrAdjusted = yItr;
  if (zItr < 0) zItrAdjusted = 0;
  else if (zItr > simulVolume.numSubCells.z -1) zItrAdjusted = simulVolume.numSubCells.z - 1;
  else zItrAdjusted = zItr;

  int currBin = xItrAdjusted + (yItrAdjusted * simulVolume.numSubCells.x) + (zItrAdjusted * simulVolume.numSubCells.x * simulVolume.numSubCells.y);

  // Based on coordinates, assign index of cell:
  mol.mySubVolIndex = currBin;
  if (VERBOSE) cout << "update_mySubVolIndex ends" << endl;
}

bool check_enter_ghosted_zone(MpiContext &mpiContext, int xBinOld, int xBinNew,
                              vector<int> &indicesEnteringMolecules,
                              SimulVolume &simulVolume) {
  if (mpiContext.rank && xBinOld == 1 && xBinNew == 0) {  // if not first rank
    return true;
  }
  if (mpiContext.rank < mpiContext.nprocs - 1 &&
      xBinOld == simulVolume.numSubCells.x - 2 &&
      xBinNew == simulVolume.numSubCells.x - 1) {  // if not last rank
    return true;
  }
  return false;
}

bool check_exit_ghosted_zone(MpiContext &mpiContext, int xBinOld, int xBinNew,
                             vector<int> &indicesEnteringMolecules,
                             SimulVolume &simulVolume) {
  if (mpiContext.rank && xBinOld == 0 && xBinNew == 1) {  // if not first rank
    return true;
  }
  if (mpiContext.rank < mpiContext.nprocs - 1 &&
      xBinOld == simulVolume.numSubCells.x - 1 &&
      xBinNew == simulVolume.numSubCells.x - 2) {  // if not last rank
    return true;
  }
  return false;
}
