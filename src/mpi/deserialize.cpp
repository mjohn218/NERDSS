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
  if (VERBOSE) cout << "nComplexes = " << nComplexes << endl;
  // Deserialize received complexes:
  for (size_t i = 0; i < nComplexes; i++) {
    if (VERBOSE) cout << "i = " << i << endl;
    // Deserialize i-th received complex:
    Complex c;
    c.deserialize(arrayRank, nArrayRank);

    c.receivedFromNeighborRank = true;
    // cout << "deserialized c.id = " << c.id << endl;

    int complexIndex = find_complex(complexList, c.id);

    // Remap received IDs into moleculeList indices on this rank:
    vector<int> memberList{};
    if (VERBOSE) cout << "members: ";
    for (auto &it : c.memberList) {  // looping over IDs
      if (VERBOSE) cout << it;
      int molIndex = find_molecule(moleculeList, it);
      if (DEBUG && (molIndex == -1))
        error(mpiContext, "5: complex member mol not found");
      memberList.push_back(molIndex);
      if (complexIndex == -1) {
        if (moleculeList[molIndex].justBoundThisStep) {
          complexIndex = moleculeList[molIndex].myComIndex;
        }
      }
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
      // if (complexList[complexIndex].trajStatus == TrajStatus::propagated) {
      //   if (VERBOSE)
      //     printf("Skip complex (id=%d) because it is already propagated.\n",
      //            complexList[complexIndex].id);
      //   complexList[complexIndex].receivedFromNeighborRank = true;
      //   continue;
      // }

      auto vecOld = complexList[complexIndex].memberList;
      if (VERBOSE)
        cout << "Old receivedFromNeighborRank = "
             << complexList[complexIndex].receivedFromNeighborRank << endl;
      complexList[complexIndex] = c;  // copy all deserialized properties
      if (VERBOSE)
        cout << "New receivedFromNeighborRank = "
             << complexList[complexIndex].receivedFromNeighborRank << endl;
      complexList[complexIndex].index = complexIndex;
      auto &vecNew = complexList[complexIndex].memberList;

      // Update member molecules to match also old molecules:
      // for (auto &it : vecOld) {
      //   // cout << "add old..." << endl;
      //   //  If old memberList element isn't found in vecNew, add it:
      //   if (std::find(vecNew.begin(), vecNew.end(), it) == vecNew.end()) {
      //     // Push only if a member molecule isn't deleted:
      //     if ((it < moleculeList.size()) &&
      //         (moleculeList[it].myComIndex != -1)) {
      //       // Push "it" only if not in a shared zone:
      //       int xBin = get_x_bin(mpiContext, moleculeList[it]);
      //       if (((!mpiContext.rank) || (xBin > 1)) &&
      //           ((mpiContext.rank == mpiContext.nprocs - 1) ||
      //            (xBin < (*(mpiContext.simulVolume)).numSubCells.x - 2))) {
      //         vecNew.push_back(it);
      //       }
      //     }
      //   }
      // }
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
      if (complexList[complexIndex].trajStatus == TrajStatus::propagated) {
        if (VERBOSE)
          printf("Skip complex (id=%d) because it is already propagated.\n",
                 complexList[complexIndex].id);
        complexList[complexIndex].receivedFromNeighborRank = true;
        continue;
      }

      auto vecOld = complexList[complexIndex].memberList;
      if (VERBOSE)
        cout << "Old receivedFromNeighborRank = "
             << complexList[complexIndex].receivedFromNeighborRank << endl;
      complexList[complexIndex] = c;  // copy all deserialized properties
      if (VERBOSE)
        cout << "New receivedFromNeighborRank = "
             << complexList[complexIndex].receivedFromNeighborRank << endl;
      complexList[complexIndex].index = complexIndex;
      auto &vecNew = complexList[complexIndex].memberList;

      // Update member molecules to match also old molecules:
      // for (auto &it : vecOld) {
      //   // cout << "add old..." << endl;
      //   //  If old memberList element isn't found in vecNew, add it:
      //   if (std::find(vecNew.begin(), vecNew.end(), it) == vecNew.end()) {
      //     // Push only if a member molecule isn't deleted:
      //     if ((it < moleculeList.size()) &&
      //         (moleculeList[it].myComIndex != -1)) {
      //       // Push "it" only if not in a shared zone:
      //       int xBin = get_x_bin(mpiContext, moleculeList[it]);
      //       if (((!mpiContext.rank) || (xBin > 1)) &&
      //           ((mpiContext.rank == mpiContext.nprocs - 1) ||
      //            (xBin < (*(mpiContext.simulVolume)).numSubCells.x - 2))) {
      //         vecNew.push_back(it);
      //       }
      //     }
      //   }
      // }
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
  if (VERBOSE)
    printf("deserializing (%d) molecules from the left rank...\n", nMols);
  // Deserialize received molecules and update structures:
  for (int iDeserialize = 0; iDeserialize < nMols; iDeserialize++) {
    // Deserialize each received molecule:
    Molecule mol;

    mol.deserialize(arrayRank, nArrayRank);

    if (VERBOSE) printf("deserializing mol (id=%d) ...\n", mol.id);

    update_mySubVolIndex(mpiContext, mol, membraneObject, simulVolume);

    mol.receivedFromNeighborRank = true;

    mol.isGhosted = false;
    if (mpiContext.rank) {  // if not first rank
      int xBin = get_x_bin(mpiContext, mol);
      if (xBin <= 0) {
        mol.isGhosted = true;
      }
    }

    int molIndex = find_molecule(moleculeList, mol.id);

    // Check whether the molecule was already in ghosted zone,
    // or it has just arrived (index == -1):
    if (molIndex != -1) {  // if found
      if (VERBOSE) printf("mol (id=%d) is already here.\n", mol.id);

      if (DEBUG && (moleculeList[molIndex].id != mol.id))
        error(mpiContext, "7 // deserialize_molecules()");
      mol.index = molIndex;

      // check whether disassociated TODO: there might be an issue when
      // dissociation included, so we need to represent the dissociation status
      // with a different TrajStatus value
      // if (moleculeList[molIndex].trajStatus == TrajStatus::propagated) {
      //   // std::cout << "skip deserialize molecule id
      //   // "<<moleculeList[molIndex].id <<std::endl;
      //   if (VERBOSE)
      //     printf(
      //         "mol (id=%d) has been propagated, so skip the
      //         deserialization.\n", mol.id);
      //   moleculeList[molIndex].receivedFromNeighborRank = true;
      //   continue;
      // }

      // Update mySubVolIndex if needed
      // (received molecule might have moved by another rank
      // since this rank serialized it)
      auto &oldCellIndex = moleculeList[molIndex].mySubVolIndex;

      // check if the mol cross the ghosted boundary
      // int xBinOld = get_x_bin(mpiContext, moleculeList[molIndex]);
      // int xBinNew = get_x_bin(mpiContext, mol);
      // if (check_enter_ghosted_zone(mpiContext, xBinOld, xBinNew,
      //                              indicesEnteringMolecules, simulVolume)) {
      //   if (VERBOSE)
      //     printf(
      //         "mol (id=%d) enter the ghosted zone; old xBin (%d), new xBin "
      //         "(%d).\n",
      //         mol.id, xBinOld, xBinNew);
      //   indicesEnteringMolecules.push_back(molIndex);
      // }
      // if (check_exit_ghosted_zone(mpiContext, xBinOld, xBinNew,
      //                             indicesEnteringMolecules, simulVolume)) {
      //   if (VERBOSE)
      //     printf(
      //         "mol (id=%d) exit the ghosted zone; old xBin (%d), new xBin "
      //         "(%d).\n",
      //         mol.id, xBinOld, xBinNew);
      //   indicesExitingMolecules.push_back(molIndex);
      // }

      moleculeList[molIndex] = mol;

      auto &newCellIndex = moleculeList[molIndex].mySubVolIndex;
      // Check if cell index has changed:
      if (newCellIndex != oldCellIndex) {
        // Update simulVolume cells memberMolList of both mySubVolIndex:
        if (oldCellIndex < simulVolume.subCellList.size()) {
          auto &vec = simulVolume.subCellList[oldCellIndex].memberMolList;
          vec.erase(remove(vec.begin(), vec.end(), molIndex), vec.end());
        }
        // Deserialized molecule will be at the same index as existing one:
        if (newCellIndex < simulVolume.subCellList.size()) {
          simulVolume.subCellList[newCellIndex].memberMolList.push_back(
              molIndex);
        }
      }
    } else {  // if not found
      if (VERBOSE) printf("mol (id=%d) is new here.\n", mol.id);
      //  The molecule doesn't exist in moleculeList. Add it:
      mol.index = moleculeList.size();
      molIndex = mol.index;
      moleculeList.push_back(mol);
    }
    indices.push_back(molIndex);
    // if(DEBUG) DEBUG_MOL("1DESERIALIZed A MOLECULE");
  }

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

    update_mySubVolIndex(mpiContext, mol, membraneObject, simulVolume);

    mol.receivedFromNeighborRank = true;

    mol.isGhosted = false;
    // If the deserialized stripe is ghost stripe, mark molecule as ghosted:
    if (mpiContext.rank < mpiContext.nprocs - 1) {  // if not last rank
      int xBin = get_x_bin(mpiContext, mol);
      if (xBin >= simulVolume.numSubCells.x - 1) {
        mol.isGhosted = true;
      }
    }

    int molIndex = find_molecule(moleculeList, mol.id);

    // Check whether the molecule was already in ghosted zone,
    // or it has just arrived (index == -1):
    if (molIndex != -1) {  // if found
      if (VERBOSE) printf("mol (id=%d) is already here.\n", mol.id);

      if (DEBUG && (moleculeList[molIndex].id != mol.id))
        error(mpiContext, "7 // deserialize_molecules()");
      mol.index = molIndex;

      // check whether disassociated TODO: there might be an issue when
      // dissociation included, so we need to represent the dissociation status
      // with a different TrajStatus value
      if (moleculeList[molIndex].trajStatus == TrajStatus::propagated) {
        // std::cout << "skip deserialize molecule id
        // "<<moleculeList[molIndex].id <<std::endl;
        if (VERBOSE)
          printf(
              "mol (id=%d) has been propagated, so skip the deserialization.\n",
              mol.id);
        moleculeList[molIndex].receivedFromNeighborRank = true;
        continue;
      }

      // Update mySubVolIndex if needed
      // (received molecule might have moved by another rank
      // since this rank serialized it)
      auto oldCellIndex = moleculeList[molIndex].mySubVolIndex;

      // check if the mol cross the ghosted boundary
      // int xBinOld = get_x_bin(mpiContext, moleculeList[molIndex]);
      // int xBinNew = get_x_bin(mpiContext, mol);
      // if (check_enter_ghosted_zone(mpiContext, xBinOld, xBinNew,
      //                              indicesEnteringMolecules, simulVolume)) {
      //   if (VERBOSE)
      //     printf(
      //         "mol (id=%d) enter the ghosted zone; old xBin (%d), new xBin "
      //         "(%d).\n",
      //         mol.id, xBinOld, xBinNew);
      //   indicesEnteringMolecules.push_back(molIndex);
      // }
      // if (check_exit_ghosted_zone(mpiContext, xBinOld, xBinNew,
      //                             indicesEnteringMolecules, simulVolume)) {
      //   if (VERBOSE)
      //     printf(
      //         "mol (id=%d) exit the ghosted zone; old xBin (%d), new xBin "
      //         "(%d).\n",
      //         mol.id, xBinOld, xBinNew);
      //   indicesExitingMolecules.push_back(molIndex);
      // }

      Molecule oldMol = moleculeList[molIndex];
      moleculeList[molIndex] = mol;

      auto newCellIndex = moleculeList[molIndex].mySubVolIndex;
      // Check if cell index has changed:
      if (newCellIndex != oldCellIndex) {
        // Update simulVolume cells memberMolList of both mySubVolIndex:
        if (oldCellIndex < simulVolume.subCellList.size()) {
          auto &vec = simulVolume.subCellList[oldCellIndex].memberMolList;
          vec.erase(remove(vec.begin(), vec.end(), molIndex), vec.end());
        }
        // Deserialized molecule will be at the same index as existing one:
        if (newCellIndex < simulVolume.subCellList.size()) {
          simulVolume.subCellList[newCellIndex].memberMolList.push_back(
              molIndex);
        }
      }
    } else {  // if not found
      if (VERBOSE) printf("mol (id=%d) is new here.\n", mol.id);
      // The molecule doesn't exist in moleculeList. Add it:
      mol.index = moleculeList.size();
      molIndex = mol.index;
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

  int currBin = xItr + (yItr * simulVolume.numSubCells.x) +
                (zItr * simulVolume.numSubCells.x * simulVolume.numSubCells.y);
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
