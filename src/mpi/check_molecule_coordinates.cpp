#include "debug/debug.hpp"
#include "mpi/mpi_function.hpp"

using namespace std;

bool is_ghosted(Molecule &mol, MpiContext &mpiContext,
                SimulVolume &simulVolume) {
  bool ghosted = false;
  int xBin = get_x_bin(mpiContext, mol);
  if (mpiContext.rank) {  // if not first rank
    if (xBin <= 0) {
      ghosted = true;
    }
  }
  if (mpiContext.rank < mpiContext.nprocs - 1) {  // if not last rank
    if (xBin >= simulVolume.numSubCells.x - 1) {
      ghosted = true;
    }
  }
  return ghosted;
}

bool is_owned_by_processor(Molecule &mol, MpiContext &mpiContext,
                           SimulVolume &simulVolume) {
  bool ghosted = false;
  int xBin = get_x_bin(mpiContext, mol);
  if (mpiContext.rank) {  // if not first rank
    if (xBin <= 0) {
      ghosted = true;
    }
  }
  if (mpiContext.rank < mpiContext.nprocs - 1) {  // if not last rank
    if (xBin >= simulVolume.numSubCells.x - 1) {
      ghosted = true;
    }
  }
  return !ghosted;
}

void check_molecule_coordinates(
    MpiContext &mpiContext, vector<Molecule> &moleculeList,
    vector<int> &molIndexList, vector<Complex> &complexList,
    std::vector<MolTemplate> &molTemplateList, SimulVolume &simulVolume,
    Membrane &membraneObject, copyCounters &counterArrays) {
  for (auto &mol_index : molIndexList) {
    auto &mol = moleculeList[mol_index];
    if (mol.myComIndex == -1 || mol.isImplicitLipid == true)
      continue;  // no checking for deleted molecules, implicit lipid
    int xBin = get_x_bin(mpiContext, mol);
    bool shared = false;
    bool ghosted = false;
    if (mpiContext.rank) {  // if not first rank
      if (xBin == 1)
        shared = true;
      else if (xBin <= 0) {
        ghosted = true;
      } else {
        // loop over all molecules within the same complex, if any of them is
        // shared, then this molecule is shared
        auto &complex = complexList[mol.myComIndex];
        for (auto &member : complex.memberList) {
          int xBinMember = get_x_bin(mpiContext, moleculeList[member]);
          if (xBinMember == 1) {
            shared = true;
            break;
          }
        }
      }
    }
    if (mpiContext.rank < mpiContext.nprocs - 1) {  // if not last rank
      if (xBin == simulVolume.numSubCells.x - 2)
        shared = true;
      else if (xBin >= simulVolume.numSubCells.x - 1) {
        ghosted = true;
      } else {
        // loop over all molecules within the same complex, if any of them is
        // shared, then this molecule is shared
        auto &complex = complexList[mol.myComIndex];
        for (auto &member : complex.memberList) {
          int xBinMember = get_x_bin(mpiContext, moleculeList[member]);
          if (xBinMember == simulVolume.numSubCells.x - 2) {
            shared = true;
            break;
          }
        }
      }
    }

  }
}

