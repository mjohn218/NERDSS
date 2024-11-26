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

    mol.isGhosted = ghosted;
    mol.isShared = shared;
  }
}

// TODO: update membraneObject.numberOfFreeLipidsEachState[0] when a implicit
// lipid-bound mol move cross ranks
void exit_ghosted_zone(MpiContext &mpiContext, Molecule &mol,
                       std::vector<MolTemplate> &molTemplateList,
                       copyCounters &counterArrays,
                       std::vector<Molecule> &moleculeList,
                       SimulVolume &simulVolume) {
  mol.isGhosted = false;

  update_copyCounters_exit_ghosted_zone(mpiContext, mol, molTemplateList,
                                        counterArrays, moleculeList,
                                        simulVolume);
}

void enter_ghosted_zone(MpiContext &mpiContext, Molecule &mol,
                        std::vector<MolTemplate> &molTemplateList,
                        copyCounters &counterArrays,
                        std::vector<Molecule> &moleculeList,
                        SimulVolume &simulVolume) {
  mol.isGhosted = true;

  update_copyCounters_enter_ghosted_zone(mpiContext, mol, molTemplateList,
                                         counterArrays, moleculeList,
                                         simulVolume);
}

void update_copyCounters_enter_ghosted_zone(
    MpiContext &mpiContext, Molecule &mol, vector<MolTemplate> &molTemplateList,
    copyCounters &counterArrays, vector<Molecule> &moleculeList,
    SimulVolume &simulVolume) {
  MolTemplate &oneTemp{molTemplateList[mol.molTypeIndex]};

  if (mol.bndpartner.empty()) {
    Molecule::numberOfMolecules--;
  } else {
    bool updateNumberOfComplexes{true};
    for (auto iface : mol.interfaceList) {
      if (iface.isBound) {
        int partnerIndex = iface.interaction.partnerIndex;
        bool partnerIsGhosted{false};
        partnerIsGhosted =
            is_ghosted(moleculeList[partnerIndex], mpiContext, simulVolume);

        // Update numberOfComplexes
        if (!moleculeList[partnerIndex].isGhosted) {
          updateNumberOfComplexes = false;
        }
      }
    }
    if (updateNumberOfComplexes) --Complex::numberOfComplexes;
  }
}

void update_copyCounters_exit_ghosted_zone(MpiContext &mpiContext,
                                           Molecule &mol,
                                           vector<MolTemplate> &molTemplateList,
                                           copyCounters &counterArrays,
                                           vector<Molecule> &moleculeList,
                                           SimulVolume &simulVolume) {
  MolTemplate &oneTemp{molTemplateList[mol.molTypeIndex]};

  if (mol.bndpartner.empty()) {
    Molecule::numberOfMolecules++;
  } else {
    bool updateNumberOfComplexes{true};
    for (auto &iface : mol.interfaceList) {
      if (iface.isBound) {
        int partnerIndex = iface.interaction.partnerIndex;
        bool partnerIsGhosted{false};
        partnerIsGhosted =
            is_ghosted(moleculeList[partnerIndex], mpiContext, simulVolume);

        // Update numberOfComplexes
        if (!moleculeList[partnerIndex].isGhosted) {
          updateNumberOfComplexes = false;
        }
      }
    }
    if (updateNumberOfComplexes) ++Complex::numberOfComplexes;
  }
}

void move_complexes_based_on_propagation(MpiContext &mpiContext,
                                         vector<Molecule> &moleculeList,
                                         SimulVolume &simulVolume,
                                         Membrane &membraneObject) {}
