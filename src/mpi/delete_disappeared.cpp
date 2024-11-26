#include <iostream>

#include "debug/debug.hpp"
#include "error/error.hpp"
#include "macro.hpp"
#include "mpi/mpi_function.hpp"

using namespace std;

void delete_disappeared_complexes(MpiContext &mpiContext,
                                  vector<Molecule> &moleculeList,
                                  vector<Complex> &complexList) {
  if (VERBOSE) cout << "delete_disappeared_complexes begins" << endl;
  int count{0};
  for (auto &com : complexList) {
    if (!com.receivedFromNeighborRank) {
      if (!com.deleteIfNotReceivedBack) {
        // TODO:ASK: (complex management):
        //  update myComIndex of com.memberList molecules to match new molecule.

      } else {
        int zone_index = get_x_bin(mpiContext, moleculeList[com.memberList[0]]);
        if (!com.isEmpty) {
          count++;
          if (VERBOSE)
            cout << "### Deleting complex id=" << com.id
                 << ", complex index=" << com.index << endl
                 << "        zone index=" << zone_index;
          com.memberList.clear();
          com.destroy(moleculeList, complexList);
        }
      }
    }
  }

  if (VERBOSE)
    cout << "delete_disappeared_complexes ends" << endl
         << count << " complexes are deleted" << endl;
}

void delete_disappeared_complexes_partial(MpiContext &mpiContext,
                                          vector<Molecule> &moleculeList,
                                          vector<Complex> &complexList,
                                          bool left) {
  if (VERBOSE) cout << "delete_disappeared_complexes_partial begins" << endl;
  int count{0};
  for (auto &com : complexList) {
    if (com.isEmpty) continue;

    int zone_index;
    int mid_index;

    if (1)
    // if(com.memberList.size() > 0)
    {
      zone_index = get_x_bin(mpiContext, moleculeList[com.memberList[0]]);
      if (mpiContext.rank == 0) {
        mid_index = (mpiContext.simulVolume->numSubCells.x - 1) / 2;
      } else if (mpiContext.rank == mpiContext.nprocs - 1) {
        mid_index = (mpiContext.simulVolume->numSubCells.x + 1) / 2;
      } else {
        mid_index = (mpiContext.simulVolume->numSubCells.x) / 2;
      }

      if (left) {
        if (zone_index >= mid_index) continue;
      } else {
        if (zone_index < mid_index) continue;
      }
    }

    if (!com.receivedFromNeighborRank) {
      if (!com.deleteIfNotReceivedBack) {
        // TODO:ASK: (complex management):
        //  update myComIndex of com.memberList molecules to match new molecule.

      } else {
        if (!com.isEmpty) {
          count++;
          if (VERBOSE)
            cout << "### Deleting complex id=" << com.id
                 << ", complex index=" << com.index << endl
                 << "        zone index=" << zone_index
                 << "        mid index=" << mid_index << endl;
          com.memberList.clear();
          com.destroy(moleculeList, complexList);
        }
      }
    }
  }

  if (VERBOSE)
    cout << "delete_disappeared_complexes_partial ends" << endl
         << count << " complexes are deleted" << endl;
}

void delete_disappeared_molecules(
    MpiContext &mpiContext, SimulVolume &simulVolume, Membrane &membraneObject,
    vector<Molecule> &moleculeList, vector<Complex> &complexList,
    vector<int> &region, vector<MolTemplate> &molTemplateList) {
  if (VERBOSE) cout << "delete_disappeared_molecules begins" << endl;
  int count{0};
  for (auto &molIndex : region) {
    Molecule &mol = moleculeList[molIndex];
    if (mol.myComIndex == -1 || mol.isImplicitLipid == true || mol.isEmpty == true)
      continue;  // no checking for deleted molecules, implicit lipid

    if (mol.isGhosted || mol.isShared) {
      if (mol.receivedFromNeighborRank == false) {
        // Remove the molecule from this rank
        // and remove it from complex as well
        mol.MPI_remove_from_one_rank(moleculeList, complexList);
        count++;
      }
    }
  }
  if (VERBOSE)
    cout << "delete_disappeared_molecules ends" << endl
         << count << " molecules are deleted" << endl;
  // cout << "delete_disappeared_molecules ends" << endl << count << " molecules
  // are deleted" << endl;
}

void disconnect_molecule_partners(unsigned &targMolIndex, Molecule &mol,
                                  vector<Molecule> &moleculeList) {
  if (VERBOSE) cout << "disconnect_molecule_partners begins" << endl;
  // Remove partnerIndex from mol.bndpartner
  for (auto &partnerIndex : mol.bndpartner) {
    if (partnerIndex != -1) {
      Molecule &partner = moleculeList[partnerIndex];

      for (int i = 0; i < partner.bndpartner.size(); i++) {
        if (partner.bndpartner[i] == targMolIndex) {
          partner.bndpartner.erase(partner.bndpartner.begin() + i);
          partner.bndlist.erase(partner.bndlist.begin() + i);
        }
      }

      // Set partnerIndex to -1 in mol.iFace[].interaction.partnerIndex
      for (auto &iFace : partner.interfaceList) {
        if (iFace.interaction.partnerIndex == targMolIndex) {
          iFace.interaction.partnerIndex = -1;
        }
      }
    }
  }
  if (VERBOSE) cout << "disconnect_molecule_partners ends" << endl;
}
