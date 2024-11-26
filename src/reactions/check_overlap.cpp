#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>

#include "boundary_conditions/reflect_functions.hpp"
#include "debug/debug.hpp"
#include "error/error.hpp"
#include "io/io.hpp"
#include "macro.hpp"
#include "math/constants.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "mpi.h"
#include "mpi/mpi_function.hpp"
#include "parser/parser_functions.hpp"
#include "reactions/association/association.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "split.cpp"
#include "system_setup/system_setup.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

void check_overlap(std::vector<int>& region, bool isLeft, unsigned simItr,
                   Parameters& params, std::vector<Molecule>& moleculeList,
                   std::vector<Complex>& complexList, SimulVolume& simulVolume,
                   std::vector<ForwardRxn>& forwardRxns,
                   std::vector<BackRxn>& backRxns,
                   std::vector<CreateDestructRxn>& createDestructRxns,
                   std::vector<MolTemplate>& molTemplateList,
                   std::map<std::string, int>& observablesList,
                   copyCounters& counterArrays, Membrane& membraneObject,
                   MpiContext& mpiContext) {
  for (auto& mol_index : region) {
    auto& mol{moleculeList[mol_index]};
    // Now track each complex (ncrosscom), and test for overlap of all
    // proteins in that complex before performing final position updates.
    if (mol.isEmpty || mol.isImplicitLipid ||mol.trajStatus == TrajStatus::propagated)
      continue;

    // if it is a left ghost mol, it can not diffuse because it can still possibly react with others in the left neighbor processor
    // if it is a left edge mol and its parent complex is left ghost, it can not diffuse because it's members in the same complex 
    // can still possibly react with others in the left neighbor processor
    // if it is a right ghost mol, it can not diffuse
    if (mol.isLeftGhost) {
      continue;
    }
    if (mol.isLeftEdge && complexList[mol.myComIndex].isLeftGhost) {
      continue;
    }
    if (mol.isRightGhost) {
      continue;
    }

    if (isLeft) {
      // if any part of this complex is at right half, it can not diffuse
      bool isRightHalf = false;
      for (auto& memMol : complexList[mol.myComIndex].memberList) {
        // check if memMol is in the vector region
        if (moleculeList[memMol].isRightHalf) {
          isRightHalf = true;
          break;
        }
      }
      if (isRightHalf) {
        continue;
      }
    }

    //  determine RS3Dinput
    double RS3Dinput{0.0};
    for (int RS3Dindex = 0; RS3Dindex < 100; RS3Dindex++) {
      if (std::abs(membraneObject.RS3Dvect[RS3Dindex + 400] -
                   mol.molTypeIndex) < 1E-2) {
        RS3Dinput = membraneObject.RS3Dvect[RS3Dindex + 300];
        break;
      }
    }

    if (complexList[mol.myComIndex].ncross > 0) {
      // This molecule can react, but it didnt
      if (mol.trajStatus == TrajStatus::none || mol.trajStatus == TrajStatus::canBeResampled) {
        // For each molecule that overlapped and did not react, check whether
        // it overlaps with its partners, do all molecules in the same complex
        // at the same time. Also, if both molecules are stuck to membrane,
        // only do xy displacement, ignore z
        // TODO: Maybe do a boundary sphere overlap check first?

        if (std::abs(complexList[mol.myComIndex].D.z) < 1E-10) {
          if (params.clusterOverlapCheck == false) {
            sweep_separation_complex_rot_memtest(
                simItr, mol.index, params, moleculeList, complexList,
                forwardRxns, molTemplateList, membraneObject);
          } else {
            sweep_separation_complex_rot_memtest_cluster(
                simItr, mol.index, params, moleculeList, complexList,
                forwardRxns, molTemplateList, membraneObject);
          }
        } else {
          // TODO:ASK: Delete the following comment?

          // if(mpiContext.simItr >= 1177)
          //                 DEBUG_FIND_MOL("6.122");
          sweep_separation_complex_rot(simItr, mol.index, params, moleculeList,
                                       complexList, forwardRxns,
                                       molTemplateList, membraneObject);
        }

        // TODO:ASK: Delete the following comment?

        // if(mpiContext.simItr >= 1177)
        //                 DEBUG_FIND_MOL("6.123");
        if (membraneObject.isSphere == true)
          reflect_complex_rad_rot(membraneObject, complexList[mol.myComIndex],
                                  moleculeList, RS3Dinput);
      }
    } else {
      if (mol.trajStatus == TrajStatus::none || mol.trajStatus == TrajStatus::canBeResampled) {
        // For proteins with ncross=0, they either moved independently, or
        // their displacements were selected based on the complex they were
        // part of, and they may not yet been moved.
        if (membraneObject.isSphere == true) {
          if (mol.trajStatus == TrajStatus::none) {
            create_complex_propagation_vectors(
                params, complexList[mol.myComIndex], moleculeList, complexList,
                molTemplateList, membraneObject);
            for (auto& memMol : complexList[mol.myComIndex].memberList)
              moleculeList[memMol].trajStatus = TrajStatus::canBeResampled;
          }
          complexList[mol.myComIndex].propagate(moleculeList, membraneObject,
                                                molTemplateList);
          reflect_complex_rad_rot(membraneObject, complexList[mol.myComIndex],
                                  moleculeList, RS3Dinput);
        } else {
          // reflect_traj_complex_rad_rot(params, moleculeList,
          // complexList[mol.myComIndex], membraneObject, RS3Dinput);
          if (mol.trajStatus == TrajStatus::none) {
            create_complex_propagation_vectors(
                params, complexList[mol.myComIndex], moleculeList, complexList,
                molTemplateList, membraneObject);
            for (auto& memMol : complexList[mol.myComIndex].memberList) {
              moleculeList[memMol].trajStatus = TrajStatus::canBeResampled;
            }
          }

          // TODO:ASK: TEMP. DISABLE SHARED COMPLEXES PROPAGATION

          if (DEBUG && (mol.myComIndex >= complexList.size())) {
            cout << "Rank: " << mpiContext.rank
                 << ", mol.myComIndex = " << mol.myComIndex
                 << ", complexList.size() = " << complexList.size() << endl;
            error("2: mol.myComIndex >= complexList.size()");
          }
          complexList[mol.myComIndex].propagate(moleculeList, membraneObject, molTemplateList);
        }
      }
    }
  }  // for (auto& mol : moleculeList)
}