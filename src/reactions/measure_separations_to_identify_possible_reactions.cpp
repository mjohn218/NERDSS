#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
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

void measure_separations_to_identify_possible_reactions(
    unsigned simItr, Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, SimulVolume& simulVolume,
    std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays,
    Membrane& membraneObject, std::vector<double>& IL2DbindingVec,
    std::vector<double>& IL2DUnbindingVec, std::vector<double>& ILTableIDs,
    std::vector<gsl_matrix*>& normMatrices,
    std::vector<gsl_matrix*>& survMatrices,
    std::vector<gsl_matrix*>& pirMatrices, int implicitlipidIndex,
    double* tableIDs, unsigned& DDTableIndex) {
  for (unsigned cellItr{0}; cellItr < simulVolume.subCellList.size();
       ++cellItr) {
    for (unsigned memItr{0};
         memItr < simulVolume.subCellList[cellItr].memberMolList.size();
         ++memItr) {
      int targMolIndex{simulVolume.subCellList[cellItr].memberMolList[memItr]};
      if (moleculeList[targMolIndex].isImplicitLipid) continue;

      // Test bimolecular reactions, and binding to implicit-lipids
      if (moleculeList[targMolIndex].freelist.size() > 0 ||
          molTemplateList[moleculeList[targMolIndex].molTypeIndex]
                  .excludeVolumeBound == true) {
        // first, check for implicit-lipid binding
        int protype = moleculeList[targMolIndex].molTypeIndex;
        if (molTemplateList[protype].bindToSurface == true) {
          check_implicit_reactions(
              targMolIndex, implicitlipidIndex, simItr, params, moleculeList,
              complexList, molTemplateList, forwardRxns, backRxns,
              counterArrays, membraneObject, IL2DbindingVec, IL2DUnbindingVec,
              ILTableIDs);
        }

        // secondly, loop over proteins in your same cell.
        for (unsigned memItr2{memItr + 1};
             memItr2 < simulVolume.subCellList[cellItr].memberMolList.size();
             ++memItr2) {
          int partMolIndex{
              simulVolume.subCellList[cellItr].memberMolList[memItr2]};

          check_bimolecular_reactions(
              targMolIndex, partMolIndex, simItr, tableIDs, DDTableIndex,
              params, normMatrices, survMatrices, pirMatrices, moleculeList,
              complexList, molTemplateList, forwardRxns, backRxns,
              counterArrays, membraneObject);
        }  // loop over protein partners in your same cell

        // thirdly, loop over all neighboring cells, and all proteins in those
        // cells. for PBC, all cells have maxnbor neighbor cells. For
        // reflecting, edge have fewer.
        for (auto& neighCellItr :
             simulVolume.subCellList[cellItr].neighborList) {
          if (DEBUG && (neighCellItr > simulVolume.subCellList.size()))
            error("neighCellItr is higher than subCellList size!");
          for (unsigned memItr2{0};
               memItr2 <
               simulVolume.subCellList[neighCellItr].memberMolList.size();
               ++memItr2) {
            int partMolIndex{
                simulVolume.subCellList[neighCellItr].memberMolList[memItr2]};
            if (DEBUG && (partMolIndex > moleculeList.size()))
              error("partMolIndex is higher than moleculeList size!");
            check_bimolecular_reactions(
                targMolIndex, partMolIndex, simItr, tableIDs, DDTableIndex,
                params, normMatrices, survMatrices, pirMatrices, moleculeList,
                complexList, molTemplateList, forwardRxns, backRxns,
                counterArrays, membraneObject);
          }  // loop over all proteins in this neighbor cell
        }    // loop over all neighbor cells
      }      // if protein i is free to bind
    }        // loop over all proteins in initial cell
  }          // End looping over all cells.
}
