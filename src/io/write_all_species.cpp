/*! \file write_all_species.cpp
 *
 * \brief
 *
 * ### Created on 2019-06-05 by Matthew Varga
 */
#include "debug/debug.hpp"
#include "io/io.hpp"
#include "mpi/mpi_function.hpp"
#include "tracing.hpp"

void write_all_species(double simTime, std::vector<Molecule>& moleculeList,
                       std::ofstream& speciesFile, copyCounters& counterArrays,
                       const Membrane& membraneObject, MpiContext& mpiContext,
                       SimulVolume& simulVolume) {
  // FOR MONOMER, EACH PROCESSOR ONLY OUTPUT THE NON-GHOST MOLECULES; FOR PAIR,
  // EACH PROCESSOR OUTPUT THE PAIR IF THE LARGER ID OF THE PAIR IS NON-GHOST

  int i, j;

  int index;
  int p1, p2;
  for (i = 0; i < counterArrays.copyNumSpecies.size(); i++)
    counterArrays.copyNumSpecies[i] = 0;

  for (p1 = 0; p1 < moleculeList.size(); p1++) {
    int numIfaces = moleculeList[p1].interfaceList.size();
    bool isGhosted = (moleculeList[p1].isLeftGhost || moleculeList[p1].isRightGhost);
    for (j = 0; j < numIfaces; j++) {
      bool isPatnerGhosted = (moleculeList[moleculeList[p1].interfaceList[j].interaction.partnerIndex].isLeftGhost || 
                              moleculeList[moleculeList[p1].interfaceList[j].interaction.partnerIndex].isRightGhost);
      if (moleculeList[p1].interfaceList[j].interaction.partnerIndex == -1) {
        // MONOMER
        if (isGhosted == true) continue;
      } else {
        // PAIR
        if (moleculeList[p1].id >
            moleculeList[p1].interfaceList[j].interaction.partnerId) {
          if (isGhosted == true) continue;
        } else {
          if (isPatnerGhosted == true) continue;
        }
      }
      // find out which state each interface on the molecule is in, and
      // increment copyNumSpecies array.
      index = moleculeList[p1].interfaceList[j].index;
      if (moleculeList[p1].isImplicitLipid == false) {
        counterArrays.copyNumSpecies[index]++;
      } else {
        // For implicit lipid, set copy numbers based on read in template.
        // int molTypeIndex = moleculeList[p1].molTypeIndex;
        for (int tmpStateIndex = 0; tmpStateIndex < membraneObject.nStates;
             tmpStateIndex++) {
          counterArrays.copyNumSpecies[index + tmpStateIndex] =
              membraneObject.numberOfFreeLipidsEachState
                  [tmpStateIndex];  // IL mol must be the first place
        }
      }
    }  // end all interfaces
  }    // end all current molecules

  speciesFile << simTime;
  for (auto i = 0; i < counterArrays.copyNumSpecies.size(); i++) {
    if (counterArrays.singleDouble[i] == 2 &&
        counterArrays.implicitDouble[i] == false) {
      // product state, contains two proteins, so will be double counted above.
      speciesFile << ',' << counterArrays.copyNumSpecies[i] * 0.5;
    } else {
      speciesFile << ',' << counterArrays.copyNumSpecies[i];
    }
  }
  speciesFile << std::endl;
}
