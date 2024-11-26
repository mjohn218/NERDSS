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
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "split.cpp"
#include "system_setup/system_setup.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

/**
 * Parses input files and initializes simulation parameters for a new
 * simulation.
 *
 * @param paramFile Name of the parameter file.
 * @param params Simulation parameters.
 * @param observablesList Map of observables to track during simulation.
 * @param forwardRxns List of forward reactions.
 * @param backRxns List of backward reactions.
 * @param createDestructRxns List of create/destruct reactions.
 * @param molTemplateList List of molecular templates.
 * @param membraneObject Membrane object.
 * @param moleculeList List of individual molecules.
 * @param complexList List of complexes.
 * @param simulVolume Simulation volume object.
 * @param observablesFileName Name of the file to write observables data.
 * @param implicitlipidIndex Index of the implicit lipid.
 * @param simItr Current simulation iteration.
 * @param mpiContext MPI context.
 * @param trajFileName Name of the file to write trajectory data.
 * @param transitionFileName Name of the file to write transition matrix data.
 */
void parse_input_for_a_new_simulation(
    std::string paramFile, Parameters& params,
    std::map<std::string, int>& observablesList,
    std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<MolTemplate>& molTemplateList, Membrane& membraneObject,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    SimulVolume& simulVolume, std::string observablesFileName,
    int implicitlipidIndex, long long int simItr, MpiContext& mpiContext,
    std::string trajFileName, std::string transitionFileName) {
  std::cout << "This is a new simulation with input file: " << paramFile
            << '\n';

  // Parse the input files
  parse_input(paramFile, params, observablesList, forwardRxns, backRxns,
              createDestructRxns, molTemplateList, membraneObject);

  // Set the number of states of implicit lipid
  for (auto& molTemplateTmp : molTemplateList) {
    if (molTemplateTmp.isImplicitLipid == true) {
      membraneObject.nStates =
          static_cast<int>(molTemplateTmp.interfaceList[0].stateList.size());
      break;
    }
  }

  // Initialize numberOfFreeLipidsEachState & numberOfProteinEachState
  for (int tmpStateIndex = 0; tmpStateIndex < membraneObject.nStates;
       tmpStateIndex++) {
    membraneObject.numberOfFreeLipidsEachState.emplace_back(0);
    membraneObject.numberOfProteinEachState.emplace_back(0);
  }

  // Verify the implicit lipid is the first
  for (auto& tempMolTemplate : molTemplateList) {
    if (tempMolTemplate.isImplicitLipid == true &&
        tempMolTemplate.molTypeIndex != 0) {
      error("implicit Lipid must be the first molecule type, exiting.");
    }
  }

  MolTemplate::numMolTypes = molTemplateList.size();

  std::cout << "NUMBER OF MOLECULE TYPES: " << params.numMolTypes
            << "NUMBER OF INTERFACES PLUS STATES, including PRODUCTS: "
            << params.numTotalSpecies << '\n';

  // Write the Observables file header and initial values
  std::ofstream observablesFile{observablesFileName};
  if (observablesList.size() == 1) {
    observablesFile << "Time (s)," << observablesList.begin()->first << '\n';
    observablesFile << "0,0\n";
  } else if (observablesList.size() > 1) {
    observablesFile << "Time (s)";
    for (auto obsItr = observablesList.begin(); obsItr != observablesList.end();
         ++obsItr) {
      observablesFile << ',' << obsItr->first;
    }
    std::cout << '\n' << "0";
    for (auto obsItr = observablesList.begin(); obsItr != observablesList.end();
         ++obsItr) {
      std::cout << ',' << obsItr->second;
    }
    observablesFile << '\n' << std::flush;
  }
  observablesFile.close();

  unsigned long reservation{};
  for (auto& molTemp : molTemplateList) {
    reservation += molTemp.copies;
  }
  moleculeList.reserve(reservation);
  complexList.reserve(reservation);

  // Create water box for sphere boundary
  if (membraneObject.isSphere) {
    membraneObject.create_water_box();
    membraneObject.sphereVol =
        (4.0 * M_PI * pow(membraneObject.sphereR, 3.0)) / 3.0;
  }

  // Generate the system coordinates, write out coordinate and topology files
  generate_coordinates(params, moleculeList, complexList, molTemplateList,
                       forwardRxns, membraneObject);
  write_psf(params, moleculeList, molTemplateList);

  // Set up some important parameters for implicit-lipid model
  initialize_paramters_for_implicitlipid_model(
      implicitlipidIndex, params, forwardRxns, backRxns, moleculeList,
      molTemplateList, complexList, membraneObject);

  // Initialize the starting copy number for each state
  initialize_states(moleculeList, molTemplateList, membraneObject);

  // Create simulation box cells
  std::cout << '\n' << "Partitioning simulation box into sub-boxes..." << '\n';
  set_rMaxLimit(params, molTemplateList, forwardRxns, 0, 0);

  simulVolume.create_simulation_volume(params, membraneObject);
  simulVolume.update_memberMolLists(params, moleculeList, complexList,
                                    molTemplateList, membraneObject, simItr,
                                    mpiContext);
  simulVolume.display();

  // Write beginning of trajectory
  std::ofstream trajFile{trajFileName};
  write_traj(0, trajFile, params, moleculeList, molTemplateList,
             membraneObject);
  trajFile.close();

  // Initialize transition matrix for each molType
  for (auto& molTemp : molTemplateList) {
    if (molTemp.countTransition == true) {
      molTemp.transitionMatrix.resize(molTemp.transitionMatrixSize);
      molTemp.lifeTime.resize(molTemp.transitionMatrixSize);
      for (int indexOne = 0; indexOne < molTemp.transitionMatrixSize;
           ++indexOne) {
        molTemp.transitionMatrix[indexOne].resize(molTemp.transitionMatrixSize);
      }
    }
  }

  // Write beginning of transition matrix
  std::ofstream transitionFile{transitionFileName};
  write_transition(0, transitionFile, molTemplateList);
  transitionFile.close();
}