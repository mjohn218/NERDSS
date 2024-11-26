/**
 * @file
 * nerdss_mpi.cpp
 *
 * @brief
 * Main function for simulation with parallel nerdss.
 *
 */
#include <stdio.h>
#include <sys/stat.h>
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

#ifdef ENABLE_PROFILING
#include <gperftools/profiler.h>
#endif

using namespace std;

/**
 * @var gsl_rng* r
 * @brief Global random number generator.
 */
gsl_rng* r;

/**
 * @typedef MDTimer
 * @brief Alias for std::chrono::system_clock used for timing.
 */
using MDTimer = std::chrono::system_clock;

/**
 * @typedef timeDuration
 * @brief Alias for std::chrono::duration used for timing.
 */
using timeDuration = std::chrono::duration<double, std::chrono::seconds>;

/**
 * @var randNum
 * @brief Global variable to store a random number.
 */
long long randNum = 0;

/**
 * @var totMatches
 * @brief Global variable to store the total number of matches.
 */
unsigned long totMatches = 0;

/**
 * @brief Main function that serves as an entrypoint to the program.
 *
 * @param argc The number of command line arguments passed to the program.
 * @param argv An array of C-strings holding the command line arguments.
 * @return An integer indicating the exit status of the program.
 */
int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  double setup_start, setup_end;
  double comp_start, comp_end;
  double comm_start, comm_end;
  double total_setup_time, total_comp_time, total_comm_time;

  total_setup_time = 0.0;
  total_comp_time = 0.0;
  total_comm_time = 0.0;

  setup_start = MPI_Wtime();

  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MpiContext mpiContext(size, NEIGHBOR_BUFFER_SIZE);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiContext.rank);

#ifdef ENABLE_PROFILING
  std::string profileFileName{"profile_output_" +
                              std::to_string(mpiContext.rank) + ".prof"};
  ProfilerStart(profileFileName.c_str());
#endif

  mkdir("PDB", 0777);
  mkdir("RESTARTS", 0777);
  mkdir("mergePDB", 0777);
  mkdir("mergeOUT", 0777);

  // Initialize random_device object to generate seed for the generator
  std::random_device rd{};
  unsigned seed{rd()};

  // Start timer
  MDTimer::time_point totalTimeStart = MDTimer::now();

  // This will store the list of molecule templates
  std::vector<MolTemplate> molTemplateList{};

  // These will store the list of reactions
  std::vector<ForwardRxn> forwardRxns{};
  std::vector<BackRxn> backRxns{};
  std::vector<CreateDestructRxn> createDestructRxns{};

  // Pre-allocate memory for the reactions vector to avoid reallocation during
  // runtime The `reserve()` function pre-allocates memory for the vector to
  // hold 10 elements.
  forwardRxns.reserve(10);
  backRxns.reserve(10);
  createDestructRxns.reserve(10);

  // This will store the system parameters
  Parameters params{};

  // Set filenames using string interpolation
  std::string observablesFileName{"observables_time_" +
                                  std::to_string(mpiContext.rank) + ".dat"};
  std::string trajFileName{"trajectory_" + std::to_string(mpiContext.rank) +
                           ".xyz"};
  std::string transitionFileName{"transition_matrix_time_" +
                                 std::to_string(mpiContext.rank) + ".dat"};
  std::string restartFileName{"restart_" + std::to_string(mpiContext.rank) +
                              ".dat"};
  std::string addFileNameInput{};
  std::string paramFile{};
  std::string restartFileNameInput;
  std::string debugFileName{"debug_output_" + std::to_string(mpiContext.rank) +
                            ".dat"};

  std::string outFileName{"output_" + std::to_string(mpiContext.rank) + ".txt"};
  if (freopen(outFileName.c_str(), "a", stdout) == nullptr) {
    std::cerr << "Failed to open the output file.\n";
    return 1;
  }

  // Create and write into debug output file
  std::ofstream debugFile{debugFileName};
  debugFile << "time, dissociation, association\n";

  // Command line flag parser
  parse_command(argc, argv, params, paramFile, restartFileNameInput,
                addFileNameInput, seed);

  if (DEBUG){
    // use specific seed for debugging
    if(mpiContext.rank == 0){
      seed = 482364341;
    }
    if(mpiContext.rank == 1){
      seed = 689863503;
    }
  }

  auto startTime = MDTimer::to_time_t(totalTimeStart);
  char charTime[24];
  std::cout << "\nStart date: ";
  if (0 <
      strftime(charTime, sizeof(charTime), "%F %T", std::localtime(&startTime)))
    std::cout << charTime << '\n';
  std::cout << "RNG Seed: " << seed << std::endl;

  // Random number generator
  const gsl_rng_type* T;
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed);

  // 2D reaction probability tables, used to calculate reaction probability
  std::vector<gsl_matrix*> survMatrices;
  std::vector<gsl_matrix*> normMatrices;
  std::vector<gsl_matrix*> pirMatrices;
  double* tableIDs = new double[params.max2DRxns * 2];

  // list of the observables
  std::map<std::string, int> observablesList;

  // simulVolume object inlcuding the information of the division of the system
  // into subVolumes
  SimulVolume simulVolume{};

  // list of all molecules in the system
  std::vector<Molecule> moleculeList{};

  // list of all complexes in the system
  std::vector<Complex> complexList{};

  // class structure that contains boundary conditions, implicit lipid model.
  Membrane membraneObject;

  // contains arrays tracking bound molecule pairs, and species copy nums.
  copyCounters counterArrays;

  // implicit-lipid index, which is also stored in
  // membraneObject.implicitlipidIndex.
  int implicitlipidIndex{0};

  // simulation iteration count
  long long int simItr{0};

  int totalSpeciesNum{0};  ///< total species num after add

  // initialize event counters to zero. Restart will update any non-zero values.
  init_association_events(counterArrays);

  std::cout << "\nParsing Input: " << std::endl;
  if (!params.fromRestart && paramFile != "") {
    parse_input_for_a_new_simulation(
        paramFile, params, observablesList, forwardRxns, backRxns,
        createDestructRxns, molTemplateList, membraneObject, moleculeList,
        complexList, simulVolume, observablesFileName, implicitlipidIndex,
        simItr, mpiContext, trajFileName, transitionFileName);
  } else if (params.fromRestart) {
    parse_input_for_a_restart_simulation(
        restartFileNameInput, addFileNameInput, params, observablesList,
        forwardRxns, backRxns, createDestructRxns, molTemplateList,
        membraneObject, moleculeList, complexList, simulVolume,
        observablesFileName, implicitlipidIndex, simItr, mpiContext,
        trajFileName, transitionFileName, seed, counterArrays);
  } else {
    error(
        "Please provide a parameter and/or restart file. Parameter file Syntax "
        "is : ./nerdss -f parameterfile.inp or ./nerdss -r restartfile.dat \n");
  }

  params.checkUnimoleculeReactionPopulation = false;
  mpiContext.checkUnimoleculeReactionPopulation =
      params.checkUnimoleculeReactionPopulation;
  params.hasRankCommunicationForLargeComplex = false;
  mpiContext.hasRankCommunicationForLargeComplex =
      params.hasRankCommunicationForLargeComplex;

  if (membraneObject.implicitLipid == true) params.implicitLipid = true;
  for (auto& oneReaction : forwardRxns) {
    if (oneReaction.rxnType == ReactionType::uniMolStateChange) {
      params.hasUniMolStateChange = true;
      break;
    }
  }
  if (createDestructRxns.empty() == false) {
    params.hasCreationDestruction = true;
    params.isNonEQ = true;

    // Set molTemp.canDestroy = ture
    for (auto& oneReaction : createDestructRxns) {
      if (oneReaction.rxnType == ReactionType::destruction) {
        molTemplateList[oneReaction.reactantMolList.at(0).molTypeIndex]
            .canDestroy = true;
      }
    }
  }

  // Setup output files
  char fnameProXYZ[100];

  std::string fileName{};
  fileName =
      "histogram_complexes_time_" + std::to_string(mpiContext.rank) + ".dat";
  std::ofstream assemblyfile{fileName};
  fileName = "mono_dimer_time_" + std::to_string(mpiContext.rank) + ".dat";
  std::ofstream dimerfile{fileName};
  fileName = "event_counters_time_" + std::to_string(mpiContext.rank) + ".dat";
  std::ofstream eventFile{fileName};
  fileName = "bound_pair_time_" + std::to_string(mpiContext.rank) + ".dat";
  std::ofstream pairOutfile{fileName};
  fileName = "copy_numbers_time_" + std::to_string(mpiContext.rank) + ".dat";
  std::ofstream speciesFile1{fileName};

  // Set the interval to write checkpoint and transition matrix file if it is
  // not specified
  if (params.checkPoint == -1) {
    params.checkPoint = params.nItr / 10;
  }
  if (params.transitionWrite == -1) {
    params.transitionWrite = params.nItr / 10;
  }

  set_exclude_volume_bound(forwardRxns, molTemplateList);

  for (auto& oneComplex : complexList) {
    oneComplex.update_properties(moleculeList, molTemplateList);
  }

  Parameters::dt = params.timeStep;
  Parameters::lastUpdateTransition.resize(molTemplateList.size());

  // Print out system information
  std::cout << "\nSimulation Parameters\n";
  params.display();
  membraneObject.display();
  std::cout << "\nMolecule Information\n";
  display_all_MolTemplates(molTemplateList);
  std::cout << "\nReactions\n";
  display_all_reactions(forwardRxns, backRxns, createDestructRxns);

  if (mpiContext.nprocs) {  // if parallel execution
    if (!params.fromRestart) {
      prepare_data_structures_for_parallel_execution(
          moleculeList, simulVolume, membraneObject, molTemplateList, params,
          forwardRxns, backRxns, createDestructRxns, counterArrays, mpiContext,
          complexList, pairOutfile);
    }
  }

  int meanComplexSize{0};

  totalSpeciesNum = init_speciesFile(speciesFile1, counterArrays,
                                     molTemplateList, forwardRxns, params);
  init_counterCopyNums(counterArrays, moleculeList, complexList,
                       molTemplateList, membraneObject, totalSpeciesNum,
                       params);

  write_all_species((simItr - params.itrRestartFrom) * params.timeStep *
                            Constants::usToSeconds +
                        params.timeRestartFrom,
                    moleculeList, speciesFile1, counterArrays, membraneObject,
                    mpiContext, simulVolume);

  init_print_dimers(dimerfile, params, molTemplateList);
  init_NboundPairs(
      counterArrays, pairOutfile, params, molTemplateList,
      moleculeList);  // initializes to zero, re-calculated for a restart!!
  write_NboundPairs(counterArrays, pairOutfile, simItr, params, moleculeList);
  print_dimers(complexList, dimerfile, simItr, params, molTemplateList);
  print_association_events(counterArrays, eventFile, simItr, params);

  int number_of_lipids = 0;  ///< sum of all states of IL
  for (size_t i = 0; i < membraneObject.numberOfFreeLipidsEachState.size();
       i++) {
    number_of_lipids += membraneObject.numberOfFreeLipidsEachState[i];
  }
  meanComplexSize = print_complex_hist(
      complexList, assemblyfile, simItr, params, molTemplateList, moleculeList,
      number_of_lipids, mpiContext, simulVolume);

  simulVolume.update_memberMolLists(params, moleculeList, complexList,
                                    molTemplateList, membraneObject, simItr,
                                    mpiContext);
  simulVolume.update_region_flags(moleculeList, complexList, mpiContext);

  std::cout << "*************** BEGIN SIMULATION **************** "
            << std::endl;

  //  begin the timer
  MDTimer::time_point simulTimeStart = MDTimer::now();
  int numSavedDurations{1000};
  std::vector<std::chrono::duration<double>> durationList(numSavedDurations);

  std::fill(durationList.begin(), durationList.end(),
            std::chrono::duration<double>(simulTimeStart - totalTimeStart));

  unsigned DDTableIndex{0};

  // Vectors to store binding probabilities for implicit lipids in 2D.
  std::vector<double> IL2DbindingVec{};
  std::vector<double> IL2DUnbindingVec{};
  std::vector<double> ILTableIDs{};

  if (params.fromRestart == true) {
    try {
      read_rng_state(mpiContext.rank);  // read the current RNG state
    } catch (std::exception& e) {
      std::cout << "Read rng_state failed. Set the rng_state with the seed: "
                << seed << std::endl;
      gsl_rng_set(r, seed);
    }

    std::ofstream restartFile{restartFileName,
                              std::ios::out};  // to show different from append

    write_rng_state(mpiContext.rank);  // write the current RNG state
    write_restart(simItr, restartFile, params, simulVolume, moleculeList,
                  complexList, molTemplateList, forwardRxns, backRxns,
                  createDestructRxns, observablesList, membraneObject,
                  counterArrays);
    restartFile.close();
  }

  long long int startSimItr = simItr + 1;
  long long int stopSimItr = params.nItr - 1;

  // write the pdb of the initial system
  write_pdb(simItr, simItr, params, moleculeList, molTemplateList,
            membraneObject, mpiContext.rank, mpiContext, simulVolume);

  setup_end = MPI_Wtime();
  total_setup_time = setup_end - setup_start;

  for (simItr += 1; simItr <= params.nItr; ++simItr) {
    // write the debug output
    debugFile << simItr * params.timeStep * Constants::usToSeconds << ", ";

    comp_start = MPI_Wtime();

    if (VERBOSE) {
      printf("simItr (%lld) start on rank (%d)...\n", simItr, mpiContext.rank);
    }

    mpiContext.simItr = simItr;

    if (DEBUG) {
      DEBUG_FIND_MOL("simItr BEGINS");
      DEBUG_FIND_COMPLEX("simItr BEGINS");
      debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                       "simItr BEGINS");
    }

    propCalled = 0;
    MDTimer::time_point startStep = MDTimer::now();

    if (VERBOSE) {
      printf(
          "Check and perform 0th and 1st order reactions for non-ghost "
          "molecules...\n");
    }
    check_perform_zeroth_first_order_reactions(
        simItr, params, moleculeList, complexList, simulVolume, forwardRxns,
        backRxns, createDestructRxns, molTemplateList, observablesList,
        counterArrays, membraneObject, IL2DbindingVec, IL2DUnbindingVec,
        ILTableIDs, mpiContext, debugFile);

    // Measure separations between proteins in neighboring cells to identify all
    // possible reactions.
    if (VERBOSE) {
      printf("Measure separations between two reactants...\n");
    }

    measure_separations_to_identify_possible_reactions(
        simItr, params, moleculeList, complexList, simulVolume, forwardRxns,
        backRxns, createDestructRxns, molTemplateList, observablesList,
        counterArrays, membraneObject, IL2DbindingVec, IL2DUnbindingVec,
        ILTableIDs, normMatrices, survMatrices, pirMatrices, implicitlipidIndex,
        tableIDs, DDTableIndex);

    // Now that separations and reaction probabilities are calculated, decide
    // whether to perform reactions for each protein.

    // Divide the subVolume into two parts and execute each part
    // synchronously with all processes
    if (VERBOSE) {
      printf("Assgin molecules to left and right parts...\n");
    }
    std::vector<int> left{};
    std::vector<int> right{};
    for (int i = 0; i < moleculeList.size(); i++) {
      int x_bin = get_x_bin(mpiContext, moleculeList[i]);
      if (mpiContext.rank == 0) {
        if (x_bin < (simulVolume.numSubCells.x - 1) / 2)
          {left.push_back(i);
          moleculeList[i].isLeftHalf = true;}
        else
          {right.push_back(i);
          moleculeList[i].isRightHalf = true;}
      } else if (mpiContext.rank == mpiContext.nprocs - 1) {
        if (x_bin < (simulVolume.numSubCells.x + 1) / 2)
          {left.push_back(i);
          moleculeList[i].isLeftHalf = true;}
        else
          {right.push_back(i);
          moleculeList[i].isRightHalf = true;}
      } else {
        if (x_bin < (simulVolume.numSubCells.x) / 2)
          {left.push_back(i);
          moleculeList[i].isLeftHalf = true;}
        else
          {right.push_back(i);
          moleculeList[i].isRightHalf = true;}
      }
    }

    if (VERBOSE) {
      cout << "Rank " << mpiContext.rank << " has " << left.size()
           << " molecules in the left part and " << right.size()
           << " molecules in the right part." << endl;
    }

    // do the left part, send to left and recieve from right
    if (VERBOSE) {
      printf("Perform reactions of molecules in the left part...\n");
    }

    int associationCount{0};

    perform_bimolecular_reactions(
        simItr, params, moleculeList, complexList, simulVolume, forwardRxns,
        backRxns, createDestructRxns, molTemplateList, observablesList,
        counterArrays, membraneObject, left, true, debugFile, associationCount);

    if (VERBOSE) {
      printf("Checking overlap and propagate molecules in the left part...\n");
    }
    check_overlap(left, true, simItr, params, moleculeList, complexList,
                  simulVolume, forwardRxns, backRxns, createDestructRxns,
                  molTemplateList, observablesList, counterArrays,
                  membraneObject, mpiContext);

    // update simVolume memberMolLists and update region flags of molecules and
    // complexes
    simulVolume.update_memberMolLists(params, moleculeList, complexList,
                                      molTemplateList, membraneObject, simItr,
                                      mpiContext);
    simulVolume.update_region_flags(moleculeList, complexList, mpiContext);

    if (DEBUG) {
      DEBUG_FIND_MOL("After Perform left half");
      DEBUG_FIND_COMPLEX("After Perform left half");
      debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                       "After Perform left half");
    }

    set<int> moleculeIndexesCheckReceive{};
    set<int> complexIndexesCheckReceive{};

    if (mpiContext.rank < mpiContext.nprocs - 1) {
      // Mark all the right ghost and right edge molecules and complexes to not
      // received

      for (int yItr{0}; yItr < simulVolume.numSubCells.y; ++yItr) {
        for (int zItr{0}; zItr < simulVolume.numSubCells.z; ++zItr) {
          int rightBoxIndex =
              (simulVolume.numSubCells.x - 1) +
              yItr * simulVolume.numSubCells.x +
              zItr * simulVolume.numSubCells.x * simulVolume.numSubCells.y;
          for (auto& molIndex :
               simulVolume.subCellList[rightBoxIndex].memberMolList) {
            moleculeList[molIndex].receivedFromNeighborRank = false;
            moleculeIndexesCheckReceive.insert(molIndex);
            complexList[moleculeList[molIndex].myComIndex]
                .receivedFromNeighborRank = false;
            complexIndexesCheckReceive.insert(moleculeList[molIndex].myComIndex);
          }
          rightBoxIndex -= 1;
          for (auto& molIndex :
               simulVolume.subCellList[rightBoxIndex].memberMolList) {
            moleculeList[molIndex].receivedFromNeighborRank = false;
            moleculeIndexesCheckReceive.insert(molIndex);
            complexList[moleculeList[molIndex].myComIndex]
                .receivedFromNeighborRank = false;
            complexIndexesCheckReceive.insert(moleculeList[molIndex].myComIndex);
            for (auto& i :
                 complexList[moleculeList[molIndex].myComIndex].memberList) {
              moleculeList[i].receivedFromNeighborRank = false;
              moleculeIndexesCheckReceive.insert(i);
            }
          }
        }
      }
      // cout << "At simItr " << simItr << endl;
      // cout << "Rank " << mpiContext.rank << " has " << moleculeIndexesCheckReceive.size()
      //    << " molecules and " << complexIndexesCheckReceive.size()
      //    << " complexes to check receive from Right side." << endl;
    }

    comp_end = MPI_Wtime();
    total_comp_time += comp_end - comp_start;
    comm_start = MPI_Wtime();

    // Send all the left ghost and left edge molecules and complexes to Left
    // neighbor processor
    if (mpiContext.rank > 0) {
      if (VERBOSE) {
        printf("Sending data from rank (%d) to rank (%d)...\n", mpiContext.rank,
               mpiContext.rank - 1);
      }
      send_data_to_left_neighboring_ranks(
          mpiContext, simItr, params, simulVolume, left, moleculeList,
          complexList, molTemplateList, membraneObject);
      if (VERBOSE) {
        printf("Sending data from rank (%d) to rank (%d) done.\n",
               mpiContext.rank, mpiContext.rank - 1);
      }
    }

    // recieve from right
    if (mpiContext.rank < mpiContext.nprocs - 1) {
      if (VERBOSE) {
        printf("Receiving ids from rank (%d) to rank (%d)...\n",
               mpiContext.rank + 1, mpiContext.rank);
      }
      receive_right_neighborhood_zones(
          mpiContext, simItr, simulVolume, right, moleculeList, complexList,
          molTemplateList, membraneObject, counterArrays);
      if (VERBOSE) {
        printf("Receiving ids from rank (%d) to rank (%d) done\n",
               mpiContext.rank + 1, mpiContext.rank);
      }
    }

    comm_end = MPI_Wtime();
    total_comm_time += comm_end - comm_start;
    comp_start = MPI_Wtime();

    // delete the molecules and complexes that not received from the right
    int delMol = 0;
    int delCom = 0;
    if (mpiContext.rank < mpiContext.nprocs - 1) {
      for (auto i : moleculeIndexesCheckReceive) {
        if (moleculeList[i].receivedFromNeighborRank == false) {
          // cout << "Rank " << mpiContext.rank << " remove molecule id " << moleculeList[i].id
          //      << " bc not recv from right." << endl;
          moleculeList[i].MPI_remove_from_one_rank(moleculeList, complexList);
          delMol++;
        }
      }

      for (auto i : complexIndexesCheckReceive) {
        if (complexList[i].receivedFromNeighborRank == false) {
          // cout << "Rank " << mpiContext.rank << " remove complex id " << complexList[i].id
          //      << " bc not recv from right." << endl;
          complexList[i].memberList.clear();
          complexList[i].destroy(moleculeList, complexList);
          delCom++;
        }
      }
      // cout << "At simItr " << simItr << endl;
      // cout << "Rank " << mpiContext.rank << " deleted " << delMol
      //    << " molecules and " << delCom << " complexes that not received from Right side."
      //    << endl;
    }

    if (DEBUG) {
      DEBUG_FIND_MOL("After Receive from right");
      DEBUG_FIND_COMPLEX("After Receive from right");
      debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                       "After Receive from right");
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////

    // do the right part, send to right and recieve from left
    if (VERBOSE) {
      printf("Update the states of molecules in the right part...\n");
    }
    perform_bimolecular_reactions(
        simItr, params, moleculeList, complexList, simulVolume, forwardRxns,
        backRxns, createDestructRxns, molTemplateList, observablesList,
        counterArrays, membraneObject, right, false, debugFile,
        associationCount);

    if (VERBOSE) {
      printf("Checking overlap and propagate molecules in the right part...\n");
    }
    check_overlap(right, false, simItr, params, moleculeList, complexList,
                  simulVolume, forwardRxns, backRxns, createDestructRxns,
                  molTemplateList, observablesList, counterArrays,
                  membraneObject, mpiContext);

    // update simVolume memberMolLists and update region flags of molecules and
    // complexes
    simulVolume.update_memberMolLists(params, moleculeList, complexList,
                                      molTemplateList, membraneObject, simItr,
                                      mpiContext);
    simulVolume.update_region_flags(moleculeList, complexList, mpiContext);

    if (DEBUG) {
      DEBUG_FIND_MOL("After Perform right half");
      DEBUG_FIND_COMPLEX("After Perform right half");
      debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                        "After Perform right half");
    }

    set<int> moleculeIndexesCheckReceive2{};
    set<int> complexIndexesCheckReceive2{};

    // Mark all the left ghost and left edge molecules and complexes to not
    // received
    if (mpiContext.rank > 0) {
      for (int yItr{0}; yItr < simulVolume.numSubCells.y; ++yItr) {
        for (int zItr{0}; zItr < simulVolume.numSubCells.z; ++zItr) {
          int leftBoxIndex =
              yItr * simulVolume.numSubCells.x +
              zItr * simulVolume.numSubCells.x * simulVolume.numSubCells.y;
          for (auto& molIndex :
               simulVolume.subCellList[leftBoxIndex].memberMolList) {
            moleculeList[molIndex].receivedFromNeighborRank = false;
            moleculeIndexesCheckReceive2.insert(molIndex);
            complexList[moleculeList[molIndex].myComIndex]
                .receivedFromNeighborRank = false;
            complexIndexesCheckReceive2.insert(moleculeList[molIndex].myComIndex);
          }
          leftBoxIndex += 1;
          for (auto& molIndex :
               simulVolume.subCellList[leftBoxIndex].memberMolList) {
            moleculeList[molIndex].receivedFromNeighborRank = false;
            moleculeIndexesCheckReceive2.insert(molIndex);
            complexList[moleculeList[molIndex].myComIndex]
                .receivedFromNeighborRank = false;
            complexIndexesCheckReceive2.insert(moleculeList[molIndex].myComIndex);
            for (auto& i :
                 complexList[moleculeList[molIndex].myComIndex].memberList) {
              moleculeList[i].receivedFromNeighborRank = false;
              moleculeIndexesCheckReceive2.insert(i);
            }
          }
        }
      }
      // cout << "At simItr " << simItr << endl;
      // cout << "Rank " << mpiContext.rank << " has " << moleculeIndexesCheckReceive2.size()
      //    << " molecules and " << complexIndexesCheckReceive2.size()
      //    << " complexes to check receive from Left side." << endl;
    }

    comp_end = MPI_Wtime();
    total_comp_time += comp_end - comp_start;
    comm_start = MPI_Wtime();

    // send to right
    if (mpiContext.rank < mpiContext.nprocs - 1) {
      if (VERBOSE) {
        printf("Sending data from rank (%d) to rank (%d)...\n", mpiContext.rank,
               mpiContext.rank + 1);
      }
      send_data_to_right_neighboring_ranks(
          mpiContext, simItr, params, simulVolume, right, moleculeList,
          complexList, molTemplateList, membraneObject);
      if (VERBOSE) {
        printf("Sending data from rank (%d) to rank (%d) done.\n",
               mpiContext.rank, mpiContext.rank + 1);
      }
    }

    // recieve from left
    if (mpiContext.rank > 0) {
      if (VERBOSE) {
        printf("Receiving data from rank (%d) to rank (%d)...\n",
               mpiContext.rank - 1, mpiContext.rank);
      }
      receive_left_neighborhood_zones(
          mpiContext, simItr, simulVolume, left, moleculeList, complexList,
          molTemplateList, membraneObject, counterArrays);
      if (VERBOSE) {
        printf("Receiving data from rank (%d) to rank (%d) done\n",
               mpiContext.rank - 1, mpiContext.rank);
      }
    }

    comm_end = MPI_Wtime();
    total_comm_time += comm_end - comm_start;
    comp_start = MPI_Wtime();

    int delMol2 = 0;
    int delCom2 = 0;
    // delete the molecules and complexes that not received from the left
    if (mpiContext.rank > 0) {
      for (auto i : moleculeIndexesCheckReceive2) {
        if (moleculeList[i].receivedFromNeighborRank == false) {
          moleculeList[i].MPI_remove_from_one_rank(moleculeList, complexList);
          delMol2++;
        }
      }

      for (auto i : complexIndexesCheckReceive2) {
        if (complexList[i].receivedFromNeighborRank == false) {
          complexList[i].memberList.clear();
          complexList[i].destroy(moleculeList, complexList);
          delCom2++;
        }
      }
      // cout << "At simItr " << simItr << endl;
      // cout << "Rank " << mpiContext.rank << " deleted " << delMol2
      //    << " molecules and " << delCom2 << " complexes that not received from Left side."
      //    << endl;
    }

    if (DEBUG) {
      DEBUG_FIND_MOL("After Receive from left");
      DEBUG_FIND_COMPLEX("After Receive from left");
      debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                       "After Receive from left");
    }

    remove_empty_slots(simItr, params, moleculeList, complexList, simulVolume,
                       forwardRxns, backRxns, createDestructRxns,
                       molTemplateList, observablesList, counterArrays,
                       membraneObject, mpiContext);

    for (auto& oneMol : moleculeList) {
      if (oneMol.isImplicitLipid) continue;

      if (oneMol.isEmpty) {
        // no empty are allowed here
        cout << "Error: Empty molecule at simItr " << simItr << endl;
        oneMol.print(mpiContext);
        exit(1);
      }

      clear_reweight_vecs(oneMol);
      oneMol.trajStatus = TrajStatus::none;
      oneMol.isDissociated = false;
      oneMol.isAssociated = false;
      oneMol.isLeftHalf = false;
      oneMol.isRightHalf = false;
      oneMol.receivedFromNeighborRank = true;
      oneMol.crossbase.clear();
      oneMol.mycrossint.clear();
      oneMol.crossrxn.clear();
      oneMol.probvec.clear();

      // update complexes. if the complex has more than one member, it'll be
      // updated more than once, but this is probably more efficient than
      // iterating over the complexes afterwards
      if (complexList[oneMol.myComIndex].isEmpty) {
        cout << "Error: Empty complex at simItr " << simItr << endl;
        complexList[oneMol.myComIndex].print(mpiContext);
        exit(1);
      }

      complexList[oneMol.myComIndex].ncross = 0;
      complexList[oneMol.myComIndex].trajStatus = TrajStatus::none;
      complexList[oneMol.myComIndex].receivedFromNeighborRank = true;
    }

    if (DEBUG) {
      set<int> comIndexes{};
      for (auto& oneMol : moleculeList) {
        comIndexes.insert(oneMol.myComIndex);
      }
      // assert that comIndexes are the same as the indexes of the complexes, and there is no difference
      for (auto& oneCom : complexList) {
        if (comIndexes.find(oneCom.index) == comIndexes.end()) {
          cout << "Error: Complex index " << oneCom.index
                    << " is not in the comIndexes set.\n";
          exit(1);
        }
      }
      if (comIndexes.size() != complexList.size()) {
        cout << "Error: comIndexes size " << comIndexes.size()
                  << " is not equal to complexList size " << complexList.size()
                  << ".\n";
        exit(1);
      }
    }

    simulVolume.update_memberMolLists(params, moleculeList, complexList,
                                      molTemplateList, membraneObject, simItr,
                                      mpiContext);
    simulVolume.update_region_flags(moleculeList, complexList, mpiContext);

    if (DEBUG) {
      DEBUG_FIND_MOL("Iteration finished");
      DEBUG_FIND_COMPLEX("Iteration finished");
      debug_molecule_complex_missmatch(mpiContext, moleculeList, complexList,
                                       "Iteration finished");
    }

    //------------------------------------------------------------------------------------

    write_output(simItr, params, trajFileName, moleculeList, molTemplateList,
                 membraneObject, mpiContext, transitionFileName, counterArrays,
                 speciesFile1, debugFile, debugFileName, restartFileName,
                 complexList, simulVolume, forwardRxns, backRxns,
                 createDestructRxns, durationList, startStep, totalTimeStart,
                 observablesList, assemblyfile);
    
    comp_end = MPI_Wtime();
    total_comp_time += comp_end - comp_start;

    if (simItr % params.timeWrite == 0) {
      std::cout << "total computation time: " << total_comp_time
                << " seconds\n";
      std::cout << "total communication time: " << total_comm_time
                << " seconds\n";

      // MERGE THE OUTPUTS FROM ALL THE PROCESSORS
      // if (mpiContext.rank == 0) {
      //   merge_outputs(mpiContext.nprocs, molTemplateList.size());
      // }
    }
  }  // end iterating over time steps

  // Write out final result
  std::cout << llinebreak << "End simulation\n";
  auto endTime = MDTimer::now();
  auto endTimeFormat = MDTimer::to_time_t(endTime);
  std::cout << "End date: ";
  //<< std::put_time(std::localtime(&endTimeFormat), "%F %T") << '\n';
  if (0 < strftime(charTime, sizeof(charTime), "%F %T",
                   std::localtime(&endTimeFormat)))
    std::cout << charTime << '\n';
  std::chrono::duration<double> wallTime = endTime - totalTimeStart;
  std::cout << "\tWall Time: ";
  std::cout << wallTime.count() << " seconds\n";
  std::cout << "total setup time: " << total_setup_time
            << " seconds\n";
  std::cout << "total computation time: " << total_comp_time
            << " seconds\n";
  std::cout << "total communication time: " << total_comm_time
            << " seconds\n";

  // MERGE THE OUTPUTS FROM ALL THE PROCESSORS
  if (mpiContext.rank == 0) {
    merge_outputs(mpiContext.nprocs, molTemplateList.size());
  }

  debugFile.close();

  delete[] tableIDs;
  gsl_rng_free(r);

#ifdef ENABLE_PROFILING
  ProfilerStop();
#endif

  MPI_Finalize();
  return 0;
}  // end main
