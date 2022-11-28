/* \file nerdss.cpp
 * \brief Main function for simulation.
 *
 * ## TODO:
 *  - Change freelist from unordered list of relative iface indices to ordered list of absolute indices.
 *    really speed up evaluate_binding_pair
 *  - Write species tracker (tracks connectivity)
 *  - Compress reflect_traj_complex_rad_rot, reflect_traj_check_span, reflect_traj_rad_rot_nocheck
 */
#include "boundary_conditions/reflect_functions.hpp"
#include "io/io.hpp"
#include "math/constants.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "mpi.h"
#include "parser/parser_functions.hpp"
#include "reactions/association/association.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "system_setup/system_setup.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

#include <chrono>
#include <cstring>
#include <iomanip>
#include <random>
#include <sstream>

using namespace std;

gsl_rng* r; /* global generator */

// easier to read timer
using MDTimer = std::chrono::system_clock;
using timeDuration = std::chrono::duration<double, std::chrono::seconds>;

/* INITIALIZE GLOBALS */
long long randNum = 0;
unsigned long totMatches = 0;

int main(int argc, char* argv[])
{
    MPI::Init(argc, argv);
    int rank, nprocs;

    rank = MPI::COMM_WORLD.Get_rank();
    nprocs = MPI::COMM_WORLD.Get_size();
    /* SIMULATION SETUP */
    // Get seed for random number generation
    // use random_device instead of time so that multiple jobs started at the same time have more unique seeds
    std::random_device rd {};
    unsigned seed { rd() };

    // get time point for total simulation time
    MDTimer::time_point totalTimeStart = MDTimer::now();

    /* SET UP SIMULATION LISTS */
    // Simulation species lists
    std::vector<MolTemplate> molTemplateList {}; // list of provided molecule templates

    // Reaction lists
    std::vector<ForwardRxn> forwardRxns {}; // list of forward reactions
    std::vector<BackRxn> backRxns {}; // list of back reactions (corresponding to forward reactions)
    std::vector<CreateDestructRxn> createDestructRxns {}; // list of creation and destruction reactions
    forwardRxns.reserve(10);
    backRxns.reserve(10);

    /* PARSE INPUT */
    Parameters params {};
    std::string paramFile {};

    // set up some output files
    // TODO: change these to open and close as needed
    std::string observablesFileName { "observables_time.dat" };
    std::string trajFileName { "trajectory.xyz" };
    std::string transitionFileName { "transition_matrix_time.dat" };
    std::string restartFileName { "restart.dat" };
    std::string addFileNameInput {}; // this is for restart with changed params or adding molecules and reactions
    std::string restartFileNameInput; //if you read it in, allow it to have its own name.
    params.rank = -1; //for serial jobs, this impacts the name of the restart file.

    // command line flag parser
    parse_command(argc, argv, params, paramFile, restartFileNameInput, addFileNameInput, seed);

    auto startTime = MDTimer::to_time_t(totalTimeStart);
    char charTime[24];
    std::cout << "\nStart date: ";
    if (0 < strftime(charTime, sizeof(charTime), "%F %T", std::localtime(&startTime)))
        std::cout << charTime << '\n';
    std::cout << "RNG Seed: " << seed << std::endl;

    //random generator
    const gsl_rng_type* T;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);

    /* SET UP SOME IMPORTANT VARIABLES */
    // 2D reaction probability tables
    std::vector<gsl_matrix*> survMatrices; // used in evaluate_binding_pair_com
    std::vector<gsl_matrix*> normMatrices; // idem
    std::vector<gsl_matrix*> pirMatrices; // idem
    double* tableIDs = new double[params.max2DRxns * 2]; // TODO: Change this?

    /* SET UP SYSTEM */
    std::map<std::string, int> observablesList;
    SimulVolume simulVolume {};
    std::vector<Molecule> moleculeList {}; // list of all molecules in the system
    std::vector<Complex> complexList {}; // list of all complexes in the system
    Membrane membraneObject; //class structure that contains boundary conditions, implicit lipid model.
    copyCounters counterArrays; // contains arrays tracking bound molecule pairs, and species copy nums.
    int implicitlipidIndex { 0 }; // implicit-lipid index, which is also stored in membraneObject.implicitlipidIndex.

    long long int simItr { 0 };

    // some variables used for parsing add.inp

    int numMolTemplateBeforeAdd { 0 }; // number of molTemp before add
    int numDoubleBeforeAdd { 0 }; // number of double complex before add
    int numForwardRxnBdeforeAdd { 0 }; //number of Forw React before add
    int numBackRxnBdeforeAdd { 0 }; //number of Back React before add
    int numCreatDestructRxnBdeforeAdd { 0 }; //num of creat and destruct react before add
    int tempLastStateIndexBeforeAdd { 0 }; // the last state index before add
    int tempLastStateIndexAfterAdd { 0 }; // the last state index after add
    int numStateAdd { 0 }; // num of added states
    int totalSpeciesNum { 0 }; //total species num after add
    init_association_events(counterArrays); //initialize event counters to zero. Restart will update any non-zero values.

    std::cout << "\nParsing Input: " << std::endl;
    if (!params.fromRestart && paramFile != "") {
        std::cout << "This is a new simulation with input file: " << paramFile << std::endl;
        // Parse the input files
        parse_input(paramFile, params, observablesList, forwardRxns, backRxns, createDestructRxns, molTemplateList, membraneObject);

        // set the number of states of implicit lipid
        for (auto& molTemplateTmp : molTemplateList) {
            if (molTemplateTmp.isImplicitLipid == true) {
                membraneObject.nStates = static_cast<int>(molTemplateTmp.interfaceList[0].stateList.size());
                break;
            }
        }
        // initialize numberOfFreeLipidsEachState & numberOfProteinEachState
        for (int tmpStateIndex = 0; tmpStateIndex < membraneObject.nStates; tmpStateIndex++) {
            membraneObject.numberOfFreeLipidsEachState.emplace_back(0);
            membraneObject.numberOfProteinEachState.emplace_back(0);
        }

        // verify the implicit lipid is the first
        for (auto& tempMolTemplate : molTemplateList) {
            if (tempMolTemplate.isImplicitLipid == true && tempMolTemplate.molTypeIndex != 0) {
                std::cerr << "Error: implicit Lipid must be the first molecule type, exiting.\n";
                exit(1);
            }
        }

        MolTemplate::numMolTypes = molTemplateList.size();
        std::cout << "NUMBER OF MOLECULE TYPES: " << params.numMolTypes
                  << "NUMBER OF INTERFACES PLUS STATES, including PRODUCTS: " << params.numTotalSpecies << std::endl;

        // write the Observables file header and initial values
        std::ofstream observablesFile { observablesFileName };
        if (observablesList.size() == 1) {
            observablesFile << "Time (s)," << observablesList.begin()->first << '\n';
            observablesFile << "0,0\n";
        } else if (observablesList.size() > 1) {
            observablesFile << "Time (s)";
            for (auto obsItr = observablesList.begin(); obsItr != observablesList.end(); ++obsItr)
                observablesFile << ',' << obsItr->first;
            std::cout << "\n0";
            for (auto obsItr = observablesList.begin(); obsItr != observablesList.end(); ++obsItr)
                std::cout << ',' << obsItr->second;
            observablesFile << '\n'
                            << std::flush;
        }

        observablesFile.close();

        unsigned long reservation {};
        for (auto& molTemp : molTemplateList)
            reservation += molTemp.copies;
        moleculeList.reserve(reservation);
        complexList.reserve(reservation);

        //create water box for sphere boundary
        if (membraneObject.isSphere) {
            membraneObject.create_water_box();
            membraneObject.sphereVol = (4.0 * M_PI * pow(membraneObject.sphereR, 3.0)) / 3.0;
        }

        // generate the system coordinates, write out coordinate and topology files
        generate_coordinates(params, moleculeList, complexList, molTemplateList, forwardRxns, membraneObject);
        write_psf(params, moleculeList, molTemplateList);

        // set up some important parameters for implicit-lipid model;
        initialize_paramters_for_implicitlipid_model(implicitlipidIndex, params, forwardRxns, backRxns,
            moleculeList, molTemplateList, complexList, membraneObject);

        // initialize the starting copy number for each state
        initialize_states(moleculeList, molTemplateList, membraneObject);

        /* CREATE SIMULATION BOX CELLS */
        std::cout << "\nPartitioning simulation box into sub-boxes..." << std::endl;
        set_rMaxLimit(params, molTemplateList, forwardRxns, 0, 0);
        simulVolume.create_simulation_volume(params, membraneObject);
        simulVolume.update_memberMolLists(params, moleculeList, complexList, molTemplateList, membraneObject, simItr);
        simulVolume.display();

        // write beginning of trajectory
        std::ofstream trajFile { trajFileName };
        write_traj(0, trajFile, params, moleculeList, molTemplateList, membraneObject);
        trajFile.close();

        // initialize transition matrix for each molType
        for (auto& molTemp : molTemplateList) {
            if (molTemp.countTransition == true) {
                molTemp.transitionMatrix.resize(molTemp.transitionMatrixSize);
                molTemp.lifeTime.resize(molTemp.transitionMatrixSize);
                for (int indexOne = 0; indexOne < molTemp.transitionMatrixSize; ++indexOne) {
                    molTemp.transitionMatrix[indexOne].resize(molTemp.transitionMatrixSize);
                }
            }

            // print transition matrix
            // if(molTemp.countTransition == true){
            //     std::cout << molTemp.molName << std::endl;
            //     for (int indexOne = 0; indexOne < molTemp.transitionMatrixSize; ++indexOne) {
            //         for (int indexTwo = 0; indexTwo < molTemp.transitionMatrixSize; ++indexTwo) {
            //             std::cout <<' '<< molTemp.transitionMatrix[indexOne][indexTwo];
            //         }
            //         std::cout << std::endl;
            //     }
            // }
        }

        // write beginning of transition matrix
        std::ofstream transitionFile { transitionFileName };
        write_transition(0, transitionFile, molTemplateList);
        transitionFile.close();
    } else if (params.fromRestart) { // && paramFile.empty()) {
        std::cout << "This is a restart simulation with restart file: " << restartFileNameInput << std::endl;
        read_rng_state(); // read the current RNG state
        std::ifstream restartFileInput { restartFileNameInput };
        if (!restartFileInput) {
            std::cerr << "Error, could not find restart file, exiting...\n";
            exit(1);
        }

        std::cout << "Reading restart file..." << std::endl;
        read_restart(simItr, restartFileInput, params, simulVolume, moleculeList, complexList, molTemplateList, forwardRxns,
            backRxns, createDestructRxns, observablesList, membraneObject, counterArrays);
        restartFileInput.close();

        // initialize numberOfProteinEachState
        for (int tmpStateIndex = 0; tmpStateIndex < membraneObject.nStates; tmpStateIndex++) {
            membraneObject.numberOfProteinEachState.emplace_back(0);
        }

        //add moldecules and reactions, modify parms according to add.inp
        if (addFileNameInput != "") {
            std::cout << "This is a restart simulation with add file: " << addFileNameInput << std::endl;
            numMolTemplateBeforeAdd = molTemplateList.size();
            numForwardRxnBdeforeAdd = static_cast<int>(forwardRxns.size());
            numBackRxnBdeforeAdd = static_cast<int>(backRxns.size());
            numCreatDestructRxnBdeforeAdd = static_cast<int>(createDestructRxns.size());

            if (numMolTemplateBeforeAdd > 0) {
                for (int forwardRxnIndex { 0 }; forwardRxnIndex < numForwardRxnBdeforeAdd; forwardRxnIndex++) {
                    ForwardRxn oneRxn {};
                    oneRxn = forwardRxns[forwardRxnIndex];
                    if (oneRxn.rxnType == ReactionType::bimolecular) {
                        numDoubleBeforeAdd++;
                    }
                }
            }

            tempLastStateIndexBeforeAdd = molTemplateList.back().interfaceList.back().stateList.back().index;

            parse_input_for_add(addFileNameInput, params, observablesList, forwardRxns, backRxns, createDestructRxns, molTemplateList, membraneObject, numDoubleBeforeAdd);
            //move the implicit lipid to the first, and unpdate mol.molTypeIndex
            for (auto& tempMolTemplate : molTemplateList) {
                if (tempMolTemplate.isImplicitLipid == true && tempMolTemplate.molTypeIndex != 0) {
                    std::cout << "Implicit Lipid must be the first molecule type!" << std::endl;
                    exit(1);
                }
            }

            tempLastStateIndexAfterAdd = molTemplateList.back().interfaceList.back().stateList.back().index;

            numStateAdd = tempLastStateIndexAfterAdd - tempLastStateIndexBeforeAdd;

            MolTemplate::numMolTypes = molTemplateList.size();

            unsigned long reservation {};
            for (auto& molTemp : molTemplateList)
                reservation += molTemp.copies;
            moleculeList.reserve(reservation);
            complexList.reserve(reservation);

            //update the iface.index for all molecules that the index is larger than tempLastStateIndexBeforeAdd
            for (auto& tempMol : moleculeList) {
                for (auto& tempIface : tempMol.interfaceList) {
                    if (tempIface.index > tempLastStateIndexBeforeAdd)
                        tempIface.index += numStateAdd;
                }
            }

            //update the index of reacts and products of forward reactions
            for (int forwardReactIndex { 0 }; forwardReactIndex < forwardRxns.size(); forwardReactIndex++) {
                if (forwardReactIndex < numForwardRxnBdeforeAdd) {
                    if (forwardRxns[forwardReactIndex].rxnType == ReactionType::bimolecular) {
                        //this is a old forward biomolecule react, update the products' index by + numStateAdd
                        for (auto& tempProduct : forwardRxns[forwardReactIndex].productListNew) {
                            tempProduct.absIfaceIndex += numStateAdd;
                        }
                    }
                }
                if (forwardReactIndex >= numForwardRxnBdeforeAdd) {
                    if (forwardRxns[forwardReactIndex].rxnType == ReactionType::bimolecular) {
                        //this is a new forward biomolecule react, update the products' index by + numDoubleBeforeAdd
                        for (auto& tempProduct : forwardRxns[forwardReactIndex].productListNew) {
                            tempProduct.absIfaceIndex += numDoubleBeforeAdd;
                        }
                    }
                }
            }

            //update the index of reacts and products of back reactions
            for (int backReactIndex { 0 }; backReactIndex < backRxns.size(); backReactIndex++) {
                if (backReactIndex < numBackRxnBdeforeAdd) {
                    //this is a old back react, update the reactants' index by + numStateAdd
                    for (auto& tempReactant : backRxns[backReactIndex].reactantListNew) {
                        tempReactant.absIfaceIndex += numStateAdd;
                    }
                } else {
                    //this is a new back react, update the reactants' index by + numDoubleBeforeAdd
                    for (auto& tempReactant : backRxns[backReactIndex].reactantListNew) {
                        tempReactant.absIfaceIndex += numDoubleBeforeAdd;
                    }
                }
            }

            // generate the coordinates, write out coordinate and topology files for added molecules
            generate_coordinates_for_restart(params, moleculeList, complexList, molTemplateList, forwardRxns, membraneObject, numMolTemplateBeforeAdd, numForwardRxnBdeforeAdd);
            for (auto& tmpComplex : complexList) {
                tmpComplex.numEachMol.clear();
                tmpComplex.numEachMol.resize(molTemplateList.size());
                for (auto& memMol : tmpComplex.memberList)
                    ++tmpComplex.numEachMol[moleculeList[memMol].molTypeIndex];

                tmpComplex.lastNumberUpdateItrEachMol.resize(molTemplateList.size());
            }

            write_psf(params, moleculeList, molTemplateList);
        }

        //create water box for sphere boundary
        if (membraneObject.isSphere) {
            membraneObject.create_water_box();
            membraneObject.sphereVol = (4.0 * M_PI * pow(membraneObject.sphereR, 3.0)) / 3.0;
        }

        std::cout << " Total number of states (reactant and product)  in the system " << RxnBase::totRxnSpecies << std::endl;
        params.numTotalSpecies = RxnBase::totRxnSpecies;
        std::cout << " Total number of molecules: " << Molecule::numberOfMolecules << " Size of molecule list : " << moleculeList.size() << std::endl;
        std::cout << "Total number of complexes: " << Complex::numberOfComplexes << " size of list: " << complexList.size() << std::endl;

        // set up some important parameters for implicit-lipid model;
        initialize_paramters_for_implicitlipid_model(implicitlipidIndex, params, forwardRxns, backRxns,
            moleculeList, molTemplateList, complexList, membraneObject);

        // initialize the starting copy number for each state
        initialize_states(moleculeList, molTemplateList, membraneObject);

        /* CREATE SIMULATION BOX CELLS */
        std::cout << "Partitioning simulation box into sub-boxes..." << std::endl;
        set_rMaxLimit(params, molTemplateList, forwardRxns, numDoubleBeforeAdd, numMolTemplateBeforeAdd);
        simulVolume.create_simulation_volume(params, membraneObject);
        simulVolume.update_memberMolLists(params, moleculeList, complexList, molTemplateList, membraneObject, simItr);
        simulVolume.display();

        // Check to make sure the trajectory length matches the restart file
        std::cout << " params.trajFile: " << params.trajFile << std::endl;
        std::ifstream trajFile { params.trajFile };
        long long int trajItr { -1 };
        if (trajFile) {
            std::string line;
            while (getline(trajFile, line)) {
                auto headerItr = line.find(':');
                if (headerItr != std::string::npos) {
                    trajItr = std::stoi(line.substr(headerItr + 1, std::string::npos)); // + 1 to ignore the colon
                }
            }
            if (trajItr == simItr) {
                std::cout << "Trajectory length matches provided restart file. Continuing...\n";
            } else {
                std::cerr << "ERROR: Trajectory length doesn't match provided restart file. Exiting...\n";
                exit(1);
            }
            trajFile.close();
        } else {
            std::cout << "WARNING: No trajectory found, writing new trajectory.\n";
        }
    } else {
        std::cerr << "Please provide a parameter and/or restart file. Parameter file Syntax is : ./rd_executable.exe "
                     "-f parameterfile.inp \n";
        exit(1);
    }

    if (membraneObject.implicitLipid == true)
        params.implicitLipid = true; //Created this parameter for convenience.
    for (auto& oneReaction : forwardRxns) {
        if (oneReaction.rxnType == ReactionType::uniMolStateChange) {
            params.hasUniMolStateChange = true;
            break;
        }
    }
    if (createDestructRxns.empty() == false) {
        params.hasCreationDestruction = true;
        params.isNonEQ = true;

        // set molTemp.canDestroy = ture
        for (auto& oneReaction : createDestructRxns) {
            if (oneReaction.rxnType == ReactionType::destruction) {
                molTemplateList[oneReaction.reactantMolList.at(0).molTypeIndex].canDestroy = true;
            }
        }
        for (auto& oneTemp : molTemplateList) {
            if (oneTemp.canDestroy == false) {
                oneTemp.monomerList.clear();
            }
        }
    }

    /* SETUP OUTPUT FILES */
    /*output files reporting bound pairs, and histogram of complex components*/

    // the variable params.numIfaces does not seem to be used anywhere for anything. It is also not correctly
    // initialized anywhere. during a restart, it is read in, but previous sims would not have set to proper value.

    char fnameProXYZ[100];
    sprintf(fnameProXYZ, "histogram_complexes_time.dat");
    std::ofstream assemblyfile(fnameProXYZ);
    sprintf(fnameProXYZ, "mono_dimer_time.dat");
    std::ofstream dimerfile(fnameProXYZ);
    sprintf(fnameProXYZ, "event_counters_time.dat");
    std::ofstream eventFile(fnameProXYZ);
    sprintf(fnameProXYZ, "bound_pair_time.dat");
    std::ofstream pairOutfile(fnameProXYZ);
    sprintf(fnameProXYZ, "copy_numbers_time.dat");
    std::ofstream speciesFile1(fnameProXYZ);

    int meanComplexSize { 0 };

    totalSpeciesNum = init_speciesFile(speciesFile1, counterArrays, molTemplateList, forwardRxns, params);
    init_counterCopyNums(counterArrays, moleculeList, complexList, molTemplateList, membraneObject, totalSpeciesNum, params); // works for default and restart

    write_all_species((simItr - params.itrRestartFrom) * params.timeStep * Constants::usToSeconds + params.timeRestartFrom, speciesFile1, counterArrays);

    init_print_dimers(dimerfile, params, molTemplateList); // works for default and restart
    init_NboundPairs(counterArrays, pairOutfile, params, molTemplateList, moleculeList); // initializes to zero, re-calculated for a restart!!
    write_NboundPairs(counterArrays, pairOutfile, simItr, params, moleculeList);
    print_dimers(complexList, dimerfile, simItr, params, molTemplateList);
    print_association_events(counterArrays, eventFile, simItr, params);
    //this will be wrong if there are no implicit lipids.
    //const int ILcopyIndex = moleculeList[implicitlipidIndex].interfaceList[0].index;

    int number_of_lipids = 0; //sum of all states of IL
    for (int i = 0; i < membraneObject.numberOfFreeLipidsEachState.size(); i++) {
        number_of_lipids += membraneObject.numberOfFreeLipidsEachState[i];
    }
    meanComplexSize = print_complex_hist(complexList, assemblyfile, simItr, params, molTemplateList, number_of_lipids);

    //set some parameters
    if (params.checkPoint == -1) {
        params.checkPoint = params.nItr / 10;
    }
    if (params.transitionWrite == -1) {
        params.transitionWrite = params.nItr / 10;
    }

    // set the excludeVolumeBoundList for each molTemplate according to the reaction list
    for (auto& oneReaction : forwardRxns) {
        if (oneReaction.excludeVolumeBound == true) {
            // this reaction declare that excludeVolumeBound, now add the partner's index to the excludeVolumeBoundList
            int react1MolIndex = oneReaction.reactantListNew[0].molTypeIndex;
            int react1InfIndex = oneReaction.reactantListNew[0].relIfaceIndex;
            int react1InfAbsIndex = oneReaction.reactantListNew[0].absIfaceIndex;
            int react2MolIndex = oneReaction.reactantListNew[1].molTypeIndex;
            int react2InfIndex = oneReaction.reactantListNew[1].relIfaceIndex;
            int react2InfAbsIndex = oneReaction.reactantListNew[1].absIfaceIndex;
            molTemplateList[react1MolIndex].interfaceList[react1InfIndex].excludeVolumeBoundList.push_back(react2MolIndex);
            molTemplateList[react2MolIndex].interfaceList[react2InfIndex].excludeVolumeBoundList.push_back(react1MolIndex);
            molTemplateList[react1MolIndex].interfaceList[react1InfIndex].excludeVolumeBoundIfaceList.push_back(react2InfIndex);
            molTemplateList[react2MolIndex].interfaceList[react2InfIndex].excludeVolumeBoundIfaceList.push_back(react1InfIndex);
            molTemplateList[react1MolIndex].interfaceList[react1InfIndex].excludeRadiusList.push_back(oneReaction.bindRadius);
            molTemplateList[react2MolIndex].interfaceList[react2InfIndex].excludeRadiusList.push_back(oneReaction.bindRadius);
            molTemplateList[react1MolIndex].interfaceList[react1InfIndex].excludeVolumeBoundReactList.push_back(oneReaction.relRxnIndex);
            molTemplateList[react2MolIndex].interfaceList[react2InfIndex].excludeVolumeBoundReactList.push_back(oneReaction.relRxnIndex);
            molTemplateList[react1MolIndex].excludeVolumeBound = true;
            molTemplateList[react2MolIndex].excludeVolumeBound = true;
        }
    }

    for (auto& oneComplex : complexList) {
        oneComplex.update_properties(moleculeList, molTemplateList);
    }

    Parameters::dt = params.timeStep;
    Parameters::lastUpdateTransition.resize(molTemplateList.size());

    /*Print out system information*/
    std::cout << "\nSimulation Parameters\n";
    params.display();
    membraneObject.display();
    std::cout << "\nMolecule Information\n";
    display_all_MolTemplates(molTemplateList);
    std::cout << "\nReactions\n";
    display_all_reactions(forwardRxns, backRxns, createDestructRxns);

    std::cout << "*************** BEGIN SIMULATION **************** " << std::endl;

    // begin the timer
    MDTimer::time_point simulTimeStart = MDTimer::now();
    int numSavedDurations { 1000 };
    std::vector<std::chrono::duration<double>> durationList(numSavedDurations);
    std::fill(durationList.begin(), durationList.end(), std::chrono::duration<double>(simulTimeStart - totalTimeStart));

    unsigned DDTableIndex { 0 };
    /*Vectors to store binding probabilities for implicit lipids in 2D.*/
    std::vector<double> IL2DbindingVec {};
    std::vector<double> IL2DUnbindingVec {};
    std::vector<double> ILTableIDs {};

    if (params.fromRestart == true) {
        read_rng_state();
        // auto endTime = MDTimer::now();
        // auto endTimeFormat = MDTimer::to_time_t(endTime);
        std::ofstream restartFile { restartFileName, std::ios::out }; // to show different from append
        // std::cout << "Writing restart file at iteration " << simItr << " ";
        // if (0 < strftime(charTime, sizeof(charTime), "%F %T", std::localtime(&endTimeFormat)))
        // std::cout << charTime << '\n';
        write_rng_state(); // write the current RNG state
        write_restart(simItr, restartFile, params, simulVolume, moleculeList, complexList, molTemplateList,
            forwardRxns, backRxns, createDestructRxns, observablesList, membraneObject, counterArrays);
        restartFile.close();
    }

    for (simItr += 1; simItr < params.nItr; ++simItr) {
        // std::cout << "simItr: " << simItr << std::endl;
        propCalled = 0;
        MDTimer::time_point startStep = MDTimer::now();

        // destruct, unimol create, and dissociation (explicit) based on population
        check_for_unimolecular_reactions_population(simItr, params, moleculeList, complexList,
            simulVolume, forwardRxns, backRxns, createDestructRxns, molTemplateList, observablesList, counterArrays,
            membraneObject);
        // NOTE: We do not need to parallize this routine. In the above routine, Destruction reaction will change the complexList, moleculeList,
        //       counterArrays.copyNumSpecies, simulVolume.subCellList[].memberMolList,
        //       obervablesList, molTemplateList[].monomerList, membraneObject.numberOfFreeLipidsEachState;
        //       uniMolCreation (A->A+B) reaction will change Molecule::emptyMolList, Complex::emptyComList,
        //       moleculeList, complexList, Complexx::numberOfComplexes, molTemplateList[].monomerList, Molecule::numberOfMolecules,
        //       params.numTotalUnits, MolTemplate::numEachMolType, molTemplateList[].monomerList,
        //       counterArrays.copyNumSpecies, obervablesList
        //       Explicit bimolecular dissociation reaction (AB->A+B) will change counterArrays.bindPairList,
        //       Complex::emptyComList, complexList, moleculeList, molTemplateList[].transitionMatrix,
        //       Parameters::lastUpdateTransition, Complex::numberOfComplexes, molTemplateList[].monomerList,
        //       counterArrays.nLoops, counterArrays.copyNumSpecies, obervablesList

        // Update member lists after creation and destruction
        simulVolume.update_memberMolLists(params, moleculeList, complexList, molTemplateList, membraneObject, simItr);
        // NOTE: moleculeList[].mySubVolIndex and simulVolume.subCellList[].memberMolList are changed

        // Zeroth order reactions (creation)
        check_for_zeroth_order_creation(simItr, params, simulVolume, forwardRxns,
            createDestructRxns, moleculeList, complexList, molTemplateList, observablesList, counterArrays, membraneObject);
        // NOTE: We do not need to parallize the above routine

        // Update member lists after creation and destruction
        simulVolume.update_memberMolLists(params, moleculeList, complexList, molTemplateList, membraneObject, simItr);
        // NOTE: moleculeList[].mySubVolIndex and simulVolume.subCellList[].memberMolList are changed

        // check for unimol state change reactions
        if (params.hasUniMolStateChange == true) {
            check_for_unimolstatechange_reactions(simItr, params, moleculeList, complexList,
                simulVolume, forwardRxns, backRxns, createDestructRxns, molTemplateList, observablesList, counterArrays, membraneObject);
            // NOTE: there is a loop over the moleculeList in the above routine. We can divide the moleculeList to pieces to parallize
        }

        /*Skip this entire loop if the system has no implicit lipids. */
        if (params.implicitLipid == true) {
            // check dissociation (implicit)
            for (unsigned molItr { 0 }; molItr < moleculeList.size(); ++molItr) {
                // only do checks if the Molecule exists
                if (moleculeList[molItr].isEmpty || moleculeList[molItr].isImplicitLipid == true || complexList[moleculeList[molItr].myComIndex].OnSurface == false || params.implicitLipid == false)
                    continue;
                check_dissociation_implicitlipid(simItr, params, simulVolume, molTemplateList, observablesList, molItr, moleculeList, complexList, backRxns, forwardRxns, createDestructRxns, counterArrays, membraneObject, IL2DbindingVec, IL2DUnbindingVec, ILTableIDs);
            }
            // NOTE: We can divide the moleculeList in the above for loop into pieces to parallize;
            //       The check_dissociation_implicitlipid routine will change moleculeList, complexList,
            //       counterArrays.copyNumSpecies, membraneObject.numberOfFreeLipidEachState, observablesList,
            //       counterArrays.nBoundPairs, molTemplateList[].monomerList
        }

        //*********************************Start Measure separations between proteins in neighboring cells to identify all possible reactions*********************************************************************************************************************
        // NOTE: this part can be parallized based on the division of the simulVolume.subCellList
        //       All the input variables in the check_*_reactions routine and the simulVolume are needed to send to the ranks
        //       The changed variables are ILTableIDs, IL2DbindingVec, moleculeList, complexList, tableIDs, pirMatrices, survMatrices, normMatrices
        for (unsigned cellItr { 0 }; cellItr < simulVolume.subCellList.size(); ++cellItr) {
            for (unsigned memItr { 0 }; memItr < simulVolume.subCellList[cellItr].memberMolList.size(); ++memItr) {
                int targMolIndex { simulVolume.subCellList[cellItr].memberMolList[memItr] };
                if (moleculeList[targMolIndex].isImplicitLipid)
                    continue;

                //Test bimolecular reactions, and binding to implicit-lipids
                if (moleculeList[targMolIndex].freelist.size() > 0 || molTemplateList[moleculeList[targMolIndex].molTypeIndex].excludeVolumeBound == true) {
                    // first, check for implicit-lipid binding
                    int protype = moleculeList[targMolIndex].molTypeIndex;
                    if (molTemplateList[protype].bindToSurface == true) {
                        check_implicit_reactions(targMolIndex, implicitlipidIndex, simItr, params, moleculeList, complexList, molTemplateList,
                            forwardRxns, backRxns, counterArrays, membraneObject, IL2DbindingVec, IL2DUnbindingVec, ILTableIDs);
                    }
                    // secondly, loop over proteins in your same cell.
                    for (unsigned memItr2 { memItr + 1 }; memItr2 < simulVolume.subCellList[cellItr].memberMolList.size(); ++memItr2) {
                        int partMolIndex { simulVolume.subCellList[cellItr].memberMolList[memItr2] };
                        check_bimolecular_reactions(targMolIndex, partMolIndex, simItr, tableIDs, DDTableIndex, params,
                            normMatrices, survMatrices, pirMatrices, moleculeList, complexList, molTemplateList,
                            forwardRxns, backRxns, counterArrays, membraneObject);
                    } // loop over protein partners in your same cell
                    // thirdly, loop over all neighboring cells, and all proteins in those cells.
                    // for PBC, all cells have maxnbor neighbor cells. For reflecting, edge have fewer.
                    for (auto& neighCellItr : simulVolume.subCellList[cellItr].neighborList) {
                        for (unsigned memItr2 { 0 }; memItr2 < simulVolume.subCellList[neighCellItr].memberMolList.size(); ++memItr2) {
                            int partMolIndex { simulVolume.subCellList[neighCellItr].memberMolList[memItr2] };
                            check_bimolecular_reactions(targMolIndex, partMolIndex, simItr, tableIDs, DDTableIndex, params,
                                normMatrices, survMatrices, pirMatrices, moleculeList, complexList, molTemplateList,
                                forwardRxns, backRxns, counterArrays, membraneObject);
                        } // loop over all proteins in this neighbor cell
                    } // loop over all neighbor cells
                } // if protein i is free to bind
            } // loop over all proteins in initial cell
        } // End looping over all cells.
        // NOTE: this part can be parallized based on the division of the simulVolume.subCellList
        //*********************************End Measure separations between proteins in neighboring cells to identify all possible reactions*********************************************************************************************************************

        //**********************Start decide whether to perform reactions for each protein****************************
        // NOTE: This part can be divided based on the moleculeList to parallize
        //       moleculeList[].trajStatus need to updated immediately for all ranks when it is changed
        //       moleculeList, complexList, counterArrays, molTemplateList, membraneObject are needed to send from other ranks to rank 0
        /*Now that separations and reaction probabilities are calculated, decide whether to perform reactions for each protein.*/
        for (int molItr { 0 }; molItr < moleculeList.size(); ++molItr) {
            // only continue if the molecule actually exists, and isn't implicit-lipid
            if (moleculeList[molItr].isEmpty || moleculeList[molItr].isImplicitLipid)
                continue;

            //Skip any proteins that just dissociated during this time step
            if (moleculeList[molItr].crossbase.size() > 0) {
                /* Evaluate whether to perform a reaction with protein i, and with whom. Flag=1 means
                 * reaction is performed. Returns correct ci1 and ci2 for this rxn.
                 * Loop over all reactions individually, instead of summing probabilities
                 */
                // these are indices in crossbase/mycrossint/crossrxn of the reaction for molecules 1 and 2,
                // should it occur
                int crossIndex1 { 0 };
                int crossIndex2 { 0 };
                bool willReact { determine_if_reaction_occurs(crossIndex1, crossIndex2, Constants::iRandMax, moleculeList[molItr], moleculeList, forwardRxns) };
                if (params.debugParams.forceAssoc) {
                    willReact = false; // we chose an association reaction
                    crossIndex1 = 0;
                    int molItr2 = moleculeList[molItr].crossbase[crossIndex1]; // crosspart[p1][ci1];
                    double pmatch = moleculeList[molItr].probvec[crossIndex1];
                    if (pmatch > 0)
                        willReact = true;
                    if (!moleculeList[molItr].isImplicitLipid) {
                        for (unsigned j = 0; j < moleculeList[molItr2].crossbase.size(); ++j) {
                            if (moleculeList[molItr2].probvec[j] == pmatch) {
                                crossIndex2 = j;
                                if (moleculeList[molItr].crossrxn[crossIndex1]
                                    == moleculeList[molItr2].crossrxn[crossIndex2])
                                    break;
                            }
                        }
                    } else {
                        crossIndex2 = 0;
                    }
                }
                if (willReact) {
                    /* This molecule will perform a bimolecular reaction
                     * either physically associate two molecules into a complex (A+B->AB)
                     * or change state of one (or both) reactants (A+B->A+B')
                     */
                    int molItr2 { moleculeList[molItr].crossbase[crossIndex1] };
                    int ifaceIndex1 { moleculeList[molItr].mycrossint[crossIndex1] };
                    int ifaceIndex2;
                    if (moleculeList[molItr2].isImplicitLipid == false) {
                        ifaceIndex2 = moleculeList[molItr2].mycrossint[crossIndex2];
                    } else {
                        ifaceIndex2 = 0;
                    }
                    std::array<int, 3> rxnIndex = moleculeList[molItr].crossrxn[crossIndex1];

                    /*First if statement is to determine if reactants are physically associating*/
                    if (forwardRxns[rxnIndex[0]].rxnType == ReactionType::bimolecular) {
                        if (moleculeList[molItr2].isImplicitLipid) {
                            // std::cout << "Performing binding of molecules " << molItr << " to the membrane surface "
                            //           << " ["
                            //           << molTemplateList[moleculeList[molItr].molTypeIndex].molName << "("
                            //           << molTemplateList[moleculeList[molItr].molTypeIndex].interfaceList[ifaceIndex1].name
                            //           << ")] with [reaction, rate] = [" << rxnIndex[0] << ',' << rxnIndex[1]
                            //           << "] Rate: " << forwardRxns[rxnIndex[0]].rateList[rxnIndex[1]].rate << " at iteration " << simItr << ".\n";
                            // std::cout << " Complex 1 size: " << complexList[moleculeList[molItr].myComIndex].memberList.size() << "\n";
                            // std::cout << "Implicit lipid, my Com Index: " << moleculeList[molItr2].myComIndex << " size of memberlist: " << complexList[moleculeList[molItr2].myComIndex].memberList.size() << std::endl;

                            if (molTemplateList[forwardRxns[rxnIndex[0]].reactantListNew[0].molTypeIndex].isImplicitLipid == false) { //IL is listed second as the reactant.
                                associate_implicitlipid(ifaceIndex1, ifaceIndex2, moleculeList[molItr], moleculeList[molItr2],
                                    complexList[moleculeList[molItr].myComIndex], complexList[moleculeList[molItr2].myComIndex], params, forwardRxns[rxnIndex[0]],
                                    moleculeList, molTemplateList, observablesList, counterArrays, complexList, membraneObject, forwardRxns, backRxns);
                            } else { //IL is listed first as the reactant.
                                associate_implicitlipid(ifaceIndex2, ifaceIndex1, moleculeList[molItr2], moleculeList[molItr],
                                    complexList[moleculeList[molItr2].myComIndex], complexList[moleculeList[molItr].myComIndex], params, forwardRxns[rxnIndex[0]],
                                    moleculeList, molTemplateList, observablesList, counterArrays, complexList, membraneObject, forwardRxns, backRxns);
                            }
                        }
                        if (moleculeList[molItr2].isImplicitLipid == false) {
                            // std::cout << "Performing association between molecules " << molItr << " and " << molItr2 << " ["
                            //           << molTemplateList[moleculeList[molItr].molTypeIndex].molName << "("
                            //           << molTemplateList[moleculeList[molItr].molTypeIndex].interfaceList[ifaceIndex1].name
                            //           << ") and " << molTemplateList[moleculeList[molItr2].molTypeIndex].molName << "("
                            //           << molTemplateList[moleculeList[molItr2].molTypeIndex].interfaceList[ifaceIndex2].name
                            //           << ")] with [reaction, rate] = [" << rxnIndex[0] << ',' << rxnIndex[1]
                            //           << "] Rate: " << forwardRxns[rxnIndex[0]].rateList[rxnIndex[1]].rate << " at iteration " << simItr << ".\n";
                            // std::cout << " Complex 1 size: " << complexList[moleculeList[molItr].myComIndex].memberList.size() << "\n";
                            // std::cout << " Complex 2 size: " << complexList[moleculeList[molItr2].myComIndex].memberList.size() << "\n";

                            // For association, molecules must be read in in the order used to define the reaction parameters.
                            if (moleculeList[molItr].interfaceList[ifaceIndex1].index
                                == forwardRxns[rxnIndex[0]].reactantListNew[0].absIfaceIndex) {
                                associate(simItr, ifaceIndex1, ifaceIndex2, moleculeList[molItr], moleculeList[molItr2],
                                    complexList[moleculeList[molItr].myComIndex], complexList[moleculeList[molItr2].myComIndex], params, forwardRxns[rxnIndex[0]],
                                    moleculeList, molTemplateList, observablesList,
                                    counterArrays, complexList, membraneObject, forwardRxns, backRxns);
                            } else {
                                associate(simItr, ifaceIndex2, ifaceIndex1, moleculeList[molItr2], moleculeList[molItr],
                                    complexList[moleculeList[molItr2].myComIndex], complexList[moleculeList[molItr].myComIndex], params, forwardRxns[rxnIndex[0]],
                                    moleculeList, molTemplateList, observablesList,
                                    counterArrays, complexList, membraneObject, forwardRxns, backRxns);
                            }
                        }
                    } else if (forwardRxns[rxnIndex[0]].rxnType == ReactionType::biMolStateChange) {
                        if (moleculeList[molItr2].isImplicitLipid) {
                            //In this case, one implicit lipid changes state.
                            int facilMolIndex { molItr };
                            int facilIfaceIndex { ifaceIndex1 };
                            int facilComIndex { moleculeList[molItr].myComIndex };
                            int stateMolIndex { molItr2 };
                            int stateComIndex { moleculeList[molItr2].myComIndex };
                            int stateIfaceIndex { ifaceIndex2 };
                            // figure out if the current molecule is the facilitator molecule or the IL which has its
                            // interface state changed
                            if (moleculeList[molItr].molTypeIndex
                                != forwardRxns[rxnIndex[0]].reactantListNew[0].molTypeIndex) {
                                facilMolIndex = molItr2;
                                facilIfaceIndex = ifaceIndex2;
                                facilComIndex = moleculeList[molItr2].myComIndex;
                                stateMolIndex = molItr;
                                stateComIndex = moleculeList[molItr].myComIndex;
                                stateIfaceIndex = ifaceIndex1;
                            }
                            // std::cout << "Performing bimolecular state change on molecule (IL) " << stateMolIndex
                            //           << " as facilitated by molecule " << facilMolIndex << " at iteration " << simItr
                            //           << '\n';
                            perform_implicitlipid_state_change(stateIfaceIndex, facilIfaceIndex, rxnIndex,
                                moleculeList[stateMolIndex], moleculeList[facilMolIndex], complexList[stateComIndex],
                                complexList[facilComIndex], counterArrays, params, forwardRxns, backRxns, moleculeList,
                                complexList, molTemplateList, observablesList, membraneObject);
                        }
                        if (moleculeList[molItr2].isImplicitLipid == false) {
                            //In this case, after two molecules collide, at least one of them changes state, rather than forming a complex.
                            int facilMolIndex { molItr };
                            int facilIfaceIndex { ifaceIndex1 };
                            int facilComIndex { moleculeList[molItr].myComIndex };
                            int stateMolIndex { molItr2 };
                            int stateComIndex { moleculeList[molItr2].myComIndex };
                            int stateIfaceIndex { ifaceIndex2 };
                            // figure out if the current molecule is the facilitator molecule or the molecule which has its
                            // interface state changed
                            if (moleculeList[molItr].molTypeIndex
                                != forwardRxns[rxnIndex[0]].reactantListNew[0].molTypeIndex) {
                                facilMolIndex = molItr2;
                                facilIfaceIndex = ifaceIndex2;
                                facilComIndex = moleculeList[molItr2].myComIndex;
                                stateMolIndex = molItr;
                                stateComIndex = moleculeList[molItr].myComIndex;
                                stateIfaceIndex = ifaceIndex1;
                            }

                            // std::cout << "Performing bimolecular state change on molecule " << stateMolIndex
                            //           << " as facilitated by molecule " << facilMolIndex << " at iteration " << simItr
                            //           << '\n';
                            perform_bimolecular_state_change(stateIfaceIndex, facilIfaceIndex, rxnIndex,
                                moleculeList[stateMolIndex], moleculeList[facilMolIndex], complexList[stateComIndex],
                                complexList[facilComIndex], counterArrays, params, forwardRxns, backRxns, moleculeList,
                                complexList, molTemplateList, observablesList, membraneObject);
                        }
                    } else {
                        std::cerr << "ERROR: Attemping bimolecular reaction which has no reaction type. Exiting..\n";
                        exit(1);
                    }
                } else {
                    /*No reaction was chosen for this molecule*/
                    if (moleculeList[molItr].trajStatus == TrajStatus::none) {
                        // Propagate complex with random translational and rotational motion
                        create_complex_propagation_vectors(params, complexList[moleculeList[molItr].myComIndex], moleculeList,
                            complexList, molTemplateList, membraneObject);
                        for (auto& memMol : complexList[moleculeList[molItr].myComIndex].memberList)
                            moleculeList[memMol].trajStatus = TrajStatus::canBeResampled;
                    }
                    // Set probability of this protein to zero in all reactions so it doesn't try to
                    // react again but the partners still will avoid overlapping.
                    for (unsigned crossItr { 0 }; crossItr < moleculeList[molItr].crossbase.size(); ++crossItr) {
                        int skipMol { moleculeList[molItr].crossbase[crossItr] };
                        for (unsigned crossItr2 { 0 }; crossItr2 < moleculeList[skipMol].crossbase.size();
                             ++crossItr2) {
                            if (moleculeList[skipMol].crossbase[crossItr2] == moleculeList[molItr].index)
                                moleculeList[skipMol].probvec[crossItr2] = 0;
                        }
                    }
                }

            } else if (moleculeList[molItr].crossbase.size() == 0) {
                /* this protein has ncross=0,
                 * meaning it neither dissociated nor tried to associate.
                 * however, it could have movestat=2 if it is part of a multi-protein
                 * complex that already displaced.
                 */
                if (moleculeList[molItr].trajStatus == TrajStatus::none) {
                    create_complex_propagation_vectors(params, complexList[moleculeList[molItr].myComIndex], moleculeList,
                        complexList, molTemplateList, membraneObject);
                    for (auto& memMol : complexList[moleculeList[molItr].myComIndex].memberList)
                        moleculeList[memMol].trajStatus = TrajStatus::canBeResampled;
                }
            }
        } // done testing all molecules for bimolecular reactions
        //**********************End decide whether to perform reactions for each protein****************************

        //**********************Start check for overlap**************************************************************
        // NOTE: The following for loop can be divided based on the moleCuleList to parallize
        // Now we have to check for overlap!!!
        for (auto& mol : moleculeList) {
            //Now track each complex (ncrosscom), and test for overlap of all proteins in that complex before
            //performing final position updates.
            if (mol.isEmpty || mol.isImplicitLipid || mol.trajStatus == TrajStatus::propagated)
                continue;

            // determine RS3Dinput
            double RS3Dinput { 0.0 };
            for (int RS3Dindex = 0; RS3Dindex < 100; RS3Dindex++) {
                if (std::abs(membraneObject.RS3Dvect[RS3Dindex + 400] - mol.molTypeIndex) < 1E-2) {
                    RS3Dinput = membraneObject.RS3Dvect[RS3Dindex + 300];
                    break;
                }
            }

            if (complexList[mol.myComIndex].ncross > 0) {
                if (mol.trajStatus == TrajStatus::none || mol.trajStatus == TrajStatus::canBeResampled) {
                    // For any protein that overlapped and did not react, check whether it overlaps with its partners,
                    // do all proteins in the same complex at the same time.
                    // Also, if both proteins are stuck to membrane, only do xy displacement, ignore z
                    // TODO: Maybe do a boundary sphere overlap check first?

                    if (std::abs(complexList[mol.myComIndex].D.z) < 1E-10) {
                        if (params.clusterOverlapCheck == false) {
                            sweep_separation_complex_rot_memtest(
                                simItr, mol.index, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);
                        } else {
                            sweep_separation_complex_rot_memtest_cluster(
                                simItr, mol.index, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);
                        }
                    } else {
                        sweep_separation_complex_rot(
                            simItr, mol.index, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);
                    }
                    if (membraneObject.isSphere == true)
                        reflect_complex_rad_rot(membraneObject, complexList[mol.myComIndex], moleculeList, RS3Dinput);
                }
            } else {
                if (mol.trajStatus == TrajStatus::none || mol.trajStatus == TrajStatus::canBeResampled) {
                    // For proteins with ncross=0, they either moved independently, or their displacements
                    // were selected based on the complex they were part of, and they may not yet been moved.
                    if (membraneObject.isSphere == true) {
                        if (mol.trajStatus == TrajStatus::none) {
                            create_complex_propagation_vectors(params, complexList[mol.myComIndex], moleculeList,
                                complexList, molTemplateList, membraneObject);
                            for (auto& memMol : complexList[mol.myComIndex].memberList)
                                moleculeList[memMol].trajStatus = TrajStatus::canBeResampled;
                        }
                        complexList[mol.myComIndex].propagate(moleculeList, membraneObject, molTemplateList);
                        reflect_complex_rad_rot(membraneObject, complexList[mol.myComIndex], moleculeList, RS3Dinput);
                    } else {
                        // reflect_traj_complex_rad_rot(params, moleculeList, complexList[mol.myComIndex], membraneObject, RS3Dinput);
                        if (mol.trajStatus == TrajStatus::none) {
                            create_complex_propagation_vectors(params, complexList[mol.myComIndex], moleculeList,
                                complexList, molTemplateList, membraneObject);
                            for (auto& memMol : complexList[mol.myComIndex].memberList)
                                moleculeList[memMol].trajStatus = TrajStatus::canBeResampled;
                        }
                        complexList[mol.myComIndex].propagate(moleculeList, membraneObject, molTemplateList);
                    }
                }
            }
        }
        //**********************End check for overlap**************************************************************

        if (simItr % params.trajWrite == 0) {
            // std::cout << "Writing trajectory...\n";
            std::ofstream trajFile { trajFileName, std::ios::app }; // for append
            write_traj(simItr, trajFile, params, moleculeList, molTemplateList, membraneObject);
            trajFile.close();
        }

        if (params.pdbWrite != -1) {
            if (simItr % params.pdbWrite == 0) {
                // std::cout << "Writing PDB file for current frame...\n";
                write_pdb(simItr, simItr, params, moleculeList, molTemplateList, membraneObject);
            }
        }

        if (params.transitionWrite != -1) {
            if (simItr % params.transitionWrite == 0) {
                // std::cout << "Writing transition matrix...\n";
                std::ofstream transitionFile { transitionFileName, std::ios::app }; // for append
                write_transition((simItr - params.itrRestartFrom) * params.timeStep * Constants::usToSeconds + params.timeRestartFrom, transitionFile, molTemplateList);
                transitionFile.close();
            }
        }

        //remove the empty complexes
        //first remove the empty complexes in the tail
        while (complexList.back().isEmpty == true) {
            int tempIndex { complexList.back().index }; // the removed complex's index
            complexList.pop_back();
            //update Complex::emptyComList
            for (auto& tempEmpty : Complex::emptyComList) {
                if (tempEmpty == tempIndex) {
                    tempEmpty = Complex::emptyComList.back();
                    Complex::emptyComList.pop_back();
                    break;
                }
            }
        }

        // put the last non-empty complex in the list to the non-last empty slot
        while (Complex::emptyComList.empty() == false) {
            int slotIndex { Complex::emptyComList.back() };
            int previousIndex { complexList.back().index };

            complexList[slotIndex] = complexList.back();
            complexList[slotIndex].index = slotIndex;
            complexList.pop_back();

            // change the mol.myComIndex with previousIndex to slotIndex
            for (auto mp : complexList[slotIndex].memberList) {
                moleculeList[mp].myComIndex = slotIndex;
            }

            //update Complex::emptyComList
            Complex::emptyComList.pop_back();

            //remove the empty complexes in the tail
            while (complexList.back().isEmpty == true) {
                int tempIndex { complexList.back().index }; // the removed complex's index
                complexList.pop_back();
                //update Complex::emptyComList
                for (auto& tempEmpty : Complex::emptyComList) {
                    if (tempEmpty == tempIndex) {
                        tempEmpty = Complex::emptyComList.back();
                        Complex::emptyComList.pop_back();
                        break;
                    }
                }
            }
        }

        //------------------------------------------------------------------------------------
        //remove the empty molecules
        //first remove the empty molecule in the tail
        while (moleculeList.back().isEmpty == true) {
            int tempIndex { moleculeList.back().index }; // the removed molecule's index
            moleculeList.pop_back();
            //update Molecule::emptyMolList
            for (auto& tempEmpty : Molecule::emptyMolList) {
                if (tempEmpty == tempIndex) {
                    tempEmpty = Molecule::emptyMolList.back();
                    Molecule::emptyMolList.pop_back();
                    break;
                }
            }
        }

        // put the last non-empty molecule in the list to the non-last empty slot
        while (Molecule::emptyMolList.empty() == false) {
            int slotIndex { Molecule::emptyMolList.back() };
            int previousIndex { moleculeList.back().index };

            moleculeList[slotIndex] = moleculeList.back();
            moleculeList[slotIndex].index = slotIndex;
            moleculeList.pop_back();

            // change the mol.index with previousIndex to slotIndex, include complex.memberlist; interface.interaction.partnerIndex;mol.bndpartner
            int tmpComIndex { moleculeList[slotIndex].myComIndex };
            for (auto& tmpMember : complexList[tmpComIndex].memberList) {
                if (tmpMember == previousIndex)
                    tmpMember = slotIndex;
            }

            for (auto& tmpPartner : moleculeList[slotIndex].bndpartner) {
                for (auto& partner : moleculeList[tmpPartner].bndpartner) {
                    if (partner == previousIndex)
                        partner = slotIndex;
                }
                for (auto& tmpIface : moleculeList[tmpPartner].interfaceList) {
                    if (tmpIface.interaction.partnerIndex == previousIndex)
                        tmpIface.interaction.partnerIndex = slotIndex;
                }
            }

            // update the oneTemp.monomerList if it is necessary
            MolTemplate& oneTemp { molTemplateList[moleculeList[slotIndex].molTypeIndex] };
            if (oneTemp.canDestroy == true) {
                std::vector<int>& oneList { oneTemp.monomerList };
                std::vector<int>::iterator result { std::find(std::begin(oneList), std::end(oneList), previousIndex) }; //check whether previousIndex in oneTemp.monomerList
                if (result != std::end(oneList)) {
                    oneList.erase(result);
                    oneList.emplace_back(slotIndex);
                }
            }

            // update the countArrays.bindPairList
            // check whether previousIndex in countArrays.bindPairList
            for (auto& oneSpecie : counterArrays.bindPairList) {
                for (auto& oneIndex : oneSpecie) {
                    if (oneIndex == previousIndex) {
                        oneIndex = slotIndex;
                    }
                }
            }
            // for (auto& oneSpecie : counterArrays.bindPairListIL2D) {
            //     for (auto& oneIndex : oneSpecie) {
            //         if (oneIndex == previousIndex) {
            //             oneIndex = slotIndex;
            //         }
            //     }
            // }
            // for (auto& oneSpecie : counterArrays.bindPairListIL3D) {
            //     for (auto& oneIndex : oneSpecie) {
            //         if (oneIndex == previousIndex) {
            //             oneIndex = slotIndex;
            //         }
            //     }
            // }

            //update Molecule::emptyMolList
            Molecule::emptyMolList.pop_back();

            //remove the empty molecules in the tail
            while (moleculeList.back().isEmpty == true) {
                int tempIndex { moleculeList.back().index }; // the removed molecule's index
                moleculeList.pop_back();
                //update Molecule::emptyMolList
                for (auto& tempEmpty : Molecule::emptyMolList) {
                    if (tempEmpty == tempIndex) {
                        tempEmpty = Molecule::emptyMolList.back();
                        Molecule::emptyMolList.pop_back();
                        break;
                    }
                }
            }
        }
        //------------------------------------------------------------------------------------

        // Clear lists used for reweighting and encounter tracking
        for (auto& oneMol : moleculeList) {
            if (oneMol.isEmpty || oneMol.isImplicitLipid)
                continue;

            clear_reweight_vecs(oneMol);
            oneMol.trajStatus = TrajStatus::none;
            oneMol.crossbase.clear();
            oneMol.mycrossint.clear();
            oneMol.crossrxn.clear();
            oneMol.probvec.clear();

            // update complexes. if the complex has more than one member, it'll be updated
            // more than once, but this is probably more efficient than iterating over the complexes afterwards
            complexList[oneMol.myComIndex].ncross = 0;
            complexList[oneMol.myComIndex].trajStatus = TrajStatus::none;
        }

        //write restart
        if (simItr % params.restartWrite == 0) {
            auto endTime = MDTimer::now();
            auto endTimeFormat = MDTimer::to_time_t(endTime);
            std::ofstream restartFile { restartFileName, std::ios::out }; // to show different from append
            // std::cout << "Writing restart file at iteration " << simItr;
            //                      << ", system time: " << std::put_time(std::localtime(&endTimeFormat), "%F %T") << '\n';
            // if (0 < strftime(charTime, sizeof(charTime), "%F %T", std::localtime(&endTimeFormat)))
            // std::cout << charTime << '\n';
            write_rng_state(); // write the current RNG state
            write_restart(simItr, restartFile, params, simulVolume, moleculeList, complexList, molTemplateList,
                forwardRxns, backRxns, createDestructRxns, observablesList, membraneObject, counterArrays);
            restartFile.close();
        }

        //write check point
        if (simItr % params.checkPoint == 0) {
            sprintf(fnameProXYZ, "restart%lld.dat", simItr);
            std::ofstream restartFile(fnameProXYZ);
            write_rng_state_simItr(simItr); // write the current RNG state
            write_restart(simItr, restartFile, params, simulVolume, moleculeList, complexList, molTemplateList,
                forwardRxns, backRxns, createDestructRxns, observablesList, membraneObject, counterArrays);
            restartFile.close();
        }

        using duration = std::chrono::duration<double>;
        durationList.erase(durationList.begin());
        durationList.emplace_back(MDTimer::now() - startStep);
        if (simItr % params.timeWrite == 0) {
            double timeSimulated { (simItr - params.itrRestartFrom) * params.timeStep * Constants::usToSeconds + params.timeRestartFrom };
            std::cout << linebreak;
            std::cout << "End iteration: " << simItr << ", simulation time: ";
            std::cout << std::scientific << timeSimulated << " seconds.\n";
            // Write out N bound pairs, histogram of complex compositions, and monomer/dimer counts.
            write_NboundPairs(counterArrays, pairOutfile, simItr, params, moleculeList);
            print_dimers(complexList, dimerfile, simItr, params, molTemplateList);
            print_association_events(counterArrays, eventFile, simItr, params);

            int number_of_lipids = 0; //sum of all states of IL
            for (int i = 0; i < membraneObject.numberOfFreeLipidsEachState.size(); i++) {
                number_of_lipids += membraneObject.numberOfFreeLipidsEachState[i];
            }
            meanComplexSize = print_complex_hist(complexList, assemblyfile, simItr, params, molTemplateList, number_of_lipids);
            auto endTime = MDTimer::now();
            auto endTimeFormat = MDTimer::to_time_t(endTime);
            std::cout << "System time: ";
            if (0 < strftime(charTime, sizeof(charTime), "%F %T", std::localtime(&endTimeFormat)))
                std::cout << charTime << '\n';
            std::cout << "Elapsed time: "
                      << std::chrono::duration_cast<std::chrono::minutes>(MDTimer::now() - totalTimeStart).count()
                      << " minutes\n";

            std::cout << "Number of molecules: " << Molecule::numberOfMolecules << '\n';
            std::cout << "Number of complexes: " << Complex::numberOfComplexes << '\n';
            std::cout << "Total reaction matches: " << totMatches << '\n';
            if (params.debugParams.printSystemInfo) {
                std::cout << "Printing full system information...\n";
                std::ofstream systemInfoFile { "system_information.dat", std::ios::app };
                print_system_information(simItr, systemInfoFile, moleculeList, complexList, molTemplateList);
                systemInfoFile.close();
            }
            // write observables
            if (!observablesList.empty()) {
                // std::cout << "Writing observables to file...\n";
                std::ofstream observablesFile { observablesFileName, std::ios::app };
                write_observables((simItr - params.itrRestartFrom) * params.timeStep * Constants::usToSeconds + params.timeRestartFrom, observablesFile, observablesList);
                observablesFile.close();
            }
            // write all species

            //std::ofstream speciesFile{ speciesFileName, std::ios::app };
            write_all_species((simItr - params.itrRestartFrom) * params.timeStep * Constants::usToSeconds + params.timeRestartFrom, speciesFile1, counterArrays);
            //speciesFile.close();

            // Estimate time remaining
            duration avgTimeStepDuration
                = std::accumulate(durationList.begin(), durationList.end(), duration { 0 }) / numSavedDurations;
            duration timeLeft = (params.nItr - simItr) * avgTimeStepDuration;
            std::cout << "Avg timestep duration: " << avgTimeStepDuration.count()
                      << ", iterations remaining: " << params.nItr - simItr
                      << ", Time left: " << std::chrono::duration_cast<std::chrono::minutes>(timeLeft).count()
                      << " minutes\n";
            auto estTimeLeft = std::chrono::time_point_cast<std::chrono::seconds>(MDTimer::now() + timeLeft);
            auto estTimeEnd = std::chrono::system_clock::to_time_t(estTimeLeft);
            std::cout << "Estimated end time: ";
            //<< std::put_time(std::localtime(&estTimeEnd), "%F %T") << '\n';
            if (0 < strftime(charTime, sizeof(charTime), "%F %T", std::localtime(&estTimeEnd)))
                std::cout << charTime << '\n';
            std::cout << llinebreak;
        }
    } // end iterating over time steps

    // Write files at last timestep
    {
        simItr--;
        // std::cout << "Writing restart file at final iteration\n.";
        std::ofstream restartFile { restartFileName, std::ios::out }; // to show different from append
        write_rng_state(); // write the current RNG state
        write_restart(simItr, restartFile, params, simulVolume, moleculeList, complexList, molTemplateList, forwardRxns,
            backRxns, createDestructRxns, observablesList, membraneObject, counterArrays);
        restartFile.close();

        // std::cout << "Writing trajectory..." << '\n';
        std::ofstream trajFile { trajFileName, std::ios::app }; // for append
        write_traj(simItr, trajFile, params, moleculeList, molTemplateList, membraneObject);
        trajFile.close();

        // std::cout << "Writing final configuration...\n";
        write_xyz("final_coords.xyz", params, moleculeList, molTemplateList);

        if (params.pdbWrite != -1) {
            // std::cout << "Writing PDB file for current frame.\n";
            write_pdb(simItr, simItr, params, moleculeList, molTemplateList, membraneObject);
        }

        if (params.transitionWrite != -1) {
            // std::cout << "Writing transition matrix...\n";
            std::ofstream transitionFile { transitionFileName, std::ios::app }; // for append
            write_transition((simItr - params.itrRestartFrom) * params.timeStep * Constants::usToSeconds + params.timeRestartFrom, transitionFile, molTemplateList);
            transitionFile.close();
        }

        if (params.debugParams.printSystemInfo) {
            // std::cout << "Printing full system information...\n";
            std::ofstream systemInfoFile { "system_information.dat", std::ios::app };
            print_system_information(simItr, systemInfoFile, moleculeList, complexList, molTemplateList);
            systemInfoFile.close();
        }

        // write observables
        if (!observablesList.empty()) {
            // std::cout << "Writing observables to file...\n";
            std::ofstream observablesFile { observablesFileName, std::ios::app };
            write_observables((simItr - params.itrRestartFrom) * params.timeStep * Constants::usToSeconds + params.timeRestartFrom, observablesFile, observablesList);
            observablesFile.close();
        }

        // write all species
        write_all_species((simItr - params.itrRestartFrom) * params.timeStep * Constants::usToSeconds + params.timeRestartFrom, speciesFile1, counterArrays);
        // Write out N bound pairs, histogram of complex compositions, and monomer/dimer counts.
        write_NboundPairs(counterArrays, pairOutfile, simItr, params, moleculeList);
        print_dimers(complexList, dimerfile, simItr, params, molTemplateList);
        print_association_events(counterArrays, eventFile, simItr, params);

        int number_of_lipids = 0; //sum of all states of IL
        for (int i = 0; i < membraneObject.numberOfFreeLipidsEachState.size(); i++) {
            number_of_lipids += membraneObject.numberOfFreeLipidsEachState[i];
        }
        meanComplexSize = print_complex_hist(complexList, assemblyfile, simItr, params, molTemplateList, number_of_lipids);
    }

    /* debug output */
    // std::cout << "propCalled: " << propCalled << std::endl;
    /* end of debug output */

    /*Write out final result*/
    std::cout << llinebreak << "End simulation\n";
    auto endTime = MDTimer::now();
    auto endTimeFormat = MDTimer::to_time_t(endTime);
    std::cout << "End date: ";
    //<< std::put_time(std::localtime(&endTimeFormat), "%F %T") << '\n';
    if (0 < strftime(charTime, sizeof(charTime), "%F %T", std::localtime(&endTimeFormat)))
        std::cout << charTime << '\n';
    std::chrono::duration<double> wallTime = endTime - totalTimeStart;
    std::cout << "\tWall Time: ";
    std::cout << wallTime.count() << " seconds\n";

    delete[] tableIDs;
    gsl_rng_free(r);
    return 0;
} // end main
