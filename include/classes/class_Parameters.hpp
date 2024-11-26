/*! \file class_parameter
 * \ingroup Parser
 * \brief Parameter and related classes to hold simulation parameters provided by the user
 * ### Created on 5/25/18 by Matthew Varga
 * ### Purpose
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 * - Simulate until $X seconds parameter (automatically determines iterations based on timestep)
 */
#pragma once

#include <fstream>
#include <map>
#include <vector>

#include "classes/mpi_functions.hpp"

/*! \enum DebugKeywords
 * \ingroup Parser
 * \brief Debug parameters read in from the command line
 */
enum class DebugKeyword : int {
    forceAssoc = 0, //!< force the probability of association to be 1
    forceDissoc = 1, //!< force the probability of dissociation to be 1
    verbosity = 2, //!< the verbosity of output to the log file
};

/*! \enum ParamKeyword
 * \ingroup Parser
 * \brief Parameter file keywords. For use with class Parameters
 */
enum class ParamKeyword : int {
    numMolTypes = 0, //!< total number of molecule types in the system
    numTotalSpecies = 1, //!< total number of interfaces in the system, including states and product states
    nItr = 2, //!< number of iterations requested
    fromRestart = 3, //!< boolean for whether or not the simulation is starting from a restart file
    timeWrite = 4, //!< interval to write out timestep information
    trajWrite = 5, //!< interval to write to coordinates file (.xyz)
    timeStep = 6, //!< timestep
    numTotalComplex = 7, //!< total number of complexes in the system
    mass = 8, //!< total mass of all molecules
    //waterBox = 9, //!< water box dimensions
    restartWrite = 10, //!< interval to write to restart file
    pdbWrite = 11, //!< interval to write pdb for non-equilibrium trajectories
    overlapSepLimit = 12, //!< distance limit between centers of mass. if less than this for any proteins in a complex,
    //!< association is canceled
    name = 13, //!< name of the simulation
    checkPoint = 14, //!< interval to write checkpoint
    scaleMaxDisplace = 15, //!< scalar of average displacement that is acceptable upon association.
    transitionWrite = 16, //!< interval to write to transition matrix file
    clusterOverlapCheck = 17, //!< is overlap checked by cluster
    assocDissocWrite = 18, //!< write association and dissociation to a file
};

/*! \enum MolKeyword
 * \ingroup Parser
 * \brief Mol file keywords. For use with class Parameters.
 */
enum class MolKeyword : int {
    name = 0, //!< molecule name
    copies = 1, //!< molecule copy numbers
    isRod = 2, //!< is the molecule a rod (one dimensional)
    isLipid = 3, //!< is the molecule part of the bilayer)
    d = 4, //!< array of translational diffusion components [x,y,z]
    dr = 5, //!< array of rotational diffusion components [x,y,z]
    com = 6, //!< molecule center of mass (should be [0,0,0]
    state = 7, //!< interface state definitions
    mass = 8, //!< mass of the molecule
    checkOverlap = 9, //!< is overlap checked for during association
    bonds = 10, //!< optional bond definitions for PSF
    isPoint = 11, //!< is the Molecule a point (i.e. do all its interfaces overlap with the COM)
    isImplicitLipid = 12, //!< Is it an implicit lipid.
    countTransition = 13, //!< is transition counted during whole simulation
    transitionMatrixSize = 14, //!< size of the transition matrix
    outsideCompartment = 15, //!< outside the compartment
    insideCompartment = 16, //!< inside the compartment
    isPromoter = 17, //!< is the molecule a promoter for transcription initiation
};

/*! \enum RxnKeyword
 * \ingroup Parser
 * \brief Reaction block keywords. For use with class RxnBase and its derived classes.
 */
enum class RxnKeyword : int {
    onRate3Dka = 0, //!< on rate of the reaction, micro
    onRate3DMacro = 1, //!< on rate of the reaction, macro
    offRatekb = 2, //!< off rate of the reaction, if reversible, micro
    offRateMacro = 3, //!< off rate of the reaction, if reversible, macro
    norm1 = 4, //!< norm of the first reactant, for use in the phi angle, if applicable
    norm2 = 5, //!< norm of the second reactant, for use in the phi angle, if applicable
    sigma = 6, //!< sigma vector (will eventually be changed to an integer)
    assocAngles = 7, //!< angles for association [$\theta_1$,$\theta_2$,$\phi_1$,$\phi_2$,$\omega$]. See \ref association and Angles
    onMem = 8, //!< boolean for if the reaction occurs only on the membrane (eventually change to detecting a lipid in
    //!< the reaction
    rate = 9,
    isCoupled = 10, //!< is the reaction coupled to another reaction?
    isObserved = 11, //!< not used, but set by the parser
    observeLabel = 12, //!< label for the observable which tracks the product of the reaction
    bindRadSameCom = 13, //!< binding separaton to force association within same complex
    irrevRingClosure = 14, //!< is ring closure irreversible, i.e. binding within complex has probability unity
    creationRadius = 15, //!< the radius of the sphere around the parent molecule where a molecule is created
    loopCoopFactor = 16, //!< the factor to multiple the loop closure rate, i.e. loopCoopFactor=exp(-dG/kT)
    length3Dto2D = 17, //!< nm, length factor converts 3D to 2D rate
    rxnLabel = 18, //!< label of the reaction, used for coupled
    coupledRxnLabel = 19, //!< lable of the coupled reaction
    kcat = 20, //!<rate for a Michaelis-Menten reaction.
    excludeVolumeBound = 21, //!< once two sites bound, still exclude their partners if this is set true
    area3Dto1D = 22, //!< nm^2, area factor converts 3D to 1D rate
};

struct Parameters {
    /*! \class Parameters
     * \ingroup Parser
     * \brief Class to hold simulation parameters
     */

    struct Debug {
        bool forceAssoc { false };
        bool forceDissoc { false };
        bool printSystemInfo { false };
        int verbosity { 0 };
        /*
        Function serialize serializes the Debug info into array of bytes.
        */
        void serialize(unsigned char *arrayRank, int &nArrayRank) {
        PUSH(forceAssoc);
        PUSH(forceDissoc);
        PUSH(printSystemInfo);
        PUSH(verbosity);
        }
        /*
        Function deserialize deserializes the Debug info from array of bytes.
        */
        void deserialize(unsigned char *arrayRank, int &nArrayRank) {
        POP(forceAssoc);
        POP(forceDissoc);
        POP(printSystemInfo);
        POP(verbosity);
        }
    };

    // parameter values
    int rank;
    int numMolTypes { 0 }; //!< number of MolTemplates. used to be Nprotypes
    int numTotalSpecies { 0 }; //!< total number of interfaces and states (including products!) possible in the system.
    long long int nItr { 0 }; //!< number of timesteps requested by user. used to be int Nit.
    int numTotalComplex { 0 }; //!< number of complexes in the system at start
    unsigned numTotalUnits { 0 }; //!< number of total molecules + interfaces in the system at start
    double timeStep { 0 }; //!< timestep, in microseconds.  used to be deltat
    double mass { 0 }; //!< total mass of the system (do we need this?)
    int max2DRxns { 10000 }; //!< maximum number of allowed unique 2D reactions
    int maxUniqueSpecies { 1000 }; //!< maximum number of allowed unique species
    long long int itrRestartFrom { 0 }; //!< number of timesteps from when the simlation restart
    double timeRestartFrom { 0.0 }; //!< time from when the simulation restart (unit: s)

    static double dt; //!< timestep, in microseconds.
    static std::vector<long long int> lastUpdateTransition; //!< steps that last update transition matrix

    int numLipids { 0 }; //!< total number of lipids in the system

    double overlapSepLimit { 0.1 }; //!< in nm. COM-COM distance less than this is cancelled, for checkOverlap molecules
    bool implicitLipid { false };
    bool hasUniMolStateChange { false };
    bool hasCreationDestruction { false };

    bool checkUnimoleculeReactionPopulation{false};
    bool hasRankCommunicationForLargeComplex{false};

    double scaleMaxDisplace { 100.0 }; //!< in nm. defaults to a large number, which allows most moves even of large magnitude upon association.
    std::string name {}; //!< name of the simulation
    Debug debugParams;

    // file names, for restart
    std::string trajFile { "trajectory.xyz" };
    std::string restartFile { "restart.dat" };

    // TODO: TEMPORARY
    bool isNonEQ { false };

    // set in code itself
    double rMaxLimit { 0 };
    double rMaxRadius { 0 };

    // IO information. Iterators need to be long long because they can exceed 2^32
    bool fromRestart { false }; //!< is this simulation initialized from a restart file. used to be int restart
    bool assocDissocWrite { false }; //!< is association and dissociation written to a file
    long long int timeWrite { 10 }; //!< timestep interval to print timestep. used to be statwrite
    long long int trajWrite { 10 }; //!< timestep interval to write coordinates file. used to be configwrite
    long long int restartWrite { 10 }; //!< timestep interval to write a restart file
    long long int pdbWrite { -1 }; //!< interval to write pdb
    long long int checkPoint { -1 }; //!< interval to write checkpoint
    long long int transitionWrite { -1 }; //!< timestep interval to write transition matrix
    bool clusterOverlapCheck { false }; //!< is overlap checked by cluster 

    void display();
    void parse_paramFile(std::ifstream& paramFile);

    Parameters() = default;
    void set_value(std::string value, ParamKeyword keywords);

    /*
    Function serialize serializes the Parameters into array of bytes.
    */
    void serialize(unsigned char *arrayRank, int &nArrayRank) {
        PUSH(rank);
        PUSH(numMolTypes);
        PUSH(numTotalSpecies);
        PUSH(nItr);
        PUSH(numTotalComplex);
        PUSH(numTotalUnits);
        PUSH(timeStep);
        PUSH(mass);
        PUSH(max2DRxns);
        PUSH(maxUniqueSpecies);
        PUSH(itrRestartFrom);
        PUSH(timeRestartFrom);
        PUSH(dt);
        serialize_primitive_vector<long long int>(lastUpdateTransition, arrayRank, nArrayRank);
        PUSH(numLipids);
        PUSH(overlapSepLimit);
        PUSH(implicitLipid);
        PUSH(hasUniMolStateChange);
        PUSH(hasCreationDestruction);
        PUSH(scaleMaxDisplace);
        serialize_string(name, arrayRank, nArrayRank);
        debugParams.serialize(arrayRank, nArrayRank);
        serialize_string(trajFile, arrayRank, nArrayRank);
        serialize_string(restartFile, arrayRank, nArrayRank);
        PUSH(isNonEQ);
        PUSH(rMaxLimit);
        PUSH(rMaxRadius);
        PUSH(fromRestart);
        PUSH(timeWrite);
        PUSH(trajWrite);
        PUSH(restartWrite);
        PUSH(pdbWrite);
        PUSH(checkPoint);
        PUSH(transitionWrite);
        PUSH(clusterOverlapCheck);
        PUSH(checkUnimoleculeReactionPopulation);
        PUSH(hasRankCommunicationForLargeComplex);
    }
    /*
    Function deserialize deserializes the Parameters from arrayRank of bytes.
    */
    void deserialize(unsigned char *arrayRank, int &nArrayRank) {
        POP(rank);
        POP(numMolTypes);
        POP(numTotalSpecies);
        POP(nItr);
        POP(numTotalComplex);
        POP(numTotalUnits);
        POP(timeStep);
        POP(mass);
        POP(max2DRxns);
        POP(maxUniqueSpecies);
        POP(itrRestartFrom);
        POP(timeRestartFrom);
        POP(dt);
        deserialize_primitive_vector<long long int>(lastUpdateTransition, arrayRank, nArrayRank);
        POP(numLipids);
        POP(overlapSepLimit);
        POP(implicitLipid);
        POP(hasUniMolStateChange);
        POP(hasCreationDestruction);
        POP(scaleMaxDisplace);
        deserialize_string(name, arrayRank, nArrayRank);
        debugParams.deserialize(arrayRank, nArrayRank);
        deserialize_string(trajFile, arrayRank, nArrayRank);
        deserialize_string(restartFile, arrayRank, nArrayRank);
        POP(isNonEQ);
        POP(rMaxLimit);
        POP(rMaxRadius);
        POP(fromRestart);
        POP(timeWrite);
        POP(trajWrite);
        POP(restartWrite);
        POP(pdbWrite);
        POP(checkPoint);
        POP(transitionWrite);
        POP(clusterOverlapCheck);
        POP(checkUnimoleculeReactionPopulation);
        POP(hasRankCommunicationForLargeComplex);
    }
};
