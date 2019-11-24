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
};

/*! \enum RxnKeyword
 * \ingroup Parser
 * \brief Reaction block keywords. For use with class RxnBase and its derived classes.
 */
enum class RxnKeyword : int {
    onRate = 0, //!< on rate of the reaction
    offRate = 1, //!< off rate of the reaction, if reversible
    norm1 = 2, //!< norm of the first reactant, for use in the phi angle, if applicable
    norm2 = 3, //!< norm of the second reactant, for use in the phi angle, if applicable
    sigma = 4, //!< sigma vector (will eventually be changed to an integer)
    assocAngles
    = 5, //!< angles for association [$\theta_1$,$\theta_2$,$\phi_1$,$\phi_2$,$\omega$]. See \ref association and Angles
    onMem = 6, //!< boolean for if the reaction occurs only on the membrane (eventually change to detecting a lipid in
               //!< the reaction
    rate = 7,
    isCoupled = 8, //!< is the reaction coupled to another reaction?
    isObserved = 9, //!< not used, but set by the parser
    observeLabel = 10, //!< label for the observable which tracks the product of the reaction
    bindRadSameCom = 11, //!< binding separaton to force association within same complex
    irrevRingClosure = 12, //!< is ring closure irreversible, i.e. binding within complex has probability unity
    creationRadius = 13, //!< the radius of the sphere around the parent molecule where a molecule is created
    loopCoopFactor = 14, //!< the factor to multiple the loop closure rate, i.e. loopCoopFactor=exp(-dG/kT)
    length3Dto2D = 15, //!< nm, length factor converts 3D to 2D rate
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
    };
    
  // parameter values
  int rank;
  int numMolTypes { 0 }; //!< number of MolTemplates. used to be Nprotypes
    //   int numIfaces { 0 }; //!< total number of interfaces on all MolTemplates. used to be nifaces
    int numTotalSpecies { 0 }; //!< total number of interfaces and states (including products!) possible in the system.
    unsigned nItr { 0 }; //!< number of timesteps requested by user. used to be int Nit.
    int numTotalComplex { 0 }; //!< number of complexes in the system at start
    unsigned numTotalUnits { 0 }; //!< number of total molecules + interfaces in the system at start
    double timeStep { 0 }; //!< timestep, in microseconds.  used to be deltat
    double mass { 0 }; //!< total mass of the system (do we need this?)
    int max2DRxns { 1000 }; //!< maximum number of allowed unique 2D reactions
    int maxUniqueSpecies { 1000 }; //!< maximum number of allowed unique species

    int numLipids { 0 }; //!< total number of lipids in the system
    
    double overlapSepLimit { 1.5 }; //!< in nm. COM-COM distance less than this is cancelled, for checkOverlap molecules

    std::string name {}; //!< name of the simulation
    Debug debugParams;

    // file names, for restart
    std::string trajFile { "trajectory.xyz" };
    std::string restartFile { "restart.dat" };

    // TODO TEMPORARY
    bool isNonEQ { false };

    // set in code itself
    double rMaxLimit { 0 };
    double rMaxRadius { 0 };

    // IO information
    bool fromRestart { false }; //!< is this simulation initialized from a restart file. used to be int restart
    int timeWrite { 10 }; //!< timestep interval to print timestep. used to be statwrite
    int trajWrite { 10 }; //!< timestep interval to write coordinates file. used to be configwrite
    int restartWrite { 10 }; //!< timestep interval to write a restart file
    int pdbWrite { -1 }; //!< interval to write pdb

    void display();
    void parse_paramFile(std::ifstream& paramFile);

    Parameters() = default;
    void set_value(std::string value, ParamKeyword keywords);
};
