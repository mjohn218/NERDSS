/*! @file classes.hpp
 * \brief Class declarations for general use in the program
 *
 * ### Created on 5/20/18 by Matthew Varga
 * ### Purpose
 * ***
 * Class declarations for all molecule and complex associated classes
 *
 * ### Notes
 * ***
 * CHANGED 6/7/18 -- moved stateChangeList on RxnBase -> RxnBase::RateState
 *
 * ### TODO List
 * ***
 *  - TODO PRIORITY MED: split into their different classes (maybe keep derived classes together?)
 *  - TODO PRIORITY LOW: I kind of want to make all data members of MolTemplate, Interface, Angles, and the reaction
 * classes private
 *  - TODO PRIORITY LOW: change all `int` indices to `size_t`
 *  - TODO PRIORITY LOW: rename class names to follow Webkit standards (no abbreviations, camel case)
 */
#pragma once

#include "classes/class_MolTemplate.hpp"
#include "classes/class_Vector.hpp"
#include "classes/class_Membrane.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <limits>

/*! \defgroup SimulClasses
 * \brief Classes actually used for simulation objects
 */

extern int propCalled;

enum class TrajStatus : int{
    none = 0,
    dissociated = 1,
    associated = 2,
    propagated = 3,
    canBeResampled = 4,
    empty = 5,
};

/*!
 * \ingroup SimulClasses
 * \brief Contains information for one specie (protein, lipid, etc.) in the system
 */
struct Molecule {
    /*!
     * \brief Storage for information on Molecules this Molecule encounters during a timestep
     */
    struct Encounter {
        size_t theirIndex{0}; //!< the encountered Molecule's index in moleculeList
        size_t myIface{ 0 }; //!< this Molecule's interface which encountered the Molecule
        size_t theirIface{ 0 }; //!< the encountered Molecule's interface which encountered this Molecule
        size_t rxnItr {0};
        double probability {0};

        Encounter() = default;
        Encounter(size_t theirIndex, size_t myIface, size_t theirIface, size_t rxnItr, double probability)
            : theirIndex(theirIndex)
            , myIface(myIface)
            , theirIface(myIface)
            , rxnItr(rxnItr)
            , probability(probability)
        {
        }
    };

    /*!
     * \brief Holds information on the interactions this Molecule is currently participarting in
     */
    struct Interaction {

        // TODO: change to int and make default to -1
        int partnerIndex{ -1 }; //!< bound partner index in moleculeList
        int partnerIfaceIndex{ -1 }; //!< interface index of the Molecule's partner
        int conjBackRxn{ -1 }; //!< back reaction of the forward reaction that formed this Interaction

        /*!
         * \brief Sets all the Interaction indices to 0. Used in break_interaction().
         *
         * TODO: Change this to -1.
         */
        void clear();

        Interaction() = default;
        Interaction(int partnerIndex, int partnerIfaceIndex)
            : partnerIndex(partnerIndex)
            , partnerIfaceIndex(partnerIfaceIndex)
        {
        }
        Interaction(int partnerIndex, int partnerIfaceIndex, int conjBackRxn)
            : partnerIndex(partnerIndex)
            , partnerIfaceIndex(partnerIfaceIndex)
            , conjBackRxn(conjBackRxn)
        {
        }
    };

    struct Iface {
        /*! \struct Molecule::Iface
         * \brief Holds the interface coordinate and the interface's state
         */

        Coord coord{ 0, 0, 0 }; //!< Coordinate of the interface (absolute, not relative)
        char stateIden{ '\0' }; //!< current state of the interface.
        int stateIndex {-1}; //!< index of the current state in MolTemplate::Interface::stateList
        int index{ -1 }; //!< this interface's absolute index (i.e. the index of its current state)
        int relIndex { -1 }; //!< this interface's relative index (index in Molecule::interfaceList)
        int molTypeIndex{ -1 }; //!< index of this interface's parent Molecule's MolTemplate (makes comparisons easier)
        bool isBound{ false }; //!< is this interface bound

        Interaction interaction;

        /*!
         * \brief This simply changes the state of the interface.
         * \param[in] newRelIndex Index of the new State in Interface::stateList
         * \param[in] newAbsIndex Constant reference to the new State in Interface::stateList.
         */
        void change_state(int newRelIndex, int newAbsIndex, char newIden);

        Iface() = default;
        explicit Iface(const Coord& coord)
            : coord(coord)
        {
        }
        Iface(const Coord& coord, char state)
            : coord(coord)
            , stateIden(state)
        {
        }
        Iface(char state, int index, int molTypeIndex, bool isBound)
            : stateIden(state)
              , index(index)
              , molTypeIndex(molTypeIndex)
              , isBound(isBound)
        {
        }
        Iface(const Coord& coord, char state, int index, int molTypeIndex, bool isBound)
            : coord(coord)
            , stateIden(state)
            , index(index)
            , molTypeIndex(molTypeIndex)
            , isBound(isBound)
        {
        }
    };
    
    int myComIndex{ -1 }; //!< which complex does the molecule belong to
    int molTypeIndex{ -1 }; //!< index of the Molecule's MolTemplate in molTemplateList
    int mySubVolIndex{ -1 };
    int index{ -1 }; //!< index of the Molecule in moleculeList
    double mass{ -1 }; //!< mass of this molecule
    bool isLipid{ false }; //!< is the molecule a lipid 
    Coord comCoord; //!< center of mass coordinate
    std::vector<Iface> interfaceList; //!< interface coordinates
    bool isEmpty{ false }; //!< true if the molecule has been destroyed and is void
    TrajStatus trajStatus{TrajStatus::none}; //!< Status of the molecule in that timestep
    
    bool isImplicitLipid = false;
    int linksToSurface {0};//!<store each proteins links to surface, to ease updating complex.
    // static variables
    static int numberOfMolecules; //!< counter for the number of molecules in the system
    static std::vector<int> emptyMolList; //!< list of indices to empty Molecules in moleculeList

    // association variables
    // temporary positions
    Coord tmpComCoord{}; //!< temporary center of mass coordinates for association
    std::vector<Coord> tmpICoords{}; //!< temporary interface coordinates for association

    // New Encounter lists
//    std::vector<Interaction> interactionList; //!< list of interactions the Molecule is participating in
//    std::vector<Encounter> encounterList; //!< list of Molecule Encounters during a timestep
//    std::vector<Encounter> prevEncounterList; //!< list of Molecule Encounters during a timestep

    // Legacy encounter lists
    /*Vectors for possible association reactions*/
    std::vector<int> freelist; // legacy, may be replaced
    std::vector<int> assoclist; // These species are capable of binding.
    std::vector<int> bndlist; // These species are capable of dissociation
    std::vector<int> bndpartner; // It if is bound, who is it bound to? Make this have the same numbering !!
    std::vector<int> bndRxnList;
    std::vector<int> bndiface; // If it is bound, though which interface!
//    int ncross = 0;
//    int movestat = 0;
    std::vector<double> probvec;
    std::vector<int> crossbase; // proteins base encountered
    std::vector<int> mycrossint; // interfaces base encountered other proteins with
//    std::vector<int> crossrxn;
    std::vector<std::array<int, 3>> crossrxn;
    
    std::vector<double> probvec_dissociate;

    /*Vectors for reweighting!*/
    std::vector<int> prevlist;
    std::vector<int> currlist;
    std::vector<int> prevmyface;
    std::vector<int> currmyface;
    std::vector<int> prevpface;
    std::vector<int> currpface;

    std::vector<double> prevnorm;
    std::vector<double> currprevnorm;
    std::vector<double> ps_prev;
    std::vector<double> currps_prev;
    std::vector<double> prevsep;
    std::vector<double> currprevsep;

    void write_crd_file(std::ofstream& os) const;
    void write_crd_file_cout() const;
    friend std::ostream& operator<<(std::ostream& os, const Molecule& mol);

    // association member functions
    void display_assoc_icoords(const std::string& name);
    void update_association_coords(const Vector& vec);
    void set_tmp_association_coords();
    void clear_tmp_association_coords();

    // other reaction member functions
    void create_random_coords(const MolTemplate& molTemplate, const Membrane &membraneObject);
    void destroy(std::vector<int>& emptyMolList);

    void display(const MolTemplate& molTemplate) const;
  void display_all() const;
    void display_my_coords(const std::string& name);

    bool operator==(const Molecule& rhs) const;
    bool operator!=(const Molecule& rhs) const;

    Molecule() {}
    Molecule(int _mycomplex, Coord _comcoords)
        : myComIndex(_mycomplex)
        , comCoord(_comcoords)
    {
    }
};

struct Complex {
    /*! \struct Complex
     * \ingroup SimulClasses
     * \brief Contains information on each complex, i.e. bound Molecules
     */

public:
    Coord comCoord; //!< Complex's center of mass coordinate
    int index{ 0 }; //!< index of this Complex in complexList
    double radius{}; //!< radius of the Complex's bounding sphere
    double mass {};
    std::vector<int> memberList{}; //!< list of member Molecule's indices
    std::vector<int> numEachMol{}; //!< list of the number of each Molecules in this complex
    Coord D{ 0, 0, 0 }; //!< Complex's translational diffusion constants
    Coord Dr{ 0, 0, 0 }; //!< Complex's rotational diffusion constants
    bool isEmpty{ false }; //!< true if the complex has been destroyed and is a void
    bool OnSurface{ false }; // to check whether on the implicit-lipid membrane.

    // static variables
    static int numberOfComplexes; //!< total number of complexes in the system. starts out equal to numberOfMolecules
    static int currNumberComTypes;
    static int currNumberMolTypes;
    static std::vector<int> emptyComList; //!< list of indices to empty Complexes in complexList

    // TODO: TEMPORARY
    static std::vector<int> obs; //!< TEMPORARY observables vector
  //std::vector<int> NofEach;//!< number of each protein type in this complex
    int linksToSurface = 0; //!< for an adsorbing surface, number of bonds/links formed between this complex and the surface.
    int iLipidIndex = 0;//!< If you need to look up the implicit lipid, this is its molecule index
    // Simulation status parameters
    int ncross{ 0 };
//    int movestat{ 0 };
    TrajStatus trajStatus {TrajStatus::none};
    Vector trajTrans;
    Coord trajRot;
    Coord tmpComCoord;

    friend std::ostream& operator<<(std::ostream& os, const Molecule& mol);

    void update_properties(const std::vector<Molecule>& moleculeList, const std::vector<MolTemplate>& molTemplateList);
    void display();
    void display(const std::string& name);
    Complex create(const Molecule& mol, const MolTemplate& molTemp);
    void destroy(std::vector<Molecule>& moleculeList, std::vector<int>& emptyMolList,
                     std::vector<int>& emptyComList);
    void put_back_into_SimulVolume(
        int& itr, Molecule& errantMol, const Membrane &membraneObject, std::vector<Molecule>& moleculeList);
    void translate(Vector transVec, std::vector<Molecule>& moleculeList);
    void propagate(std::vector<Molecule>& moleculeList);

    Complex() = default;
//    Complex(Molecule mol, Coord D, Coord Dr);
    Complex(const Molecule& mol, const MolTemplate& oneTemp);
    Complex(int _index, const Molecule& _memMol, const MolTemplate& _molTemp);
    Complex(Coord comCoord, Coord D, Coord Dr);
};

