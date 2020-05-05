/*! @file bngl_parser_classes.hpp
 * \brief Functions related to association
 *
 * ### Created on 5/2/18 by Matthew Varga
 * ### Purpose
 * ***
 *  - class definitions for use only in reaction file parsing
 * Contains functions specific to association of two Molecules
 *
 * ### TODO List
 * ***
 *  - Write comments
 *  - split into different class headers and function definitions
 *  - PRIORITY SUPER LOW: Make ParsedMol just use RxnIface, not the nested class
 */
#pragma once

#include "classes/class_Rxns.hpp"
#include <array>
#include <iostream>
#include <limits>
#include <map>

/*! \enum Involvement
 * \ingroup Parser
 * \brief Enumeration of reactant involvements in a parsed reaction.
 */
enum struct Involvement : int {
    none = -1, //!< default for initialization
    ancillary = 0, //!< interface is totally ancillary to the reaction, but needs to be present
    possible = 1, //!< nterface may or may not be involved. assigned during parsing of the BNGL
    intermInvolved = 2, //!< intermediate involvement, for assignment when determining indices. will be changed later
    //!< (not implemented)
    interactionChange = 3, //!< involved and only the interaction changes
    stateChange = 4, //!< interface only changes state
    facilitator = 5, //!< interface which facilitates a bimolecular state change
    bothStateAndInteractionChange = 6, //!< the interface changes both interaction and state (not implemented)
};
std::ostream& operator<<(std::ostream& os, const Involvement& involve);

// Struct to hold information for each molecule parsed from file
struct ParsedMol {
    /*! \struct ParsedMol
     * \ingroup Parser
     * \brief Holds MolTemplate information temporarily during input file parsing
     */
    struct IfaceInfo {
        /*! \struct ParsedMol::IfaceInfo
         * \ingroup Parser
         * \brief Holds interface information temporarily during input file parsing
         */

        int molTypeIndex { -1 }; //!< index of the MolTemplate this interface's parent ParsedMol (Molecule) belongs to
        std::string ifaceName {};
        int absIndex { -1 }; //!< absolute index of this interface
        int relIndex { -1 };
        char state { '\0' }; //!< the state identifier of the interface
        bool isBound { false }; //!< is the interface bound to anything
        bool indexFound { false }; //!< boolean for determine_iface_indices()
        int bondIndex { -1 }; // -1 if not bound, or if reactant side, 0 if wildcard
        std::string ifaceAndStateName; //!< full name of the interface + state (if it has one)
        Involvement ifaceRxnStatus { Involvement::none };
        int speciesIndex { -1 }; //!< which species is the interface on (only used for bimolecular reactions)

        friend std::ostream& operator<<(std::ostream& os, const ParsedMol::IfaceInfo& oneInfo);
        bool operator==(const ParsedMol::IfaceInfo& iface) const;
        bool operator!=(const ParsedMol::IfaceInfo& iface) const;
        bool operator==(const RxnIface& rxnIface) const;
        bool operator<(const IfaceInfo& rhs) const;
        bool operator>(const IfaceInfo& rhs) const;

        // change the status of the reaction, for products
        void change_ifaceRxnStatus(int newIndex, Involvement newRxnStatus);

        void display()
        {
            std::cout << absIndex << ", ";
            if (state != '\0')
                std::cout << state << ", ";
            else
                std::cout << "NO STATE, ";
            std::cout << std::boolalpha << isBound << ", " << bondIndex << ", " << ifaceRxnStatus;
        }

        IfaceInfo() = default;
        IfaceInfo(std::string _ifaceName, char _state, bool _isBound, Involvement _ifaceRxnStatus, int speciesIndex);
        IfaceInfo(std::string _ifaceName, char _state, bool _isBound, int _bondIndex, Involvement _ifaceRxnStatus);
        IfaceInfo(std::string _ifaceName, char _state, bool _isBound, int _bondIndex, Involvement _ifaceRxnStatus,
            int speciesIndex);

        // TODO: Temporary
        explicit IfaceInfo(const Interface& _iface)
            : ifaceName(_iface.name)
            , absIndex(_iface.stateList.at(0).index)
            , relIndex(_iface.index)
        {
        }
    };

    int molTypeIndex { -1 }; //!< index of the corresponding MolTemplate
    int specieIndex { -1 }; //!< The index of the specie in the reaction
    std::string molName; //!< name of the corresponding MolTemplate (makes it easier for matching)
    std::vector<IfaceInfo> interfaceList;

    bool operator==(const ParsedMol& mol1) const;
    std::ostream& display_full_name(std::ostream& os) const;

    // member functions
    //    void set_ifaceAndStateNames();
    void set_molTypeIndex(const std::vector<MolTemplate>& molTemplateList);
    //    bool has_no_Z_iface();

    void display() const;

    ParsedMol() = default;
    explicit ParsedMol(const MolTemplate& oneTemp);
};

// I am so sick of writing this out, so here's a type alias
using parsedIface = std::pair<const std::string, ParsedMol::IfaceInfo>;
struct ParsedRxn : public ForwardRxn {
    /*!\struct ParsedRxn
     * \ingroup Parser
     * \brief Contains all the information parsed from the user provided reaction file.
     * Later converted into a ForwardRxn or CreateDestructRxn, based on the reaction type.
     */

    // flags created during parsing for use later on
    bool willBeMultipleRxns { false }; //!< does an interface have no explicit state
    std::multimap<int, ParsedMol::IfaceInfo> noStateList; //!< list of ifaces which have states but do not declare them.
    //!< Form of the multimap is {index of ParsedMol the interface is on, index of the interface in ifaceList}
    bool includesIntCoords { false }; //!< does the file include internal coordinates
    std::vector<Vector> norms;
    double creationRadius { 1.0 }; //!< see CreateDesructRxn

    double onRate3Dka { std::numeric_limits<double>::quiet_NaN() }; //!< the forward rate of the reaction, micro
    double onRate3DMacro { std::numeric_limits<double>::quiet_NaN() }; //!< the forward rate of the reaction, macro
    double offRatekb { std::numeric_limits<double>::quiet_NaN() }; //!< the rate of a reversible reaction's back reaction, micro
    double offRateMacro { std::numeric_limits<double>::quiet_NaN() }; //!< the rate of a reversible reaction's back reaction, macro
    double kcat {std::numeric_limits<double>::quiet_NaN() }; //!< the rate for a Michaelis-Menten reaction catalysis.
    std::vector<std::vector<RxnIface>> otherIfaceLists; //!< ancillary interfaces which must be present
    std::pair<RxnIface, RxnIface> stateChangeIface; //!< interfaces which don't change interaction but change state

    // TODO: Rename productList and reactantList, to allow RxnBase to have those names
    std::vector<ParsedMol> productList; //!< list of parsed product molecules (ParsedMol)
    std::vector<ParsedMol> reactantList; //!< list of parsed reactant molecules (ParsedMol)
    std::vector<ParsedMol::IfaceInfo> rxnReactants; //!< reactant interfaces with Involvement::interactionChange
    std::vector<ParsedMol::IfaceInfo> rxnProducts; //!< reactant interfaces with Involvement::interactionChange

    /* REACTION CREATION FUNCTIONS */
    /*!
     * \brief Checks previously parsed ForwardRxns for the reactant of a unimolecular state change reaction with a
     * bound reactant.
     *
     * If the reactants are bound in a unimolecular state change reaction, check through all previous forward
     * reactions for the reactants. If none of the reactions already have the state change product as a reaction
     * forming the state change reactant exists, create a new absolute interface index.
     */
    void check_previous_bound_states(int& totSpecies, const std::vector<ForwardRxn>& forwardRxns, const std::vector<MolTemplate>& molTemplateList);

    /*!\ingroup Parser
     * \brief The purpose of this function is to determine which of the parsed interfaces are reactants.
     *
     * These are the interfaces which change their interaction status and will go into the rxnReactants and rxnProducts
     * lists
     *   - I want to do this by looking at the products which have the ifaceRxnStatus =
     *   interactionChange/stateChange/bothStateAndInteractionChange, and comparing them to the products which
     *   have possible/stateChange (since the other involved Involvments are not determine for the products yet to the
     *   reactants which have the ifaceRxnStatus = possible
     */
    void determine_reactants();

    /*! \ingroup Parser
     * \brief Need to determine creation products separately as they are not determined elsewhere
     *
     * This is because there are no reactants to compare to the products (as in check_for_valid_states() and
     * check_for_state_change())
     *
     * @param molTemplateList list of MolTemplates to compare to product molName
     */
    void determine_creation_products(const std::vector<MolTemplate>& molTemplateList);

    /*! \ingroup Parser
     * \brief Function to determine non-interacting interfaces' involvement in the reaction
     * This function is meant to go through the reactantList of a parsed reaction and determine which ifaces are:
     *   1. ancillary - required ifaces which don't change => otherIfaceList
     *   2. change state - ifaces which don't change interaction but to change state => stateChangeIface
     */
    void create_other_iface_lists(const std::vector<MolTemplate>& molTemplateList);

    /*! \ingroup Parser
     * \brief  This function is meant to check through each of the existing ForwardRxns and determine if the currently
     * parsed reaction is a duplicate
     *
     * - if same products && reactants, different rate || different isOnMem || different ancillary ifaces
     *              => make a conditional rate in the existing ForwardRxn
     * - if same products && reactants && rate && isOnMem && ancillary ifaces
     *              => this reaction is a duplicate and is discarded
     *
     * TODO: Sigma and the association angles should not be changed, maybe check if they are and kick back
     */
    bool check_for_conditional_rates(
        int& totSpecies, std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns);

    /*! \ingroup Parser
     * \brief This function is meant to make on individual reaction when splitting a reaction into multiples
     *
     * This only makes a split reaction for situations in which there is only one interface with no explicit state
     *
     */
    ParsedRxn make_one_split_reaction(
        int& totSpecies, int parsedMolIndex, ParsedMol::IfaceInfo& reactIface, const Interface::State& tempState);

    /*! \ingroup Parser
     * \brief This function splits a ParsedRxn containing interfaces which don't declare an explicit state and splits
     * the reaction ito multiple ParsedRxns.
     *
     * If an interface, which has a state as per its MolTemplate, is included in the reactants and products but does not
     * declare a specific state, that means the interface must be present for the reaction to occur, but that its state
     * does not matter. So, we split the one ParsedRxn into multiple ParsedRxns, wherein the interface is included with
     * each of its available states
     *
     * \param[in] forwardRxns list of all current ForwardRxns
     * \param[in] molTemplateList list of all MolTemplates, as provided by the user
     * \param[out] std::vector<ParsedRxn> list of new ParsedRxns for each state detected
     */
    std::vector<ParsedRxn> split_into_multiple_reactions(
        int& totSpecies, std::vector<ForwardRxn>& forwardRxns, const std::vector<MolTemplate>& molTemplateList);

    /*!\ingroup Parser
     * \brief Takes the information from the ParsedRxn (reaction parsed from the reactions block) and creates the
     * relevant reactions for use in the simulations
     */
    void assemble_reactions(std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
        std::vector<CreateDestructRxn>& createDestructRxns, const std::vector<MolTemplate>& molTemplateList);

    // Parsing function
    void set_value(std::string& line, RxnKeyword rxnKeyword);
    std::pair<bool, std::string> isComplete(const std::vector<MolTemplate>& molTemplateList);

    // display functions
    void display_angles() const;
    void display() const override;

    ParsedRxn() = default;
};

bool operator==(const std::vector<RxnIface>& rxnIfaceList, const std::vector<ParsedMol::IfaceInfo>& parsedIfaceList);

// struct IntCoordCont {
//    /*! \struct IntCoordCont
//     * \ingroup Parser
//     * \brief This class simply holds internal coordinates to calculate association angles (and sigma) if user
//     provided
//     */
//    struct Mol {
//        Coord comCoord;
//        int proType;
//        int reactantIface;
//        std::vector<Coord> ifaceCoords;
//    };
//
//    Mol reactant1;
//    Mol reactant2;
//
//    void display() const;
//    void calc_angles_from_crds(ParsedRxn& parsedRxn, const std::vector<MolTemplate>& molTemplateList);
//
//    IntCoordCont() = default;
//};
