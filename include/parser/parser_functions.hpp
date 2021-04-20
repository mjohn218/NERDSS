/*! \file parser_math_functions.hpp
 * \ingroup Parser
 * Created on 6/7/18 by Matthew Varga
 * Notes:
 * ### Purpose
 * ***
 * Main function to take in formatted reaction files and convert them into ForwardRxn, BackRxn, and CreateDestructRxn
 * for use in the simulation
 *
 * ### Notes
 * ***
 * #### Adding Keywords
 * #### Limitations
 *  - Implementing the ability for an interface to change both interaction and state at the same time is too complicated
 *    right now.
 *  - Cannot have any unbinding events in reaction file, e.g. A(a!1).B(b!1) <-> A(a) + B(b)
 *  - an iface having a wildcard bond in the product but is a reactant has to be handled at association time
 *  - Limitation: Kinases (or other molecules which facilitate bimolecular state change reactions) can only have one
 *    interface listed in the reaction, e.g. K(k) + A(a~x) <-> K(k) + A(a~y),
 *    not K(k1, k2) + A(a~x) <-> K(k1, k2) + A(a~y)
 *  - Each reaction can only have two molecules. Things can be bound, but can't have three things, e.g.
 *    A(a) + B(b) -> whatever and A(a!*).A(a!*) + B(b) -> whatever are acceptable whereas A(a!*) + A(a!*) + B(b) ->
 *    whatever is not
 *
 * #### Programming Notes
 *  - `std::find_if()` has the same complexity as a for-loop based on the same predicate (at most first-last actions of
 *    the predicate). Follows this format:
 *    `std:find_if(wholep.begin(), wholep.end(), [](const MolTemplate& oneTemp)->bool {return oneTemp.name ==
 *    whever.name;})`
 *   - the last field in the above implementation of `find_if()` is a lambda expression (basically a one-off function).
 *     Follows this format: `[captured variables](function arguments)->return type {function itself}`
 *  - the `?` denotes a ternary (conditional) operator, follows this format:
 *    `(condition) ? : if true : if false`
 *    so it evaluates the condition and does the if true statement if so and the if false statement if not
 *  - You'll find a lot of try-catch blocks and if-else statements in the actual reading in of the values. This is on
 *    purpose: I wanted to write it so that it is robust enough to deal with the reaction file not being perfectly
 *    formatted.
 *  - All ifaces are read in and converted to lowercase, and all states are converted to caps, to make sure comparisons
 *    work.
 *
 * ### TODO List
 * ***
 *  - TODO PRIORITY HIGH: Need a check for whether or not the two sides have the same molecules
 *  - TODO PRIORITY HIGH: need a check for unbalanced iface numbers
 *  - PRIORITY MED: add a check for initialized, but unset, variables (in all three blocks)
 *  - TODO PRIORITY LOW: when determining reversibility, (iterator +2/3) is hardcoded
 *  - TODO PRIORITY LOW: add calculation of total interfaces, etc., if they aren't provided
 *  - TODO: Currently no way to tell program which interface state to use. Just uses first listed state as default state
 *  - TODO PERFORMANCE: A lot of this would probably be easier if I made stateChange and interactionChange vectors, so i
 *  - TODO: As of yet, creation and destruction reactions need to define at least one interface. It's overall pretty
 * unwieldly for creation an destruction in general. don't have to loop
 */

/*! \defgroup Parser
 * \brief Functions specifically used in the parsing of input files
 *
 * NOTES:
 * #### Adding Keywords
 * To add keywords, add the desired keyword to the following places for the following keyword types
 * - Parameter file keywords
 *   - ParamKeyword in classes/class_Parameters.hpp
 *   - paramKeywords and in the switch statement within Parameters::set_value() in classes/class_Parameters.cpp
 * - Molecule Information File keywords
 *   - MolKeyword in classes/class_Parameters.hpp
 *   - molKeywords in parse_molFile() in parser/parser_functions.cpp
 *   - in the switch statement in MolTemplate::set_value() in classes/class_MolTemplate.cpp
 * - Reaction block keywords
 *   - RxnKeyword in classes/class_Parameters.hpp
 *   - rxnKeywords in parse_reaction() in parser/parser_functions.cpp
 *   - in the switch statement in ParsedRxn::set_value() in parser/class_bngl_parser_functions.cpp
 */

#pragma once

#include "classes/class_Membrane.hpp"
#include "classes/class_Observable.hpp"
#include "classes/class_Quat.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_bngl_parser.hpp"
#include <algorithm>

/* MAIN */
/*! \ingroup Parser
 * \brief Main input file parsing function.
 *
 * Delegates to parse_reactionFile(), parse_paramFile(), and parse_molFile().
 */
void parse_input(std::string& fileName, Parameters& params, std::map<std::string, int>& observableList,
    std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<CreateDestructRxn>& createDestructRxns, std::vector<MolTemplate>& molTemplateList, Membrane& membraneObject);

/*!\ingroup Parser
 * \brief This function parses input for restart with add.inp.
 *
 * TODO: 
 */
void parse_input_for_add(std::string& fileName, Parameters& params, std::map<std::string, int>& observableList,
    std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<CreateDestructRxn>& createDestructRxns, std::vector<MolTemplate>& molTemplateList, Membrane& membraneObject, int numDoubleBeforeAdd);

/*!\ingroup Parser
 * \brief This function parses command line flag.
 *
 * TODO: 
 */
void parse_command(int argc, char* argv[], Parameters& params, std::string& paramFileName, std::string& restartFileName, std::string& addFileName, unsigned int& seed);

/* INDEX DETERMINATION FUNCTIONS */

/*!\ingroup Parser
 * \brief This function takes a target ParsedMol and makes sure that all of the states in its ifaceList are valid.
 *
 * If it finds a reactant with a state, it then checks the products for a matching interface to see if its state
 * changed. If so, its Involvement gets set to stateChanged
 *
 * @param[in] targetMol ParsedMol to check for valid Interface::States
 * @param[in] parsedRxn& the reaction currently being parsed
 * @param[in] molTemplateList list of all MolTemplates
 *
 * TODO: now that ParsedMol has a molTypeIndex, take out the search for the MolTemplate and just pass the correct
 * one to this function
 */
void check_for_valid_states(
    size_t parsedMolIndex, ParsedMol& targMol, ParsedRxn& parsedRxn, const std::vector<MolTemplate>& molTemplateList);

/*!\ingroup Parser
 * \brief The purpose of this function is to take some target iface, which has a state, of the reactants and compare
 * it to the product ifaces to see if the state changed
 *
 * @param[in] targetIface interface to check for a state change
 * @param[in] the molecule on which the targetIface is located
 */
void check_for_state_change(ParsedMol::IfaceInfo& targetIface, ParsedMol& targetMol, ParsedRxn& parsedRxn);

/*!\ingroup Parser
 * \brief This function determines the iface indices of a molecule based on its name
 * @param[in] totSpecies total number of unique reactants and products in the list of reactions
 * @param[in] targMol target molecule for which we want to find the iface indices
 */
void determine_iface_indices(int specieIndex, int& totSpecies, ParsedMol& targMol, ParsedRxn& parsedRxn,
    const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList);

/*!\ingroup Parser
 * \brief Takes the string Interface names of two bound Interfaces and determines their Interface ifaceIndex
 *
 * Does this by searching through the forwardRxns list and checking each reactantList productList for their
 * corresponding interfaces. If no reactions are found with those interfaces, a new product is formed and the total
 * number of reactions is iterated upon. NOTES:
 *  - if the code has gotten to this point, it means that the reactant/product being parsed is not an original
 *  late interface and must be a product of a provided reaction
 *  - several ways I can think of implementing this:
 *    - find what the iface is bound to and search the forwardRxns list for a reaction whose reactants are those
 *      interfaces
 *    - find what the iface is bound to and search the backRxns for the product which corresponds to that
 *      defined interface
 *    - the former is probably easier. Is there a situation in which it wouldn't work?
 */
void determine_bound_iface_index(int& totSpecies, ParsedMol::IfaceInfo& targIface,
    ParsedRxn& parsedRxn, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList);

bool areCorrespondingMolecules(std::pair<std::string, ParsedMol::IfaceInfo>& prodMolIface,
    std::pair<std::string, ParsedMol::IfaceInfo>& reactMolIface);

/*! \ingroup Parser
 * \brief returns true if the two ifaces are identical but for their states
 */
bool areSameExceptState(const ParsedMol::IfaceInfo& iface1, const ParsedMol::IfaceInfo& iface2);
/*************************/

/* PARSING FUNCTIONS */
bool skipLine(std::string line);

/*!\ingroup Parser
 * \brief Takes a molecule from a reaction and parses its molecule (MolTemplate) type, interfaces, interactions, and
 * states.
 *
 * @param[in] totSpecies total species (Interfaces + States) in the system so far
 * @param[in] speciesIndex which of the two reactant species is this molecule on
 * @param[in] isProductSide is this molecule on the product side
 * @param[in] oneMol the molecule to be parsed
 * @param[out] ParsedMol temporary molecule containing the parsed information
 *
 * - Format of each molecule in the reaction will be mol(interface~state!bondIndex), with $state and $bondIndex
 * being optional.
 */
ParsedMol parse_molecule_bngl(int& totSpecies, bool isProductSide, std::pair<std::string, int> oneMol);

/*!\ingroup Parser
 * \brief Takes starting copy numbers of states from a input line
 *
 * @param[in] oneLine the line to be parsed
 * @param[out] ParsedMolNumState struct containing the parsed information
 *
 * - Format of input line  will be 100 (interface~state1), 100 (interface~state2)
 */
ParsedMolNumState parse_number_bngl(std::string oneLine);

/*!\ingroup Parser
 * \brief This function determines the reaction type and reversibility, and parses the reaction file accordingly
 *
 * @param[in] reactionFile user provided reaction file
 * @param[in] totSpecies total number of species (reactants and products)
 * @param[in] molTemplateList vector of MolTemplates
 * @param[in] forwardRxns vector of ForwardRxns
 * @param[in] backRxns vector of BackRxns, inverse reactions of a corresponding reversible ForwardRxn
 * @param[in] createDestructRxns vector of CreateDestructRxns
 */
void parse_reaction(std::ifstream& reactionFile, int& totSpecies, int& numProvidedRxns,
    std::vector<MolTemplate>& molTemplateList, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<CreateDestructRxn>& createDestructRxns,
    std::map<std::string, int>& observablesList, Membrane& membraneObject);

bool read_boolean(std::string fileLine);

/*! \ingroup Parser
 * \brief Reads a formatted block of internal coordinates starting with the center of mass coordinate
 *
 * Format:
 *
 *     COM   x y z
 *     Iface x y z
 */
void read_internal_coordinates(std::ifstream& molFile, MolTemplate& molTemplate);

/*!\ingroup Parser
 * \brief just a simple function to remove comments, if they exist
 */
void remove_comment(std::string& line);

/*!\ingroup Parser
 * \brief This function reads an unformatted (other than the internal coordinates) molecule information file.
 *
 * @param[in] mol molecule to whom the molecule information file belongs
 * @param[out] completed MolTemplate
 */
MolTemplate parse_molFile(std::string& mol);

/*! \ingroup Parser
 * \brief Reads the diffusion constant arrays from the parameters block
 */
std::vector<double> parse_input_array(std::string& line);

/*!\ingroup Parser
 * \brief Parses the state lines from the molecule information files and adds them to the MolTemplate
 */
void parse_states(std::string& line, MolTemplate& molTemplate);

/*!\ingroup Parser
 * \brief Reads a boolean in either numeric or alphabetical format
 *
 * @param[in] fileLine line from the input file containing a boolean
 * @param[out] bool parsed boolean
 */
void read_bonds(int numBonds, std::ifstream& molFile, MolTemplate& molTemplate);
/*************************/

/*************************/

/* LIST POPULATION */
size_t find_molTypeIndex_from_ifaceIndex(const int& targIfaceIndex, const std::vector<MolTemplate>& molTemplateList);

/*! \ingroup Parser
 * \brief Just sets the conjugate reaction iterators to the iterator of the last element in forwardRxns
 *  and backRxns;
 */
void create_conjugate_reaction_itrs(std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns);

/*!\ingroup Parser
 * \brief This function populates the vectors of the Interface's pertinent reactions for each State
 */
void populate_reaction_lists(const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    const std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<MolTemplate>& molTemplateList);
/*******************/
void populate_reaction_lists_for_add(const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    const std::vector<CreateDestructRxn>& createDestructRxns,
    std::vector<MolTemplate>& molTemplateList, int addForwardRxnNum, int addBackRxnNum, int addCreateDestructRxnNum);

/* DISPLAY */
std::string write_mol_iface(std::string mol, std::string iface);

/*! \ingroup Parser
 * \brief Just displays all reactions as parsed from the input file
 */
void display_all_reactions(const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    const std::vector<CreateDestructRxn>& createDestructRxns);
void display_all_MolTemplates(const std::vector<MolTemplate>& molTemplates);
/*******************/
