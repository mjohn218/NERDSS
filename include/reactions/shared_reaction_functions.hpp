/*! \file shared_reaction_functions.hpp

 * ### Created on 11/6/18 by Matthew Varga
 * ### Purpose
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */
#pragma once

#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"
#include "gsl/gsl_matrix.h"

#include <algorithm>

/*!
 * \brief Determines if the Interface::State of the reactant is equivalent to the reactant as contained in the
 * reaction
 *
 * \param reactant Interface::State of the reacting interface
 * \param tempReactant a reactant of the ForwardRxn/BackRxn/CreateDestructRxn being evaluated
 */
bool isReactant(const Molecule::Iface& reactIface, const Molecule& reactMol, const RxnIface& tempReactant);

/*!
 * \brief Determines if the Molecule is a reactant for a Creation/Destruction reaction.
 */
bool isReactant(const Molecule& currMol, const Complex& currCom, const CreateDestructRxn& currRxn,
    const std::vector<Molecule>& moleculeList);

/*!
 * \brief Determines if the associating Molecules have the ancillary interfaces required by this reaction.
 */
bool hasIntangibles(int reactIndex1, int reactIndex2, const Molecule& reactMol1, const Molecule& reactMol2,
    const RxnBase::RateState& currRxnState);

/*!
 * \brief Determines if the Molecule in a unimolecular reaction has he ancillary interaces required
 */
bool hasIntangibles(int reactantIndex, const Molecule& reactMol, const RxnBase::RateState& currRxnState);

/* FUNCTIONS TO DETERMINE WHICH REACTION TO PERFORM */

/*!
 * \brief Determines which RxnState to use for the association/dissociation reaction.
 */
int find_reaction_rate_state(int simItr, int relIfaceIndex1, int relIfaceIndex2, const Molecule& reactMol1,
    const Molecule& reactMol2, const BackRxn& backRxn,
    const std::vector<MolTemplate>& molTemplateList);

size_t find_reaction_rate_state(const Molecule& reactMol1, const Molecule& reactMol2, const ForwardRxn& forwardRxn,
    const std::vector<MolTemplate>& molTemplateList);

/*!
 * \brief This function determines which reaction to use based on the identities of the reacting Molecules and
 * interfaces, and returns its reaction index and rate index
 *
 * \param[out] std::array<size_t, 2> array of two indices, [index of ForwardRxn, index of RxnBase::RateState]
 */
void find_which_reaction(int ifaceIndex1, int ifaceIndex2, int& rxnIndex, int& rateIndex, bool& isStateChangeBackRxn,
    const Interface::State& currState, const Molecule& reactMol1, const Molecule& reactMol2,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    const std::vector<MolTemplate>& molTemplateList);

/*!
 * \brief Determines which state change reaction to use based on the identity of the current Interface::State of the
 * target Interface.
 *
 * \params[out] std::array<int, 3> [rxnIndex, rateIndex, rxnType], where rxnType = 0 if forward, 1 if back
 */
void find_which_state_change_reaction(int ifaceIndex, int& rxnIndex, int& rateIndex, bool& isStateChangeBackRxn,
    const Molecule& reactMol, const Interface::State& currState,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns);

std::array<int, 2> find_which_reaction(int ifaceIndex, const Molecule& reactMol, const Interface::State& currState,
    const std::vector<BackRxn>& forwardRxns);

/*!
 * \brief
 */
// double passocF(double r0, double tCurr, double Dtot, double bindRadius, double alpha, double cof);

extern int evalBindNum;
/*!
 * \brief Main function for evaluating the potential interactions between two Molecules.
 *
 * \param[in] pro1Index Index of first protein in moleculeList.
 * \param[in] pro2Index Index of second protein in moleculeList.
 * \param[in] simItr Current simulation iteration.
 * \param[in] tableIDs C-style array containing the rates and their total diffusion constant.
 * \param[in] DDTableIndex Current index in tableIDs
 * \param[in] params User-provided simulation Parameters.
 * \param[in] normMatrices List of all previously calculated normMatrix
 * \param[in] survMatrices List of all previously calculated survMatrix
 * \param[in] pirMatrices List of all previously calculated pirMatrix
 * \param[in] moleculeList List of all Molecules in the system.
 * \param[in] complexList List of all Complexes in the system.
 * \param[in] molTemplateList List of all user-provided MolTemplates.
 * \param[in] forwardRxns List of all user-provided ForwardRxns.
 *
 * If 2D, also calculates/looks up values for 2D reaction tables (normMatrix, survMatrix, and pirMatrix)
 */
void check_bimolecular_reactions(int pro1Index, int pro2Index, int simItr, double* tableIDs, unsigned& DDTableIndex,
    const Parameters& params, std::vector<gsl_matrix*>& normMatrices, std::vector<gsl_matrix*>& survMatrices,
    std::vector<gsl_matrix*>& pirMatrices, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, copyCounters& counterArrays, Membrane& membraneObject);

/*!
 * \brief Determines if binding of two molecules within the same complex can occur.
 *
 * If the two molecules are within Rmax, probability of reaction is set to 1.0. If not, it is set to 0.
 */

void evaluate_binding_within_complex(int pro1Index, int pro2Index, int iface1Index, int iface2Index, int rxnIndex,
    int rateIndex, bool isBiMolStateChange, const Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const ForwardRxn& oneRxn,
    const std::vector<BackRxn>& backRxns, Membrane& membraneObject, copyCounters& counterArrays);

bool determine_if_reaction_occurs(int& crossIndex1, int& crossIndex2, const double maxRandInt, Molecule& mol,
    std::vector<Molecule>& moleculeList, const std::vector<ForwardRxn>& forwardRxns);

void update_Nboundpairs(int ptype1, int ptype2, int chg, const Parameters& params, copyCounters& counterArrays);

void check_implicit_reactions(int pro1Index, int pro2Index, int simItr,
    const Parameters& params, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, copyCounters& counterArrays, Membrane& membraneObject, std::vector<double>& IL2DbindingVec, std::vector<double>& IL2DUnbindingVec, std::vector<double>& ILTableIDs);

