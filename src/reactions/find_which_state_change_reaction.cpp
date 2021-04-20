#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"

void find_which_state_change_reaction(int ifaceIndex, int& rxnIndex, int& rateIndex, bool& isStateChangeBackRxn,
    const Molecule& reactMol, const Interface::State& currState, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns)
{
    // TRACE();
    for (auto rxnItr : currState.stateChangeRxns) {
        const ForwardRxn& oneRxn = forwardRxns[rxnItr.first]; // rxnItr.first is the reaction's index in forwardRxns
        if (oneRxn.rxnType == ReactionType::uniMolStateChange) {
            if (reactMol.interfaceList[ifaceIndex].index == oneRxn.reactantListNew[0].absIfaceIndex) {
                int matches { 0 };
                for (const auto& oneRate : oneRxn.rateList) {
                    if (hasIntangibles(0, reactMol, oneRate)) {
                        ++matches;
                        rateIndex = static_cast<int>(&oneRate - &oneRxn.rateList[0]);
                    }
                }

                if (matches == 0) {
                    // TODO: either an error or the molecule doesn't match all the requirements. do what?
                } else if (matches > 1) {
                    // std::cout << "Multiple rates found, using last found rate.\n";
                    rxnIndex = oneRxn.relRxnIndex;
                } else {
                    rxnIndex = oneRxn.relRxnIndex;
                }
            } else if (reactMol.interfaceList[ifaceIndex].index == oneRxn.productListNew[0].absIfaceIndex) {
                int matches { 0 };
                for (const auto& oneRate : backRxns[oneRxn.conjBackRxnIndex].rateList) {
                    if (hasIntangibles(0, reactMol, oneRate)) {
                        ++matches;
                        rateIndex = static_cast<int>(&oneRate - &backRxns[oneRxn.conjBackRxnIndex].rateList[0]);
                    }
                }

                if (matches == 0) {
                    // TODO: either an error or the molecule doesn't match all the requirements. do what?
                } else if (matches > 1) {
                    // std::cout << "Multiple rates found, using last found rate.\n";
                    rateIndex = backRxns[oneRxn.conjBackRxnIndex].relRxnIndex;
                    isStateChangeBackRxn = true;
                } else {
                    rateIndex = backRxns[oneRxn.conjBackRxnIndex].relRxnIndex;
                    isStateChangeBackRxn = true;
                }
            } else {
                continue;
            }
        }
    }
}
