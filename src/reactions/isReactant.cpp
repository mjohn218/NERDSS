#include "reactions/shared_reaction_functions.hpp"

bool isReactant(const Molecule::Iface& reactIface, const Molecule& reactMol, const RxnIface& tempReactant)
{
    if ((tempReactant.requiresInteraction && !reactIface.isBound)
        || (!tempReactant.requiresInteraction && reactIface.isBound))
        return false;

    return (reactMol.molTypeIndex == tempReactant.molTypeIndex) && (tempReactant.absIfaceIndex == reactIface.index)
        && (reactIface.stateIden == tempReactant.requiresState);
}

bool isReactant(const Molecule& currMol, const Complex& currCom, const CreateDestructRxn& currRxn,
    const std::vector<Molecule>& moleculeList)
{
    if (currRxn.reactantMolList.size() == 1) {
        if (currMol.molTypeIndex != currRxn.reactantMolList[0].molTypeIndex
            || currMol.interfaceList.size() != currRxn.reactantMolList[0].interfaceList.size()) {
            return false;
        }

        for (unsigned ifaceItr { 0 }; ifaceItr < currMol.interfaceList.size(); ++ifaceItr) {
            if (!isReactant(
                    currMol.interfaceList[ifaceItr], currMol, currRxn.reactantMolList[0].interfaceList[ifaceItr])) {
                // if the current Molecule's interfaces don't match up, it's not a reactant Molecule
                return false;
            }
        }
        return true;
    } else if (currRxn.reactantMolList.size() > 1 && currCom.memberList.size() > 1) {
        // the only way this happens is if the reactant is in a bound state
        bool isReactOne { currMol.molTypeIndex == currRxn.reactantMolList[0].molTypeIndex };
        bool isReactTwo { currMol.molTypeIndex == currRxn.reactantMolList[1].molTypeIndex };

        int currIndex { 0 };
        int otherIndex { 0 };
        // if the reactant is two of the same Molecule types bound together, just use the first reactant
        if (isReactOne || (isReactOne && isReactTwo)) { // TODO: check the condition here
            currIndex = 0;
            otherIndex = 1;
        } else if (isReactTwo)
            currIndex = 1;
        else
            return false;

        const CreateDestructRxn::CreateDestructMol& currReact = currRxn.reactantMolList[currIndex];
        const CreateDestructRxn::CreateDestructMol& otherReact = currRxn.reactantMolList[otherIndex];
        std::vector<unsigned> reactMatchList {}; // list of indices other other molecules to check for bonds
        for (unsigned ifaceItr { 0 }; ifaceItr < currMol.interfaceList.size(); ++ifaceItr) {
            if (currMol.interfaceList[ifaceItr].stateIden != currReact.interfaceList[ifaceItr].requiresState
                || currMol.interfaceList[ifaceItr].isBound != currReact.interfaceList[ifaceItr].requiresInteraction) {
                // if the current Molecule's interfaces don't match up, it's not a reactant Molecule
                return false;
            }

            // if the interface is bound, add the partner interface to list to check later
            // only add to list if the index of the partner is larger (to avoid double checking)
            if (currMol.interfaceList[ifaceItr].isBound
                && (currMol.index < currMol.interfaceList[ifaceItr].interaction.partnerIndex)
                && moleculeList[currMol.interfaceList[ifaceItr].interaction.partnerIndex].molTypeIndex
                    == otherReact.molTypeIndex) {
                reactMatchList.push_back(currMol.interfaceList[ifaceItr].interaction.partnerIndex);
            }
        }

        // now we need to check the other reactant's interfaces to see if they match
        for (auto matchIndex : reactMatchList) {
            unsigned matchingIfaces { 0 };
            for (unsigned ifaceItr { 0 }; ifaceItr < moleculeList[matchIndex].interfaceList.size(); ++ifaceItr) {
                if (moleculeList[matchIndex].interfaceList[ifaceItr].stateIden
                        == otherReact.interfaceList[ifaceItr].requiresState
                    && moleculeList[matchIndex].interfaceList[ifaceItr].isBound
                        == otherReact.interfaceList[ifaceItr].requiresInteraction) {
                    ++matchingIfaces;
                }
            }

            if (matchingIfaces == moleculeList[matchIndex].interfaceList.size())
                return true; // if there's a single match, break out and return true
        }
    }

    return false; // default state is true
}
