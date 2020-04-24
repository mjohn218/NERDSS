#include "reactions/shared_reaction_functions.hpp"

bool hasIntangibles(int reactantIndex, const Molecule& reactMol, const RxnBase::RateState& currRxnState)
{
    int matches { 0 }; // number of matching ancillary interfaces found

    for (const auto& anccIface : currRxnState.otherIfaceLists[reactantIndex]) {
        auto mol1Itr = find_if(reactMol.interfaceList.begin(), reactMol.interfaceList.end(),
            [&](const Molecule::Iface& oneIface) -> bool { return anccIface == oneIface; });

        // if the ancillary interface was found in either of the two reactants, iterate the number of matches
        // TODO: Does it matter if there are matches in both reactants?
        if (mol1Itr != reactMol.interfaceList.end())
            ++matches;
    }

    return matches == currRxnState.otherIfaceLists[reactantIndex].size();
}

bool hasIntangibles(int reactIndex1, int reactIndex2, const Molecule& reactMol1, const Molecule& reactMol2,
    const RxnBase::RateState& currRxnState)
{
    // if there are no ancillary ifaces, just return true
    unsigned long totalAnccIfaces { currRxnState.otherIfaceLists[0].size() + currRxnState.otherIfaceLists[1].size() };
    if (totalAnccIfaces == 0)
        return true;

    int matches { 0 }; // number of matching ancillary interfaces found
    for (const auto& anccIface : currRxnState.otherIfaceLists[reactIndex1]) {
        auto molItr = find_if(reactMol1.interfaceList.begin(), reactMol1.interfaceList.end(),
            [&](const Molecule::Iface& oneIface) -> bool { return anccIface == oneIface; });
        if (molItr != reactMol1.interfaceList.end())
            ++matches;
    }
    for (const auto& anccIface : currRxnState.otherIfaceLists[reactIndex2]) {
        auto molItr = find_if(reactMol2.interfaceList.begin(), reactMol2.interfaceList.end(),
            [&](const Molecule::Iface& oneIface) -> bool { return anccIface == oneIface; });
        if (molItr != reactMol2.interfaceList.end())
            ++matches;
    }

    return matches == totalAnccIfaces;
}
