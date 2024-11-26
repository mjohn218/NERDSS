#include "reactions/shared_reaction_functions.hpp"

bool hasIntangibles(int reactantIndex, const Molecule& reactMol,
                    const RxnBase::RateState& currRxnState) {
  for (const auto& anccIface : currRxnState.otherIfaceLists[reactantIndex]) {
    // auto mol1Itr = find_if(reactMol.interfaceList.begin(),
    // reactMol.interfaceList.end(),
    //     [&](const Molecule::Iface& oneIface) -> bool { return anccIface ==
    //     oneIface; });

    // if the ancillary interface was found in either of the two reactants,
    // iterate the number of matches
    // TODO: Does it matter if there are matches in both reactants?
    // if (mol1Itr == reactMol.interfaceList.end())
    //     return false;

    bool foundFlag = false;
    for (const auto& oneIface : reactMol.interfaceList) {
      if (anccIface.molTypeIndex == oneIface.molTypeIndex &&
          anccIface.relIfaceIndex == oneIface.relIndex &&
          anccIface.requiresInteraction == oneIface.isBound &&
          anccIface.requiresState == oneIface.stateIden) {
        foundFlag = true;
        break;
      }
    }
    if (foundFlag == false) {
      return false;
    }
  }
  return true;
}

bool hasIntangibles(int reactIndex1, int reactIndex2, const Molecule& reactMol1,
                    const Molecule& reactMol2,
                    const RxnBase::RateState& currRxnState) {
  if (currRxnState.otherIfaceLists[0].empty() &&
      currRxnState.otherIfaceLists[1].empty()) {
    return true;
  }

  for (const auto& anccIface : currRxnState.otherIfaceLists[reactIndex1]) {
    // auto molItr = find_if(reactMol1.interfaceList.begin(),
    // reactMol1.interfaceList.end(),
    //     [&](const Molecule::Iface& oneIface) -> bool { return anccIface ==
    //     oneIface; });
    // if (molItr == reactMol1.interfaceList.end())
    //     return false;

    bool foundFlag = false;
    for (const auto& oneIface : reactMol1.interfaceList) {
      if (anccIface.molTypeIndex == oneIface.molTypeIndex &&
          anccIface.relIfaceIndex == oneIface.relIndex &&
          anccIface.requiresInteraction == oneIface.isBound &&
          anccIface.requiresState == oneIface.stateIden) {
        foundFlag = true;
        break;
      }
    }
    if (foundFlag == false) {
      return false;
    }
  }
  for (const auto& anccIface : currRxnState.otherIfaceLists[reactIndex2]) {
    // auto molItr = find_if(reactMol2.interfaceList.begin(),
    // reactMol2.interfaceList.end(),
    //     [&](const Molecule::Iface& oneIface) -> bool { return anccIface ==
    //     oneIface; });
    // if (molItr == reactMol2.interfaceList.end())
    //     return false;

    bool foundFlag = false;
    for (const auto& oneIface : reactMol2.interfaceList) {
      if (anccIface.molTypeIndex == oneIface.molTypeIndex &&
          anccIface.relIfaceIndex == oneIface.relIndex &&
          anccIface.requiresInteraction == oneIface.isBound &&
          anccIface.requiresState == oneIface.stateIden) {
        foundFlag = true;
        break;
      }
    }
    if (foundFlag == false) {
      return false;
    }
  }

  return true;
}
