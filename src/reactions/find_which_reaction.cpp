#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"

void find_which_reaction(int ifaceIndex1, int ifaceIndex2, int& rxnIndex, int& rateIndex, bool& isStateChangeBackRxn,
    const Interface::State& currState, const Molecule& reactMol1, const Molecule& reactMol2,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
    const std::vector<MolTemplate>& molTemplateList)
{
    // TRACE();
    for (auto rxnItr : currState.myForwardRxns) {
        const ForwardRxn& oneRxn = forwardRxns[rxnItr];
        // see if we can find both of the reactants in the reaction's reactantList
        // using iterators because I want to easily see if the reactMol{1,2} isn't in the reactantList
        int reactIndex1 { -1 };
        int reactIndex2 { -1 };
        for (int reactItr { 0 }; reactItr < oneRxn.reactantListNew.size(); ++reactItr) {
            if (reactMol1.isImplicitLipid && molTemplateList[oneRxn.reactantListNew[reactItr].molTypeIndex].isImplicitLipid) {
                if (reactIndex1 == -1) {
                    reactIndex1 = reactItr;
                    continue;
                }
            } else {
                if (reactMol1.interfaceList[ifaceIndex1].index == oneRxn.reactantListNew[reactItr].absIfaceIndex) {
                    if (reactIndex1 == -1) {
                        reactIndex1 = reactItr;
                        continue;
                    }
                }
            }
            if (reactMol2.isImplicitLipid && molTemplateList[oneRxn.reactantListNew[reactItr].molTypeIndex].isImplicitLipid) {
                if (reactIndex2 == -1) {
                    reactIndex2 = reactItr;
                    continue;
                }
            } else {
                if (reactMol2.interfaceList[ifaceIndex2].index == oneRxn.reactantListNew[reactItr].absIfaceIndex) {
                    if (reactIndex2 == -1) {
                        reactIndex2 = reactItr;
                        continue;
                    }
                }
            }
        }

        if (reactMol2.isImplicitLipid) {
            reactIndex2 = 1;
        }
        if (reactMol1.isImplicitLipid) {
            reactIndex1 = 1;
        }
	//std::cout <<" curr rxn: "<<oneRxn.relRxnIndex<<"reactIndex1 and 2: "<<reactIndex1<<' '<<reactIndex2<<std::endl;
	if(oneRxn.conjBackRxnIndex > 0){
	  if (oneRxn.rxnType == ReactionType::biMolStateChange && (reactIndex1 == -1 || reactIndex2 == -1)) {
            std::vector<std::array<int, 2>> matchList;
            /*Why is this allowed here, using the products?*/
	    for (int prodItr { 0 }; prodItr < oneRxn.productListNew.size(); ++prodItr) {
	      if (reactMol1.interfaceList[ifaceIndex1].index == oneRxn.productListNew[prodItr].absIfaceIndex) {
		reactIndex1 = prodItr;
	      }
	      if (reactMol2.interfaceList[ifaceIndex2].index == oneRxn.productListNew[prodItr].absIfaceIndex) {
		reactIndex2 = prodItr;
	      }
	      //std::cout <<"test bimolecular state change using: reactIndex 1 and 2: "<<reactIndex1 <<' '<<reactIndex2<<" print backRxnIndex: "<<oneRxn.conjBackRxnIndex<<std::endl;
	      if (reactIndex1 != -1 && reactIndex2 != -1 ){
		int matches { 0 };
		
		
		for (auto& oneRate : backRxns[oneRxn.conjBackRxnIndex].rateList) {
		  if (hasIntangibles(reactIndex1, reactIndex2, reactMol1, reactMol2, oneRate)) {
		    //                            matchList.emplace_back(static_cast<int>(&oneRate -
		    //                            &oneRxn.rateList[0]),
		    //                                (oneRate.otherIfaceLists[0].size() +
		    //                                oneRate.otherIfaceLists[1].size()));
		    std::array<int, 2> tmpArr { { static_cast<int>(&oneRate - &oneRxn.rateList[0]),
			  static_cast<int>(
					   oneRate.otherIfaceLists[0].size() + oneRate.otherIfaceLists[1].size()) } };
		    matchList.push_back(tmpArr);
		    ++matches;
		  }
		}
		if (matches == 1) {
		  rxnIndex = oneRxn.relRxnIndex;
		  isStateChangeBackRxn = true;
		  return;
		} else if (matches > 1) {
		  return;
		} else {
		  // if there are multiple matching rates, use the one with the most required ancillary interfaces
		  int bestFitIndex { 0 };
		  int numAncIfaces { 0 };
		  for (auto& match : matchList) {
		    if (match[1] > numAncIfaces) {
		      bestFitIndex = &match - &matchList[0];
		      numAncIfaces = match[1];
		    }
		  }
		  rateIndex = matchList[bestFitIndex][0];
		  rxnIndex = oneRxn.relRxnIndex;
		  isStateChangeBackRxn = true;
		  return;
		}
	      }
	    }
	  }
	}//This loop should not be attempted if there is no conjBackRxn
	//	std::cout <<"Did not assign in bimolecular state change loop "<<std::endl;
        if (reactIndex1 != -1 && reactIndex2 != -1) {
            int matches { 0 };
            std::vector<std::array<int, 2>> matchList;
            for (auto& oneRate : oneRxn.rateList) {
                if (hasIntangibles(reactIndex1, reactIndex2, reactMol1, reactMol2, oneRate)) {
                    std::array<int, 2> tmpArr { { static_cast<int>(&oneRate - &oneRxn.rateList[0]),
                        static_cast<int>(oneRate.otherIfaceLists[0].size() + oneRate.otherIfaceLists[1].size()) } };
                    matchList.push_back(tmpArr);
                    ++matches;
                }
            }
            if (matches == 1) {
                rateIndex = matchList[0][0];
                rxnIndex = oneRxn.relRxnIndex;
                ++totMatches;
                return;
            } else if (matches == 0) {
                // if there are no matching rates and the reaction symmetric (i.e. interface 1 binding to interface 1,
                // swap the reactants and check their ancillary interfaces
                if (oneRxn.rxnType == ReactionType::bimolecular
                    && (oneRxn.intReactantList[0] == oneRxn.intReactantList[1])) {
                    for (auto& oneRate : oneRxn.rateList) {
                        if (hasIntangibles(reactIndex2, reactIndex1, reactMol1, reactMol2, oneRate)) {
                            std::array<int, 2> tmpArr { { static_cast<int>(&oneRate - &oneRxn.rateList[0]),
                                static_cast<int>(
                                    oneRate.otherIfaceLists[0].size() + oneRate.otherIfaceLists[1].size()) } };
                            matchList.push_back(tmpArr);
                            ++matches;
                        }
                    }
                    if (matches == 1) {
                        rateIndex = matchList[0][0];
                        rxnIndex = oneRxn.relRxnIndex;
                        ++totMatches;
                        return;
                    } else if (matches == 0) {
                        return;
                    } else {
                        int bestFitIndex { 0 };
                        int numAncIfaces { 0 };
                        for (auto& match : matchList) {
                            if (match[1] > numAncIfaces) {
                                bestFitIndex = &match - &matchList[0];
                                numAncIfaces = match[1];
                            }
                        }
                        rateIndex = matchList[bestFitIndex][0];
                        rxnIndex = oneRxn.relRxnIndex;
                        ++totMatches;
                        return;
                    }
                }
                return;
            } else {
                // if there are multiple matching rates, use the one with the most required ancillary interfaces
                int bestFitIndex { 0 };
                int numAncIfaces { 0 };
                for (auto& match : matchList) {
                    if (match[1] > numAncIfaces) {
                        bestFitIndex = &match - &matchList[0];
                        numAncIfaces = match[1];
                    }
                }
                rateIndex = matchList[bestFitIndex][0];
                rxnIndex = oneRxn.relRxnIndex;
                ++totMatches;
                return;
            }
        }
    }
}
