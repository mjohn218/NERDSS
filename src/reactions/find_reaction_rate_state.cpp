#include "io/io.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"

int find_reaction_rate_state(int simItr, int relIfaceIndex1, int relIfaceIndex2, const Molecule& reactMol1,
    const Molecule& reactMol2, const BackRxn& backRxn, const std::vector<MolTemplate>& molTemplateList)
{
    // TRACE();
    int matches { 0 };
    std::vector<std::array<int, 2>> matchList {};
    int reactIndex1 { -1 };
    int reactIndex2 { -1 };
    for (int reactItr { 0 }; reactItr < backRxn.reactantListNew.size(); ++reactItr) {
        if (reactMol1.molTypeIndex == backRxn.reactantListNew[reactItr].molTypeIndex
            && relIfaceIndex1 == backRxn.reactantListNew[reactItr].relIfaceIndex) {
            if (reactIndex1 == -1) {
                reactIndex1 = reactItr;
                continue;
            }
        }
        if (reactMol2.molTypeIndex == backRxn.reactantListNew[reactItr].molTypeIndex
            && relIfaceIndex2 == backRxn.reactantListNew[reactItr].relIfaceIndex) {
            if (reactIndex2 == -1) {
                reactIndex2 = reactItr;
                continue;
            }
        }
    }

    if (reactIndex1 == -1 || reactIndex2 == -1) {
        std::cerr << llinebreak;
        std::cerr << "ERROR: Association product has incorrect dissociation reaction.\n";
        reactMol1.display(molTemplateList[reactMol1.molTypeIndex]);
        std::cerr << linebreak;
        reactMol2.display(molTemplateList[reactMol2.molTypeIndex]);
        std::cerr << llinebreak;
        backRxn.display();
        std::cout << std::endl;
        exit(1);
    }

    for (auto& oneRate : backRxn.rateList) {
        if (hasIntangibles(reactIndex1, reactIndex2, reactMol1, reactMol2, oneRate)) {
            // increment matches and add the rate state's index to the list of matches
            ++matches;
            //            matchList.emplace_back(static_cast<int>(&oneRate - &backRxn.rateList[0]),
            //            (oneRate.otherIfaceLists[0].size() + oneRate.otherIfaceLists[1].size()));
            std::array<int, 2> tmpArr { { static_cast<int>(&oneRate - &backRxn.rateList[0]),
                static_cast<int>(oneRate.otherIfaceLists[0].size() + oneRate.otherIfaceLists[1].size()) } };
            matchList.push_back(tmpArr);
        }
    }

    if (matches == 1) {
        return matchList.front()[0];
        //        std::cerr << "Error, no matching dissociation rate state was found for product " <<
        //        exit(1);
    } else if (matches == 0)
        return -1;
    else {
        // if there are multiple matches, see which fits best, e.g. the rate with the most required ancillary interfaces
        int bestFitIndex { 0 };
        int numAnccIfaces { 0 };
        for (auto& match : matchList) {
            if (match[1] > numAnccIfaces) {
                bestFitIndex = &match - &matchList[0];
                numAnccIfaces = match[1];
            }
        }
        return matchList[bestFitIndex][0];
    }

    return 0; // shouldn't ever get to this point
}
