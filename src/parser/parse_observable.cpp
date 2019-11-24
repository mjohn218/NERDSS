#include "parser/parser_functions.hpp"

#include <sstream>

SpeciesTracker::Observable parse_observable(const std::string& line, const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns, const std::vector<CreateDestructRxn>& createDestructRxns)
{
    std::map<const std::string, SpeciesTracker::ObservableType> observableTypes
        = { { "molecule", SpeciesTracker::ObservableType::molecule },
              { "complex", SpeciesTracker::ObservableType::complex } };
    // set up some temporary variables
    SpeciesTracker::Observable tmpObs;
    std::string specie;
    {
        std::string obsTypeStr;
        std::stringstream lineStream { line }; // to read in a formatted manner

        // read Observable information
        lineStream >> obsTypeStr >> tmpObs.name >> specie;

        // determine the Observable type
        std::transform(obsTypeStr.begin(), obsTypeStr.end(), obsTypeStr.begin(), ::tolower);
        auto obsTypeItr = observableTypes.find(obsTypeStr);
        if (obsTypeItr == observableTypes.end()) {
            std::cerr << "FATAL ERROR: Observable type " << obsTypeStr << " unknown. Exiting.\n";
            exit(1);
        } else {
            tmpObs.observableType = obsTypeItr->second;
        }
    }

    // parse Observable BNGL
    std::vector<ParsedMol> parsedMols;
    // split into Molecules based on period (.)
    {
        size_t position { 0 };
        std::string tmpSpecie { specie }; // so i don't erase the real thing
        std::vector<std::string> molecules;
        while ((position = tmpSpecie.find('.')) != std::string::npos) {
            molecules.emplace_back(tmpSpecie.substr(0, position));
            tmpSpecie.erase(0, position + 1);
        }
        molecules.emplace_back(tmpSpecie.substr(0, std::string::npos)); // so it actually includes the last molecule

        for (auto& oneMol : molecules) {
            std::pair<std::string, int> tmpMol { oneMol, 0 };
            int tmpI { 1 }; // just so i don't have to make another function
            parsedMols.emplace_back(parse_molecule_bngl(tmpI, false, tmpMol));
        }
    }

    bool molFound { false };
    for (auto& parsedObs : parsedMols) {
        for (const auto& oneTemp : molTemplateList) {
            if (oneTemp.molName == parsedObs.molName) {
                parsedObs.molTypeIndex = oneTemp.molTypeIndex;

                bool ifaceFound { false };
                for (auto& obsIface : parsedObs.interfaceList) {
                    for (const auto& tempIface : oneTemp.interfaceList) {
                        if (tempIface.name == obsIface.ifaceName) {
                            obsIface.relIndex = tempIface.index;
                            for (auto& tempState : tempIface.stateList) {
                                if (tempState.iden == obsIface.state) {
                                    obsIface.absIndex = tempState.index;
                                    break;
                                }
                            }
                            ifaceFound = true;
                        }
                        if (ifaceFound)
                            break;
                    }
                    ifaceFound = false;
                }
                molFound = true;
            }
            if (molFound)
                break;
        }
        molFound = false;
    }

    if (tmpObs.observableType == SpeciesTracker::ObservableType::complex) {
        for (auto& oneRxn : forwardRxns) {
            if (oneRxn.rxnType == ReactionType::bimolecular) {
                // if the molTypeIndex of the constituents and products match, go on
                bool firstObsfirstProd { parsedMols[0].molTypeIndex == oneRxn.productListNew[0].molTypeIndex
                    && parsedMols[0].interfaceList[0].relIndex == oneRxn.productListNew[0].relIfaceIndex
                    && parsedMols[0].interfaceList[0].isBound == oneRxn.productListNew[0].requiresInteraction };
                bool firstObsSecondProd { parsedMols[0].molTypeIndex == oneRxn.productListNew[1].molTypeIndex
                    && parsedMols[0].interfaceList[0].relIndex == oneRxn.productListNew[1].relIfaceIndex
                    && parsedMols[0].interfaceList[0].isBound == oneRxn.productListNew[1].requiresInteraction };
                bool secondObsfirstProd { parsedMols[1].molTypeIndex == oneRxn.productListNew[0].molTypeIndex
                    && parsedMols[1].interfaceList[0].relIndex == oneRxn.productListNew[0].relIfaceIndex
                    && parsedMols[1].interfaceList[0].isBound == oneRxn.productListNew[0].requiresInteraction };
                bool secondObsSecondProd { parsedMols[1].molTypeIndex == oneRxn.productListNew[1].molTypeIndex
                    && parsedMols[1].interfaceList[0].relIndex == oneRxn.productListNew[1].relIfaceIndex
                    && parsedMols[1].interfaceList[0].isBound == oneRxn.productListNew[1].requiresInteraction };

                if (oneRxn.isSymmetric && firstObsfirstProd && firstObsSecondProd && secondObsfirstProd
                    && secondObsSecondProd) {
                }
            }
        }
    }

    for (auto& parsedMol : parsedMols) {
        tmpObs.constituentList.emplace_back();
        SpeciesTracker::Observable::Constituent& newConst = tmpObs.constituentList.back();
        newConst.molTypeIndex = parsedMol.molTypeIndex;
        for (auto& parsedIface : parsedMol.interfaceList)
            newConst.interfaceList.emplace_back(parsedIface.relIndex, parsedIface.absIndex, parsedIface.state,
                parsedIface.isBound, parsedIface.bondIndex);
    }
    return tmpObs;
}
