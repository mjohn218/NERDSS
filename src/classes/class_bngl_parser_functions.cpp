/*! @file bngl_parser_member_functions.cpp
 * \ingroup Parser
 * \brief Class member functions for reaction file parsing
 *
 * ### Created on 5/15/18 by Matthew Varga
 * ### Purpose
 * ***
 * Contains class member functions for classes ParsedMol and ParsedRxn found in bngl_parsing_class.hpp
 *
 * ### TODO List
 * ***
 *  - TODO (PRIORITY MED): Function to do final checking of ParsedRxn before conversion into
 * ForwardRxn(s)/CreateDestructRxn
 */

#include "classes/class_Rxns.hpp"
#include "io/io.hpp"
#include "parser/parser_functions.hpp"

#include <iomanip>

std::ostream& operator<<(std::ostream& os, const Involvement& involve)
{
    int iType { static_cast<std::underlying_type<Involvement>::type>(involve) };
    switch (iType) {
    case 0: {
        return os << "Ancillary";
    }
    case 1: {
        return os << "Possibly involved";
    }
    case 2: {
        return os << "Interim involvement";
    }
    case 3: {
        return os << "Only the Interaction changes";
    }
    case 4: {
        return os << "Only the State changes";
    }
    case 5: {
        return os << "Both the Interaction and State changes";
    }
    default:
        //        involvement_err(__func__, __LINE__);
        exit(1);
    }
}

std::ostream& operator<<(std::ostream& os, const ReactionType& rxnType)
{
    // scoped enumerations are not implicilty convertible. must be explicitly converted
    int rType { static_cast<std::underlying_type<ReactionType>::type>(rxnType) };
    try {
        switch (rType) {
        case 0: {
            return os << "Unimolecular state change";
        }
        case 1: {
            return os << "Bimolecular association";
        }
        case 2: {
            return os << "Bimolecular state change";
        }
        case 3: {
            return os << "Creation from concentration";
        }
        case 4: {
            return os << "Destruction";
        }
        case 5: {
            return os << "Creation from molecule";
        }
        default:
            throw std::invalid_argument("Invalid reaction type");
        }
    } catch (const std::invalid_argument& e) {
        std::cout << e.what() << '\n';
        exit(1);
    }
}

bool operator==(const std::vector<RxnIface>& rxnIfaceList, const std::vector<ParsedMol::IfaceInfo>& parsedIfaceList)
{
    if (parsedIfaceList.size() != rxnIfaceList.size())
        return false;

    unsigned numMatches { 0 };
    for (unsigned rxnItr { 0 }; rxnItr < rxnIfaceList.size(); ++rxnItr) {
        for (unsigned parsedItr { 0 }; parsedItr < parsedIfaceList.size(); ++parsedItr) {
            if (rxnIfaceList[rxnItr].requiresInteraction == parsedIfaceList[parsedItr].isBound
                && rxnIfaceList[rxnItr].requiresState == parsedIfaceList[parsedItr].state
                && rxnIfaceList[rxnItr].molTypeIndex == parsedIfaceList[parsedItr].molTypeIndex) {
                ++numMatches;
                break;
            }
        }
    }

    return (numMatches == parsedIfaceList.size());
}

// I am so sick of writing this out, so here's a type alias
using parsedIface = std::pair<const std::string, ParsedMol::IfaceInfo>;

/* CONSTRUCTORS */
ParsedMol::IfaceInfo::IfaceInfo(
    std::string _ifaceName, char _state, bool _isBound, Involvement _ifaceRxnStatus, int speciesIndex)
    : ifaceName(_ifaceName)
    , state(_state)
    , isBound(_isBound)
    , ifaceRxnStatus(_ifaceRxnStatus)
    , speciesIndex(speciesIndex)
{
    if (_state != '\0')
        ifaceAndStateName = _ifaceName + "~" + _state;
    else
        ifaceAndStateName = ifaceName;
}

ParsedMol::IfaceInfo::IfaceInfo(
    std::string _ifaceName, char _state, bool _isBound, int _bondIndex, Involvement _ifaceRxnStatus)
    : ifaceName(_ifaceName)
    , state(_state)
    , isBound(_isBound)
    , bondIndex(_bondIndex)
    , ifaceRxnStatus(_ifaceRxnStatus)
{
    if (_state != '\0')
        ifaceAndStateName = ifaceName + "~" + _state;
    else
        ifaceAndStateName = ifaceName;
}

ParsedMol::IfaceInfo::IfaceInfo(
    std::string _ifaceName, char _state, bool _isBound, int _bondIndex, Involvement _ifaceRxnStatus, int speciesIndex)
    : ifaceName(_ifaceName)
    , state(_state)
    , isBound(_isBound)
    , bondIndex(_bondIndex)
    , ifaceRxnStatus(_ifaceRxnStatus)
    , speciesIndex(speciesIndex)
{
    if (_state != '\0')
        ifaceAndStateName = _ifaceName + "~" + _state;
    else
        ifaceAndStateName = ifaceName;
}

bool ParsedMol::IfaceInfo::operator!=(const ParsedMol::IfaceInfo& iface) const
{
    return ((bool { this->absIndex != iface.absIndex }) || (bool { this->state != iface.state })
        || (bool { this->isBound != iface.isBound }) || (bool { this->bondIndex != iface.bondIndex }));
}

bool ParsedMol::IfaceInfo::operator==(const ParsedMol::IfaceInfo& iface) const
{
    return ((bool { this->absIndex == iface.absIndex }) && (bool { this->state == iface.state })
        && (bool { this->isBound == iface.isBound }) && (bool { this->bondIndex == iface.bondIndex }));
}

bool ParsedMol::IfaceInfo::operator==(const RxnIface& rxnIface) const
{
    return (bool { this->absIndex == rxnIface.absIfaceIndex } && bool { this->state == rxnIface.requiresState }
        && bool { this->isBound == rxnIface.requiresInteraction });
}

std::ostream& operator<<(std::ostream& os, const ParsedMol::IfaceInfo& oneInfo)
{
    os << oneInfo.absIndex << ", ";
    if (oneInfo.state != '\0')
        os << oneInfo.state << ", ";
    else
        os << "NO STATE, ";
    os << std::boolalpha << oneInfo.isBound << ", " << oneInfo.bondIndex << ", " << oneInfo.ifaceRxnStatus;
    return os;
}

bool ParsedMol::IfaceInfo::operator<(const ParsedMol::IfaceInfo& rhs) const { return speciesIndex < rhs.speciesIndex; }

bool ParsedMol::IfaceInfo::operator>(const ParsedMol::IfaceInfo& rhs) const { return rhs < *this; }

void ParsedMol::IfaceInfo::change_ifaceRxnStatus(int newIndex, Involvement newRxnStatus)
{
    this->indexFound = true;
    this->absIndex = newIndex;
    this->ifaceRxnStatus = newRxnStatus;
}

/* ParsedMol */
// Constructors
ParsedMol::ParsedMol(const MolTemplate& oneTemp)
{
    molTypeIndex = oneTemp.molTypeIndex;
    molName = oneTemp.molName;
    for (const auto& tempIface : oneTemp.interfaceList)
        interfaceList.emplace_back(tempIface);
}

bool ParsedMol::operator==(const ParsedMol& mol1) const
{
    if ((this->molName != mol1.molName) || (this->interfaceList.size() != mol1.interfaceList.size()))
        return false;
    for (auto& thisTargetIface : this->interfaceList) {
        // check each iface on this mol to see if it exists on mol1
        //        auto ifacePos { mol1.interfaceList.find(thisTargetIface.first) };
        auto ifacePos = std::find_if(mol1.interfaceList.begin(), mol1.interfaceList.end(),
            [&](const ParsedMol::IfaceInfo& iface) { return iface.ifaceName == thisTargetIface.ifaceName; });
        if (ifacePos != mol1.interfaceList.end())
            return thisTargetIface == *ifacePos;
    }
    return false;
}

std::ostream& ParsedMol::display_full_name(std::ostream& os) const
{
    os << this->molName << '(';
    size_t ifaceSizeIterator { 0 };
    for (auto& iface : interfaceList) {
        if (ifaceSizeIterator != interfaceList.size() - 1)
            os << iface.ifaceName << ", ";
        else
            os << iface.ifaceName << ']';
    }
    return os;
}

void ParsedMol::display() const
{
    std::cout << "Parsed Mol: " << molName << '\n';
    for (auto& iface : interfaceList)
        std::cout << '[' << iface.ifaceName << ":{" << iface << "}]\n";
}

void ParsedMol::set_molTypeIndex(const std::vector<MolTemplate>& molTemplateList)
{
    for (auto& oneTemp : molTemplateList) {
        if (oneTemp.molName == this->molName) {
            this->molTypeIndex = static_cast<int>(oneTemp.molTypeIndex);
            break;
        }
    }

    if (this->molTypeIndex == -1) {
        std::cerr << "Error, molecule " << this->molName << " not found. Exiting.\n";
        exit(1);
    }
}

void ParsedRxn::assemble_reactions(std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<CreateDestructRxn>& createDestructRxns, const std::vector<MolTemplate>& molTemplateList)
{

    if (this->rxnType == ReactionType::bimolecular || this->rxnType == ReactionType::uniMolStateChange
        || this->rxnType == ReactionType::biMolStateChange) {
        forwardRxns.emplace_back(*this, molTemplateList);
        forwardRxns.back().relRxnIndex = forwardRxns.size() - 1;
        // if the reaction is reversible, create the conjugate BackRxn
        if (forwardRxns.back().isReversible) {
            backRxns.emplace_back(this->offRatekb, forwardRxns.back());
            backRxns.back().relRxnIndex = backRxns.size() - 1;
            create_conjugate_reaction_itrs(forwardRxns, backRxns);
        } else {
            forwardRxns.back().conjBackRxnIndex = -1;
        }
    } else {
        createDestructRxns.emplace_back(*this, molTemplateList);
        createDestructRxns.back().relRxnIndex = createDestructRxns.size() - 1;
    }
}

/* ParsedRxn */
std::pair<bool, std::string> ParsedRxn::isComplete(const std::vector<MolTemplate>& molTemplateList)
{
    if (std::isnan(onRate3Dka))
        return std::make_pair(false, "onRate3Dka is not set for reaction");

    if (rxnType == ReactionType::bimolecular) {
        if (isReversible && std::isnan(offRatekb))
            return std::make_pair(false, "offRatekb is not set for reaction");
        if (std::isnan(bindRadius))
            return std::make_pair(false, "Sigma is not defined");
        if (!molTemplateList[reactantList[0].molTypeIndex].isPoint && !molTemplateList[reactantList[1].molTypeIndex].isPoint) {
            if (!molTemplateList.at(reactantList[0].molTypeIndex).isRod
                && (norm1.magnitude == 0 || std::isnan(norm1.magnitude))) {
                return std::make_pair(false, "Norm1 is not set or is not a vector for reaction");
            }
            if (!molTemplateList.at(reactantList[1].molTypeIndex).isRod
                && (norm2.magnitude == 0 || std::isnan(norm2.magnitude))) {
                return std::make_pair(false, "Norm2 is not set or is not a vector for reaction");
            }
        }
    }

    return std::make_pair(true, "Parameter set complete for reaction");
}

void ParsedRxn::display_angles() const
{
    std::cout << std::setw(10) << std::left << "Theta 1" << std::setw(7) << std::right << assocAngles.theta1 << '\n';
    std::cout << std::setw(10) << std::left << "Theta 2" << std::setw(7) << std::right << assocAngles.theta2 << '\n';
    std::cout << std::setw(10) << std::left << "Phi 1" << std::setw(7) << std::right << assocAngles.phi1 << '\n';
    std::cout << std::setw(10) << std::left << "Phi 2" << std::setw(7) << std::right << assocAngles.phi2 << '\n';
    std::cout << std::setw(10) << std::left << "Omega" << std::setw(7) << std::right << assocAngles.omega << '\n';
}

void ParsedRxn::set_value(std::string& line, RxnKeyword rxnKeyword)
{
    int key = static_cast<std::underlying_type<RxnKeyword>::type>(rxnKeyword);
    try {
        switch (key) {
        case 0: {
            onRate3Dka = std::stod(line);
            std::cout << "Read in value of onRate3Dka: " << onRate3Dka << "nm^3us^-1\n";
            break;
        }
        case 1: {
            onRate3DMacro = std::stod(line);
            std::cout << "Read in value of onRate3DMacro: " << onRate3DMacro << "uM^-1us^-1\n";
            break;
        }
        case 2: {
            offRatekb = std::stod(line);
            std::cout << "Read in value of offRatekb: " << offRatekb << "s^-1\n";
            break;
        }
        case 3: {
            offRateMacro = std::stod(line);
            std::cout << "Read in value of offRateMacro: " << offRateMacro << "s^-1\n";
            break;
        }
        case 4: {
            norm1 = Vector { parse_input_array(line) };
            norm1.calc_magnitude();
            std::cout << "Read in value of norm1: " << norm1 << '\n';
            break;
        }
        case 5: {
            norm2 = Vector { parse_input_array(line) };
            norm2.calc_magnitude();
            std::cout << "Read in value of norm2: " << norm2 << '\n';
            break;
        }
        case 6: {
            bindRadius = std::stod(line);
            std::cout << "Read in value of sigma: " << bindRadius << "nm\n";
            break;
        }
        case 7: {
            assocAngles = Angles { parse_input_array(line) };
            std::cout << "Read in value of assocAngles: [theta1: " << assocAngles.theta1 << ", theta2: " << assocAngles.theta2 << ", phi1: " << assocAngles.phi1 << ", phi2: " << assocAngles.phi2 << ", omega: " << assocAngles.omega << "]" << '\n';
            break;
        }
        case 8: {
            isOnMem = read_boolean(line);
            break;
        }
        case 9: {
            onRate3Dka = std::stod(line);
            std::cout << "Read in value of rate: " << onRate3Dka << '\n';
            break;
        }
        case 10: {
            isCoupled = true;
            coupledRxn = CoupledRxn { std::stoi(line) };
            std::cout << "Read in value of coupledRxn: absRxnIndex, " << coupledRxn.absRxnIndex << '\n';
            break;
        }
        case 11: {
            isObserved = true;
            break;
        }
        case 12: {
            isObserved = true;
            observeLabel = line;
            std::cout << "Read in value of observeLabel: " << observeLabel << '\n';
            break;
        }
        case 13: {
            bindRadSameCom = std::stod(line);
            std::cout << "Read in value of bindRadSameCom: " << bindRadSameCom << '\n';
            break;
        }
        case 14: {
            irrevRingClosure = read_boolean(line);
            break;
        }
        case 15: {
            creationRadius = std::stod(line);
            break;
        }
        case 16: {
            loopCoopFactor = std::stod(line);
            std::cout << "Read in value of loopCoopFactor: " << loopCoopFactor << '\n';
            break;
        }
        case 17: {
            length3Dto2D = std::stod(line);
            std::cout << "Read in value of length3Dto2D: " << length3Dto2D << "nm\n";
            break;
        }
        case 18: {
            rxnLabel = line;
            std::cout << "Read in value of rxnLabel: " << rxnLabel << '\n';
            break;
        }
        case 19: {
            isCoupled = true;
            coupledRxn = CoupledRxn { line };
            std::cout << "Read in value of coupledRxn: label, " << coupledRxn.label << '\n';
            break;
        }
        case 20: {
            kcat = std::stod(line);
            std::cout << "Read in value of kcat: " << kcat << '\n';
            break;
        }
        case 21: {
            excludeVolumeBound = read_boolean(line);
            break;
        }
        default: {
            throw std::invalid_argument("Not a valid keyword.");
        }
        }
    } catch (std::invalid_argument& e) {
        std::cout << e.what() << '\n';
        exit(1);
    }
}

void ParsedRxn::display() const
{
    /*THIS ROUTINE IS DUPLICATED IN CLASS_RXNS>CPP, FORWARDRXN::DISPLAY, and OTHER RXN TYPES, THIS IS NOT CALLED!*/
    std::cout << "reactantList:\n"
              << std::setw(10) << std::setfill('-') << ' ' << std::setfill(' ') << '\n';
    for (auto& oneReactant : reactantList) {
        oneReactant.display();
    }
    std::cout << "productList:\n"
              << std::setw(10) << std::setfill('-') << ' ' << std::setfill(' ') << '\n';
    for (auto& oneProduct : productList) {
        oneProduct.display();
    }

    if (hasStateChange) {
        std::cout << "Interface " << stateChangeIface.first.ifaceName << " on molecule "
                  << stateChangeIface.first.molTypeIndex << " changes state from "
                  << stateChangeIface.first.requiresState << " to " << stateChangeIface.second.requiresState << '\n';
    }

    if (rxnType == ReactionType::bimolecular) {
        std::cout << std::setw(10) << std::setfill('-') << ' ' << std::setfill(' ') << "\nAssociation Angles:\n";
        display_angles();
        std::cout << "Association sigma vector:\n";
        std::cout << "Sigma: " << bindRadius << '\n';
        std::cout << "Reactant 1 normal: " << norm1 << '\n';
        std::cout << "Reactant 2 normal: " << norm2 << '\n';
    }
    std::cout << "\nOn membrane? " << std::boolalpha << isOnMem << '\n';
    std::cout << "bindRadSameCom " << bindRadSameCom << '\n';
    std::cout << "loopCoopFactor " << loopCoopFactor << '\n';
    std::cout << "length3Dto2D " << length3Dto2D << '\n';
    std::cout << "isCoupled? " << isCoupled << '\n';
    if (isCoupled)
        std::cout << " coupledRxn Number: " << coupledRxn.absRxnIndex << " type: " << coupledRxn.rxnType << '\n';
    std::cout << "microRate3D: " << onRate3Dka << '\n';
    std::cout << "macroRate3D: " << onRate3DMacro << '\n';
    if (isReversible)
        std::cout << "micro Off Rate: " << offRatekb << '\n';
    std::cout << "macro Off Rate: " << offRateMacro << '\n';
}

void ParsedRxn::check_previous_bound_states(int& totSpecies, const std::vector<ForwardRxn>& forwardRxns, const std::vector<MolTemplate>& molTemplateList)
{
    for (auto& oneRxn : forwardRxns) {
        if (oneRxn.rxnType == ReactionType::bimolecular) {
            for (auto& oneProduct : oneRxn.productListNew) {
                if (oneProduct.requiresState == this->rxnReactants[0].state) {
                    // add species for the biomolecular Product with other states of oneProduct
                    for (auto& oneState : molTemplateList[this->rxnReactants[0].molTypeIndex].interfaceList[this->rxnReactants[0].relIndex].stateList) {
                        if (oneState.iden != this->rxnReactants[0].state) {
                            ++totSpecies;
                            this->rxnProducts.emplace_back();
                            this->rxnProducts.back().absIndex = totSpecies;
                        }
                    }
                }
            }
        }
    }

    std::cout << "Warning: must declare association reaction before declaring a state change of its product.\n";
}

void ParsedRxn::determine_reactants()
{
    // TODO: My god, I need to fix this. It's so bad.

    if (rxnType == ReactionType::uniMolStateChange) {
        // - if it is stateChange, we already know the reactants. this could be done earlier, but it's better for code
        // legibility to have all the reactant determination done in one place, I think
        int numberOfMatches { 0 }; // index to keep track of matches. should only be 1 at the end
        // find the iface which changes state
        for (auto& reactant : reactantList) {
            for (auto& reactIface : reactant.interfaceList) {
                reactIface.molTypeIndex = reactant.molTypeIndex;
                if (reactIface.ifaceRxnStatus == Involvement::stateChange) {
                    intReactantList.push_back(reactIface.absIndex);
                    rxnReactants.emplace_back(reactIface);
                    ++numberOfMatches;
                } else if (reactIface.ifaceRxnStatus == Involvement::possible) {
                    // if the reactant is still listed as possible, look for a matching possible iface in the product.
                    // if they're equal, they're ancillary
                    for (auto& product : productList) {
                        auto prodIfaceItr = std::find_if(product.interfaceList.begin(), product.interfaceList.end(),
                            [&](const ParsedMol::IfaceInfo& prodIface) -> bool { return reactIface == prodIface; });
                        if (prodIfaceItr != product.interfaceList.end()) {
                            reactIface.ifaceRxnStatus = Involvement::ancillary;
                            prodIfaceItr->ifaceRxnStatus = Involvement::ancillary;
                            break;
                        }
                        // TODO: add handling if they're a match isn't found
                    }
                }
            }
            if (numberOfMatches == 1) {
                for (auto& product : productList) {
                    for (auto& prodIface : product.interfaceList) {
                        if (prodIface.molTypeIndex == reactant.molTypeIndex
                            && prodIface.ifaceRxnStatus == Involvement::stateChange) {
                            intProductList.emplace_back(prodIface.absIndex);
                            rxnProducts.emplace_back(prodIface);
                        }
                    }
                }
            } else if (numberOfMatches > 1) {
                // TODO: Error handling
                exit(1);
            } else {
                continue;
            }
        }
    } else {
        std::vector<int> speciesUsed; // keep track of which reaction species you've found a reactant on
        for (auto& product : productList) {
            for (auto& prodIface : product.interfaceList) {
                prodIface.molTypeIndex = product.molTypeIndex;
                if (prodIface.ifaceRxnStatus == Involvement::interactionChange) {
                    // This interface only changes binding status.

                    // make sure there's only one match. if two, throw to something to see if it matters
                    int numberOfMatches { 0 };
                    for (auto& reactant : reactantList) {
                        // only check if the interface lists are the same size
                        if (reactant.specieIndex == product.specieIndex && reactant.molTypeIndex == product.molTypeIndex
                            && reactant.interfaceList.size() == product.interfaceList.size()) {
                            // find the interface with the same name and only on a different species than we've already
                            // gotten a reactant from
                            auto reactIfaceItr = std::find_if(reactant.interfaceList.begin(),
                                reactant.interfaceList.end(), [&](const ParsedMol::IfaceInfo& reactIface) -> bool {
                                    return reactIface.ifaceName == prodIface.ifaceName
                                        && std::find(speciesUsed.begin(), speciesUsed.end(), reactIface.speciesIndex)
                                        == speciesUsed.end();
                                });

                            if (reactIfaceItr != reactant.interfaceList.end()) {
                                std::cout << "Found corresponding interfaces: React: " << reactant.molName << '('
                                          << reactIfaceItr->ifaceName << ')' << " -- Prod: " << product.molName << '('
                                          << prodIface.ifaceName << ')' << '\n';

                                // we found a matching interface in the reactants. Make it a reactant
                                reactIfaceItr->molTypeIndex = prodIface.molTypeIndex;
                                reactIfaceItr->ifaceRxnStatus = Involvement::interactionChange;
                                intReactantList.push_back(reactIfaceItr->absIndex);
                                rxnReactants.emplace_back(*reactIfaceItr);

                                // iterate the number of matches
                                ++numberOfMatches;
                                speciesUsed.emplace_back(reactIfaceItr->speciesIndex);
                                break;
                            }
                        }
                    }
                } else if (prodIface.ifaceRxnStatus == Involvement::possible) {
                    for (auto& reactant : reactantList) {
                        if (reactant.specieIndex == product.specieIndex && reactant.molTypeIndex == product.molTypeIndex
                            && reactant.interfaceList.size() == product.interfaceList.size()) {
                            for (auto& reactIface : reactant.interfaceList) {
                                reactIface.molTypeIndex = reactant.molTypeIndex;
                                int numOfMatches { 0 };
                                if (reactIface == prodIface) {
                                    reactIface.ifaceRxnStatus = Involvement::ancillary;
                                    prodIface.ifaceRxnStatus = Involvement::ancillary;
                                    break;
                                }

                                if (numOfMatches > 1)
                                    exit(1);
                            }
                        }
                    }
                } else if (prodIface.ifaceRxnStatus == Involvement::stateChange) {
                    for (auto& reactant : reactantList) {
                        if ((prodIface.speciesIndex == reactant.specieIndex)
                            && (reactant.molTypeIndex == product.molTypeIndex)
                            && (reactant.interfaceList.size() == product.interfaceList.size())) {
                            auto reactIfaceItr = std::find_if(reactant.interfaceList.begin(),
                                reactant.interfaceList.end(), [&](const ParsedMol::IfaceInfo& reactIface) {
                                    return reactIface.ifaceName == prodIface.ifaceName
                                        && std::find(speciesUsed.begin(), speciesUsed.end(), reactIface.speciesIndex)
                                        == speciesUsed.end();
                                });
                            if (reactIfaceItr != reactant.interfaceList.end()) {
                                std::cout << "Found corresponding interfaces: React: " << reactant.molName << '('
                                          << reactIfaceItr->ifaceName << ") -- Prod: " << product.molName << '('
                                          << prodIface.ifaceName << ')' << '\n';

                                this->rxnType = ReactionType::biMolStateChange;
                                reactIfaceItr->molTypeIndex = prodIface.molTypeIndex;
                                reactIfaceItr->ifaceRxnStatus = Involvement::stateChange;
                                intReactantList.push_back(reactIfaceItr->absIndex);
                                rxnReactants.emplace_back(*reactIfaceItr);
                                intProductList.push_back(prodIface.absIndex);
                                rxnProducts.emplace_back(prodIface);
                            }
                        }
                    }
                }
            }
        }

        if (this->rxnType == ReactionType::biMolStateChange) {
            int facilIndex { (rxnReactants[0].speciesIndex == 0) ? 1 : 0 };
            ParsedMol::IfaceInfo& facilIface = reactantList[facilIndex].interfaceList.at(0);
            facilIface.ifaceRxnStatus = Involvement::facilitator;

            // make the facilitator the first thing in the reactant and product list
            intReactantList.emplace(intReactantList.begin(), facilIface.absIndex);
            rxnReactants.emplace(rxnReactants.begin(), facilIface);
            intProductList.emplace(intProductList.begin(), facilIface.absIndex);
            rxnProducts.emplace(rxnProducts.begin(), facilIface);
        }

        // if the two reactants have the same ifaceIndex, it's a symmetric reaction
        if (rxnReactants.size() != 2 && this->rxnType == ReactionType::bimolecular) {
            std::cerr << "Error determining reactants, check reaction file for missing interfaces." << '\n';
            for (auto& reactant : rxnReactants) {
                std::cerr << reactant << ' ';
            }
            exit(1);
        }
        this->isSymmetric = (rxnReactants[0] == rxnReactants[1]);
    }
}

void ParsedRxn::determine_creation_products(const std::vector<MolTemplate>& molTemplateList)
{
    for (auto& reactant : reactantList) {
        try {
            auto tempNameItr = std::find_if(molTemplateList.begin(), molTemplateList.end(),
                [&](const MolTemplate& oneTemp) -> bool { return oneTemp.molName == reactant.molName; });
            reactant.molTypeIndex = static_cast<int>(tempNameItr->molTypeIndex);
            reactant.specieIndex = &reactant - &reactantList[0];

            for (auto& iface : reactant.interfaceList) {
                iface.molTypeIndex = reactant.molTypeIndex;
                intReactantList.emplace_back(iface.absIndex);
                rxnReactants.emplace_back(iface);
            }
        } catch (std::out_of_range& e) {
            // TODO: write this
            std::cout << e.what() << '\n';
            exit(1);
        }
    }

    for (auto& product : productList) {
        try {
            auto tempNameItr = std::find_if(molTemplateList.begin(), molTemplateList.end(),
                [&](const MolTemplate& oneTemp) -> bool { return oneTemp.molName == product.molName; });
            product.molTypeIndex = static_cast<int>(tempNameItr->molTypeIndex);
            product.specieIndex = &product - &productList[0];

            for (auto& iface : product.interfaceList) {
                iface.molTypeIndex = product.molTypeIndex;
                intProductList.emplace_back(iface.absIndex);
                rxnProducts.emplace_back(iface);
            }
        } catch (std::out_of_range& e) {
            // TODO: write this
            std::cout << e.what() << '\n';
            exit(1);
        }
    }
}

void ParsedRxn::create_other_iface_lists(const std::vector<MolTemplate>& molTemplateList)
{
    for (auto& reactant : reactantList) {
        std::vector<RxnIface> tmpIfaceVec;
        for (auto& reactIface : reactant.interfaceList) {
            // if the iface's status is ancillary, create a new element in the vector holding those ifaces
            if (reactIface.ifaceRxnStatus == Involvement::ancillary) {
                int relIface { molTemplateList[reactant.molTypeIndex].find_relIndex_from_absIndex(
                    reactIface.absIndex) };
                tmpIfaceVec.emplace_back(reactIface.ifaceName, reactant.molTypeIndex, reactIface.absIndex, relIface,
                    (reactIface.state == '\0') ? '\0' : reactIface.state, reactIface.isBound);
            }
        }
        otherIfaceLists.emplace_back(tmpIfaceVec);
    }
}

bool ParsedRxn::check_for_conditional_rates(
    int& totSpecies, std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns)
{
    // sort the reactants in order according to their iface index, to make comparing reactions easier
    // don't do this for bimolecular state changes, since we put the facilitator first, regardless of its index
    //    if (rxnType != ReactionType::biMolStateChange) {
    //        if (intReactantList.size() == 2 && intReactantList[0] > intReactantList[1]) {
    //            std::swap(intReactantList[0], intReactantList[1]);
    //            std::swap(rxnReactants[0], rxnReactants[1]);
    //        }
    //    } else {
    if (rxnType == ReactionType::biMolStateChange) {
        if (intProductList[0] != intReactantList[0]) {
            std::swap(intProductList[0], intProductList[1]);
            std::swap(rxnProducts[0], rxnProducts[1]);
        }
    }
    //    }

    std::vector<int> tmpVec1 { intReactantList };
    std::sort(tmpVec1.begin(), tmpVec1.end());
    for (auto& oneRxn : forwardRxns) {
        std::vector<int> tmpVec2 { oneRxn.intReactantList };
        std::sort(tmpVec2.begin(), tmpVec2.end());
        // if the reactants are the same...
        if (tmpVec1 == tmpVec2 && rxnType == oneRxn.rxnType) {
            // if the values had to be swapped, swap everything so the order is the same as the parent reaction
            if ((intReactantList != tmpVec1 && oneRxn.intReactantList == tmpVec2)
                || (intReactantList == tmpVec1 && oneRxn.intReactantList != tmpVec2)) {
                std::swap(intReactantList[0], intReactantList[1]);
                std::swap(rxnReactants[0], rxnReactants[1]);
                std::swap(rxnProducts[0], rxnProducts[1]);
                std::swap(otherIfaceLists[0], otherIfaceLists[1]);
            }
            // ...give that ForwardRxn this on rate and list of ancillary interfaces...
            oneRxn.rateList.emplace_back(onRate3Dka, otherIfaceLists);
            // ...and the conjugate BackRxn the off rate and list of ancillary interfaces. and interfaces which change
            // state (swapped, of course)
            if (oneRxn.isReversible) {
                backRxns[oneRxn.conjBackRxnIndex].rateList.emplace_back(this->offRatekb, otherIfaceLists);
            }
            --totSpecies; // since it's not a new reaction, reduce the number of total species...
            std::cout << "Forward Reaction " << &oneRxn - &forwardRxns[0]
                      << " has been updated with a new rate:\nRate:" << onRate3Dka << '\n';
            std::cout << "Reactant 1 requires interfaces:\n";
            for (auto& iface : otherIfaceLists[0])
                std::cout << iface << '\n';
            if (rxnType == ReactionType::bimolecular || rxnType == ReactionType::biMolStateChange) {
                std::cout << "Reactant 2 requires interfaces:\n";
                for (auto& iface : otherIfaceLists[1])
                    std::cout << iface << '\n';
            }
            return true;
        }
    }
    return false;
}

ParsedRxn ParsedRxn::make_one_split_reaction(
    int& totSpecies, int parsedMolIndex, ParsedMol::IfaceInfo& reactIface, const Interface::State& tempState)
{
    ++totSpecies;
    ParsedRxn newRxn(*this);

    auto newIfaceList = newRxn.reactantList[parsedMolIndex].interfaceList;
    // look through the list of interfaces for the interface which we need to make the multiple
    // reactions for
    for (auto& newIface : newIfaceList) {
        if (newIface.ifaceName == reactIface.ifaceName) {
            newIface.state = tempState.iden;
            // if the interface isn't bound, give it the index of the MolTemplate interface state
            if (!newIface.isBound)
                newIface.absIndex = tempState.index;

            // update the reactant lists
            for (auto& oneReact : newRxn.rxnReactants) {
                // only need to update the IfaceInfo, the name should stay the same (also it doesn't compile if you
                // try to change the map key)
                if (oneReact.ifaceName == newIface.ifaceName)
                    oneReact = newIface;
            }
            // change the Reactant index to the updated, correct index
            for (auto& oneReact : newRxn.intReactantList) {
                if (oneReact == -1) {
                    oneReact = newIface.absIndex;
                    break;
                }
            }

            // update the product lists
            newRxn.intProductList = std::vector<int>(newRxn.intProductList.size(), totSpecies);
            for (auto& oneProd : newRxn.rxnProducts) {
                oneProd.absIndex = totSpecies;
                if (oneProd.ifaceName == newIface.ifaceName)
                    oneProd.state = newIface.state;
            }
        }
    }
    newRxn.reactantList[parsedMolIndex].interfaceList = newIfaceList;
    return newRxn;
}

std::vector<ParsedRxn> ParsedRxn::split_into_multiple_reactions(
    int& totSpecies, std::vector<ForwardRxn>& forwardRxns, const std::vector<MolTemplate>& molTemplateList)
{
    struct StateCont {
        std::string molName;
        parsedIface ifaceInfo;
        char state;
    };

    std::vector<ParsedRxn> newRxns;
    // explicit states (molName, ifaceName, iterator to the iface)
    for (auto& reactIface : noStateList) {
        --totSpecies; // need to remove a species for each interface we split into multiple reactions
        ParsedMol thisMol { this->reactantList[reactIface.first] };

        auto tempNameItr = std::find_if(molTemplateList.begin(), molTemplateList.end(),
            [&](const MolTemplate& oneTemp) -> bool { return oneTemp.molName == thisMol.molName; });
        // reactIface.first looks bad, but reactIface.second = the iface element of noStateList,
        // reactIface.first is that iface's name
        auto tempIfaceItr = std::find_if(tempNameItr->interfaceList.begin(), tempNameItr->interfaceList.end(),
            [&](const Interface& oneIface) -> bool { return oneIface.name == reactIface.second.ifaceName; });

        // TODO: error for out of range

        // if only one interface has no explicit state, we can just make the new reactions immediately
        for (auto& tempState : tempIfaceItr->stateList)
            newRxns.emplace_back(make_one_split_reaction(totSpecies, reactIface.first, reactIface.second, tempState));
    }

    return newRxns;
}
