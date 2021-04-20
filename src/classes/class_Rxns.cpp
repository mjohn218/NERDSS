/*! \file class_rxn_functions

 * ### Created on 10/18/18 by Matthew Varga
 * ### Purpose
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */

#include "classes/class_Rxns.hpp"
#include "classes/class_bngl_parser.hpp"

#include <iomanip>

// Defined some static variables
unsigned RxnBase::numberOfRxns = 0;
int RxnBase::totRxnSpecies = 0;

/* RXNIFACE */
RxnIface::RxnIface(std::string ifaceName, int molTypeIndex, int absIfaceIndex, int relIfaceIndex, char requiresState,
    bool requiresInteraction)
    : ifaceName(std::move(ifaceName))
    , molTypeIndex(molTypeIndex)
    , absIfaceIndex(absIfaceIndex)
    , relIfaceIndex(relIfaceIndex)
    , requiresState(requiresState)
    , requiresInteraction(requiresInteraction)
{
}

std::ostream& operator<<(std::ostream& os, const RxnIface& rxnIface)
{
    os << "[iface name: " << rxnIface.ifaceName << ", iface index: " << rxnIface.absIfaceIndex;
    if (rxnIface.requiresState != '\0')
        os << ", requires state " << rxnIface.requiresState;
    if (rxnIface.requiresInteraction)
        os << ", requires interaction";
    else
        os << ", required Free";
    os << ']';
    return os;
}

bool RxnIface::operator==(const Molecule::Iface& molIface) const
{
    //    return this->molTypeIndex == molIface.molTypeIndex && this->absIfaceIndex == molIface.index
    //        && this->requiresState == molIface.stateIden && this->requiresInteraction == molIface.isBound;
    // this can't check absIfaceIndex, since bound interfaces in the reactants can't have indexed bonds, yet
    return this->molTypeIndex == molIface.molTypeIndex && this->relIfaceIndex == molIface.relIndex
        && this->requiresState == molIface.stateIden && this->requiresInteraction == molIface.isBound;
}

/* FORWARDRXN */
ForwardRxn::ForwardRxn(ParsedRxn& parsedRxn, const std::vector<MolTemplate>& molTemplateList)
{
    // set booleans
    isReversible = parsedRxn.isReversible;
    isOnMem = parsedRxn.isOnMem;
    isSymmetric = parsedRxn.isSymmetric;
    hasStateChange = parsedRxn.hasStateChange;
    isCoupled = parsedRxn.isCoupled;
    coupledRxn = parsedRxn.coupledRxn;
    isObserved = parsedRxn.isObserved;
    observeLabel = parsedRxn.observeLabel;
    bindRadSameCom = parsedRxn.bindRadSameCom;
    loopCoopFactor = parsedRxn.loopCoopFactor;
    irrevRingClosure = parsedRxn.irrevRingClosure;
    length3Dto2D = parsedRxn.length3Dto2D;
    rxnLabel = parsedRxn.rxnLabel;
    // set the product name for species output file headers
    productName = parsedRxn.productName;
    //    fullProductName = parsedRxn.fullProductName;

    // set reaction type
    rxnType = parsedRxn.rxnType;

    absRxnIndex = RxnBase::numberOfRxns;
    ++RxnBase::numberOfRxns;

    assocAngles = parsedRxn.assocAngles;
    bindRadius = parsedRxn.bindRadius;

    // set index lists
    intProductList = parsedRxn.intProductList;
    intReactantList = parsedRxn.intReactantList;

    // set norms
    norm1 = parsedRxn.norm1;
    norm2 = parsedRxn.norm2;

    // set excludeVolumeBound
    excludeVolumeBound = parsedRxn.excludeVolumeBound;

    // create reactant and product lists
    rateList.emplace_back();
    rateList.back().rate = parsedRxn.onRate3Dka;
    rateList.back().otherIfaceLists = parsedRxn.otherIfaceLists;
    // reactants (and while we're in the loop, create the RateState otherIfaceList
    for (auto& mol : parsedRxn.reactantList) {
        for (auto& iface : mol.interfaceList) {
            if (iface.ifaceRxnStatus == Involvement::ancillary) {
                //                rateList.back().otherIfaceLists.emplace_back(
                //                    iface.ifaceName, mol.molTypeIndex, iface.absIndex, iface.relIndex, iface.state,
                //                    iface.isBound);
            } else if (iface.ifaceRxnStatus == Involvement::stateChange) {
                stateChangeIface.first = RxnIface(
                    iface.ifaceName, mol.molTypeIndex, iface.absIndex, iface.relIndex, iface.state, iface.isBound);
            }
        }
    }

    // products
    for (auto& mol : parsedRxn.productList) {
        for (auto& iface : mol.interfaceList) {
            if (iface.ifaceRxnStatus == Involvement::stateChange) {
                stateChangeIface.second = RxnIface(
                    iface.ifaceName, mol.molTypeIndex, iface.absIndex, iface.relIndex, iface.state, iface.isBound);
            }
        }
    }

    // set up the ForwardRxn reactant and product lists
    // TODO: replace these with constructors
    for (auto& reactant : parsedRxn.rxnReactants) {
        reactantListNew.emplace_back(reactant.ifaceName, reactant.molTypeIndex, reactant.absIndex, reactant.relIndex,
            reactant.state, reactant.isBound);
    }

    for (auto& product : parsedRxn.rxnProducts) {
        productListNew.emplace_back(product.ifaceName, product.molTypeIndex, product.absIndex, product.relIndex,
            product.state, product.isBound);
    }

    // TODO: this is a temporary fix for a problem elsewhere in the parser
    for (unsigned prodItr { 0 }; prodItr < productListNew.size(); ++prodItr) {
        if (productListNew[prodItr].relIfaceIndex == -1)
            productListNew[prodItr].relIfaceIndex = reactantListNew[prodItr].relIfaceIndex;
    }
}

void ForwardRxn::display() const
{
    std::cout << "Absolute index: " << absRxnIndex << std::endl;
    std::cout << "Type: " << rxnType << std::endl;

    if (rxnType == ReactionType::bimolecular) {
        std::cout << "Reactants:\n";
        for (auto& reactant : reactantListNew)
            std::cout << ' ' << reactant << std::endl;
        std::cout << std::endl;
        std::cout << "Products:\n";
        for (auto& product : productListNew)
            std::cout << ' ' << product << std::endl;
        std::cout << std::endl;
    } else if (rxnType == ReactionType::biMolStateChange) {
        std::cout << "Facilitator: ";
        for (auto& reactant : reactantListNew) {
            if (reactant.absIfaceIndex != stateChangeIface.first.absIfaceIndex)
                std::cout << reactant << std::endl;
        }
    }

    if (hasStateChange) {
        std::cout << "State Change Reactant: " << stateChangeIface.first << std::endl;
        std::cout << "State Change Product: " << stateChangeIface.second << std::endl;
    }

    std::cout << "Rate(s):" << std::endl;
    if (rxnType != ReactionType::uniMolStateChange) {
        for (auto& rate : rateList) {
            std::cout << "Rate " << &rate - &rateList[0] << ": " << rate.rate << std::endl;
            if (!rate.otherIfaceLists.empty()) {
                std::cout << "Reactant 1 requires interfaces:" << std::endl;
                for (auto& iface : rate.otherIfaceLists[0]) {
                    std::cout << ' ' << iface << std::endl;
                }
                std::cout << "Reactant 2 requires interfaces:" << std::endl;
                for (auto& iface : rate.otherIfaceLists[1]) {
                    std::cout << ' ' << iface << std::endl;
                }
            }
        }
    } else {
        for (auto& rate : rateList) {
            std::cout << "Rate " << &rate - &rateList[0] << ": " << rate.rate << std::endl;
        }
    }

    if (rxnType == ReactionType::bimolecular) {
        std::cout << "Sigma: " << bindRadius;
        std::cout << std::endl;
        assocAngles.display();
        std::cout << std::endl;
    }

    std::cout << "label: " << rxnLabel << std::endl;

    if (rxnType != ReactionType::uniMolStateChange) {
        std::cout << "Is On Membrane: " << std::boolalpha << isOnMem << std::endl;
        std::cout << "bindRadSameCom " << bindRadSameCom << std::endl;
        std::cout << "loopCoopFactor " << loopCoopFactor << std::endl;
        std::cout << "length3Dto2D " << length3Dto2D << std::endl;
        std::cout << "isCoupled? " << isCoupled << std::endl;
        if (isCoupled)
            std::cout << " coupledRxn Number: " << coupledRxn.absRxnIndex << " type: " << coupledRxn.rxnType << " prob to perform coupled: " << coupledRxn.probCoupled << std::endl;
    }
}

void ForwardRxn::assoc_display(const std::vector<MolTemplate>& molTemplateList) const
{
    std::cout << "\nReaction data:\n";
    // TODO: break this off into the constructor, and maybe put it onto RateState, because each conditional rate will
    // have different otherIfaceList and/or stateChangeIface.
    std::cout << "Product: ";
    std::cout << molTemplateList[productListNew.front().molTypeIndex].molName << '(' << productListNew.front().ifaceName
              << "!1)";
    for (const auto& product : productListNew)
        std::cout << '.' << molTemplateList[product.molTypeIndex].molName << '(' << product.ifaceName << "!1)";
    std::cout << "\nAngles:\n";
    assocAngles.display();
    std::cout << std::endl;
}

void ForwardRxn::Angles::display() const
{
    std::cout << "Association angles:\n";
    std::cout << std::setw(10) << std::left << "Theta 1" << std::setw(7) << std::right << theta1 << '\n';
    std::cout << std::setw(10) << std::left << "Theta 2" << std::setw(7) << std::right << theta2 << '\n';
    std::cout << std::setw(10) << std::left << "Phi 1" << std::setw(7) << std::right << phi1 << '\n';
    std::cout << std::setw(10) << std::left << "Phi 2" << std::setw(7) << std::right << phi2 << '\n';
    std::cout << std::setw(10) << std::left << "Omega" << std::setw(7) << std::right << omega << '\n';
}

ForwardRxn ForwardRxn::bngl_copy_rxn()
{
    ForwardRxn tmp {};
    tmp.intProductList = intProductList;
    tmp.intReactantList = intReactantList;
    tmp.rxnType = rxnType;

    return tmp;
}

/* BACKRXN */
BackRxn::BackRxn(double offRatekb, ForwardRxn& forwardRxn)
{
    // set booleans
    isOnMem = forwardRxn.isOnMem;
    isSymmetric = forwardRxn.isSymmetric;
    hasStateChange = forwardRxn.hasStateChange;
    isCoupled = forwardRxn.isCoupled;
    coupledRxn = forwardRxn.coupledRxn;
    isObserved = forwardRxn.isObserved;
    observeLabel = forwardRxn.observeLabel;
    bindRadSameCom = forwardRxn.bindRadSameCom;
    loopCoopFactor = forwardRxn.loopCoopFactor;
    length3Dto2D = forwardRxn.length3Dto2D;
    rxnLabel = forwardRxn.rxnLabel;
    // set reaction type
    rxnType = forwardRxn.rxnType;

    absRxnIndex = RxnBase::numberOfRxns;
    ++RxnBase::numberOfRxns;

    // swap reactants and products
    stateChangeIface
        = std::pair<RxnIface, RxnIface> { forwardRxn.stateChangeIface.second, forwardRxn.stateChangeIface.first };

    reactantListNew = forwardRxn.productListNew;
    productListNew = forwardRxn.reactantListNew;

    intReactantList = forwardRxn.intProductList;
    intProductList = forwardRxn.intReactantList;

    // sort them, so they're the same as the forward reaction
    //        if (intReactantList.size() > 1 && intReactantList[0] > intReactantList[1]) {
    //            std::swap(intReactantList[0], intReactantList[1]);
    //            std::swap(reactantListNew[0], reactantListNew[1]);
    //        }

    rateList.emplace_back(offRatekb, forwardRxn.rateList.back().otherIfaceLists);
}

void BackRxn::display() const
{
    std::cout << "Absolute index: " << absRxnIndex << '\n';
    std::cout << "Type: " << rxnType << '\n';
    if (!hasStateChange && rxnType == ReactionType::bimolecular) {
        std::cout << "Reactants:\n";
        for (auto& reactant : reactantListNew)
            std::cout << ' ' << reactant << '\n';
        std::cout << std::endl;
        std::cout << "Products:\n";
        for (auto& product : productListNew)
            std::cout << ' ' << product << '\n';
        std::cout << std::endl;
    } else if (rxnType == ReactionType::biMolStateChange) {
        std::cout << "Facilitator: ";
        for (auto& reactant : reactantListNew) {
            if (reactant.absIfaceIndex != stateChangeIface.first.absIfaceIndex)
                std::cout << reactant << '\n';
        }
    }

    if (hasStateChange) {
        std::cout << "State Change Reactant: " << stateChangeIface.first << '\n';
        std::cout << "State Change Product: " << stateChangeIface.second << '\n';
    }

    std::cout << "\nRate(s):\n";
    for (auto& rate : rateList) {
        std::cout << "Rate " << &rate - &rateList[0] << ": " << rate.rate << '\n';
        if (!rate.otherIfaceLists.empty()) {
            std::cout << "Reactant 1 requires interfaces:\n";
            for (auto& iface : rate.otherIfaceLists[0]) {
                std::cout << ' ' << iface << '\n';
            }
            std::cout << "Reactant 2 requires interfaces:\n";
            for (auto& iface : rate.otherIfaceLists[1]) {
                std::cout << ' ' << iface << '\n';
            }
        }
    }
    std::cout << "On Membrane? " << std::boolalpha << isOnMem << std::endl;
}

/* CREATEDESTRUCTRXNS */

CreateDestructRxn::CreateDestructRxn(ParsedRxn& parsedRxn, const std::vector<MolTemplate>& molTemplateList)
{
    // set booleans
    isOnMem = parsedRxn.isOnMem;
    isSymmetric = parsedRxn.isSymmetric;
    hasStateChange = parsedRxn.hasStateChange;
    isObserved = parsedRxn.isObserved;
    observeLabel = parsedRxn.observeLabel;
    creationRadius = parsedRxn.creationRadius;
    rxnLabel = parsedRxn.rxnLabel;

    // set reaction type
    rxnType = parsedRxn.rxnType;

    // initialize interface lists with default interfaces from MolTemplate
    for (const auto& reactant : parsedRxn.reactantList) {
        CreateDestructMol tmpMol {};
        tmpMol.molTypeIndex = reactant.molTypeIndex;
        tmpMol.molName = reactant.molName;
        for (const auto& tempIface : molTemplateList[reactant.molTypeIndex].interfaceList) {
            RxnIface tmpIface {};
            tmpIface.ifaceName = tempIface.name;
            tmpIface.molTypeIndex = reactant.molTypeIndex;
            tmpIface.relIfaceIndex = tempIface.index;
            tmpIface.absIfaceIndex = tempIface.stateList[0].index;
            tmpMol.interfaceList.emplace_back(tmpIface);
        }
        this->reactantMolList.emplace_back(tmpMol);
    }

    for (const auto& product : parsedRxn.productList) {
        CreateDestructMol tmpMol {};
        tmpMol.molTypeIndex = product.molTypeIndex;
        tmpMol.molName = product.molName;
        for (const auto& tempIface : molTemplateList[product.molTypeIndex].interfaceList) {
            RxnIface tmpIface {};
            tmpIface.ifaceName = tempIface.name;
            tmpIface.molTypeIndex = product.molTypeIndex;
            tmpIface.relIfaceIndex = tempIface.index;
            tmpIface.absIfaceIndex = tempIface.stateList[0].index;
            tmpMol.interfaceList.emplace_back(tmpIface);
        }
        this->productMolList.emplace_back(tmpMol);
    }

    for (const auto& reactant : parsedRxn.rxnReactants) {
        reactantMolList[reactant.speciesIndex].interfaceList[reactant.relIndex].requiresState = reactant.state;
        reactantMolList[reactant.speciesIndex].interfaceList[reactant.relIndex].requiresInteraction = reactant.isBound;
    }

    for (const auto& product : parsedRxn.rxnProducts) {
        productMolList[product.speciesIndex].interfaceList[product.relIndex].requiresState = product.state;
        productMolList[product.speciesIndex].interfaceList[product.relIndex].requiresInteraction = product.isBound;
    }

    rateList.emplace_back();
    rateList.back().rate = parsedRxn.onRate3Dka;
    rateList.back().otherIfaceLists = parsedRxn.otherIfaceLists;

    absRxnIndex = numberOfRxns;
    ++numberOfRxns;

    // set reactant/product lists, depending on the reaction type
}

void CreateDestructRxn::display() const
{
    // TODO: Flesh this out so it outputs interfaces too
    std::cout << "Absolute index: " << absRxnIndex << '\n';
    std::cout << "Type: " << rxnType << '\n';
    if (rxnType != ReactionType::zerothOrderCreation) {
        std::cout << "Reactants:";
        for (auto& reactant : reactantMolList)
            std::cout << " [" << reactant.molName << ']';
        std::cout << std::endl;
    }
    if (rxnType != ReactionType::destruction) {
        std::cout << "Products:";
        for (auto& product : productMolList)
            std::cout << " [" << product.molName << ']';
        std::cout << std::endl;
    }
    std::cout << "On Membrane? " << std::boolalpha << isOnMem << '\n';
    std::cout << " Number of rates: " << rateList.size() << " first rate: " << rateList[0].rate << '\n';
    std::cout << "otherIfaceListSize: " << rateList[0].otherIfaceLists.size() << '\n';
    std::cout << "Rate(s):\n";
    for (auto& rate : rateList) {
        std::cout << "Rate " << &rate - &rateList[0] << ": " << rate.rate << '\n';
        if (!rate.otherIfaceLists.empty()) {
            std::cout << "Reactant 1 requires interfaces:\n";
            for (auto& iface : rate.otherIfaceLists[0]) {
                std::cout << ' ' << iface << '\n';
            }
        }
    }
    std::cout << "label: " << rxnLabel << std::endl;
}
