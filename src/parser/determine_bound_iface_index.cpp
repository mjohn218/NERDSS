#include "parser/parser_functions.hpp"

void determine_bound_iface_index(int& totSpecies, ParsedMol::IfaceInfo& targIface, ParsedRxn& parsedRxn,
    const std::vector<ForwardRxn>& forwardRxns, const std::vector<MolTemplate>& molTemplateList)
{
    // make sure they're actually bound (this is just a double check)
    // TODO: I think this might just accidentally work for symmetric reactions
    for (auto& product : parsedRxn.productList) {
        for (auto& prodIface : product.interfaceList) {

            // if it isn't the same iface and the two ifaces are bound
            if ((targIface != prodIface) && (targIface.bondIndex == prodIface.bondIndex)) {
                // iterate the number of species in the system
                ++totSpecies;

                // change the two interfaces to reflect that they change interaction
                targIface.change_ifaceRxnStatus(totSpecies, Involvement::interactionChange);
                prodIface.change_ifaceRxnStatus(totSpecies, Involvement::interactionChange);

                // make sure the matching product actually has a molecule index
                if (prodIface.molTypeIndex == -1)
                    prodIface.molTypeIndex = product.molTypeIndex;

                // add the products to various lists
                parsedRxn.intProductList.emplace_back(totSpecies);
                parsedRxn.rxnProducts.emplace_back(prodIface);
                parsedRxn.rxnProducts.emplace_back(targIface);
                return;
            }
        }
    }
}
