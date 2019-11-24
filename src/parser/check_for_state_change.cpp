#include "parser/parser_functions.hpp"

void check_for_state_change(ParsedMol::IfaceInfo& targetIface, ParsedMol& targetMol, ParsedRxn& parsedRxn)
{
    // find the product(s) which correspond to the target iface molecule
    for (auto& product : parsedRxn.productList) {
        if (product.molName == targetMol.molName) {
            product.molTypeIndex = targetMol.molTypeIndex; // just for some bookkeeping

            for (auto& prodIface : product.interfaceList) {
                // if the two interfaces are on the same species and are identical but for the state,
                // that means its state changed during the reaction
                if (areSameExceptState(targetIface, prodIface)) {
                    std::cout << "State change found on species " << product.specieIndex << ". Interface "
                              << prodIface.ifaceName << " changes state from " << targetIface.state << " to "
                              << prodIface.state << ".\n";
                    targetIface.ifaceRxnStatus = Involvement::stateChange;
                    prodIface.ifaceRxnStatus = Involvement::stateChange;
                    parsedRxn.hasStateChange = true;
                }
            }
        }
    }
}
