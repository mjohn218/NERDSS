#include "parser/parser_functions.hpp"

void check_for_valid_states( size_t parsedMolIndex, ParsedMol& targMol, ParsedRxn& parsedRxn, const std::vector<MolTemplate>& molTemplateList)
{

    // find the template which corresponds to the reactant molName
    std::cout << "Checking for valid states for reactant " << targMol.molName << '\n';
    auto tempNameItr = std::find_if(molTemplateList.begin(), molTemplateList.end(),
        [&](const MolTemplate& oneTemp) -> bool { return oneTemp.molName == targMol.molName; });

    if (tempNameItr == molTemplateList.end()) {
        std::cerr << "Error, target molecule's type cannot be found in list of templates.";
        exit(1);
    }

    targMol.molTypeIndex = static_cast<int>(tempNameItr->molTypeIndex);

    //    for (auto& reactIfaceItr : targMol.interfaceList) {
    for (auto ifaceItr = targMol.interfaceList.begin(); ifaceItr != targMol.interfaceList.end(); ++ifaceItr) {
        // find the template iface which corresponds to the reactant iface
        auto tempIfaceItr = std::find_if(tempNameItr->interfaceList.begin(), tempNameItr->interfaceList.end(),
            [&](const Interface& oneTempIface) -> bool { return oneTempIface.name == ifaceItr->ifaceName; });

        if (tempIfaceItr == tempNameItr->interfaceList.end()) {
            std::cout << ifaceItr->ifaceName << " is not a valid interface for molecule template "
                      << tempNameItr->molName << '\n';
            exit(120);
        }

        auto tempStateItr = tempIfaceItr->stateList.begin();

        if (tempIfaceItr->stateList.size() != 1 && ifaceItr->state == '\0') {
            // if no state is explicitly provided, but the iface has specific states, we need to create separate
            // reactions for each state
            std::cout << "No interface state provided for reactant "
                      << write_mol_iface(targMol.molName, ifaceItr->ifaceName)
                      << ". Will create separate reactions for each possible interface." << '\n';
            parsedRxn.noStateList.insert({ parsedMolIndex, *ifaceItr });
        } else if (tempIfaceItr->stateList.size() != 1 && ifaceItr->state != '\0') {
            // find the state which corresponds to the reactant iface state
            tempStateItr = std::find_if(tempIfaceItr->stateList.begin(), tempIfaceItr->stateList.end(),
                [&](const Interface::State& oneTempState) -> bool { return oneTempState.iden == ifaceItr->state; });

            if (tempStateItr == tempIfaceItr->stateList.end()) {
                // if it doesn't exist, exit
                std::cout << ifaceItr->state << " is not a valid state for interface "
                          << write_mol_iface(targMol.molName, ifaceItr->ifaceName) << '\n';
                exit(1);
            }

            // if the state does exist, look for a state change
            if (parsedRxn.rxnType != ReactionType::destruction) {
                std::cout << "Found state, looking for state change..." << '\n';
                check_for_state_change(*ifaceItr, targMol, parsedRxn);
            } else
                std::cout << "Found state for reactant " << write_mol_iface(targMol.molName, ifaceItr->ifaceName)
                          << '\n';
        }
    }
}
