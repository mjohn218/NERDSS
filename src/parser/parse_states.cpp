#include "parser/parser_functions.hpp"

void parse_states(std::string& line, MolTemplate& molTemplate)
{
    auto statePos = line.find('~');

    // get the interface name
    std::string ifaceName { line.substr(0, statePos) };
    line.erase(0, statePos + 1);

    // get each of its states
    std::vector<std::string> states;
    while ((statePos = line.find('~')) != std::string::npos) {
        states.emplace_back(line.substr(0, statePos));
        line.erase(0, statePos + 1);
    }
    states.emplace_back(line.substr(0, std::string::npos)); // the last value will be out of scope of the above loop

    // now find the interface which corresponds to this one
    auto ifaceNameItr = std::find_if(molTemplate.interfaceList.begin(), molTemplate.interfaceList.end(),
        [&](const Interface& oneIface) -> bool { return oneIface.name == ifaceName; });

    if (ifaceNameItr == molTemplate.interfaceList.end()) {
        std::cout << ifaceName << " is not a valid interface for molecule " << molTemplate.molName << ".\n";
        exit(1);
    } else {
        // TODO: replace the placeholder integer
        for (auto& oneState : states)
            ifaceNameItr->stateList.emplace_back(static_cast<char>(std::toupper(oneState[0])), -1);

        std::cout << ifaceName << " has " << states.size() << " state(s): ";
        for (auto oneState : states)
            std::cout << oneState << "\t";
        std::cout << std::endl;
    }
}
