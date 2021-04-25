#include "parser/parser_functions.hpp"

MolTemplate parse_molFile(std::string& mol)
{
    /* NOTE: need to edit both this and enum class MolKeyword if you want to add keywords  */
    std::map<const std::string, MolKeyword> molKeywords = { { "name", MolKeyword::name },
        { "copies", MolKeyword::copies }, { "isrod", MolKeyword::isRod }, { "islipid", MolKeyword::isLipid },
        { "d", MolKeyword::d }, { "dr", MolKeyword::dr }, { "com", MolKeyword::com }, { "state", MolKeyword::state },
        { "mass", MolKeyword::mass }, { "checkoverlap", MolKeyword::checkOverlap }, { "bonds", MolKeyword::bonds }, { "isimplicitlipid", MolKeyword::isImplicitLipid },
        { "ispoint", MolKeyword::isPoint } , { "counttransition", MolKeyword::countTransition }, { "transitionmatrixsize", MolKeyword::transitionMatrixSize }};

    std::cout << mol + ".mol" << '\n';
    std::ifstream molFile { mol + ".mol" };
    if (!molFile) {
        std::cout << "Cannot open molecule config file for molecule " << mol << '\n';
        exit(1);
    }

    MolTemplate tmpTemplate;

    std::string line;
    auto initialPos = molFile.tellg();
    while (getline(molFile, line)) {
        line.erase(
            std::remove_if(line.begin(), line.end(), [](unsigned char x) { return std::isspace(x); }), line.end());

        // skip if entire line is a comment, or remove the trailing comment
        if (line[0] == '#') {
            initialPos = molFile.tellg();
            continue;
        } else
            remove_comment(line);

        std::string buffer;
        for (auto lineItr = line.begin(); lineItr != line.end(); ++lineItr) {
            // if the line starts with com, it's the beginning of the coordinates block
            if (std::isdigit(*lineItr) && molKeywords.find(buffer)->second == MolKeyword::com) {
                molFile.seekg(initialPos);
                std::cout << "Coordinates: " << std::endl;
                read_internal_coordinates(molFile, tmpTemplate);
                break;
            }

            if (std::isalnum(*lineItr))
                buffer += std::tolower(static_cast<char>(*lineItr));
            else if (*lineItr == '=') {
                auto keyFind = molKeywords.find(buffer);
                if (keyFind == molKeywords.end()) {
                    std::cout << buffer + " is an invalid argument.";
                    exit(1);
                }

                line.erase(line.begin(), lineItr + 1); // + 1 removes the '=' sign. could make this erase(remove_if)

                // find the value type from the keyword and then set that parameter
                if (keyFind->second == MolKeyword::state) {
                    std::cout << "States: " << std::endl;
                    parse_states(line, tmpTemplate);
                    break;
                } else if (keyFind->second == MolKeyword::bonds) {
                    // std::cout << "Found bonds for molecule " << tmpTemplate.molName << ".\n";
                    std::cout << "Bonds: " << std::endl;
                    int numBonds = std::stoi(line);
                    read_bonds(numBonds, molFile, tmpTemplate);
                    break;
                } else {
                    tmpTemplate.set_value(line, keyFind->second);
                    break;
                }
            }
        }

        if (molFile.eof()) {
            break;
        }
        initialPos = molFile.tellg();
    }

    // set state indices
    int ifaceItr { 0 };
    for (auto& iface : tmpTemplate.interfaceList) {
        iface.index = ifaceItr;
        /*initialize values of the molecule templates interfaces, including if it can be in multiple states
         stored in molTemplate.interfaceList[].stateList[].index*/
        if (iface.stateList.empty()) {
            MolTemplate::absToRelIface.push_back(ifaceItr);
            iface.stateList.emplace_back(iface.name, Interface::State::totalNumOfStates);
            ++Interface::State::totalNumOfStates;
        } else {
            tmpTemplate.ifacesWithStates.push_back(&iface - &tmpTemplate.interfaceList[0]);
            for (auto& state : iface.stateList) {
                MolTemplate::absToRelIface.push_back(ifaceItr);
                state.index = Interface::State::totalNumOfStates;
                state.ifaceAndStateName = iface.name + "~" + state.iden;
                ++Interface::State::totalNumOfStates;
            }
        }
        ++ifaceItr;
    }

    // calculate 'radius'
    for (auto& iface : tmpTemplate.interfaceList) {
        Vector tmpVec { tmpTemplate.comCoord - iface.iCoord };
        tmpVec.calc_magnitude();
        if (tmpVec.magnitude > tmpTemplate.radius)
            tmpTemplate.radius = tmpVec.magnitude;
    }
    std::cout << "radius calculated from the ifacesCoord: " << tmpTemplate.radius << " nm" << std::endl;
    // set mass to radius
    if (tmpTemplate.mass < 0) {
        tmpTemplate.mass = tmpTemplate.radius;
        std::cout << "mass auto calculated from the radius: " << tmpTemplate.mass << std::endl;
    }

    std::cout << std::endl;
    return tmpTemplate;
}
