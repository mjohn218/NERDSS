#include "parser/parser_functions.hpp"

ParsedMol parse_molecule_bngl(int& totSpecies, bool isProductSide,
    std::pair<std::string, int> oneMol) // totSpecies is not altered in this routine
{
    std::string buffer;
    ParsedMol tmpMol;
    tmpMol.specieIndex = oneMol.second;
    for (auto molIterator = oneMol.first.begin(); molIterator != oneMol.first.end(); std::advance(molIterator, 1)) {
        if (isalnum(*molIterator))
            // if the character is alphanumeric, just append it to the buffer
            buffer += *molIterator;
        else {
            switch (*molIterator) {
            case '~': {
                // if the character is '~' it indicates a state. add an iface to the vector with
                // iface~state as the ifaceName and state as the state
                std::string iface { buffer };
                ++molIterator;
                buffer += "~";
                buffer += static_cast<char>(std::toupper(*molIterator));
                if (*(molIterator + 1) != '!') {
                    tmpMol.interfaceList.emplace_back(iface, static_cast<char>(std::toupper(*molIterator)), false,
                        Involvement::possible, oneMol.second);
                    buffer.clear(); // flush the buffer, no bond
                }
                break;
            }
            case '(': {
                /* -if the character is '(', it indicates the beginning of an iface list belong to a
                 * protein reactant set the molecule name to the buffer before the parenthesis
                 */
                tmpMol.molName = buffer;
                buffer.clear();
                break;
            }
            case '!': { // '!' indicates a bond follows, with * being a wildcard
                ++molIterator;
                std::string iface { buffer };
                buffer.clear();

                // look for the bond index or wildcard
                while ((*molIterator != ')') && (*molIterator != ',')) {
                    buffer += *molIterator;
                    ++molIterator;
                }

                bool isWildcard { (buffer.size() == 1) && (buffer[0] == '*') };
                auto stateLocation = iface.find('~');

                /* -need to check if it is the product side or reactant side for two reasons:
                 * 1) theoretically, bond indices could be double digit (though unlikely)
                 * 2) only wildcards can exist in the reactant side, no bond indices
                 */
                if (!isProductSide) {
                    /* -if product side, only look for wildcard bonds, and add to
                     * ancillaryIfaceList, since they can't create new bonds
                     * -unless it has a state then, add to possiblyInvolvedIfaceList, since the state could change
                     */
                    if (!isWildcard) {
                        // TODO: Write this
                        std::cerr
                            << "Error, no indexed interactions are allowed in the reactants. Converting to wildcard.\n";
                        exit(1);
                    }

                    // if a state existed, create an Iface with that state required
                    if (stateLocation != std::string::npos) {
                        auto state = iface[stateLocation + 1];
                        iface = iface.substr(0, stateLocation);
                        tmpMol.interfaceList.emplace_back(iface, state, true, 0, Involvement::possible, oneMol.second);
                    } else
                        tmpMol.interfaceList.emplace_back(iface, '\0', true, 0, Involvement::ancillary, oneMol.second);
                } else {
                    // if on the product side, you need to look for a bond index
                    /* -the gist of this is that if the character after the bond indicator is a
                     * wildcard, the iface belongs in the ancillaryIfaceList, since that means it
                     * couldn't have created a bond in the reaction
                     * -if it doesn't have a wildcard, check if there's a state, and create the Iface accordingly
                     */

                    if ((stateLocation != std::string::npos) && !isWildcard) {
                        auto state = iface[stateLocation + 1];
                        iface = iface.substr(0, stateLocation);
                        tmpMol.interfaceList.emplace_back(
                            iface, state, true, std::stoi(buffer), Involvement::interactionChange, oneMol.second);
                    } else if ((stateLocation != std::string::npos) && isWildcard) {
                        auto state = iface[stateLocation + 1];
                        iface = iface.substr(0, stateLocation);
                        tmpMol.interfaceList.emplace_back(iface, state, true, 0, Involvement::possible, oneMol.second);
                    } else if ((stateLocation == std::string::npos) && !isWildcard)
                        tmpMol.interfaceList.emplace_back(
                            iface, '\0', true, std::stoi(buffer), Involvement::interactionChange, oneMol.second);
                    else
                        tmpMol.interfaceList.emplace_back(iface, '\0', true, 0, Involvement::ancillary, oneMol.second);
                }
                buffer.clear();
                break;
            }
            case ',': {
                // If an unbound iface precedes this, make a new ancillary Iface with the information
                if (!buffer.empty())
                    tmpMol.interfaceList.emplace_back(buffer, '\0', false, Involvement::possible, oneMol.second);
                buffer.clear();
                break;
            }
            case ')': {
                // a ')' indicates the end of an iface list
                if (!buffer.empty()) { // if the buffer is empty, there's nothing to do
                    tmpMol.interfaceList.emplace_back(buffer, '\0', false, Involvement::possible, oneMol.second);
                    buffer.clear();
                }
                break;
            }
            default: {
                //                gen_read_err(__func__, __LINE__);
                std::cerr << "ERROR: Character " << *molIterator << " is not valid in reactions. Exiting...\n.";
                exit(1);
            }
            }
        }
    }
    return tmpMol;
}

ParsedMolNumState parse_number_bngl(std::string oneLine)
{
    std::string buffer;
    ParsedMolNumState tmpMolNum;
    for (auto molIterator = oneLine.begin(); molIterator != oneLine.end(); std::advance(molIterator, 1)) {
        if (isalnum(*molIterator))
            // if the character is alphanumeric, just append it to the buffer
            buffer += *molIterator;
        else {
            switch (*molIterator) {
            case '~': {
                // if the character is '~' it indicates a state. add an iface to the vector with
                // iface~state as the ifaceName and state as the state
                std::string iface { buffer };
                ++molIterator;
                buffer += "~";
                buffer += static_cast<char>(std::toupper(*molIterator));
                break;
            }
            case '(': {
                /* -if the character is '(', it indicates the beginning of an iface list
                 * set the copy number to the buffer before the parenthesis
                 */
                tmpMolNum.totalCopyNumbers += std::stoi(buffer);
                tmpMolNum.numberEachState.emplace_back(std::stoi(buffer));
                buffer.clear();
                break;
            }
            case ',': {
                if (!buffer.empty())
                    buffer += *molIterator;
                break;
            }
            case ')': {
                // a ')' indicates the end of an iface list
                if (!buffer.empty()) { // if the buffer is empty, there's nothing to do
                    tmpMolNum.nameEachState.emplace_back(buffer);
                    buffer.clear();
                }
                break;
            }
            default: {
                std::cerr << "ERROR: Character " << *molIterator << " is not valid in starting copy numbers for each state. Exiting...\n.";
                exit(1);
            }
            }
        }
    }
    if (!buffer.empty()) {
        tmpMolNum.totalCopyNumbers += std::stoi(buffer);
    }
    return tmpMolNum;
}
