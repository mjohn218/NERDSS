#include "system_setup/system_setup.hpp"

void initialize_states(std::vector<Molecule>& moleculeList, std::vector<MolTemplate>& molTemplateList, Membrane& membraneObject)
{
    for (auto& tmpMol : molTemplateList) {
        //std::cout << "molName: " << tmpMol.molName << std::endl;
        //std::cout << "copies: " << tmpMol.copies << std::endl;
        int tmpLength { 0 };
        tmpLength = tmpMol.startingNumState.nameEachState.size();
        for (int i = 0; i < tmpLength; i++) {
            //std::cout << tmpMol.startingNumState.nameEachState[i] << ": " << tmpMol.startingNumState.numberEachState[i] << std::endl;
        }

        if (tmpLength > 0) {
            if (tmpMol.isImplicitLipid == true) {
                // for the implicit case, update the membraneObject
                std::vector<int> tmpCopyNumbers {};
                std::vector<std::string> tmpIfaceName {};
                std::vector<char> tmpStateName {};
                for (int i = 0; i < tmpLength; i++) {
                    // parse the nameEachState, for the implicit case, only one interface can exist
                    std::string oneLine { tmpMol.startingNumState.nameEachState[i] };
                    std::string buffer;
                    tmpCopyNumbers.emplace_back(tmpMol.startingNumState.numberEachState[i]);

                    for (auto molIterator = oneLine.begin(); molIterator != oneLine.end(); std::advance(molIterator, 1)) {
                        if (isalnum(*molIterator))
                            // if the character is alphanumeric, just append it to the buffer
                            buffer += *molIterator;
                        else {
                            switch (*molIterator) {
                            case '~': {
                                // if the character is '~' it indicates a state. add an iface to the vector with
                                // iface~state as the ifaceName and state as the state
                                tmpIfaceName.emplace_back(buffer);
                                ++molIterator;
                                buffer += "~";
                                buffer += static_cast<char>(std::toupper(*molIterator));
                                tmpStateName.emplace_back(static_cast<char>(std::toupper(*molIterator)));
                                break;
                            }
                            case ',': {
                                std::cerr << "ERROR: Implicit molecule can only have one interface. Exiting...\n.";
                                exit(1);
                                break;
                            }
                            default: {
                                std::cerr << "ERROR: Character " << *molIterator << " is not valid in starting copy numbers for each state. Exiting...\n.";
                                exit(1);
                            }
                            }
                        }
                    }
                }

                // initialize copy numbers according to State
                for (int i = 0; i < tmpLength; i++) {
                    int stateNum = tmpMol.interfaceList[0].stateList.size();
                    for (int stateIndex = 0; stateIndex < stateNum; stateIndex++) {
                        if (tmpMol.interfaceList[0].stateList[stateIndex].iden == tmpStateName[i]) {
                            membraneObject.numberOfFreeLipidsEachState[stateIndex] = tmpCopyNumbers[i];
                        }
                    }
                }
            } else {
                // for the non-implicit case
                std::vector<int> tmpCopyNumbers {};
                std::vector<std::vector<std::string>> tmpIfaceName {};
                std::vector<std::vector<char>> tmpStateName {};
                for (int i = 0; i < tmpLength; i++) {
                    // parse the nameEachState, for the implicit case, only one interface can exist
                    std::string oneLine { tmpMol.startingNumState.nameEachState[i] };
                    std::string buffer;
                    tmpCopyNumbers.emplace_back(tmpMol.startingNumState.numberEachState[i]);
                    tmpIfaceName.emplace_back();
                    tmpStateName.emplace_back();

                    for (auto molIterator = oneLine.begin(); molIterator != oneLine.end(); std::advance(molIterator, 1)) {
                        if (isalnum(*molIterator))
                            // if the character is alphanumeric, just append it to the buffer
                            buffer += *molIterator;
                        else {
                            switch (*molIterator) {
                            case '~': {
                                // if the character is '~' it indicates a state. add an iface to the vector with
                                // iface~state as the ifaceName and state as the state
                                tmpIfaceName.back().emplace_back(buffer);
                                ++molIterator;
                                buffer += "~";
                                buffer += static_cast<char>(std::toupper(*molIterator));
                                tmpStateName.back().emplace_back(static_cast<char>(std::toupper(*molIterator)));
                                break;
                            }
                            case ',': {
                                buffer.clear();
                                break;
                            }
                            default: {
                                std::cerr << "ERROR: Character " << *molIterator << " is not valid in starting copy numbers for each state. Exiting...\n.";
                                exit(1);
                            }
                            }
                        }
                    }
                }

                // initialize starting copy numbers according to states
                int tmpCount { 0 };
                int tmpCopyNumbersIndex { -1 };
                int firstIndex { -1 };
                for (auto& oneMol : moleculeList) {
                    if (oneMol.molTypeIndex == tmpMol.molTypeIndex) {
                        firstIndex = oneMol.index;
                        break;
                    }
                }
                for (int i = 0; i < tmpMol.startingNumState.totalCopyNumbers; i++) {
                    tmpCount = 0;
                    tmpCopyNumbersIndex = -1;
                    while (true) {
                        if (i < tmpCount) {
                            break;
                        } else {
                            tmpCopyNumbersIndex++;
                            tmpCount += tmpCopyNumbers[tmpCopyNumbersIndex];
                        }
                    }

                    // update the index for the state
                    int molItr { -1 };
                    int ifaceItr { -1 };
                    int stateItr { -1 };
                    molItr = firstIndex + i;

                    int tmpIndex = tmpStateName[tmpCopyNumbersIndex].size();
                    for (int j = 0; j < tmpIndex; j++) {
                        std::string ifaceName = tmpIfaceName[tmpCopyNumbersIndex][j];
                        char stateName = tmpStateName[tmpCopyNumbersIndex][j];

                        for (auto& oneIface : tmpMol.interfaceList) {
                            if (oneIface.name == ifaceName) {
                                ifaceItr = oneIface.index;
                                break;
                            }
                        }

                        for (auto& oneState : tmpMol.interfaceList[ifaceItr].stateList) {
                            if (oneState.iden == stateName) {
                                stateItr = static_cast<int>(&oneState - &tmpMol.interfaceList[ifaceItr].stateList[0]);
                                break;
                            }
                        }

                        moleculeList[molItr].interfaceList[ifaceItr].index = tmpMol.interfaceList[ifaceItr].stateList[stateItr].index;
                        moleculeList[molItr].interfaceList[ifaceItr].relIndex = ifaceItr;
                        moleculeList[molItr].interfaceList[ifaceItr].stateIden = tmpMol.interfaceList[ifaceItr].stateList[stateItr].iden;
                        moleculeList[molItr].interfaceList[ifaceItr].stateIndex = stateItr;
                    }
                }
            }
        }
    }
}