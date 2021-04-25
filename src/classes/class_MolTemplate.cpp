/*! \file class_moltemplate_functions

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

#include "classes/class_MolTemplate.hpp"
#include "parser/parser_functions.hpp"

#include <iostream>

/* INTERFACE */
// Static Variables
int Interface::State::totalNumOfStates = 0;
std::vector<int> MolTemplate::absToRelIface {};
unsigned MolTemplate::numMolTypes = 0;
std::vector<int> MolTemplate::numEachMolType {};

// Constructors
Interface::State::State(int index)
    : index(index)
{
}

Interface::Interface(std::string name, const Coord& iCoord)
    : iCoord(iCoord)
    , name(std::move(name))
{
}

Interface::State::State(const std::string& ifaceAndStateName, char iden, int index)
    : ifaceAndStateName(ifaceAndStateName)
    , iden(iden)
    , index(index)
{
}
Interface::State::State(const std::string& ifaceAndStateName, int index)
    : ifaceAndStateName(ifaceAndStateName)
    , index(index)
{
}

Interface::Interface(std::string name, std::vector<Interface::State> states, Coord iCoord)
    : iCoord(iCoord)
    , name(name)
    , stateList(states)
{
    set_ifaceAndStateNames();
}

// Functions
void Interface::set_ifaceAndStateNames()
{
    // set the full name of the iface+state automatically during construction
    if (stateList.size() > 1) {
        for (auto& state : stateList) {
            if (state.iden != '\0') {
                state.ifaceAndStateName = name + "~" + state.iden;
            }
        }
    }
}
/*******************/

/* MOLTEMPLATE */

// Constructors
MolTemplate::MolTemplate(std::string molName, std::vector<Interface>& interfaceList)
    : molName(molName)
    , interfaceList(interfaceList)
{
}

MolTemplate::MolTemplate(Coord& comCoord, std::vector<Interface>& Interfaces)
    : comCoord(comCoord)
    , interfaceList(Interfaces)
{
}

// Functions
void MolTemplate::display() const
{
    std::cout << "Molecule template " << molTypeIndex << '\n';
    std::cout << "Name: " << molName << '\n';
    std::cout << "Copy number:" << copies << '\n';
    std::cout <<" Diffusion trans: "<<D.x <<' '<<D.y<<' '<<D.z<<'\n';
    std::cout <<" Diffusion Rot: "<<Dr.x <<' '<<Dr.y<<' '<<Dr.z<<'\n';
    
    if (isLipid) {
        if (isImplicitLipid)
            std::cout << "Is a implicitLipid: " << std::boolalpha << isImplicitLipid << '\n';
        else
            std::cout << "Is a lipid: " << std::boolalpha << isLipid << '\n';
    }
    std::cout << "Is a rod: " << std::boolalpha << isRod << '\n';
    std::cout << "Is a point: " << std::boolalpha << isPoint << '\n';
    std::cout << "Radius: " << radius << '\n';
    std::cout << "\nInterfaces:\n";
    for (auto& iface : interfaceList) {
        if (iface.stateList.size() == 1) {
            std::cout << "Name: " << iface.name << '\n';
            std::cout << "Relative index: " << iface.index << '\n';
            std::cout << "Absolute index: " << iface.stateList[0].index << '\n';
            if (!iface.stateList[0].myForwardRxns.empty()) {
                std::cout << "Forward Reactions: ";
                for (auto& rxn : iface.stateList[0].myForwardRxns)
                    std::cout << " (" << rxn << ")";
                std::cout << '\n';
            }
            if (!iface.stateList[0].myCreateDestructRxns.empty()) {
                std::cout << "Creation/Destruction Reactions: ";
                for (auto& rxn : iface.stateList[0].myCreateDestructRxns)
                    std::cout << " (" << rxn << ")";
                std::cout << '\n';
            }
            std::cout << '\n';
        } else {
            std::cout << "Name: " << iface.name << '\n';
            std::cout << "Relative index: " << iface.index << '\n';
            std::cout << "States:" << '\n';
            for (auto& state : iface.stateList) {
                std::cout << ""
                          << "(identity: " << state.iden << ", absolute index: " << state.index << ")\n";
                //                if (!state.myForwardRxns.empty()) {
                //                    std::cout << "Forward Reactions: ";
                //                    for (auto& rxn : state.myForwardRxns)
                //                        for (auto& rxn : iface.stateList[0].myCreateDestructRxns)
                //                    std::cout << '\n';
                //                }
                //                if (!state.stateChangeRxns.empty()) {
                //                    std::cout << "State Change Reactions: ";
                //                    for (auto& rxn : state.stateChangeRxns)
                //                        for (auto& rxn : iface.stateList[0].myCreateDestructRxns)
                //                    std::cout << '\n';
                //                }
                //                if (!state.myCreateDestructRxns.empty()) {
                //                    std::cout << "Creation/Destruction Reactions: ";
                //                    for (auto& rxn : state.myCreateDestructRxns)
                //                        for (auto& rxn : iface.stateList[0].myCreateDestructRxns)
                //                    std::cout << '\n';
                //                }
                std::cout << '\n';
            }
        }
    }
}

void MolTemplate::display(const std::string& name) const
{
    std::cout << "MolTemplate " << name << ':' << '\n';
    std::cout << "Rod? " << std::boolalpha << isRod << '\n';
    std::cout << "Radius: " << radius << '\n';
    std::cout << comCoord << '\n';
    for (auto& iface : interfaceList)
        std::cout << iface.iCoord << '\n';
    std::cout << '\n';
}

int MolTemplate::find_relIndex_from_absIndex(int targStateIndex) const
{
    /*!
     * \brief Finds the relative index of Interface given the absolute index of the State
     */

    // subroutine to find the Interface index from the State index
    int i { 0 };
    for (auto& iface : interfaceList) {
        for (auto& state : iface.stateList) {
            if (state.index == targStateIndex) {
                return i;
            }
        }
        ++i;
    }
    exit(1);
    //    not_in_cont_err(__PRETTY_FUNCTION__, __LINE__);
    return 0; // it'll terminate before it hits this, this is just so it
    // doesn't yell at me for not returning an int
}

int MolTemplate::find_absIndex_from_relIndex(int relIndex, char state) const
{
    for (auto& iface : interfaceList) {
        for (auto& tempState : iface.stateList) {
            if (state == tempState.iden && relIndex == iface.index)
                return tempState.index;
        }
    }

    std::cerr << "Absolute index for state " << state << " cannot be found on interface " << relIndex << ".\n";
    exit(1);
}

void MolTemplate::set_value(std::string& line, MolKeyword molKeyword)
{
    /*! \ingroup Parser
     * \brief sets the value of a MolTemplate's variable based on the keyword parsed from the parameter input file
     */

    auto key = static_cast<std::underlying_type<MolKeyword>::type>(molKeyword);
    switch (key) {
    case 0: {
        molName = line;
        break;
    }
    case 1: {
        copies = std::stoi(line);
        break;
    }
    case 2: {
        isRod = read_boolean(line);
        break;
    }
    case 3: {
        isLipid = read_boolean(line);
        std::cout << "Read in isLipid: " << std::boolalpha << isLipid << std::endl;
        break;
    }
    case 4: {
        D = Coord(parse_input_array(line));
        std::cout << "Read in D: [" << D.x << "um^2s^-1, " << D.y << "um^2s^-1, " << D.z << "um^2s^-1]" << std::endl;
        break;
    }
    case 5: {
        Dr = Coord(parse_input_array(line));
        std::cout << "Read in Dr: [" << Dr.x << "rad^2s^-1, " << Dr.y << "rad^2s^-1, " << Dr.z << "rad^2s^-1]" << std::endl;
        break;
    }
    case 8: {
        mass = std::stod(line);
        std::cout << "Read in mass: " << std::boolalpha << mass << std::endl;
        break;
    }
    case 9: {
        checkOverlap = read_boolean(line);
        std::cout << "Read in checkOverlap: " << std::boolalpha << checkOverlap << std::endl;
        break;
    }
    case 11: {
        isPoint = read_boolean(line);
        break;
    }
    case 12: {
        isImplicitLipid = read_boolean(line);
        if (isImplicitLipid == true) {
            isLipid = true;
        }
        std::cout << "Read in isImplicitLipid: " << std::boolalpha << isImplicitLipid << std::endl;
        break;
    }
    case 13: {
        countTransition = read_boolean(line);
        std::cout << "Read in countTransition: " << std::boolalpha << countTransition << std::endl;
        break;
    }
    case 14: {
        transitionMatrixSize = std::stoi(line);
        break;
    }
    default: {
        std::cout << "Keyword [BLANK] is not a valid keyword, ignoring." << '\n';
    }
    }
}

size_t find_molTypeIndex_from_ifaceIndex(const int& targIfaceIndex, const std::vector<MolTemplate>& molTemplateList)
{
    /* this function should find the molTypeIndex of the parent MolTemplate of am absolute ifaceIndex */
    try {
        for (auto& oneTemplate : molTemplateList) {
            for (auto& oneIface : oneTemplate.interfaceList) {
                for (auto& oneState : oneIface.stateList) {
                    if (oneState.index == targIfaceIndex)
                        return oneTemplate.molTypeIndex;
                }
            }
        }
        throw std::out_of_range("Interface::State index cannot be found in list of MolTemplates...");
    } catch (const std::out_of_range& e) {
        // TODO: write some handling here
        exit(1);
    }
}

/*******************/
