#include "math/rand_gsl.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"
#include <numeric>

Molecule initialize_molecule_after_transmission_reaction(int index, const Molecule& parentMol, Parameters& params,
    MolTemplate& molTemplate, const TransmissionRxn& currRxn, const Coord& newPos, const bool plusRand, const Membrane& membraneObject)
{
    // TRACE();
    Molecule tmp {};
    // TODO: Need to finish this, state should be set to reaction's product's state
    tmp.molTypeIndex = molTemplate.molTypeIndex;
    tmp.mass = molTemplate.mass;
    tmp.isLipid = molTemplate.isLipid;

    // Set up interface state vectors
    tmp.freelist = std::vector<int>(molTemplate.interfaceList.size());
    std::iota(tmp.freelist.begin(), tmp.freelist.end(), 0);
    tmp.interfaceList = std::vector<Molecule::Iface>(molTemplate.interfaceList.size());

    Coord tempCoord { newPos };
    double molRadius = molTemplate.radius;
    if (plusRand) {
        bool retryNewPos = true;
        // put the specie in a random position around the specie which created it
        while (retryNewPos) {
            double theta { rand_gsl() * 2 * M_PI };
            double phi { std::acos(rand_gsl() * 2 - 1) };
            double cosTheta { std::cos(theta) };
            double sinTheta { std::sin(theta) };
            double cosPhi { std::cos(phi) };
            double sinPhi { std::sin(phi) };
            Coord transVec { molRadius * cosTheta * sinPhi, molRadius * sinTheta * sinPhi,
                molRadius * cosPhi }; // use mol radius here to make sure that the mol has been moved further enough
            tempCoord = newPos + transVec;
            // get distance from origin to check if new position is outside compartment
            double distanceToOrigin = tempCoord.get_magnitude();
            if ((distanceToOrigin > membraneObject.compartmentR) && (molTemplate.outsideCompartment)) {
                retryNewPos = false;
            } else if ((distanceToOrigin < membraneObject.compartmentR) && (molTemplate.insideCompartment)) {
                retryNewPos = false;
            }
        }
    }
    // move molecule to the new coordinates "newPos"
    tmp.comCoord = tempCoord;
    for (unsigned ifaceItr { 0 }; ifaceItr < molTemplate.interfaceList.size(); ++ifaceItr) {
        tmp.interfaceList[ifaceItr].coord = molTemplate.interfaceList[ifaceItr].iCoord + tmp.comCoord;
        tmp.interfaceList[ifaceItr].index = molTemplate.interfaceList[ifaceItr].stateList[0].index;
        tmp.interfaceList[ifaceItr].relIndex = ifaceItr;
        tmp.interfaceList[ifaceItr].stateIden = molTemplate.interfaceList[ifaceItr].stateList[0].iden;
        tmp.interfaceList[ifaceItr].stateIndex = 0;
        tmp.interfaceList[ifaceItr].molTypeIndex = molTemplate.molTypeIndex;
    }

    // set the interface states to the states defined in the reaction
    if (currRxn.productMolList.back().molTypeIndex == molTemplate.molTypeIndex) {
        for (auto& rxnIface : currRxn.productMolList.back().interfaceList) {
            tmp.interfaceList[rxnIface.relIfaceIndex].stateIden = rxnIface.requiresState;
            // set the stateIndex by finding the state in the MolTemplate Interface's stateList
            for (unsigned stateItr { 0 }; stateItr < molTemplate.interfaceList[rxnIface.relIfaceIndex].stateList.size();
                 ++stateItr) {
                if (molTemplate.interfaceList[rxnIface.relIfaceIndex].stateList[stateItr].iden
                    == rxnIface.requiresState) {
                    tmp.interfaceList[rxnIface.relIfaceIndex].stateIndex = stateItr;
                    break;
                }
            }
        }
    }

    // clean up
    tmp.isEmpty = false;

    // iterate number of molecules in the system and set index
    tmp.index = index;
    ++Molecule::numberOfMolecules;
    params.numTotalUnits = params.numTotalUnits + molTemplate.interfaceList.size() + 1;

    // keep track of molecule types
    ++MolTemplate::numEachMolType[molTemplate.molTypeIndex];

    return tmp;
}
