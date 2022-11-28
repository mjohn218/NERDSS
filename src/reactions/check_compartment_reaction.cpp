#include "classes/class_Membrane.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"
#include <algorithm>

// pro2Index is a lipid's index.

void check_compartment_reaction(int pro1Index, int pro2Index, int simItr,
    const Parameters& params, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<TransmissionRxn>& transmissionRxns,
    const std::vector<BackRxn>& backRxns, copyCounters& counterArrays, Membrane& membraneObject, std::vector<double>& IL2DbindingVec, std::vector<double>& IL2DUnbindingVec, std::vector<double>& ILTableIDs)
{

    /*
    * Calculate distance between molecule and compartment.
    * Establish whether entering or exiting.  (Determined by molecule type.)
    */


    int pro1MolType = moleculeList[pro1Index].molTypeIndex;

    if (molTemplateList[pro1MolType].crossesCompartment == false) {
        return;
    }

    // Find distance between molecule's COM and compartment origin.
    bool isEntering = false;
    double distanceToOrigin = moleculeList[pro1Index].comCoord.get_magnitude();
    if (distanceToOrigin > membraneObject.compartmentR) {
        isEntering = true;
    }

    if (isEntering) {
        int rxnIndex { -1 };
        rxnIndex = molTemplateList[pro1MolType].transmissionRxnIndex;

        double Dtot = 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.x +
                                   complexList[moleculeList[pro1Index].myComIndex].D.y +
                                   complexList[moleculeList[pro1Index].myComIndex].D.z);
        double cf = cos(sqrt(4.0 * complexList[moleculeList[pro1Index].myComIndex].Dr.z * params.timeStep));
        double Dr1;
        int relIface1 {transmissionRxns[rxnIndex].reactantListNew[0].relIfaceIndex};
        double relIfaceDistance = moleculeList[pro1Index].interfaceList[relIface1].coord.get_magnitude();

        Vector ifaceVec { moleculeList[pro1Index].interfaceList[relIface1].coord
                        - complexList[moleculeList[pro1Index].myComIndex].comCoord };
        double magMol1 { ifaceVec.x * ifaceVec.x + ifaceVec.y * ifaceVec.y + ifaceVec.z * ifaceVec.z };
        Dr1 = 2.0 * magMol1 * (1.0 - cf);
        Dtot += Dr1 / (6.0 * params.timeStep);
        Dtot += membraneObject.droplet.D;


        double onRate {transmissionRxns[rxnIndex].rateList[0].rate};
        double distToCompartment{relIfaceDistance - membraneObject.compartmentR};
        //this definition of Rmax may not need bindRadius depending on the binding model.

        double Rmax { 3.0 * sqrt(6.0 * Dtot * params.timeStep) + transmissionRxns[rxnIndex].bindRadius };

        if (distToCompartment < Rmax){
            // Checking if molecule is part of a complex and if so applying boundary condition of compartment
            if (complexList[moleculeList[pro1Index].myComIndex].memberList.size() > 1) {
                moleculeList[pro1Index].transmissionProb = 0;
                moleculeList[pro1Index].enforceCompartmentBC = true;
                return;
            }
            //Entering probability
            moleculeList[pro1Index].transmissionProb = 0;
            determine_entering_compartment_probability(distToCompartment, transmissionRxns, rxnIndex, pro1Index, moleculeList, Dtot, params, membraneObject);

        } else {
          moleculeList[pro1Index].transmissionProb = -1;
        }


    } else {
        // check the exiting
        int rxnIndex{-1};
        rxnIndex = molTemplateList[pro1MolType].transmissionRxnIndex;

        double Dtot = 1.0 / 3.0 * (complexList[moleculeList[pro1Index].myComIndex].D.x + complexList[moleculeList[pro1Index].myComIndex].D.y + complexList[moleculeList[pro1Index].myComIndex].D.z);
        double cf = cos(sqrt(4.0 * complexList[moleculeList[pro1Index].myComIndex].Dr.z * params.timeStep));
        double Dr1;
        int relIface1{transmissionRxns[rxnIndex].reactantListNew[0].relIfaceIndex};
        double relIfaceDistance = moleculeList[pro1Index].interfaceList[relIface1].coord.get_magnitude();

        Vector ifaceVec{moleculeList[pro1Index].interfaceList[relIface1].coord - complexList[moleculeList[pro1Index].myComIndex].comCoord};
        double magMol1{ifaceVec.x * ifaceVec.x + ifaceVec.y * ifaceVec.y + ifaceVec.z * ifaceVec.z};
        Dr1 = 2.0 * magMol1 * (1.0 - cf);
        Dtot += Dr1 / (6.0 * params.timeStep);
        Dtot += membraneObject.droplet.D;

        double onRate{transmissionRxns[rxnIndex].rateList[0].rate};
        double distToCompartment{membraneObject.compartmentR - relIfaceDistance};
        //this definition of Rmax may not need bindRadius depending on the binding model.

        double Rmax{3.0 * sqrt(6.0 * Dtot * params.timeStep) + transmissionRxns[rxnIndex].bindRadius};

        if (distToCompartment < Rmax)
        {
            // Checking if molecule is part of a complex and if so applying boundary condition of compartment
            if (complexList[moleculeList[pro1Index].myComIndex].memberList.size() > 1)
            {
                moleculeList[pro1Index].transmissionProb = 0;
                moleculeList[pro1Index].enforceCompartmentBC = true;
                return;
            }
            //Entering probability
            moleculeList[pro1Index].transmissionProb = 0;
            // TODO: This is to be defined
            determine_exiting_compartment_probability(distToCompartment, transmissionRxns, rxnIndex, pro1Index, moleculeList, Dtot, params, membraneObject);
        }
        else
        {
            moleculeList[pro1Index].transmissionProb = -1;
        }
    }
}
