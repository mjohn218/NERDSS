#include "boundary_conditions/reflect_functions.hpp"
#include "classes/class_Rxns.hpp"
#include "io/io.hpp"
#include "reactions/association/association.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"

#include <cmath>
#include <iomanip>

void perform_transmission_reaction(int moleculeIndex, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, Membrane& membraneObject, copyCounters& counterArrays, Parameters& params,
    std::vector<TransmissionRxn>& transmissionRxns, const std::vector<ForwardRxn>& forwardRxns, SimulVolume& simulVolume) {

    /*
    * Calculate distance between molecule and compartment.
    * Establish whether entering or exiting.  (Determined by molecule type.)
    */


    int pro1MolType = moleculeList[moleculeIndex].molTypeIndex;

    // Find distance between molecule's COM and compartment origin.
    bool isEntering = false;
    double distanceToOrigin = moleculeList[moleculeIndex].comCoord.get_magnitude();
    int rxnIndex {molTemplateList[pro1MolType].transmissionRxnIndex};
    int molTypeIndex2 {transmissionRxns[rxnIndex].productMolList.back().molTypeIndex};
    double distToCompartment {distanceToOrigin - membraneObject.compartmentR};
    double displacement {2*distToCompartment};
    Coord currPos {moleculeList[moleculeIndex].comCoord};
    double scaler {(distanceToOrigin-displacement)/distanceToOrigin};
    Coord newPos {scaler*currPos};
    MolTemplate& oneTemp = molTemplateList[pro1MolType];

    //create new molecule at newPos
    int newMolIndex { 0 };
    int newComIndex { 0 };
    create_molecule_and_complex_from_transmission_rxn(moleculeIndex, newMolIndex, newComIndex,
                                                    molTemplateList[molTypeIndex2], params, transmissionRxns[rxnIndex],
                                                    simulVolume, moleculeList, complexList, molTemplateList, forwardRxns, membraneObject, newPos);

    //update the counterArray 
	for (const auto& iface : moleculeList[newMolIndex].interfaceList)
	  ++counterArrays.copyNumSpecies[iface.index];

    // forbid this molecule to react with other molecules
    for (unsigned crossItr { 0 }; crossItr < moleculeList[moleculeIndex].crossbase.size(); ++crossItr) {
        int skipMol { moleculeList[moleculeIndex].crossbase[crossItr] };
        for (unsigned crossItr2 { 0 }; crossItr2 < moleculeList[skipMol].crossbase.size(); ++crossItr2) {
            if (moleculeList[skipMol].crossbase[crossItr2] == moleculeList[moleculeIndex].index){
                moleculeList[skipMol].probvec[crossItr2] = 0.0;
            }
        }
    }
    complexList[moleculeList[moleculeIndex].myComIndex].ncross = -1;
    moleculeList[moleculeIndex].crossbase.clear();

    //destroy molecule
    for (auto& memMol : complexList[moleculeList[moleculeIndex].myComIndex].memberList) {
        for (auto& iface : moleculeList[memMol].interfaceList) {
            --counterArrays.copyNumSpecies[iface.index];
        }
    }
    
    complexList[moleculeList[moleculeIndex].myComIndex].destroy(moleculeList, complexList); //This is not finished

    // remove the molecule from the SimulVolume subsCellList
    // have this here to avoid circular header calls with SimulVolume and Molecule_Complex
    int molItr { moleculeIndex };
    for (auto itr = simulVolume.subCellList[moleculeList[molItr].mySubVolIndex].memberMolList.begin(); itr != simulVolume.subCellList[moleculeList[molItr].mySubVolIndex].memberMolList.end(); ++itr) {
        if (*itr == moleculeList[molItr].index) {
            simulVolume.subCellList[moleculeList[molItr].mySubVolIndex].memberMolList.erase(itr);
            break;
        }
    }
    moleculeList[molItr].mySubVolIndex = -1;

    //Also need to adjust the monomer list, only if this molecule can be destroyed. 
	if(oneTemp.canDestroy == true)
	  oneTemp.monomerList.erase(std::find_if(oneTemp.monomerList.begin(), oneTemp.monomerList.end(), [&](const size_t& mol) { return mol == molItr; }));

}
