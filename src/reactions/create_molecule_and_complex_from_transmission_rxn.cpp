#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "tracing.hpp"

void create_molecule_and_complex_from_transmission_rxn(int parentMolIndex, int& newMolIndex, int& newComIndex,
    MolTemplate& createdMolTemp, Parameters& params, const TransmissionRxn& currRxn, SimulVolume& simulVolume,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns, const Membrane& membraneObject,
    const Coord& newPos)
{
    newMolIndex = 0;
    newComIndex = 0;
    if (Molecule::emptyMolList.size() > 0) {
        // Check to make sure the object pointed to by the last element in the each list is actually empty
        // If it isn't delete it from the list and move on until you find one that is
        try {
            while (!moleculeList[Molecule::emptyMolList.back()].isEmpty)
                Molecule::emptyMolList.pop_back();

            // if there's an available empty Molecule spot, make the new Molecule in place
            newMolIndex = Molecule::emptyMolList.back();
            Molecule::emptyMolList.pop_back(); // remove the empty Molecule spot index from the list
        } catch (std::out_of_range& e) {
            newMolIndex = moleculeList.size();
            moleculeList.emplace_back(); // create new empty Molecule spot
        }
    } else {
        newMolIndex = moleculeList.size();
        moleculeList.emplace_back(); // create new empty Molecule spot
    }

    if (Complex::emptyComList.size() > 0) {
        // Check to make sure the object pointed to by the last element in the each list is actually empty
        // If it isn't delete it from the list and move on until you find one that is
        try {
            while (!complexList[Complex::emptyComList.back()].isEmpty)
                Complex::emptyComList.pop_back();
            // if there's an available empty Complex spot, make the new Complex in place
            newComIndex = Complex::emptyComList.back();
            Complex::emptyComList.pop_back(); // remove the empty Complex spot index from the list
        } catch (std::out_of_range) {
            newComIndex = complexList.size();
            complexList.emplace_back(); // create new empty Complex spot
        }

    } else {
        newComIndex = complexList.size();
        complexList.emplace_back(); // create new empty Complex spot
    }

    // Now create the new species
    
    bool needsResampling { true };
    int overlapCounter { 0 };
    while (needsResampling) {
        if (overlapCounter == 0) {
            moleculeList[newMolIndex] = initialize_molecule_after_transmission_reaction(
                newMolIndex, moleculeList[parentMolIndex], params, createdMolTemp, currRxn, newPos, false, membraneObject);
            needsResampling = moleculeOverlaps(params, simulVolume, moleculeList[newMolIndex], moleculeList,
            complexList, forwardRxns, molTemplateList, membraneObject);
        } else {
            moleculeList[newMolIndex] = initialize_molecule_after_transmission_reaction(
                newMolIndex, moleculeList[parentMolIndex], params, createdMolTemp, currRxn, newPos, true, membraneObject);
            needsResampling = moleculeOverlaps(params, simulVolume, moleculeList[newMolIndex], moleculeList,
                complexList, forwardRxns, molTemplateList, membraneObject);
        }
        ++overlapCounter;
    }

    moleculeList[newMolIndex].myComIndex = newComIndex;
    moleculeList[newMolIndex].trajStatus = TrajStatus::propagated;
    moleculeList[newMolIndex].isDissociated = true;
    complexList[newComIndex] = Complex { newComIndex, moleculeList.at(newMolIndex), createdMolTemp };
    complexList[newComIndex].trajStatus = TrajStatus::propagated;
    ++Complex::numberOfComplexes;

    // add to monomerList if canDestroy = true
    {
        Molecule& oneMol { moleculeList[newMolIndex] };
        MolTemplate& oneTemp { molTemplateList[oneMol.molTypeIndex] };
        if (oneTemp.canDestroy) {
            oneTemp.monomerList.emplace_back(oneMol.index);
        }
    }
}
