#include "io/io.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "system_setup/system_setup.hpp"
#include <numeric>

//this is used to generate the added molecules and complexes
void generate_coordinates_for_restart(Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, std::vector<MolTemplate>& molTemplateList,
    const std::vector<ForwardRxn>& forwardRxns, const Membrane& membraneObject, int numMolTemplateBeforeAdd, int numForwardRxnBdeforeAdd)
{
    for (int molTemplateListIndex = numMolTemplateBeforeAdd; molTemplateListIndex < molTemplateList.size(); molTemplateListIndex++) {
        MolTemplate oneTemp {};
        oneTemp = molTemplateList[molTemplateListIndex];
        for (unsigned itr { 0 }; itr < oneTemp.copies; ++itr) {
            create_molecule_and_complex_for_restart(oneTemp, params, moleculeList, complexList, molTemplateList, forwardRxns, membraneObject);
        }
    }
}

void create_molecule_and_complex_for_restart(MolTemplate& createdMolTemp, Parameters& params, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns, const Membrane& membraneObject)
{
    int newMolIndex = 0;
    int newComIndex = 0;
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
    while (needsResampling) {
        moleculeList[newMolIndex]
            = initialize_molecule_for_restart(newMolIndex, params, createdMolTemp, membraneObject);
        needsResampling = moleculeOverlapsForRestart(params, moleculeList[newMolIndex], moleculeList,
            complexList, forwardRxns, molTemplateList, membraneObject);
    }

    moleculeList[newMolIndex].myComIndex = newComIndex;
    complexList[newComIndex] = Complex { newComIndex, moleculeList.at(newMolIndex), createdMolTemp };
    ++Complex::numberOfComplexes;

    createdMolTemp.monomerList.emplace_back(newMolIndex); //add this new molecule to the monomerList
}

Molecule initialize_molecule_for_restart(
    int index, Parameters& params, MolTemplate& molTemplate, const Membrane& membraneObject)
{
    /*!
     * \brief Creates a molecule according to add.inp, assigns
     * coords.
     *
     * params[in] molTemplate MolTemplate of the molecule to be created
     */

    Molecule tmp {};
    // TODO: Need to finish this, state should be set to reaction's product's state
    tmp.molTypeIndex = molTemplate.molTypeIndex;
    tmp.mass = molTemplate.mass;
    tmp.isLipid = molTemplate.isLipid;

    // Set up interface state vectors
    tmp.freelist = std::vector<int>(molTemplate.interfaceList.size());
    std::iota(tmp.freelist.begin(), tmp.freelist.end(), 0);
    tmp.interfaceList = std::vector<Molecule::Iface>(molTemplate.interfaceList.size());

    // Create center of mass
    tmp.create_random_coords(molTemplate, membraneObject);

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

bool moleculeOverlapsForRestart(const Parameters& params, Molecule& createdMol,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject)
{
    // if it's inside the box, check if it overlaps with any molecule
    std::vector<unsigned> checkedMols {};
    for (auto memMol : moleculeList) {
        if (memMol.index == createdMol.index) {
            continue;
        }

        const Complex& oneCom = complexList[memMol.myComIndex]; // legibility

        // check bounding sphere
        Vector tmpVec { createdMol.comCoord - oneCom.comCoord };
        tmpVec.calc_magnitude();
        if (tmpVec.magnitude > (molTemplateList[createdMol.molTypeIndex].radius + oneCom.radius))
            return false;

        // TODO: this is awful
        for (auto& comMemMol : oneCom.memberList) {
            // check if the two Molecules can even react first

            const Molecule& partMol = moleculeList[comMemMol];
            for (const auto& oneRxn : forwardRxns) {
                for (unsigned iface1Itr { 0 }; iface1Itr < createdMol.interfaceList.size(); ++iface1Itr) {
                    for (unsigned iface2Itr { 0 }; iface2Itr < partMol.interfaceList.size(); ++iface2Itr) {
                        if (isReactant(createdMol.interfaceList[iface1Itr], createdMol, oneRxn.reactantListNew[0])
                            && isReactant(partMol.interfaceList[iface2Itr], partMol, oneRxn.reactantListNew[1])) {

                            // if they're reactants, check if they're within the binding radius
                            Vector ifaceVec { createdMol.interfaceList[iface1Itr].coord
                                - partMol.interfaceList[iface2Itr].coord };
                            ifaceVec.calc_magnitude();

                            if (ifaceVec.magnitude > oneRxn.bindRadius)
                                return true;
                        } else if (isReactant(
                                       createdMol.interfaceList[iface1Itr], createdMol, oneRxn.reactantListNew[1])
                            && isReactant(partMol.interfaceList[iface2Itr], partMol, oneRxn.reactantListNew[0])) {

                            // if they're reactants, check if they're within the binding radius
                            Vector ifaceVec { createdMol.interfaceList[iface1Itr].coord
                                - partMol.interfaceList[iface2Itr].coord };
                            ifaceVec.calc_magnitude();

                            if (ifaceVec.magnitude > oneRxn.bindRadius)
                                return true;
                        } else {
                            continue;
                        }
                    }
                }
            }
        }
    }

    return false;
}