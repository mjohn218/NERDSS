#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "tracing.hpp"

bool moleculeOverlaps(const Parameters& params, SimulVolume& simulVolume, Molecule& createdMol,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject)
{
    // TRACE();
    // get which box the Molecule belongs to
    int xItr { int((createdMol.comCoord.x + membraneObject.waterBox.x / 2) / simulVolume.subCellSize.x) };
    int yItr { int((createdMol.comCoord.y + membraneObject.waterBox.y / 2) / simulVolume.subCellSize.y) };
    int zItr { int(-(createdMol.comCoord.z + 1E-6 - membraneObject.waterBox.z / 2.0) / simulVolume.subCellSize.z) };
    int currBin
        = xItr + (yItr * simulVolume.numSubCells.x) + (zItr * simulVolume.numSubCells.x * simulVolume.numSubCells.y);

    if (molTemplateList[createdMol.molTypeIndex].D.z == 0
        && std::abs(createdMol.comCoord.z) - std::abs((membraneObject.waterBox.z / 2)) > 1E-6) {
        // std::cerr << "Molecule " << createdMol.index << " of type " << molTemplateList[createdMol.molTypeIndex].molName
        //           << " is off the membrane. Writing coordinates and exiting.\n";
        return true;
    }

    // Now make sure the Molecule is still inside the box in all dimensions
    if (createdMol.comCoord.z > (membraneObject.waterBox.z / 2) || createdMol.comCoord.z + 1E-6 < -(membraneObject.waterBox.z / 2)) {
        // std::cout << "Molecule " << createdMol.index
        //           << " is outside simulation volume in the z-dimension, with center of mass coordinates ["
        //           << createdMol.comCoord << "]. Attempting to fit back into box.\n";
        return true;
    } else if (createdMol.comCoord.y > (membraneObject.waterBox.y / 2)
        || createdMol.comCoord.y + 1E-6 < -(membraneObject.waterBox.y / 2)) {
        // std::cout << "Molecule " << createdMol.index
        //           << " is outside simulation volume in the y-dimension, with center of mass coordinates ["
        //           << createdMol.comCoord << "]. Attempting to fit back into box.\n";
        return true;
    } else if (createdMol.comCoord.x > (membraneObject.waterBox.x / 2)
        || createdMol.comCoord.x + 1E-6 < -(membraneObject.waterBox.x / 2)) {
        // std::cout << "Molecule " << createdMol.index
        //           << " is outside simulation volume in the x-dimension, with center of mass coordinates ["
        //           << createdMol.comCoord << "]. Attempting to fit back into box.\n";
        return true;
    } else if (currBin > (simulVolume.numSubCells.tot) || currBin < 0) {
        // std::cout << "Molecule " << createdMol.index
        //           << " is outside simulation volume with center of mass coordinates [" << createdMol.comCoord
        //           << "]. Attempting to fit back into box.\n";
        return true;
    } else {
        // if it's inside the box, check if it overlaps with any molecule
        std::vector<unsigned> checkedMols {};
        for (auto memMol : simulVolume.subCellList[currBin].memberMolList) {
            const Complex& oneCom = complexList[moleculeList[memMol].myComIndex]; // legibility

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

        createdMol.mySubVolIndex = currBin;
        simulVolume.subCellList[currBin].memberMolList.push_back(createdMol.index);
        return false;
    }
}

void create_molecule_and_complex_from_rxn(int parentMolIndex, int& newMolIndex, int& newComIndex, bool createInVicinity,
    MolTemplate& createdMolTemp, Parameters& params, const CreateDestructRxn& currRxn, SimulVolume& simulVolume,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns, const Membrane& membraneObject)
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
    if (createInVicinity) {
        bool needsResampling { true };
        while (needsResampling) {
            moleculeList[newMolIndex] = initialize_molecule_after_uni_reaction(
                newMolIndex, moleculeList[parentMolIndex], params, createdMolTemp, currRxn);
            needsResampling = moleculeOverlaps(params, simulVolume, moleculeList[newMolIndex], moleculeList,
                complexList, forwardRxns, molTemplateList, membraneObject);
        }
    } else {
        bool needsResampling { true };
        while (needsResampling) {
            moleculeList[newMolIndex]
                = initialize_molecule_after_zeroth_reaction(newMolIndex, params, createdMolTemp, currRxn, membraneObject);
            needsResampling = moleculeOverlaps(params, simulVolume, moleculeList[newMolIndex], moleculeList,
                complexList, forwardRxns, molTemplateList, membraneObject);
        }
    }

    moleculeList[newMolIndex].myComIndex = newComIndex;
    moleculeList[newMolIndex].trajStatus = TrajStatus::propagated;
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
