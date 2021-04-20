#include "io/io.hpp"
#include "system_setup/system_setup.hpp"
#include "tracing.hpp"

void generate_coordinates(const Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, std::vector<MolTemplate>& molTemplateList,
    const std::vector<ForwardRxn>& forwardRxns, const Membrane& membraneObject)
{
    // First create all the molecules and their corresponding complexes, with random center of mass coordinates
    for (auto& oneTemp : molTemplateList) {
        if (oneTemp.isImplicitLipid == false) {
            for (unsigned itr { 0 }; itr < oneTemp.copies; ++itr) {
                moleculeList.emplace_back(initialize_molecule(Complex::numberOfComplexes, params, oneTemp, membraneObject));
                complexList.emplace_back(initialize_complex(moleculeList.back(), molTemplateList[moleculeList.back().molTypeIndex]));
                oneTemp.monomerList.emplace_back(moleculeList.back().index);
            }
        } else {
            moleculeList.emplace_back(initialize_molecule(Complex::numberOfComplexes, params, oneTemp, membraneObject));
            complexList.emplace_back(initialize_complex(moleculeList.back(), molTemplateList[moleculeList.back().molTypeIndex]));
        }
    }
    std::cout << "NUMBER OF MOLECULES IN GEN COORDS: " << moleculeList.size() << std::endl;

    if (moleculeList.size() == 0 && params.numTotalUnits == 0) {
        std::cout << "No molecules present, skipping coordinate generation.\n";
        return;
    }

    std::cout << "\nFinding and fixing overlapping proteins.\n";
    int currItr { 0 };
    while (currItr < 50) {
        bool hasOverlap { false };
        ++currItr;
        int numOverlap { 0 };

        unsigned long molListSize { moleculeList.size() };
        auto tmp = moleculeList[0];
        for (unsigned long mol1Itr { 0 }; mol1Itr < molListSize; ++mol1Itr) {
            auto& mol1 = moleculeList[mol1Itr]; // get a reference just to make the code below less messy
            for (unsigned long mol2Itr { 0 }; mol2Itr < molListSize; ++mol2Itr) {
                auto& mol2 = moleculeList[mol2Itr];
                const MolTemplate& mol2Temp { molTemplateList[mol2.molTypeIndex] };

                if ((mol1Itr != mol2Itr) && areInVicinity(mol1, mol2, molTemplateList)) {
                    for (unsigned int iface1Itr { 0 }; iface1Itr < mol1.interfaceList.size(); ++iface1Itr) {
                        auto& iface1 = mol1.interfaceList[iface1Itr];
                        for (unsigned int iface2Itr { 0 }; iface2Itr < mol2.interfaceList.size(); ++iface2Itr) {
                            auto& iface2 = mol2.interfaceList[iface2Itr];
                            int rxnIndex { 0 };
                            bool theyInteract { false };
                            //                            for (auto& rxn : forwardRxns) {
                            for (auto rxnItr : molTemplateList[moleculeList[mol1Itr].molTypeIndex]
                                                   .interfaceList[iface1Itr]
                                                   .stateList[0]
                                                   .myForwardRxns) {
                                const ForwardRxn& rxn = forwardRxns[rxnItr];
                                if ((rxn.reactantListNew[0].molTypeIndex == mol1.molTypeIndex
                                        && rxn.reactantListNew[1].molTypeIndex == mol2.molTypeIndex)
                                    || (rxn.reactantListNew[0].molTypeIndex == mol2.molTypeIndex
                                        && rxn.reactantListNew[1].molTypeIndex == mol1.molTypeIndex)) {
                                    theyInteract = true;
                                    rxnIndex = &rxn - &forwardRxns[0];
                                    break;
                                }
                            }
                            if (theyInteract) {
                                Vector tmpVec { iface1.coord - iface2.coord };
                                double mag = tmpVec.x * tmpVec.z + tmpVec.y * tmpVec.y + tmpVec.z * tmpVec.z;
                                if (mag < forwardRxns[rxnIndex].bindRadius * forwardRxns[rxnIndex].bindRadius) {
                                    mol2.create_random_coords(mol2Temp, membraneObject);
                                    complexList[mol2.myComIndex].comCoord = mol2.comCoord;
                                    hasOverlap = true;
                                    ++numOverlap;
                                    iface1Itr = mol1.interfaceList.size(); // break from interface loops
                                    iface2Itr = mol2.interfaceList.size();
                                }
                                //                                }
                            }
                        }
                    }
                }
            }
        }
        std::cout << "Iteration: " << currItr << "\nOverlapping Proteins (" << std::boolalpha << hasOverlap
                  << "): " << numOverlap << '\n';
        if (!hasOverlap) {
            std::cout << "No overlapping proteins found.\n";
            break;
        }
    }
    write_xyz("initial_crds.xyz", params, moleculeList, molTemplateList);
}
