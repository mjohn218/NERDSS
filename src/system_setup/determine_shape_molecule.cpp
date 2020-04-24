#include "classes/class_Coord.hpp"
#include "system_setup/system_setup.hpp"

void determine_shape_molecule(std::vector<MolTemplate>& molTemplateList)
{
    for (auto& molTemplate : molTemplateList) {
        if (molTemplate.isLipid == true && std::abs(molTemplate.D.z - 0) > 1e-12) {
            std::cerr << "WARNING: Dz of Lipid Must Be Zero. Dz of " << molTemplate.molName << " has been set to zero." << std::endl;
            molTemplate.D.z = 0;
        }
        if (molTemplate.isImplicitLipid == true && std::abs(molTemplate.D.z - 0) > 1e-12) {
            std::cerr << "WARNING: Dz of Implicit Lipid Must Be Zero. Dz of " << molTemplate.molName << " has been set to zero." << std::endl;
            molTemplate.D.z = 0;
        }

        molTemplate.isPoint = true;
        molTemplate.isRod = true;
        // if all interfaces are on the COM, this is a point
        for (auto& molInterface : molTemplate.interfaceList) {
            if (molInterface.iCoord != molTemplate.comCoord) {
                molTemplate.isPoint = false;
            }
        }
        if (molTemplate.isPoint == true) {
            molTemplate.isRod = false;
            continue; // isRod = false
        } else {
            // not point, determine isRod
            // if only one interface, this is rod when not point
            if (molTemplate.interfaceList.size() == 1) {
                molTemplate.isRod = true;
                continue;
            } else {
                //more than one interface, each pair must through the COM if it is a Rod
                for (int i = 0; i < molTemplate.interfaceList.size() - 1; i++) {
                    for (int j = i + 1; j < molTemplate.interfaceList.size(); j++) {
                        molTemplate.isRod = is_co_linear(molTemplate.comCoord, molTemplate.interfaceList[i].iCoord, molTemplate.interfaceList[j].iCoord);
                        if (molTemplate.isRod == false) {
                            break;
                        }
                    }
                    if (molTemplate.isRod == false) {
                        break;
                    }
                }
            }
        }
    }
}
