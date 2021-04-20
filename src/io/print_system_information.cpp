/*! \file print_system_information.cpp
 *
 * \brief
 *
 * ### Created on 2019-06-12 by Matthew Varga
 */

#include "io/io.hpp"
#include "tracing.hpp"

void print_system_information(long long int simItr, std::ofstream& systemFile, const std::vector<Molecule>& moleculeList,
    const std::vector<Complex>& complexList, const std::vector<MolTemplate>& molTemplateList)
{
    // TRACE();
    systemFile << "Iteration: " << simItr << '\n';
    systemFile << llinebreak << llinebreak << "\t\t\tMOLECULES\n"
               << llinebreak << llinebreak << linebreak;
    for (const auto& mol : moleculeList) {
        systemFile << "Index: " << mol.index << '\n';
        systemFile << "Is empty: " << std::boolalpha << mol.isEmpty << '\n';
        if (!mol.isEmpty) {
            systemFile << "Type: " << molTemplateList[mol.molTypeIndex].molName << '\n';
            systemFile << "Parent complex index: " << mol.myComIndex << '\n';
            systemFile << "Sub volume index: " << mol.mySubVolIndex << '\n';
            systemFile << "Is a lipid: " << std::boolalpha << mol.isLipid << '\n';
            systemFile << "Center of mass coordinate: " << mol.comCoord << '\n';
            systemFile << "Interfaces:\n";
            for (const auto& iface : mol.interfaceList) {
                systemFile << "\tRelative index: " << iface.relIndex << '\n';
                systemFile << "\tAbsolute index: " << iface.index << '\n';
                systemFile << "\tInterface name: "
                           << molTemplateList[mol.molTypeIndex].interfaceList[iface.relIndex].name << '\n';
                systemFile << "\tCoordinate: " << iface.coord << '\n';
                if (iface.stateIden != '\0')
                    systemFile << "\tCurrent state: " << iface.stateIden << '\n';
                if (iface.isBound) {
                    systemFile << "\tInteraction:\n";
                    systemFile << "\t\tPartner index: " << iface.interaction.partnerIndex << '\n';
                    systemFile << "\t\tPartner interface index " << iface.interaction.partnerIfaceIndex << '\n';
                }
            }
        }
        systemFile << linebreak;
    }

    systemFile << llinebreak << llinebreak << "\t\t\tCOMPLEXES\n"
               << llinebreak << llinebreak << linebreak;
    for (const auto& complex : complexList) {
        systemFile << "Index: " << complex.index << '\n';
        systemFile << "Is empty: " << complex.isEmpty << '\n';
        if (!complex.isEmpty) {
            systemFile << "Mass: " << complex.mass << '\n';
            systemFile << "Radius: " << complex.radius << '\n';
            systemFile << "Center of mass coordinate: " << complex.comCoord << '\n';
            systemFile << "Translational diffusion constants:\n";
            systemFile << "\tDx = " << complex.D.x << '\n';
            systemFile << "\tDy = " << complex.D.y << '\n';
            systemFile << "\tDz = " << complex.D.z << '\n';
            systemFile << "Rotational diffusion constants:\n";
            systemFile << "\tDrx = " << complex.Dr.x << '\n';
            systemFile << "\tDry = " << complex.Dr.y << '\n';
            systemFile << "\tDrz = " << complex.Dr.z << '\n';
            systemFile << "Member molecules:";
            for (auto memMol : complex.memberList)
                systemFile << ' ' << memMol;
            systemFile << '\n';
            systemFile << "Number of each molecule type:";
            for (unsigned i { 0 }; i < complex.numEachMol.size(); ++i)
                systemFile << '\t' << molTemplateList[i].molName << ": " << complex.numEachMol[i] << '\n';
        }
        systemFile << linebreak;
    }
}
