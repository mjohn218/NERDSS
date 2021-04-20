#include "parser/parser_functions.hpp"

void read_bonds(int numBonds, std::ifstream& molFile, MolTemplate& molTemplate)
{
    for (unsigned bondItr { 0 }; bondItr < numBonds; ++bondItr) {
        std::string atom1 {};
        std::string atom2 {};
        molFile >> atom1 >> atom2;
        molFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::transform(atom1.begin(), atom1.end(), atom1.begin(), ::towlower);
        std::transform(atom2.begin(), atom2.end(), atom2.begin(), ::towlower);

        // TODO: Replace the below
        std::array<int, 2> tmpBond { { -1, -1 } };
        if (atom1 == std::string { "com" })
            tmpBond[0] = 0;
        if (atom2 == std::string { "com" })
            tmpBond[1] = 0;

        for (unsigned ifaceItr { 0 }; ifaceItr < molTemplate.interfaceList.size(); ++ifaceItr) {
            if (molTemplate.interfaceList[ifaceItr].name == atom1)
                tmpBond[0] = ifaceItr + 1;
            if (molTemplate.interfaceList[ifaceItr].name == atom2)
                tmpBond[1] = ifaceItr + 1;

            if (tmpBond[0] != -1 && tmpBond[1] != -1)
                break;
        }

        if (tmpBond[0] == -1 || tmpBond[1] == -1) {
            std::cout << "Invalid atom in " << molTemplate.molName << " mol file. Ignoring bonds...\n";
            molTemplate.bondList.clear();
            return;
        } else {
            std::sort(tmpBond.begin(), tmpBond.end());
            molTemplate.bondList.emplace_back(tmpBond);
            std::cout << atom1 << "-" << atom2 << ", ";
        }
    }
    std::cout << std::endl;
}
