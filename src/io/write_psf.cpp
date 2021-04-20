#include "io/io.hpp"
#include "tracing.hpp"
#include <chrono>
#include <ctime>
#include <iomanip>

void write_psf(const Parameters& params, const std::vector<Molecule>& moleculeList,
    const std::vector<MolTemplate>& molTemplateList)
{
    // TODO: Not currently working properly

    // std::cout << "Writing system PSF to system.psf...\n";
    std::ofstream outFile("system.psf");
    // Write PSF Header
    outFile << "PSF CMAP CHEQ\n\n";
    outFile << std::setw(8) << "2"
            << " !NTITLE\n";
    outFile << "REMARKS PSF for entire system\n";
    outFile << "REMARKS total molecules: " << Molecule::numberOfMolecules << '\n';
    outFile << "REMARKS total complexes: " << Complex::numberOfComplexes << '\n';
    /*Exclude Implicit lipids HERE*/

    //add up the number of sites to print
    int nAtom = 0;
    for (auto& mol : moleculeList) {
        const MolTemplate& molTemp { molTemplateList[mol.molTypeIndex] };
        if (molTemp.isImplicitLipid == false)
            nAtom += mol.interfaceList.size() + 1; //+1 is for COM.
    }
    outFile << '\n'
            << std::setw(8) << nAtom << " !NATOM\n";
    int numWritten { 1 };

    /*Add up the nuber of bonds to print here.*/
    int nBond { 0 };
    std::string color1 = "N";
    std::string color2 = "LPA";
    for (auto& mol : moleculeList) {
        const MolTemplate& molTemp { molTemplateList[mol.molTypeIndex] };
        if (molTemp.isImplicitLipid == true)
            continue;

        outFile << std::setw(8) << numWritten << ' ' << std::setw(4) << molTemp.molName.substr(0, 3) << ' '
                << std::setw(4) << mol.index << ' ' << std::setw(4) << "COM" << ' ' << std::setw(4) << "O" << ' '
                << std::setw(4) << "0" << ' ' << std::setw(10) << "0" << ' ' << std::setw(10) << mol.mass << ' '
                << std::setw(10) << "0\n";
        ++numWritten;
        nBond += molTemp.bondList.size();
        for (size_t ifaceItr { 0 }; ifaceItr < mol.interfaceList.size(); ++ifaceItr) {
            if (ifaceItr < 3) {
                outFile << std::setw(8) << numWritten << ' ' << std::setw(4) << molTemp.molName.substr(0, 3) << ' '
                        << std::setw(4) << mol.index << ' ' << std::setw(4)
                        << molTemp.interfaceList[ifaceItr].name.substr(0, 3) << ' ' << std::setw(4) << color1 << ' '
                        << std::setw(4) << "0" << ' ' << std::setw(10) << "0" << ' ' << std::setw(10) << mol.mass << ' '
                        << std::setw(10) << "0\n";

            } else if (ifaceItr < 6) {
                outFile << std::setw(8) << numWritten << ' ' << std::setw(4) << molTemp.molName.substr(0, 3) << ' '
                        << std::setw(4) << mol.index << ' ' << std::setw(4)
                        << molTemp.interfaceList[ifaceItr].name.substr(0, 3) << ' ' << std::setw(4) << color2 << ' '
                        << std::setw(4) << "0" << ' ' << std::setw(10) << "0" << ' ' << std::setw(10) << mol.mass << ' '
                        << std::setw(10) << "0\n";

            } else {

                outFile << std::setw(8) << numWritten << ' ' << std::setw(4) << molTemp.molName.substr(0, 3) << ' '
                        << std::setw(4) << mol.index << ' ' << std::setw(4)
                        << molTemp.interfaceList[ifaceItr].name.substr(0, 3) << ' ' << std::setw(4) << "O" << ' '
                        << std::setw(4) << "0" << ' ' << std::setw(10) << "0" << ' ' << std::setw(10) << mol.mass << ' '
                        << std::setw(10) << "0\n";
            }

            ++numWritten;
        }
    }
    int emptyItr { int(moleculeList.size()) };
    while (numWritten < params.numTotalUnits + 1) {
        outFile << std::setw(8) << numWritten << ' ' << std::setw(4) << "EMY" << ' ' << std::setw(4) << emptyItr << ' '
                << std::setw(4) << "EMY" << ' ' << std::setw(4) << "EMY" << ' ' << std::setw(4) << "0" << ' '
                << std::setw(10) << "0" << ' ' << std::setw(10) << 0 << ' ' << std::setw(10) << "0\n";
        ++emptyItr;
        ++numWritten;
    }
    // Write Bond Info
    numWritten = 1;
    /*Loop over all molecules in system and print bonds for their type. */
    outFile << "\n\n"
            << std::setw(8) << nBond << " !NBOND: bonds\n";
    int colNum { 0 };
    for (auto& mol : moleculeList) {

        //    for (auto& oneTemp : molTemplateList) {
        //nBond += (oneTemp.bondList.size() == 0) ? oneTemp.copies * oneTemp.interfaceList.size()
        //                            : oneTemp.copies * oneTemp.bondList.size();
        const MolTemplate& molTemp { molTemplateList[mol.molTypeIndex] };
        if (molTemp.isImplicitLipid == true) {
            numWritten -= 1; //because you skip this index of the molecule
            continue;
        }

        int i = mol.index;
        for (size_t bondItr { 0 }; bondItr < molTemp.bondList.size(); ++bondItr) {
            int tmp = numWritten + i + molTemp.bondList[bondItr][0];
            int tmp2 = numWritten + i + molTemp.bondList[bondItr][1];
            outFile << std::setw(8) << tmp << std::setw(8) << tmp2;
            ++colNum;
            if (colNum % 4 == 0)
                outFile << std::endl;
        }

        numWritten += molTemp.interfaceList.size();
    }

    //     for (auto& molTemp : molTemplateList) {
    //         if (molTemp.bondList.size() == 0) {
    //             for (size_t i { 0 }; i < molTemp.copies; ++i) {
    //                 for (size_t itr { 0 }; itr < molTemp.interfaceList.size(); ++itr) {
    //                     int tmp = numWritten + itr + 1;
    //                     outFile << std::setw(8) << numWritten << std::setw(8) << tmp;
    //                     ++colNum;
    //                     if (colNum % 4 == 0)
    //                         outFile << '\n';
    //                 }
    //                 numWritten += molTemp.interfaceList.size() + 1;
    //             }
    //         } else {
    //             for (size_t i { 0 }; i < molTemp.copies; ++i) {
    //                 for (size_t bondItr { 0 }; bondItr < molTemp.bondList.size(); ++bondItr) {
    //                     int tmp = numWritten + i + molTemp.bondList[bondItr][0];
    //                     int tmp2 = numWritten + i + molTemp.bondList[bondItr][1];
    //                     outFile << std::setw(8) << tmp << std::setw(8) << tmp2;
    //                     ++colNum;
    //                     if (colNum % 4 == 0)
    //                         outFile << '\n';
    //                 }
    //                 numWritten += molTemp.interfaceList.size();
    //             }
    //             numWritten += molTemp.copies;
    //         }
    //     }
}
