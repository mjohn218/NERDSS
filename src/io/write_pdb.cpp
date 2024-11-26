#include "debug/debug.hpp"
#include "io/io.hpp"
#include "mpi/mpi_function.hpp"
#include "tracing.hpp"
#include <chrono>
#include <ctime>
#include <iomanip>

void write_pdb(long long int simItr, unsigned frameNum, const Parameters& params, const std::vector<Molecule>& moleculeList,
    const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject)
{
    // TRACE();
    long long int totFrames { params.nItr / params.trajWrite };
    //    std::ofstream pdbFile { "pdb/" + std::to_string(frameNum) + ".pdb" };
    // if (!pdbFile) {
    std::ofstream pdbFile { "PDB/" + std::to_string(frameNum) + ".pdb" };
    //}
    auto printTime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    pdbFile << std::left << std::setw(6) << "TITLE" << ' ' << std::left << std::setw(70) << "PDB TIMESTEP " << simItr
            << " CREATED " << std::ctime(&printTime);
    pdbFile << std::left << std::setw(6) << "CRYST1  " << std::setw(9) << membraneObject.waterBox.x << std::setw(9)
            << membraneObject.waterBox.y << std::setw(9) << membraneObject.waterBox.z << std::setw(7) << 90 << std::setw(7) << 90
            << std::setw(7) << 90 << ' ' << 'P' << std::setw(4) << 1 << std::endl;

    int i { 0 };
    // for (const auto& oneTemp : molTemplateList) {
    //     // this is just so visualization software reads species in at the same color whether or not they're in the
    //     // system
    //     if (oneTemp.isImplicitLipid)
    //         continue;
    //     pdbFile << "ATOM  " << std::right << std::setw(5) << i << " " << std::left << std::setw(4) << "COM"
    //             << " " << std::right << std::setw(3)
    //             << molTemplateList[i].molName.substr(0, 3) << "  " << std::right << std::setw(4) << i << "    "
    //             << std::right << std::setw(8) << membraneObject.waterBox.x << std::right << std::setw(8) << membraneObject.waterBox.y << std::right << std::setw(8)
    //             << membraneObject.waterBox.z << std::setw(6) << 0.00 << std::setw(6) << 0.00 << std::left << std::setw(2)
    //             << "CL" << std::endl;
    //     ++i;
    // }
    int molCounter { 0 };
    for (const auto& mol : moleculeList) {
        if (mol.isImplicitLipid)
            continue;

        if (!mol.isEmpty) {
            const MolTemplate& oneTemp = molTemplateList[mol.molTypeIndex];
            pdbFile << "ATOM  " << std::right << std::setw(5) << i % 100000 << " " << std::left << std::setw(4) << "COM"
                    << " " << std::right
                    << std::setw(3) << oneTemp.molName.substr(0, 3) << "  " << std::right << std::setw(4) << molCounter % 10000 << "    "
                    << std::right << std::setw(8) << std::fixed << std::setprecision(1) << (mol.comCoord.x + membraneObject.waterBox.x / 2) << std::right << std::setw(8)
                    << (mol.comCoord.y + membraneObject.waterBox.y / 2) << std::right << std::setw(8)
                    << (mol.comCoord.z + membraneObject.waterBox.z / 2);
            pdbFile.unsetf(std::ios_base::fixed);
            pdbFile << std::setw(6) << 0.00 << std::setw(6) << 0.00
                    << std::left << std::setw(2) << "CL" << std::endl;
            ++i;

            for (unsigned j { 0 }; j < mol.interfaceList.size(); ++j) {
                pdbFile << "ATOM  " << std::right << std::setw(5) << i % 100000 << " " << std::left << std::setw(4) // << ' '
                        << oneTemp.interfaceList[j].name.substr(0, 4) << " " << std::right << std::setw(3)
                        << oneTemp.molName.substr(0, 3) << "  " << std::right << std::setw(4) << molCounter % 10000 << "    "
                        << std::right << std::setw(8) << std::fixed << std::setprecision(1) << (mol.interfaceList[j].coord.x + membraneObject.waterBox.x / 2) << std::right << std::setw(8)
                        << (mol.interfaceList[j].coord.y + membraneObject.waterBox.y / 2) << std::right << std::setw(8)
                        << (mol.interfaceList[j].coord.z + membraneObject.waterBox.z / 2);
                pdbFile.unsetf(std::ios_base::fixed);
                pdbFile << std::setw(6) << 0.00
                        << std::setw(6) << 0.00 << std::left << std::setw(2) << "CL" << std::endl;
                ++i;
            }
            ++molCounter;
        }
    }
}
