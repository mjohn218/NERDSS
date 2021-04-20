/*! \file class_SimulVolume.cpp
 * \ingroup SimulClasses
 * ### Created on 10/19/18 by Matthew Varga
 * ### Purpose
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */
#include "classes/class_SimulVolume.hpp"
#include "io/io.hpp"

#include <chrono>
#include <classes/class_SimulVolume.hpp>
#include <iostream>

/* SIMULBOX::SUBBOX */
// Member Functions
void SimulVolume::SubVolume::display()
{
    std::cout << "SubVolume " << absIndex << '\n';
    std::cout << "\tRel. Indices: [" << xIndex << ", " << yIndex << ", " << zIndex << "]\n";
    std::cout << "\tMolecule Members:";
    for (auto& mol : memberMolList)
        std::cout << ' ' << mol;
    std::cout << "\n\tNeighbors (abs. index):";
    for (auto& cell : neighborList)
        std::cout << ' ' << cell;
    std::cout << std::endl;
}

/* SIMULBOX::DIMENSIONS */
// Constructors
SimulVolume::Dimensions::Dimensions(const Parameters& params, const Membrane& membraneObject)
{
    //    double cellLength { params.rMaxLimit * 1.5 }; // give some buffer to the cell sizes
    double cellLength { params.rMaxLimit };
    x = int(floor(membraneObject.waterBox.x / cellLength));
    y = int(floor(membraneObject.waterBox.y / cellLength));

    int scale { 1 }; //(params.rMaxLimit / ((params.rMaxLimit - params.rMaxRadius)) > 2) ? 1 : 4 };
    z = std::max(4, int(floor(membraneObject.waterBox.z / cellLength)) / scale);

    tot = x * y * z;
}

// Member Functions
void SimulVolume::Dimensions::check_dimensions(const Parameters& params, const Membrane& membraneObject)
{
    // now check to make sure none of the dimensions are too small
    // if they are, change the scale
    //    double cellLength { params.rMaxLimit * 1.2 };
    double cellLength { params.rMaxLimit };
    if (membraneObject.waterBox.x / (x * 1.0) < cellLength)
        x = std::max(4, int(floor(membraneObject.waterBox.x / cellLength)) / 2);
    if (membraneObject.waterBox.y / (y * 1.0) < cellLength)
        y = std::max(4, int(floor(membraneObject.waterBox.y / cellLength)) / 2);
    if (membraneObject.waterBox.z / (z * 1.0) < cellLength) {
        if (membraneObject.waterBox.z > 0)
            z = std::max(4, int(floor(membraneObject.waterBox.z / cellLength)) / 2);
        else
            z = 1;
    }
    tot = x * y * z;

    // Now make sure There aren't too many pairs of cells.
    //At the same time, if the numberOfMolecules changes throughout the simulation,
    //do not want the subvolumes to be too large.
    //based on total molecules
    int totMol = Molecule::numberOfMolecules;
    double maxPairsMols { 0.5 * totMol * totMol };
    double setLowerMax = 4000; //allow this many subvolumes, even if it is larger than maxPairsMols.
    double maxPairs = std::max(setLowerMax, maxPairsMols);
    int minCells { 64 };
    double scale { 2.0 };

    // TODO: this if statement is temporary. noneq reactions starting with few reactants will always use minimum
    // number of cells, which can get VERY slow
    if (tot > maxPairs) { //&& !params.isNonEQ) {
        while (tot > maxPairs && tot > minCells) {
            std::cout << "CELL PAIR MAX EXCEEDED\n"
                      << "\tCurrent number of cells: " << tot << "\n\tMax number of cells: " << maxPairs << '\n';
            std::cout << "Scaling down number of cells.\n";

            x = std::max(2, int(floor(membraneObject.waterBox.x / cellLength) / scale));
            y = std::max(2, int(floor(membraneObject.waterBox.y / cellLength) / scale));
            z = std::max(2, int(floor(membraneObject.waterBox.z / cellLength) / (2 * scale)));
            tot = x * y * z;
            scale += 1;
        }
    }

    tot = x * y * z;
}

/* SIMULBOX */
// Member Functions

void SimulVolume::display()
{
    std::cout << "Simulation volume parameters:\n";
    std::cout << "Total sub-volumes: " << numSubCells.tot << '\n';
    std::cout << "\tDimensions: [" << numSubCells.x << ", " << numSubCells.y << ", " << numSubCells.z << "]\n";
    std::cout << "\tMaximum sub-volume neighbors: " << maxNeighbors << '\n';
    std::cout << "\tSub-volume size: [" << subCellSize.x << ", " << subCellSize.y << ", " << subCellSize.z << "]\n";
}

void SimulVolume::create_simulation_volume(const Parameters& params, const Membrane& membraneObject)
{
    // Determine the number of boxes there will be in each dimension
    numSubCells = Dimensions(params, membraneObject);
    numSubCells.check_dimensions(params, membraneObject);

    // Calculate the cells' dimensions in nanometers
    if (membraneObject.waterBox.z > 0)
        subCellSize = Coord { membraneObject.waterBox.x / (numSubCells.x * 1.0), membraneObject.waterBox.y / (numSubCells.y * 1.0),
            membraneObject.waterBox.z / (numSubCells.z * 1.0) };
    else
        subCellSize = Coord { membraneObject.waterBox.x / (numSubCells.x * 1.0), membraneObject.waterBox.y / (numSubCells.y * 1.0),
            1 };
    // Create cell neighborlists.
    subCellList = std::vector<SubVolume>(numSubCells.tot);
    create_cell_neighbor_list_cubic();
}

void SimulVolume::create_cell_neighbor_list_cubic()
{
    int cellNum { 0 };
    for (unsigned zItr { 0 }; zItr < numSubCells.z; ++zItr) {
        for (unsigned yItr { 0 }; yItr < numSubCells.y; ++yItr) {
            for (unsigned xItr { 0 }; xItr < numSubCells.x; ++xItr) {
                // set up SubVolume
                subCellList[cellNum].absIndex = cellNum;
                subCellList[cellNum].xIndex = xItr;
                subCellList[cellNum].yIndex = yItr;
                subCellList[cellNum].zIndex = zItr;

                /*For each cell figure out its 13 neighbors that are ~forward and up*/
                // This only works for cubic subvolumes
                if (xItr < (numSubCells.x - 1)) {
                    // count all the plus x boxes
                    subCellList[cellNum].neighborList.push_back(
                        (xItr + 1) + yItr * numSubCells.x + zItr * (numSubCells.x * numSubCells.y));
                    if (yItr < (numSubCells.y - 1)) {
                        subCellList[cellNum].neighborList.push_back(
                            (xItr + 1) + (yItr + 1) * numSubCells.x + zItr * (numSubCells.x * numSubCells.y));
                        if (zItr < (numSubCells.z - 1)) {
                            subCellList[cellNum].neighborList.push_back(
                                (xItr + 1) + (yItr + 1) * numSubCells.x + (zItr + 1) * (numSubCells.x * numSubCells.y));
                        }
                    }
                    if (zItr < (numSubCells.z - 1)) {
                        subCellList[cellNum].neighborList.push_back(
                            (xItr + 1) + yItr * numSubCells.x + (zItr + 1) * (numSubCells.x * numSubCells.y));
                    }
                }
                if (yItr < (numSubCells.y - 1)) {
                    subCellList[cellNum].neighborList.push_back(
                        xItr + (yItr + 1) * numSubCells.x + zItr * (numSubCells.x * numSubCells.y));
                    if (xItr > 0) {
                        subCellList[cellNum].neighborList.push_back(
                            (xItr - 1) + (yItr + 1) * numSubCells.x + zItr * (numSubCells.x * numSubCells.y));
                        if (zItr < (numSubCells.z - 1)) {
                            subCellList[cellNum].neighborList.push_back(
                                (xItr - 1) + (yItr + 1) * numSubCells.x + (zItr + 1) * (numSubCells.x * numSubCells.y));
                        }
                    }
                    if (zItr < (numSubCells.z - 1)) {
                        subCellList[cellNum].neighborList.push_back(
                            xItr + (yItr + 1) * numSubCells.x + (zItr + 1) * (numSubCells.x * numSubCells.y));
                    }
                }
                if (zItr < (numSubCells.z - 1)) {
                    subCellList[cellNum].neighborList.push_back(
                        xItr + yItr * numSubCells.x + (zItr + 1) * (numSubCells.x * numSubCells.y));
                    if (xItr > 0) {
                        subCellList[cellNum].neighborList.push_back(
                            xItr - 1 + yItr * numSubCells.x + (zItr + 1) * (numSubCells.x * numSubCells.y));
                        if (yItr > 0) {
                            subCellList[cellNum].neighborList.push_back(
                                xItr - 1 + (yItr - 1) * numSubCells.x + (zItr + 1) * (numSubCells.x * numSubCells.y));
                        }
                    }
                    if (yItr > 0) {
                        subCellList[cellNum].neighborList.push_back(
                            xItr + (yItr - 1) * numSubCells.x + (zItr + 1) * (numSubCells.x * numSubCells.y));
                        if (xItr < (numSubCells.x - 1)) {
                            subCellList[cellNum].neighborList.push_back(
                                xItr + 1 + (yItr - 1) * numSubCells.x + (zItr + 1) * (numSubCells.x * numSubCells.y));
                        }
                    }
                }
                if (subCellList[cellNum].neighborList.size() > maxNeighbors) {
                    std::cerr << "ERROR: Maximum number of neighbors exceeded for SubVolume " << cellNum
                              << ". Exiting\n";
                    exit(1);
                }
                ++cellNum;
            } // end looping over z cells
        } // end looping over y cells
    } // end looping over x cells
}

void SimulVolume::update_memberMolLists(const Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject, int simItr)
{
    // make sure the list of member molecules is empty
    for (auto& subBox : subCellList)
        subBox.memberMolList.clear();

    int itrCheck = 1000; //no need to check every step if it violates box boundaries.

    int itr { 0 };
    if (simItr % itrCheck != 0) {
        /*just assign bins, don't check bin limits/errors*/
        for (unsigned molItr { 0 }; molItr < moleculeList.size(); ++molItr) {
            Molecule& mol = moleculeList[molItr]; // just for legibility

            if (mol.isEmpty || mol.isImplicitLipid)
                continue;

            // get which box the Molecule belongs to
            int xItr { int((mol.comCoord.x + membraneObject.waterBox.x / 2) / subCellSize.x) };
            int yItr { int((mol.comCoord.y + membraneObject.waterBox.y / 2) / subCellSize.y) };
            int zItr;
            if (membraneObject.waterBox.z > 0)
                zItr = int(-(mol.comCoord.z + 1E-6 - membraneObject.waterBox.z / 2.0) / subCellSize.z);
            else
                zItr = 0;

            // allow the modecule a bit out of the box
            if (xItr == -1)
                xItr = 0;
            if (xItr == numSubCells.x)
                xItr = numSubCells.x - 1;
            if (yItr == -1)
                yItr = 0;
            if (yItr == numSubCells.y)
                yItr = numSubCells.y - 1;
            if (zItr == -1)
                zItr = 0;
            if (zItr == numSubCells.z)
                zItr = numSubCells.z - 1;

            int currBin = xItr + (yItr * numSubCells.x) + (zItr * numSubCells.x * numSubCells.y);

            mol.mySubVolIndex = currBin;
            if (currBin >= numSubCells.tot) {
                std::cerr << "Molecule " << mol.index
                          << " seems outside simulation volume, with center of mass coordinates ["
                          << mol.comCoord << "].\n";
                exit(1);
            }
            subCellList[currBin].memberMolList.push_back(mol.index);
        }
    } else {
        /*make sure proteins are within bin limits, lipids are on membrane*/
        for (unsigned molItr { 0 }; molItr < moleculeList.size(); ++molItr) {
            Molecule& mol = moleculeList[molItr]; // just for legibility

            if (mol.isEmpty || mol.isImplicitLipid)
                continue;

            // get which box the Molecule belongs to
            int xItr { int((mol.comCoord.x + membraneObject.waterBox.x / 2) / subCellSize.x) };
            int yItr { int((mol.comCoord.y + membraneObject.waterBox.y / 2) / subCellSize.y) };
            int zItr;
            if (membraneObject.waterBox.z > 0)
                zItr = int(-(mol.comCoord.z + 1E-6 - membraneObject.waterBox.z / 2.0) / subCellSize.z);
            else
                zItr = 0;

            if (xItr == -1)
                xItr = 0;
            if (xItr == numSubCells.x)
                xItr = numSubCells.x - 1;
            if (yItr == -1)
                yItr = 0;
            if (yItr == numSubCells.y)
                yItr = numSubCells.y - 1;
            if (zItr == -1)
                zItr = 0;
            if (zItr == numSubCells.z)
                zItr = numSubCells.z - 1;
            int currBin = xItr + (yItr * numSubCells.x) + (zItr * numSubCells.x * numSubCells.y);

            // Make sure the Molecule is still on the membrane if its supposed to be
            if (std::abs(molTemplateList[mol.molTypeIndex].D.z) < 1E-10) {
                //define RS3Dinput
                double RS3Dinput { 0.0 };

                for (int RS3Dindex = 0; RS3Dindex < 100; RS3Dindex++) {
                    if (std::abs(membraneObject.RS3Dvect[RS3Dindex + 400] - mol.molTypeIndex) < 1E-2) {
                        RS3Dinput = membraneObject.RS3Dvect[RS3Dindex + 300];
                        //   std::cout << mol.molTypeIndex << "\t" << membraneObject.RS3Dvect[RS3Dindex + 400] << "\t" << RS3Dindex << "\n";
                        break;
                    }
                }

                if (mol.comCoord.z - 0.1 > -membraneObject.waterBox.z * 0.5 + RS3Dinput && mol.isImplicitLipid == false) {
                    //            && std::abs(mol.comCoord.z) - std::abs((membraneObject.waterBox.z / 2)) > 1E-6)
                    std::cerr << "Molecule " << mol.index << " of type " << molTemplateList[mol.molTypeIndex].molName
                              << " is off the membrane. Writing coordinates and exiting.\n";
                    //  std::cout << mol.molTypeIndex << "\t" << -membraneObject.waterBox.z * 0.5 + RS3Dinput << "\t" << RS3Dinput << "\n";
                    write_xyz(std::string { "error_coord_dump.xyz" }, params, moleculeList, molTemplateList);
                    exit(1);
                }
            }

            // Now make sure the Molecule is still inside the box in all dimensions
            if (mol.comCoord.z > (membraneObject.waterBox.z / 2) || mol.comCoord.z + 1E-6 < -(membraneObject.waterBox.z / 2)) {
                std::cout << "Molecule " << mol.index
                          << " is outside simulation volume in the z-dimension, with center of mass coordinates ["
                          << mol.comCoord << "]. Attempting to fit back into box.\n";
                complexList[mol.myComIndex].put_back_into_SimulVolume(itr, mol, membraneObject, moleculeList, molTemplateList);
                // reset member search
                molItr = 0;
                for (auto& subBox : subCellList)
                    subBox.memberMolList.clear();
            } else if (mol.comCoord.y > (membraneObject.waterBox.y / 2) || mol.comCoord.y + 1E-6 < -(membraneObject.waterBox.y / 2)) {
                std::cout << "Molecule " << mol.index
                          << " is outside simulation volume in the y-dimension, with center of mass coordinates ["
                          << mol.comCoord << "]. Attempting to fit back into box.\n";
                complexList[mol.myComIndex].put_back_into_SimulVolume(itr, mol, membraneObject, moleculeList, molTemplateList);
                // reset member search
                molItr = 0;
                for (auto& subBox : subCellList)
                    subBox.memberMolList.clear();
            } else if (mol.comCoord.x > (membraneObject.waterBox.x / 2) || mol.comCoord.x + 1E-6 < -(membraneObject.waterBox.x / 2)) {
                std::cout << "Molecule " << mol.index
                          << " is outside simulation volume in the x-dimension, with center of mass coordinates ["
                          << mol.comCoord << "]. Attempting to fit back into box.\n";
                complexList[mol.myComIndex].put_back_into_SimulVolume(itr, mol, membraneObject, moleculeList, molTemplateList);
                // reset member search
                molItr = 0;
                for (auto& subBox : subCellList)
                    subBox.memberMolList.clear();
            } else if (currBin > (numSubCells.tot) || currBin < 0) {
                std::cout << "Molecule " << mol.index << " is outside simulation volume with center of mass coordinates ["
                          << mol.comCoord << "]. Attempting to fit back into box.\n";
                complexList[mol.myComIndex].put_back_into_SimulVolume(itr, mol, membraneObject, moleculeList, molTemplateList);
                // reset member search
                molItr = 0;
                for (auto& subBox : subCellList)
                    subBox.memberMolList.clear();
            } else {
                // The Molecule is in the simulation volume, okay to proceed
                mol.mySubVolIndex = currBin;
                subCellList[currBin].memberMolList.push_back(mol.index);
            }
        } //loop over all molecules.
    } //check all boundary limits are OK.
}
