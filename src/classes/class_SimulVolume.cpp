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

#include <chrono>
#include <classes/class_SimulVolume.hpp>
#include <iostream>

#include "error/error.hpp"
#include "io/io.hpp"

/* SIMULBOX::SUBBOX */
// Member Functions
void SimulVolume::SubVolume::display() {
  std::cout << "SubVolume " << absIndex << '\n';
  std::cout << "\tRel. Indices: [" << xIndex << ", " << yIndex << ", " << zIndex
            << "]\n";
  std::cout << "\tMolecule Members:";
  for (auto& mol : memberMolList) std::cout << ' ' << mol;
  std::cout << "\n\tNeighbors (abs. index):";
  for (auto& cell : neighborList) std::cout << ' ' << cell;
  std::cout << std::endl;
}

/* SIMULBOX::DIMENSIONS */
// Constructors
SimulVolume::Dimensions::Dimensions(const Parameters& params,
                                    const Membrane& membraneObject) {
  //    double cellLength { params.rMaxLimit * 1.5 }; // give some buffer to the
  //    cell sizes
  double cellLength{params.rMaxLimit};
  x = int(floor(membraneObject.waterBox.x / cellLength));
  y = int(floor(membraneObject.waterBox.y / cellLength));

  int scale{1};  //(params.rMaxLimit / ((params.rMaxLimit - params.rMaxRadius))
                 //> 2) ? 1 : 4 };
  z = std::max(4, int(floor(membraneObject.waterBox.z / cellLength)) / scale);

  tot = x * y * z;
}

// Member Functions
void SimulVolume::Dimensions::check_dimensions(const Parameters& params,
                                               const Membrane& membraneObject) {
  // now check to make sure none of the dimensions are too small
  // if they are, change the scale
  //    double cellLength { params.rMaxLimit * 1.2 };
  double cellLength{params.rMaxLimit};
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

  /*
  HARD CODING TO LIMIT THE NUMBER OF SUBVOLUMES, TOO MANY CAN BE VERY SLOW!
  FOR THE MPI BRANCH, WE WANT THE XBINS AS LARGE AS POSSIBLE AND THE YBINS
  AND ZBINS AS SMALL AS POSSIBLE. THEN THE COMMUNICATION VOLUME WILL BE SMALL
  AND EFFICIENT. WE WILL KEEP THE TOTAL BINS NO MORE THAN 27000
  */
  if (x * y * z > 27000) {
    y = static_cast<int>(std::sqrt(27000 / x));
    z = y;
  }

  tot = x * y * z;

  // Now make sure There aren't too many pairs of cells.
  // At the same time, if the numberOfMolecules changes throughout the
  // simulation, do not want the subvolumes to be too large. based on total
  // molecules
  int totMol = Molecule::numberOfMolecules;
  double maxPairsMols{0.5 * totMol * totMol};
  double setLowerMax = 2000000;  // allow this many subvolumes, even if it is
                              // larger than maxPairsMols.
  double maxPairs = std::max(setLowerMax, maxPairsMols);
  int minCells{64};
  double scale{2.0};

  // TODO: this if statement is temporary. noneq reactions starting with few
  // reactants will always use minimum number of cells, which can get VERY slow
  if (tot > maxPairs) {  //&& !params.isNonEQ) {
    while (tot > maxPairs && tot > minCells) {
      std::cout << "CELL PAIR MAX EXCEEDED\n"
                << "\tCurrent number of cells: " << tot
                << "\n\tMax number of cells: " << maxPairs << '\n';
      std::cout << "Scaling down number of cells.\n";

      x = std::max(2,
                   int(floor(membraneObject.waterBox.x / cellLength) / scale));
      y = std::max(2,
                   int(floor(membraneObject.waterBox.y / cellLength) / scale));
      z = std::max(
          2, int(floor(membraneObject.waterBox.z / cellLength) / (2 * scale)));
      tot = x * y * z;
      scale += 1;
    }
  }

  tot = x * y * z;
}

/* SIMULBOX */
// Member Functions

void SimulVolume::display() {
  std::cout << "Simulation volume parameters:\n";
  std::cout << "Total sub-volumes: " << numSubCells.tot << '\n';
  std::cout << "\tDimensions: [" << numSubCells.x << ", " << numSubCells.y
            << ", " << numSubCells.z << "]\n";
  std::cout << "\tMaximum sub-volume neighbors: " << maxNeighbors << '\n';
  std::cout << "\tSub-volume size: [" << subCellSize.x << ", " << subCellSize.y
            << ", " << subCellSize.z << "]\n";
}

void SimulVolume::create_simulation_volume(const Parameters& params,
                                           const Membrane& membraneObject) {
  // Determine the number of boxes there will be in each dimension
  numSubCells = Dimensions(params, membraneObject);
  numSubCells.check_dimensions(params, membraneObject);

  // Calculate the cells' dimensions in nanometers
  if (membraneObject.waterBox.z > 0)
    subCellSize = Coord{membraneObject.waterBox.x / (numSubCells.x * 1.0),
                        membraneObject.waterBox.y / (numSubCells.y * 1.0),
                        membraneObject.waterBox.z / (numSubCells.z * 1.0)};
  else
    subCellSize = Coord{membraneObject.waterBox.x / (numSubCells.x * 1.0),
                        membraneObject.waterBox.y / (numSubCells.y * 1.0), 1};
  // Create cell neighborlists.
  subCellList = std::vector<SubVolume>(numSubCells.tot);
  create_cell_neighbor_list_cubic();
}

void SimulVolume::create_cell_neighbor_list_cubic() {
  int cellNum{0};
  for (unsigned zItr{0}; zItr < numSubCells.z; ++zItr) {
    for (unsigned yItr{0}; yItr < numSubCells.y; ++yItr) {
      for (unsigned xItr{0}; xItr < numSubCells.x; ++xItr) {
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
              (xItr + 1) + yItr * numSubCells.x +
              zItr * (numSubCells.x * numSubCells.y));
          if (yItr < (numSubCells.y - 1)) {
            subCellList[cellNum].neighborList.push_back(
                (xItr + 1) + (yItr + 1) * numSubCells.x +
                zItr * (numSubCells.x * numSubCells.y));
            if (zItr < (numSubCells.z - 1)) {
              subCellList[cellNum].neighborList.push_back(
                  (xItr + 1) + (yItr + 1) * numSubCells.x +
                  (zItr + 1) * (numSubCells.x * numSubCells.y));
            }
          }
          if (zItr < (numSubCells.z - 1)) {
            subCellList[cellNum].neighborList.push_back(
                (xItr + 1) + yItr * numSubCells.x +
                (zItr + 1) * (numSubCells.x * numSubCells.y));
          }
        }
        if (yItr < (numSubCells.y - 1)) {
          subCellList[cellNum].neighborList.push_back(
              xItr + (yItr + 1) * numSubCells.x +
              zItr * (numSubCells.x * numSubCells.y));
          if (xItr > 0) {
            subCellList[cellNum].neighborList.push_back(
                (xItr - 1) + (yItr + 1) * numSubCells.x +
                zItr * (numSubCells.x * numSubCells.y));
            if (zItr < (numSubCells.z - 1)) {
              subCellList[cellNum].neighborList.push_back(
                  (xItr - 1) + (yItr + 1) * numSubCells.x +
                  (zItr + 1) * (numSubCells.x * numSubCells.y));
            }
          }
          if (zItr < (numSubCells.z - 1)) {
            subCellList[cellNum].neighborList.push_back(
                xItr + (yItr + 1) * numSubCells.x +
                (zItr + 1) * (numSubCells.x * numSubCells.y));
          }
        }
        if (zItr < (numSubCells.z - 1)) {
          subCellList[cellNum].neighborList.push_back(
              xItr + yItr * numSubCells.x +
              (zItr + 1) * (numSubCells.x * numSubCells.y));
          if (xItr > 0) {
            subCellList[cellNum].neighborList.push_back(
                xItr - 1 + yItr * numSubCells.x +
                (zItr + 1) * (numSubCells.x * numSubCells.y));
            if (yItr > 0) {
              subCellList[cellNum].neighborList.push_back(
                  xItr - 1 + (yItr - 1) * numSubCells.x +
                  (zItr + 1) * (numSubCells.x * numSubCells.y));
            }
          }
          if (yItr > 0) {
            subCellList[cellNum].neighborList.push_back(
                xItr + (yItr - 1) * numSubCells.x +
                (zItr + 1) * (numSubCells.x * numSubCells.y));
            if (xItr < (numSubCells.x - 1)) {
              subCellList[cellNum].neighborList.push_back(
                  xItr + 1 + (yItr - 1) * numSubCells.x +
                  (zItr + 1) * (numSubCells.x * numSubCells.y));
            }
          }
        }
        if (subCellList[cellNum].neighborList.size() > maxNeighbors) {
          std::cerr
              << "ERROR: Maximum number of neighbors exceeded for SubVolume "
              << cellNum << ". Exiting\n";
          exit(1);
        }
        ++cellNum;
      }  // end looping over z cells
    }    // end looping over y cells
  }      // end looping over x cells
}

void SimulVolume::update_memberMolLists(
    const Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject,
    int simItr,
    MpiContext& mpiContext)  // xItr bin should be reduced by xOffset in order
                             // to calculate the box id in the current rank
{
  // make sure the list of member molecules is empty
  for (auto& subBox : subCellList) subBox.memberMolList.clear();

  int itr{0};

  for (unsigned molItr{0}; molItr < moleculeList.size(); ++molItr) {
    Molecule& mol = moleculeList[molItr];  // just for legibility
    mol.isLeftGhost = false;
    mol.isLeftEdge = false;
    mol.isRightGhost = false;
    mol.isRightEdge = false;

    if (mol.isEmpty) continue;

    Complex& com = complexList[mol.myComIndex];
    com.isLeftGhost = false;
    com.isLeftEdge = false;
    com.isRightGhost = false;
    com.isRightEdge = false;

    if (mol.isImplicitLipid) continue;

    // Get which box the Molecule belongs to.
    // xOffset represents starting x coordinate of cell at current rank;
    // Namely, xItr was not customized for each rank, but left as it was.
    // Therefore, calculating currBin requires substracting this offset.
    int xItr{
        int((mol.comCoord.x + membraneObject.waterBox.x / 2) / subCellSize.x) -
        mpiContext.xOffset};
    int yItr{
        int((mol.comCoord.y + membraneObject.waterBox.y / 2) / subCellSize.y)};
    int zItr;
    if (membraneObject.waterBox.z > 0)
      zItr = int(-(mol.comCoord.z + 1E-6 - membraneObject.waterBox.z / 2.0) /
                 subCellSize.z);
    else
      zItr = 0;

    int xItrAdjusted, yItrAdjusted, zItrAdjusted;

    if (xItr < 0) xItrAdjusted = 0;
    else if (xItr > numSubCells.x -1) xItrAdjusted = numSubCells.x - 1;
    else xItrAdjusted = xItr;
    if (yItr < 0) yItrAdjusted = 0;
    else if (yItr > numSubCells.y -1) yItrAdjusted = numSubCells.y - 1;
    else yItrAdjusted = yItr;
    if (zItr < 0) zItrAdjusted = 0;
    else if (zItr > numSubCells.z -1) zItrAdjusted = numSubCells.z - 1;
    else zItrAdjusted = zItr;

    int currBin = xItrAdjusted + (yItrAdjusted * numSubCells.x) + (zItrAdjusted * numSubCells.x * numSubCells.y);

    mol.mySubVolIndex = currBin;
    if (currBin >= numSubCells.tot) {
      std::cerr << "Molecule " << mol.index << " (ID=" << mol.id << ")"
                << " seems outside simulation volume, with center of mass "
                    "coordinates ["
                << mol.comCoord << "].\n";
      std::cerr << "numSubCells.x = " << numSubCells.x
                << ", numSubCells.y = " << numSubCells.y
                << ", numSubCells.z = " << numSubCells.z << "\n";
      std::cerr << "xItr = " << xItr << ", yItr = " << yItr
                << ", zItr = " << zItr << "\n";
      std::cerr << "currBin = " << currBin
                << ", numSubCells.tot = " << numSubCells.tot
                << ", mpiContext.xOffset = " << mpiContext.xOffset << "\n";
      error("mol outside box.");
      exit(1);
    }
    subCellList[currBin].memberMolList.push_back(molItr);
  }    // loop over all molecules.
}

void SimulVolume::update_region_flags(std::vector<Molecule> &moleculeList,
      std::vector<Complex> &complexList,MpiContext& mpiContext) {

    // if not the first rank, locate the most left sub box along x-axis to update isLeftGhost, isLeftEdge for molecules and complexes
    if (mpiContext.rank > 0) {
      for (int yItr{0}; yItr < numSubCells.y; ++yItr) {
        for (int zItr{0}; zItr < numSubCells.z; ++zItr) {
          int leftBoxIndex = yItr * numSubCells.x + zItr * numSubCells.x * numSubCells.y;
          for (auto& molIndex : subCellList[leftBoxIndex].memberMolList) {
            moleculeList[molIndex].isLeftGhost = true;
            complexList[moleculeList[molIndex].myComIndex].isLeftGhost = true;
          }
          leftBoxIndex += 1;
          for (auto& molIndex : subCellList[leftBoxIndex].memberMolList) {
            moleculeList[molIndex].isLeftEdge = true;
            complexList[moleculeList[molIndex].myComIndex].isLeftEdge = true;
            for (auto& i : complexList[moleculeList[molIndex].myComIndex].memberList) {
              if (moleculeList[i].isLeftGhost == false)
                moleculeList[i].isLeftEdge = true;
            }
          }
        }
      }
    }
    // if not the last rank, locate the most right sub box along x-axis to update isRightGhost, isRightEdge for molecules and complexes
    if (mpiContext.rank < mpiContext.nprocs - 1) {
      for (int yItr{0}; yItr < numSubCells.y; ++yItr) {
        for (int zItr{0}; zItr < numSubCells.z; ++zItr) {
          int rightBoxIndex = (numSubCells.x - 1) + yItr * numSubCells.x + zItr * numSubCells.x * numSubCells.y;
          for (auto& molIndex : subCellList[rightBoxIndex].memberMolList) {
            moleculeList[molIndex].isRightGhost = true;
            complexList[moleculeList[molIndex].myComIndex].isRightGhost = true;
          }
          rightBoxIndex -= 1;
          for (auto& molIndex : subCellList[rightBoxIndex].memberMolList) {
            moleculeList[molIndex].isRightEdge = true;
            complexList[moleculeList[molIndex].myComIndex].isRightEdge = true;
            for (auto& i : complexList[moleculeList[molIndex].myComIndex].memberList) {
              if (moleculeList[i].isRightGhost == false)
                moleculeList[i].isRightEdge = true;
            }
          }
        }
      }
    }
}
