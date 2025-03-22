#include "classes/class_Coord.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Vector.hpp"
#include "io/io.hpp"
#include "parser/parser_functions.hpp"
#include "system_setup/system_setup.hpp"
#include "tracing.hpp"
#include <cctype>
#include <cmath>
#include <iostream>

void generate_coordinates(const Parameters &params,
                          std::vector<Molecule> &moleculeList,
                          std::vector<Complex> &complexList,
                          std::vector<MolTemplate> &molTemplateList,
                          const std::vector<ForwardRxn> &forwardRxns,
                          const Membrane &membraneObject,
                          std::string &coordinateFileName) {
  // First create all the molecules and their corresponding complexes, with
  // random center of mass coordinates
  for (auto &oneTemp : molTemplateList) {
    if (oneTemp.isImplicitLipid == false) {
      for (unsigned itr{0}; itr < oneTemp.copies; ++itr) {
        moleculeList.emplace_back(initialize_molecule(
            Complex::numberOfComplexes, params, oneTemp, membraneObject));
        complexList.emplace_back(initialize_complex(
            moleculeList.back(),
            molTemplateList[moleculeList.back().molTypeIndex]));
        oneTemp.monomerList.emplace_back(moleculeList.back().index);
        moleculeList.back().complexId = complexList.back().id;
      }
    } else {
      oneTemp.isLipid = true;
      moleculeList.emplace_back(initialize_molecule(
          Complex::numberOfComplexes, params, oneTemp, membraneObject));
      complexList.emplace_back(initialize_complex(
          moleculeList.back(),
          molTemplateList[moleculeList.back().molTypeIndex]));
      moleculeList.back().complexId = complexList.back().id;
    }
  }
  std::cout << "NUMBER OF MOLECULES IN GEN COORDS: " << moleculeList.size()
            << std::endl;

  if (moleculeList.size() == 0 && params.numTotalUnits == 0) {
    std::cout << "No molecules present, skipping coordinate generation.\n";
    return;
  }

  std::cout << "\nFinding and fixing overlapping proteins.\n";
  int currItr{0};
  while (currItr < 50) {
    bool hasOverlap{false};
    ++currItr;
    int numOverlap{0};

    unsigned long molListSize{moleculeList.size()};
    auto tmp = moleculeList[0];
    for (unsigned long mol1Itr{0}; mol1Itr < molListSize; ++mol1Itr) {
      auto &mol1 = moleculeList[mol1Itr]; // get a reference just to make the
                                          // code below less messy
      for (unsigned long mol2Itr{0}; mol2Itr < molListSize; ++mol2Itr) {
        auto &mol2 = moleculeList[mol2Itr];
        const MolTemplate &mol2Temp{molTemplateList[mol2.molTypeIndex]};

        if ((mol1Itr != mol2Itr) &&
            are_molecules_in_vicinity(mol1, mol2, molTemplateList)) {
          for (unsigned int iface1Itr{0}; iface1Itr < mol1.interfaceList.size();
               ++iface1Itr) {
            auto &iface1 = mol1.interfaceList[iface1Itr];
            for (unsigned int iface2Itr{0};
                 iface2Itr < mol2.interfaceList.size(); ++iface2Itr) {
              auto &iface2 = mol2.interfaceList[iface2Itr];
              int rxnIndex{0};
              bool theyInteract{false};
              //                            for (auto& rxn : forwardRxns) {
              for (auto rxnItr :
                   molTemplateList[moleculeList[mol1Itr].molTypeIndex]
                       .interfaceList[iface1Itr]
                       .stateList[0]
                       .myForwardRxns) {
                const ForwardRxn &rxn = forwardRxns[rxnItr];
                if ((rxn.reactantListNew[0].molTypeIndex == mol1.molTypeIndex &&
                     rxn.reactantListNew[1].molTypeIndex ==
                         mol2.molTypeIndex) ||
                    (rxn.reactantListNew[0].molTypeIndex == mol2.molTypeIndex &&
                     rxn.reactantListNew[1].molTypeIndex ==
                         mol1.molTypeIndex)) {
                  theyInteract = true;
                  rxnIndex = &rxn - &forwardRxns[0];
                  break;
                }
              }
              if (theyInteract) {
                Vector tmpVec{iface1.coord - iface2.coord};
                double mag = tmpVec.x * tmpVec.x + tmpVec.y * tmpVec.y +
                             tmpVec.z * tmpVec.z;
                if (mag < forwardRxns[rxnIndex].bindRadius *
                              forwardRxns[rxnIndex].bindRadius) {
                  mol2.create_random_coords(mol2Temp, membraneObject);
                  complexList[mol2.myComIndex].comCoord = mol2.comCoord;
                  hasOverlap = true;
                  ++numOverlap;
                  iface1Itr =
                      mol1.interfaceList.size(); // break from interface loops
                  iface2Itr = mol2.interfaceList.size();
                }
                //                                }
              }
            }
          }
        }
      }
    } // end iterating over all molecules.
    std::cout << "Iteration: " << currItr << "\nOverlapping Proteins ("
              << std::boolalpha << hasOverlap << "): " << numOverlap << '\n';
    if (!hasOverlap) {
      std::cout << "No overlapping proteins found.\n";
      break;
    }
  } // end while loop

  // replace molecules with given coordinates
  std::cout << "\ncoordinate file: " << coordinateFileName << std::endl;
  // The format is as follows:
  // newtype
  // <molecule name>
  // <number of coordinates>
  // <x1> <y1> <z1> 
  // <x2> <y2> <z2>
  // ...
  if (coordinateFileName != ""){
    std::cout << "Parsing given coordinates..." << std::endl;
    std::vector<int> changedMoleculeIndex{};
    std::ifstream inputFile{coordinateFileName};
    if (!inputFile) {
      std::cerr << "Coordinates file cannot be opened. Exiting..." << std::endl;
      exit(1);
    }
    while (!inputFile.eof()) {
      std::string line;
      getline(inputFile, line);
      std::string tmpLine{create_tmp_line(line)};

      // std::cout << tmpLine << std::endl;

      if (skipLine(tmpLine)) {
        continue;
      } else if (tmpLine == "newtype") {
        // read the type of molecule
        getline(inputFile, line);
        std::string tmpMolName = create_tmp_line(line);
        // std::cout << tmpLine << std::endl;
        int stratIndex{};
        for (int i = 0; i < moleculeList.size(); i++) {
          std::string molName{
              molTemplateList[moleculeList[i].molTypeIndex].molName};
          std::transform(molName.begin(), molName.end(), molName.begin(),
                         ::tolower);
          if (molName == tmpMolName) {
            stratIndex = i;
            break;
          }
        }

        // read the coordinate of center of mass
        getline(inputFile, line);
        tmpLine = create_tmp_line(line);
        int numCoords{std::stoi(tmpLine)};
        // std::cout << numCoords << std::endl;
        for (int i = 0; i < numCoords; i++) {
          // read the coordinate of center of mass
          double x{0};
          double y{0};
          double z{0};
          inputFile >> x >> y >> z;
          // update the coordinates
          changedMoleculeIndex.emplace_back(stratIndex + i);
          Molecule &movingMol = moleculeList[stratIndex + i];
          Vector transVec{x - movingMol.comCoord.x, y - movingMol.comCoord.y,
                          z - movingMol.comCoord.z};
          movingMol.comCoord = transVec + movingMol.comCoord;
          for (unsigned int j{0}; j < movingMol.interfaceList.size(); ++j)
            movingMol.interfaceList[j].coord =
                transVec + movingMol.interfaceList[j].coord;
        }
      }
    }
    // check overlap for changed molecules
    std::cout << "|\n";
    std::cout << "|-Finding and fixing overlapping proteins.\n";
    currItr = 0;
    while (currItr < 50) {
      bool hasOverlap{false};
      ++currItr;
      int numOverlap{0};

      unsigned long molListSize{moleculeList.size()};
      auto tmp = moleculeList[0];
      for (int mol1Itr : changedMoleculeIndex) {
        // std::cout << mol1Itr << std::endl;
        auto &mol1 = moleculeList[mol1Itr]; // get a reference just to make the
                                            // code below less messy
        for (unsigned long mol2Itr{0}; mol2Itr < molListSize; ++mol2Itr) {
          // if (std::find(v.begin(), v.end(), key) != v.end())
          // if molecule 2 is also changed?
          bool mol2Changed{false};
          if (std::find(changedMoleculeIndex.begin(),
                        changedMoleculeIndex.end(),
                        mol2Itr) != changedMoleculeIndex.end()) {
            mol2Changed = true;
          }

          auto &mol2 = moleculeList[mol2Itr];
          const MolTemplate &mol2Temp{molTemplateList[mol2.molTypeIndex]};

          if ((mol1Itr != mol2Itr) &&
              are_molecules_in_vicinity(mol1, mol2, molTemplateList)) {
            for (unsigned int iface1Itr{0};
                 iface1Itr < mol1.interfaceList.size(); ++iface1Itr) {
              auto &iface1 = mol1.interfaceList[iface1Itr];
              for (unsigned int iface2Itr{0};
                   iface2Itr < mol2.interfaceList.size(); ++iface2Itr) {
                auto &iface2 = mol2.interfaceList[iface2Itr];
                int rxnIndex{0};
                bool theyInteract{false};
                for (auto rxnItr :
                     molTemplateList[moleculeList[mol1Itr].molTypeIndex]
                         .interfaceList[iface1Itr]
                         .stateList[0]
                         .myForwardRxns) {
                  const ForwardRxn &rxn = forwardRxns[rxnItr];
                  if ((rxn.reactantListNew[0].molTypeIndex ==
                           mol1.molTypeIndex &&
                       rxn.reactantListNew[1].molTypeIndex ==
                           mol2.molTypeIndex) ||
                      (rxn.reactantListNew[0].molTypeIndex ==
                           mol2.molTypeIndex &&
                       rxn.reactantListNew[1].molTypeIndex ==
                           mol1.molTypeIndex)) {
                    theyInteract = true;
                    rxnIndex = &rxn - &forwardRxns[0];
                    break;
                  }
                }
                if (theyInteract) {
                  Vector tmpVec{iface1.coord - iface2.coord};
                  double mag = tmpVec.x * tmpVec.x + tmpVec.y * tmpVec.y +
                               tmpVec.z * tmpVec.z;
                  if (mag < forwardRxns[rxnIndex].bindRadius *
                                forwardRxns[rxnIndex].bindRadius) {
                    if (mol2Changed) {
                      std::cout
                          << "|-WARNING!!! "
                          << molTemplateList[mol1.molTypeIndex].molName
                          << " at " << mol1.comCoord.x << " " << mol1.comCoord.y
                          << " " << mol1.comCoord.z << " and "
                          << mol2Temp.molName << " at " << mol2.comCoord.x
                          << " " << mol2.comCoord.y << " " << mol2.comCoord.z
                          << " overlaps. Please check your input!" << std::endl;
                    } else {
                      mol2.create_random_coords(mol2Temp, membraneObject);
                      complexList[mol2.myComIndex].comCoord = mol2.comCoord;
                      hasOverlap = true;
                      ++numOverlap;
                      // break from interface loops
                      iface1Itr = mol1.interfaceList.size();
                      iface2Itr = mol2.interfaceList.size();
                    }
                  }
                }
              }
            }
          }
        }
      } // end iterating over changed molecules
      std::cout << "|-Iteration: " << currItr << "\n|-Overlapping Proteins ("
                << std::boolalpha << hasOverlap << "): " << numOverlap << '\n';
      if (!hasOverlap) {
        std::cout << "|-No overlapping proteins found.\n";
        break;
      }
    }
  }
  write_xyz("DATA/initial_crds.xyz", params, moleculeList, molTemplateList);
}
