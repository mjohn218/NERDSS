#include "io/io.hpp"
#include "tracing.hpp"
#include <chrono>
#include <ctime>

void read_restart(long long int& simItr, std::ifstream& restartFile, Parameters& params, SimulVolume& simulVolume,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<CreateDestructRxn>& createDestructRxns,
    std::map<std::string, int>& observablesList, Membrane& membraneObject, copyCounters& counterArrays)
{
    // TRACE();
    try {
        // Read parameters
        std::cout << "READ IN PARMATERS from restart file" << std::endl;
        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        {
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.nItr;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> simItr;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "Restarting simulation from iteration " << simItr << '\n';
            params.itrRestartFrom = simItr;

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.timeRestartFrom;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "Current simulation time (s): " << params.timeRestartFrom << '\n';

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.numMolTypes;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.numTotalSpecies;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.numTotalComplex;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.numTotalUnits;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.numLipids;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.timeStep;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.max2DRxns;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> membraneObject.waterBox.x >> membraneObject.waterBox.y >> membraneObject.waterBox.z;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            membraneObject.waterBox.volume = membraneObject.waterBox.x * membraneObject.waterBox.y * membraneObject.waterBox.z;

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> membraneObject.implicitlipidIndex >> membraneObject.nSites >> membraneObject.nStates >> membraneObject.No_free_lipids >> membraneObject.No_protein >> membraneObject.totalSA;

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            for (int i = 0; i < membraneObject.nStates; i++) {
                membraneObject.numberOfFreeLipidsEachState.emplace_back(0);
                restartFile >> membraneObject.numberOfFreeLipidsEachState[i];
            }

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> membraneObject.implicitLipid >> membraneObject.TwoD >> membraneObject.isBox >> membraneObject.isSphere >> membraneObject.sphereR;

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.overlapSepLimit;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.rMaxLimit;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.timeWrite;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.trajWrite;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.restartWrite;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.pdbWrite;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.checkPoint;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.scaleMaxDisplace;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.transitionWrite;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '=');
            restartFile >> params.clusterOverlapCheck;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            unsigned long lastUpdateTransitionSize { 0 };
            restartFile >> lastUpdateTransitionSize;
            for (unsigned itr { 0 }; itr < lastUpdateTransitionSize; ++itr) {
                int index { 0 };
                restartFile >> index;
                Parameters::lastUpdateTransition.push_back(index);
            }
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        std::cout << "restart write, pdbWrite: " << params.restartWrite << ' ' << params.pdbWrite << std::endl;
        /*	std::cout<<"READ IN SUB volume PARTITIONING from restart file"<<std::endl;
        // Read Simulation Volume
        {
            restartFile >> simulVolume.numSubCells.x >> simulVolume.numSubCells.y >> simulVolume.numSubCells.z
                >> simulVolume.numSubCells.tot;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile >> simulVolume.subCellSize.x >> simulVolume.subCellSize.y >> simulVolume.subCellSize.z;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            for (unsigned subCellItr { 0 }; subCellItr < simulVolume.numSubCells.tot; ++subCellItr) {
                SimulVolume::SubVolume subCell {};
                restartFile >> subCell.absIndex >> subCell.xIndex >> subCell.yIndex >> subCell.zIndex;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                unsigned neighborListSize { 0 };
                restartFile >> neighborListSize;
                for (unsigned itr { 0 }; itr < neighborListSize; ++itr) {
                    unsigned neighbor { 0 };
                    restartFile >> neighbor;
                    subCell.neighborList.push_back(neighbor);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                unsigned memberListSize { 0 };
                restartFile >> memberListSize;
                for (unsigned itr { 0 }; itr < memberListSize; ++itr) {
                    unsigned member { 0 };
                    restartFile >> member;
                    subCell.memberMolList.push_back(member);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                simulVolume.subCellList.push_back(subCell);
            }
        }
	*/
        std::cout << "READ IN MOL TEMPLATE from restart file" << std::endl;
        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        // Read MolTemplates
        {
            restartFile >> MolTemplate::numMolTypes;
            for (unsigned itr { 0 }; itr < MolTemplate::numMolTypes; ++itr) {
                unsigned num { 0 };
                restartFile >> num;
                MolTemplate::numEachMolType.push_back(num);
            }
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << " Num moltypes: " << MolTemplate::numMolTypes << '\n';
            unsigned absToRelIfaceSize { 0 };
            restartFile >> absToRelIfaceSize;
            for (unsigned itr { 0 }; itr < absToRelIfaceSize; ++itr) {
                unsigned iface { 0 };
                restartFile >> iface;
                MolTemplate::absToRelIface.push_back(iface);
            }
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            restartFile >> Interface::State::totalNumOfStates;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            for (unsigned itr { 0 }; itr < MolTemplate::numMolTypes; ++itr) {
                MolTemplate oneTemp {};
                restartFile >> oneTemp.molTypeIndex >> oneTemp.molName;
                std::cout << " protein index, name: " << oneTemp.molTypeIndex << ' ' << oneTemp.molName << '\n';
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> oneTemp.copies >> oneTemp.mass >> oneTemp.radius;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> oneTemp.isLipid >> oneTemp.isImplicitLipid >> oneTemp.isRod >> oneTemp.isPoint >> oneTemp.checkOverlap >> oneTemp.countTransition >> oneTemp.transitionMatrixSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> oneTemp.comCoord.x >> oneTemp.comCoord.y >> oneTemp.comCoord.z;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> oneTemp.D.x >> oneTemp.D.y >> oneTemp.D.z;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> oneTemp.Dr.x >> oneTemp.Dr.y >> oneTemp.Dr.z;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // reaction partners
                {
                    unsigned partSize { 0 };
                    restartFile >> partSize;

                    for (unsigned partItr { 0 }; partItr < partSize; ++partItr) {
                        int partner { 0 };
                        restartFile >> partner;
                        oneTemp.rxnPartners.push_back(partner);
                    }
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }

                // optional bonds
                {
                    unsigned bondSize { 0 };
                    restartFile >> bondSize;

                    for (unsigned bondItr { 0 }; bondItr < bondSize; ++bondItr) {
                        int iface1 { 0 };
                        int iface2 { 0 };
                        restartFile >> iface1 >> iface2;
                        oneTemp.bondList.push_back(std::array<int, 2> { { iface1, iface2 } });
                    }
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }

                // write interfaces
                unsigned oneTempIfaceSize { 0 };
                restartFile >> oneTempIfaceSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                for (unsigned ifaceItr { 0 }; ifaceItr < oneTempIfaceSize; ++ifaceItr) {
                    Interface tmpIface {};
                    restartFile >> tmpIface.index >> tmpIface.name;
                    std::cout << " iface index, name: " << tmpIface.index << ' ' << tmpIface.name << '\n';
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> tmpIface.iCoord.x >> tmpIface.iCoord.y >> tmpIface.iCoord.z;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                    unsigned stateListSize { 0 };
                    restartFile >> stateListSize;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    for (unsigned stateItr { 0 }; stateItr < stateListSize; ++stateItr) {
                        Interface::State tmpState {};
                        restartFile >> tmpState.index >> tmpState.iden;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                        // partner list
                        unsigned rxnPartnersSize { 0 };
                        restartFile >> rxnPartnersSize;
                        for (unsigned partItr { 0 }; partItr < rxnPartnersSize; ++partItr) {
                            unsigned partner { 0 };
                            restartFile >> partner;
                            tmpState.rxnPartners.push_back(partner);
                        }
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                        // reaction lists
                        unsigned myForwardRxnsSize { 0 };
                        restartFile >> myForwardRxnsSize;
                        for (unsigned rxnItr { 0 }; rxnItr < myForwardRxnsSize; ++rxnItr) {
                            int rxn { 0 };
                            restartFile >> rxn;
                            tmpState.myForwardRxns.push_back(rxn);
                        }
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                        unsigned myCreateDestructRxnsSize { 0 };
                        restartFile >> myCreateDestructRxnsSize;
                        for (unsigned rxnItr { 0 }; rxnItr < myCreateDestructRxnsSize; ++rxnItr) {
                            int rxn { 0 };
                            restartFile >> rxn;
                            tmpState.myCreateDestructRxns.push_back(rxn);
                        }
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                        unsigned stateChangeRxnsSize { 0 };
                        restartFile >> stateChangeRxnsSize;
                        for (unsigned rxnItr { 0 }; rxnItr < stateChangeRxnsSize; ++rxnItr) {
                            unsigned elem1 { 0 };
                            unsigned elem2 { 0 };
                            restartFile >> elem1 >> elem2;
                            tmpState.stateChangeRxns.emplace_back(elem1, elem2);
                        }
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                        // done reading state, add to state list
                        tmpIface.stateList.emplace_back(tmpState);
                    }
                    // done reading interface, add to interface list
                    oneTemp.interfaceList.emplace_back(tmpIface);
                }

                // read ifaceWithStates
                unsigned ifacesWithStatesSize { 0 };
                restartFile >> ifacesWithStatesSize;
                for (unsigned rxnItr { 0 }; rxnItr < ifacesWithStatesSize; ++rxnItr) {
                    unsigned elem { 0 };
                    restartFile >> elem;
                    oneTemp.ifacesWithStates.push_back(elem);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // read monomerList
                unsigned monomerListSize { 0 };
                restartFile >> monomerListSize;
                for (unsigned rxnItr { 0 }; rxnItr < monomerListSize; ++rxnItr) {
                    unsigned elem { 0 };
                    restartFile >> elem;
                    oneTemp.monomerList.push_back(elem);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // read lifetime
                if(oneTemp.countTransition == true) {
                    unsigned lifeTimeSize {0};
                    oneTemp.lifeTime.resize(oneTemp.transitionMatrixSize);
                    for(int indexOne = 0; indexOne < oneTemp.transitionMatrixSize; ++indexOne) {
                        restartFile >> lifeTimeSize;
                        for (unsigned rxnItr { 0 }; rxnItr < lifeTimeSize; ++rxnItr) {
                            double elem { 0.0 };
                            restartFile >> elem;
                            oneTemp.lifeTime[indexOne].push_back(elem);
                        }
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    }
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }

                // read transition matrix
                if(oneTemp.countTransition == true) {
                    oneTemp.transitionMatrix.resize(oneTemp.transitionMatrixSize);
                    for (int indexOne = 0; indexOne < oneTemp.transitionMatrixSize; ++indexOne) {
                        oneTemp.transitionMatrix[indexOne].resize(oneTemp.transitionMatrixSize);
                    }

                    for(int indexOne = 0; indexOne < oneTemp.transitionMatrixSize; ++indexOne){
                        for (int indexTwo = 0; indexTwo < oneTemp.transitionMatrixSize; ++indexTwo){
                            restartFile >> oneTemp.transitionMatrix[indexOne][indexTwo];
                        }
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    }
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                }

                // done reading template, add to template list
                molTemplateList.emplace_back(oneTemp);
            }
        }
        std::cout << "READ IN REACTIONs from restart file" << std::endl;
        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        // write Reactions
        {
            unsigned forwardRxnsSize { 0 };
            unsigned backRxnsSize { 0 };
            unsigned createDestructRxnsSize { 0 };
            restartFile >> RxnBase::numberOfRxns >> forwardRxnsSize >> backRxnsSize >> createDestructRxnsSize
                >> RxnBase::totRxnSpecies;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << " Num rxns: " << RxnBase::numberOfRxns << " forwardRxns: " << forwardRxnsSize << '\n';
            // forward reactions
            for (unsigned rxnItr { 0 }; rxnItr < forwardRxnsSize; ++rxnItr) {
                ForwardRxn tmpRxn;
                restartFile >> tmpRxn.absRxnIndex >> tmpRxn.relRxnIndex >> tmpRxn.rxnLabel;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                std::cout << " itr: " << rxnItr << " abs index, relindex: " << tmpRxn.absRxnIndex << ' ' << tmpRxn.relRxnIndex << '\n';
                int rxnType { -1 };
                restartFile >> rxnType >> tmpRxn.isSymmetric >> tmpRxn.isOnMem >> tmpRxn.hasStateChange;
                std::cout << "Rxntype: " << rxnType << '\n';
                if (rxnType != -1) {
                    tmpRxn.rxnType = static_cast<ReactionType>(rxnType); // turn the int rxnType into ReactionType
                } else {
                    std::cerr << "ERROR: Cannot parse reaction type for reaction " << tmpRxn.absRxnIndex
                              << ". Exiting.\n";
                    exit(1);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpRxn.isObserved;
                if (tmpRxn.isObserved)
                    restartFile >> tmpRxn.observeLabel;
                restartFile >> tmpRxn.productName;
                std::cout << "Product name: " << tmpRxn.productName << std::endl;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpRxn.isReversible >> tmpRxn.conjBackRxnIndex >> tmpRxn.irrevRingClosure >> tmpRxn.bindRadSameCom >> tmpRxn.loopCoopFactor >> tmpRxn.length3Dto2D;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpRxn.bindRadius;
                std::string th1, th2, ph1, ph2, omega;
                restartFile >> th1 >> th2 >> ph1 >> ph2 >> omega;

                tmpRxn.assocAngles.theta1 = std::stod(th1);
                tmpRxn.assocAngles.theta2 = std::stod(th2);
                tmpRxn.assocAngles.phi1 = std::stod(ph1);
                tmpRxn.assocAngles.phi2 = std::stod(ph2);
                tmpRxn.assocAngles.omega = std::stod(omega);
                if (th1 == "nan")
                    tmpRxn.assocAngles.theta1 = std::numeric_limits<double>::quiet_NaN();
                if (th2 == "nan")
                    tmpRxn.assocAngles.theta2 = std::numeric_limits<double>::quiet_NaN();
                if (ph1 == "nan")
                    tmpRxn.assocAngles.phi1 = std::numeric_limits<double>::quiet_NaN();
                if (ph2 == "nan")
                    tmpRxn.assocAngles.phi2 = std::numeric_limits<double>::quiet_NaN();
                if (omega == "nan")
                    tmpRxn.assocAngles.omega = std::numeric_limits<double>::quiet_NaN();
                //>> tmpRxn.assocAngles.theta1 >> tmpRxn.assocAngles.theta2
                //>> tmpRxn.assocAngles.phi1 >> tmpRxn.assocAngles.phi2 >> tmpRxn.assocAngles.omega;
                std::cout << "RXN angles " << tmpRxn.bindRadius << ' ' << tmpRxn.assocAngles.theta1 << ' ' << tmpRxn.assocAngles.theta2
                          << ' ' << tmpRxn.assocAngles.phi1 << ' ' << tmpRxn.assocAngles.phi2 << ' ' << tmpRxn.assocAngles.omega << '\n';
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> tmpRxn.norm1.x >> tmpRxn.norm1.y >> tmpRxn.norm1.z;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> tmpRxn.norm2.x >> tmpRxn.norm2.y >> tmpRxn.norm2.z;
                std::cout << " norm 2: " << tmpRxn.norm2.x << ' ' << tmpRxn.norm2.y << ' ' << tmpRxn.norm2.z << '\n';
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpRxn.excludeVolumeBound;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpRxn.isCoupled;
                std::cout << " reaction is coupled? " << tmpRxn.isCoupled << std::endl;
                if (tmpRxn.isCoupled) {
                    rxnType = -1;
                    std::cout << "did not enter iscoupled loop " << '\n';
                    restartFile >> tmpRxn.coupledRxn.absRxnIndex >> tmpRxn.coupledRxn.relRxnIndex >> rxnType >> tmpRxn.coupledRxn.label >> tmpRxn.coupledRxn.probCoupled;
                    if (rxnType != -1) {
                        tmpRxn.coupledRxn.rxnType = static_cast<ReactionType>(rxnType);
                    } else {
                        std::cerr << "ERROR: Cannot parse reaction type for reaction coupled to reaction "
                                  << tmpRxn.absRxnIndex << ". Exiting.\n";
                        exit(1);
                    }
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // integer reactants
                unsigned intReactantListSize { 0 };
                restartFile >> intReactantListSize;
                std::cout << "Nreactant first round " << intReactantListSize << '\n';
                for (unsigned itr { 0 }; itr < intReactantListSize; ++itr) {
                    int reactant { -1 };
                    restartFile >> reactant;
                    tmpRxn.intReactantList.push_back(reactant);
                    std::cout << " reactant: " << reactant << '\n';
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // integer products
                unsigned intProductListSize { 0 };
                restartFile >> intProductListSize;
                for (unsigned itr { 0 }; itr < intProductListSize; ++itr) {
                    int product { -1 };
                    restartFile >> product;
                    tmpRxn.intProductList.push_back(product);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // reactant list
                unsigned reactantListNewSize { 0 };
                restartFile >> reactantListNewSize;
                std::cout << "N reactants: " << reactantListNewSize << '\n';
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < reactantListNewSize; ++itr) {
                    RxnIface oneReact {};
                    restartFile >> oneReact.molTypeIndex;
                    std::cout << " molTypeIndex: " << oneReact.molTypeIndex << '\n';
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneReact.ifaceName >> oneReact.absIfaceIndex >> oneReact.relIfaceIndex;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneReact.requiresState >> oneReact.requiresInteraction;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    tmpRxn.reactantListNew.emplace_back(oneReact);
                    std::cout << " requiresState, requiresInteraction: " << oneReact.requiresState << ' ' << oneReact.requiresInteraction << '\n';
                }

                // product list
                unsigned productListNewSize { 0 };
                restartFile >> productListNewSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < productListNewSize; ++itr) {
                    RxnIface oneProd {};
                    restartFile >> oneProd.molTypeIndex;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneProd.ifaceName >> oneProd.absIfaceIndex >> oneProd.relIfaceIndex;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneProd.requiresState >> oneProd.requiresInteraction;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    tmpRxn.productListNew.emplace_back(oneProd);
                }

                // rate list
                unsigned rateListSize { 0 };
                restartFile >> rateListSize;
                std::cout << "Nrates: " << rateListSize << '\n';
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < rateListSize; ++itr) {
                    RxnBase::RateState oneRate {};
                    restartFile >> oneRate.rate;
                    std::cout << " rate: " << oneRate.rate << '\n';
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                    unsigned otherIfaceListsSize { 0 };
                    unsigned ifaceItr { 0 };
                    restartFile >> otherIfaceListsSize;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    for (unsigned listItr { 0 }; listItr < otherIfaceListsSize; ++listItr) {
                        unsigned oneListSize { 0 };
                        std::vector<RxnIface> tmpIfaceVec {};
                        restartFile >> oneListSize;
                        std::cout << "onelistsize: " << oneListSize << '\n';
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        for (unsigned anccIfaceItr { 0 }; anccIfaceItr < oneListSize; ++anccIfaceItr) {
                            RxnIface otherIface {};
                            restartFile >> otherIface.molTypeIndex >> otherIface.ifaceName >> otherIface.absIfaceIndex
                                >> otherIface.relIfaceIndex >> otherIface.requiresState
                                >> otherIface.requiresInteraction;
                            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                            tmpIfaceVec.emplace_back(otherIface);
                        }
                        oneRate.otherIfaceLists.push_back(tmpIfaceVec);
                    }
                    tmpRxn.rateList.emplace_back(oneRate);
                }
                tmpRxn.display();
                forwardRxns.emplace_back(tmpRxn);
            }
            std::cout << " Done with forward reactions " << '\n';
            // backRxns
            for (unsigned rxnItr { 0 }; rxnItr < backRxnsSize; ++rxnItr) {
                BackRxn tmpRxn;
                restartFile >> tmpRxn.absRxnIndex >> tmpRxn.relRxnIndex;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                int rxnType { -1 };
                restartFile >> rxnType >> tmpRxn.isSymmetric >> tmpRxn.isOnMem >> tmpRxn.hasStateChange;
                tmpRxn.rxnType = static_cast<ReactionType>(rxnType); // turn the int rxnType into ReactionType
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> tmpRxn.isObserved;
                if (tmpRxn.isObserved)
                    restartFile >> tmpRxn.observeLabel;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> tmpRxn.conjForwardRxnIndex;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> tmpRxn.isCoupled;
                std::cout << " reaction is coupled? " << tmpRxn.isCoupled << std::endl;
                if (tmpRxn.isCoupled) {
                    rxnType = -1;
                    std::cout << "did not enter iscoupled loop " << '\n';
                    restartFile >> tmpRxn.coupledRxn.absRxnIndex >> tmpRxn.coupledRxn.relRxnIndex >> rxnType >> tmpRxn.coupledRxn.label >> tmpRxn.coupledRxn.probCoupled;
                    if (rxnType != -1) {
                        tmpRxn.coupledRxn.rxnType = static_cast<ReactionType>(rxnType);
                    } else {
                        std::cerr << "ERROR: Cannot parse reaction type for reaction coupled to reaction "
                                  << tmpRxn.absRxnIndex << ". Exiting.\n";
                        exit(1);
                    }
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // integer reactants
                unsigned intReactantListSize { 0 };
                restartFile >> intReactantListSize;
                for (unsigned itr { 0 }; itr < intReactantListSize; ++itr) {
                    int reactant { -1 };
                    restartFile >> reactant;
                    tmpRxn.intReactantList.push_back(reactant);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // integer products
                unsigned intProductListSize { 0 };
                restartFile >> intProductListSize;
                for (unsigned itr { 0 }; itr < intProductListSize; ++itr) {
                    int product { -1 };
                    restartFile >> product;
                    tmpRxn.intProductList.push_back(product);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // reactant list
                unsigned reactantListNewSize { 0 };
                restartFile >> reactantListNewSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < reactantListNewSize; ++itr) {
                    RxnIface oneReact {};
                    restartFile >> oneReact.molTypeIndex;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneReact.ifaceName >> oneReact.absIfaceIndex >> oneReact.relIfaceIndex;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneReact.requiresState >> oneReact.requiresInteraction;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    tmpRxn.reactantListNew.emplace_back(oneReact);
                }

                // product list
                unsigned productListNewSize { 0 };
                restartFile >> productListNewSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < productListNewSize; ++itr) {
                    RxnIface oneProd {};
                    restartFile >> oneProd.molTypeIndex;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneProd.ifaceName >> oneProd.absIfaceIndex >> oneProd.relIfaceIndex;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> oneProd.requiresState >> oneProd.requiresInteraction;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    tmpRxn.productListNew.emplace_back(oneProd);
                }

                // rate list
                unsigned rateListSize { 0 };
                restartFile >> rateListSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < rateListSize; ++itr) {
                    RxnBase::RateState oneRate {};
                    restartFile >> oneRate.rate;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                    unsigned otherIfaceListSize { 0 };
                    restartFile >> otherIfaceListSize;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    for (unsigned listItr { 0 }; listItr < otherIfaceListSize; ++listItr) {
                        unsigned oneListSize { 0 };
                        restartFile >> oneListSize;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        std::vector<RxnIface> tmpIfaceVec;
                        for (unsigned anccIfaceItr { 0 }; anccIfaceItr < oneListSize; ++anccIfaceItr) {
                            RxnIface otherIface {};
                            restartFile >> otherIface.molTypeIndex >> otherIface.ifaceName >> otherIface.absIfaceIndex
                                >> otherIface.relIfaceIndex >> otherIface.requiresState
                                >> otherIface.requiresInteraction;
                            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                            tmpIfaceVec.emplace_back(otherIface);
                        }
                        oneRate.otherIfaceLists.push_back(tmpIfaceVec);
                    }
                    tmpRxn.rateList.emplace_back(oneRate);
                }
                backRxns.emplace_back(tmpRxn);
            }
            std::cout << "Done with back reactions " << '\n';
            // creation and destruction reactions
            std::cout << "Now creation and destruction " << '\n';
            for (unsigned rxnItr { 0 }; rxnItr < createDestructRxnsSize; ++rxnItr) {
                CreateDestructRxn tmpRxn {};
                restartFile >> tmpRxn.absRxnIndex >> tmpRxn.relRxnIndex;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                int rxnType { -1 };
                restartFile >> rxnType >> tmpRxn.isOnMem;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpRxn.isObserved; //>>tmpRxn.observeLabel;
                if (tmpRxn.isObserved)
                    restartFile >> tmpRxn.observeLabel;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                tmpRxn.rxnType = static_cast<ReactionType>(rxnType);
                restartFile >> tmpRxn.creationRadius;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // integer reactants
                unsigned intReactantListSize { 0 };
                restartFile >> intReactantListSize;
                for (unsigned itr { 0 }; itr < intReactantListSize; ++itr) {
                    int reactant { -1 };
                    restartFile >> reactant;
                    tmpRxn.intReactantList.push_back(reactant);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // integer products
                unsigned intProductListSize { 0 };
                restartFile >> intProductListSize;
                for (unsigned itr { 0 }; itr < intProductListSize; ++itr) {
                    int product { -1 };
                    restartFile >> product;
                    tmpRxn.intProductList.push_back(product);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // reactant list
                unsigned reactantMolListSize { 0 };
                restartFile >> reactantMolListSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < reactantMolListSize; ++itr) {
                    CreateDestructRxn::CreateDestructMol oneMol {};
                    unsigned interfaceListSize { 0 };
                    restartFile >> oneMol.molTypeIndex >> oneMol.molName >> interfaceListSize;
                    std::cout << " Destroy molecule: " << oneMol.molName << std::endl;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    for (unsigned ifaceItr { 0 }; ifaceItr < interfaceListSize; ++ifaceItr) {
                        RxnIface tmpIface {};
                        restartFile >> tmpIface.molTypeIndex;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        restartFile >> tmpIface.ifaceName >> tmpIface.absIfaceIndex >> tmpIface.relIfaceIndex;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        restartFile >> tmpIface.requiresState >> tmpIface.requiresInteraction;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        oneMol.interfaceList.emplace_back(tmpIface);
                    }
                    tmpRxn.reactantMolList.emplace_back(oneMol);
                }

                // product list
                unsigned productMolListSize { 0 };
                restartFile >> productMolListSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                for (unsigned itr { 0 }; itr < productMolListSize; ++itr) {
                    CreateDestructRxn::CreateDestructMol oneMol {};
                    unsigned interfaceListSize { 0 };
                    restartFile >> oneMol.molTypeIndex >> oneMol.molName >> interfaceListSize;
                    std::cout << " CREATION OF MOL: " << oneMol.molName << std::endl;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    for (unsigned ifaceItr { 0 }; ifaceItr < interfaceListSize; ++ifaceItr) {
                        RxnIface tmpIface {};
                        restartFile >> tmpIface.molTypeIndex;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        restartFile >> tmpIface.ifaceName >> tmpIface.absIfaceIndex >> tmpIface.relIfaceIndex;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        restartFile >> tmpIface.requiresState >> tmpIface.requiresInteraction;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        oneMol.interfaceList.emplace_back(tmpIface);
                    }
                    tmpRxn.productMolList.emplace_back(oneMol);
                }

                unsigned rateListSize { 0 };
                restartFile >> rateListSize;
                for (unsigned itr { 0 }; itr < rateListSize; ++itr) {
                    RxnBase::RateState tmpRate {};
                    restartFile >> tmpRate.rate;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                    unsigned otherIfaceListsSize { 0 };
                    restartFile >> otherIfaceListsSize;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                    // for (unsigned ifaceItr { 0 }; ifaceItr < otherIfaceListSize; ++ifaceItr) {
                    //     RxnIface tmpIface {};
                    //     restartFile >> tmpIface.molTypeIndex >> tmpIface.ifaceName >> tmpIface.absIfaceIndex
                    //         >> tmpIface.relIfaceIndex >> tmpIface.requiresState >> tmpIface.requiresInteraction;
                    //     restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    //     tmpRate.otherIfaceLists.emplace_back(std::vector<RxnIface> { tmpIface });
                    // }

                    // restartFile >> otherIfaceListsSize;
                    // restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    for (unsigned listItr { 0 }; listItr < otherIfaceListsSize; ++listItr) {
                        unsigned oneListSize { 0 };
                        std::vector<RxnIface> tmpIfaceVec {};
                        restartFile >> oneListSize;
                        std::cout << "onelistsize: " << oneListSize << '\n';
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        for (unsigned anccIfaceItr { 0 }; anccIfaceItr < oneListSize; ++anccIfaceItr) {
                            RxnIface otherIface {};
                            restartFile >> otherIface.molTypeIndex >> otherIface.ifaceName >> otherIface.absIfaceIndex
                                >> otherIface.relIfaceIndex >> otherIface.requiresState
                                >> otherIface.requiresInteraction;
                            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                            tmpIfaceVec.emplace_back(otherIface);
                        }
                        tmpRate.otherIfaceLists.push_back(tmpIfaceVec);
                    }
                    //copied new version up to here.
                    tmpRxn.rateList.emplace_back(tmpRate);
                }
                createDestructRxns.emplace_back(tmpRxn);
            }
        }
        std::cout << "Now read in coordinates " << std::endl;
        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        // write Molecules
        {
            int molListSize { 0 };

            restartFile >> molListSize >> Molecule::numberOfMolecules;
            std::cout << "Mol list size and molecule.numberofMolecules: " << molListSize << ' ' << Molecule::numberOfMolecules << std::endl;
            for (unsigned molItr { 0 }; molItr < molListSize; ++molItr) {
                Molecule tmpMol {};
                restartFile >> tmpMol.index >> tmpMol.isEmpty >> tmpMol.myComIndex >> tmpMol.molTypeIndex
                    >> tmpMol.mySubVolIndex;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpMol.mass >> tmpMol.isLipid >> tmpMol.isImplicitLipid >> tmpMol.linksToSurface >> tmpMol.isEmpty;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // center of mass
                restartFile >> std::fixed >> tmpMol.comCoord.x >> tmpMol.comCoord.y >> tmpMol.comCoord.z;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // interface lists
                unsigned freeListSize { 0 };
                restartFile >> freeListSize;
                for (unsigned itr { 0 }; itr < freeListSize; ++itr) {
                    int freeIface { 0 };
                    restartFile >> freeIface;
                    tmpMol.freelist.emplace_back(freeIface);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                unsigned bndListSize { 0 };
                restartFile >> bndListSize;
                for (unsigned itr { 0 }; itr < bndListSize; ++itr) {
                    int bndIface { 0 };
                    restartFile >> bndIface;
                    tmpMol.bndlist.emplace_back(bndIface);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                unsigned bndPartnerListSize { 0 };
                restartFile >> bndPartnerListSize;
                for (unsigned itr { 0 }; itr < bndPartnerListSize; ++itr) {
                    int bndPartner { 0 };
                    restartFile >> bndPartner;
                    tmpMol.bndpartner.emplace_back(bndPartner);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // interfaces
                unsigned interfaceListSize { 0 };
                restartFile >> interfaceListSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                for (unsigned itr { 0 }; itr < interfaceListSize; ++itr) {
                    Molecule::Iface tmpIface {};
                    restartFile >> tmpIface.index >> tmpIface.relIndex >> tmpIface.molTypeIndex >> tmpIface.stateIndex
                        >> tmpIface.stateIden >> tmpIface.isBound;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    restartFile >> std::fixed >> tmpIface.coord.x >> tmpIface.coord.y >> tmpIface.coord.z;
                    restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                    if (tmpIface.isBound) {
                        restartFile >> tmpIface.interaction.partnerIndex >> tmpIface.interaction.partnerIfaceIndex
                            >> tmpIface.interaction.conjBackRxn;
                        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    }
                    tmpMol.interfaceList.emplace_back(tmpIface);
                }

                // reweighting lists
                unsigned listSize { 0 };
                restartFile >> listSize;
                for (unsigned itr { 0 }; itr < listSize; ++itr) {
                    int elem { 0 };
                    restartFile >> elem;
                    tmpMol.prevlist.emplace_back(elem);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> listSize;
                for (unsigned itr { 0 }; itr < listSize; ++itr) {
                    int elem { 0 };
                    restartFile >> elem;
                    tmpMol.prevmyface.emplace_back(elem);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> listSize;
                for (unsigned itr { 0 }; itr < listSize; ++itr) {
                    int elem { 0 };
                    restartFile >> elem;
                    tmpMol.prevpface.emplace_back(elem);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> listSize;
                for (unsigned itr { 0 }; itr < listSize; ++itr) {
                    double elem { 0 };
                    restartFile >> elem;
                    tmpMol.prevnorm.emplace_back(elem);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> listSize;
                for (unsigned itr { 0 }; itr < listSize; ++itr) {
                    double elem { 0 };
                    restartFile >> elem;
                    tmpMol.ps_prev.emplace_back(elem);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> listSize;
                for (unsigned itr { 0 }; itr < listSize; ++itr) {
                    double elem { 0 };
                    restartFile >> elem;
                    tmpMol.prevsep.emplace_back(elem);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                moleculeList.emplace_back(tmpMol);
                //std::cout <<"read in : "<<tmpMol.index<<" first interface z crd: "<<tmpMol.comCoord.z<<std::endl;
            }

            unsigned long emptyMolListSize { 0 };
            restartFile >> emptyMolListSize;
            for (unsigned itr { 0 }; itr < emptyMolListSize; ++itr) {
                int index { 0 };
                restartFile >> index;
                Molecule::emptyMolList.push_back(index);
            }
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "N empty molecules: " << emptyMolListSize << std::endl;
        }
        std::cout << "Now read in complexes from RESTART " << '\n';
        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        // read Complexes
        {
            int comListSize { 0 };
            restartFile >> comListSize >> Complex::numberOfComplexes;
            std::cout << " Ncomplexes including empties: " << comListSize << " N actual complexes: " << Complex::numberOfComplexes << std::endl;
            for (unsigned comItr { 0 }; comItr < comListSize; ++comItr) {
                Complex tmpCom {};
                restartFile >> tmpCom.index >> tmpCom.isEmpty >> tmpCom.radius >> tmpCom.mass;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpCom.linksToSurface >> tmpCom.iLipidIndex >> tmpCom.OnSurface;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpCom.comCoord.x >> tmpCom.comCoord.y >> tmpCom.comCoord.z;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpCom.D.x >> tmpCom.D.y >> tmpCom.D.z;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                restartFile >> tmpCom.Dr.x >> tmpCom.Dr.y >> tmpCom.Dr.z;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // member molecule lists
                unsigned listSize { 0 };
                restartFile >> listSize;
                for (unsigned itr { 0 }; itr < listSize; ++itr) {
                    int memMol { -1 };
                    restartFile >> memMol;
                    tmpCom.memberList.emplace_back(memMol);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                // numEachMol list
                restartFile >> listSize;
                for (unsigned itr { 0 }; itr < listSize; ++itr) {
                    int memMol { -1 };
                    restartFile >> memMol;
                    tmpCom.numEachMol.emplace_back(memMol);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

                restartFile >> listSize;
                for (unsigned itr { 0 }; itr < listSize; ++itr) {
                    int memMol { -1 };
                    restartFile >> memMol;
                    tmpCom.lastNumberUpdateItrEachMol.emplace_back(memMol);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                complexList.emplace_back(tmpCom);
            }

            unsigned long emptyComListSize { 0 };
            restartFile >> emptyComListSize;
            std::cout << "N empty complexes " << emptyComListSize << '\t';
            for (unsigned itr { 0 }; itr < emptyComListSize; ++itr) {
                int index { 0 };
                restartFile >> index;
                Complex::emptyComList.push_back(index);
            }
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        } //done reading complexes
        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        // read observables
        {
            int numObs { 0 };
            restartFile >> numObs;
            std::cout << "N observables " << numObs << '\t';
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            for (unsigned itr { 0 }; itr < numObs; ++itr) {
                std::string obsName;
                int obsCount;
                restartFile >> obsName >> obsCount;
                observablesList.emplace(obsName, obsCount);
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
        }
        restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        // read counterArrays
        {
            restartFile >> counterArrays.nLoops >> counterArrays.nCancelOverlapPartner >> counterArrays.nCancelOverlapSystem >> counterArrays.nCancelDisplace2D >> counterArrays.nCancelDisplace3D >> counterArrays.nCancelDisplace3Dto2D >> counterArrays.nCancelSpanBox >> counterArrays.nAssocSuccess >> counterArrays.eventArraySize;
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            //read events3D, 3Dto2D, and 2D
            for (int i = 0; i < counterArrays.eventArraySize; i++)
                restartFile >> counterArrays.events3D[i];
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            for (int i = 0; i < counterArrays.eventArraySize; i++)
                restartFile >> counterArrays.events3Dto2D[i];
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            for (int i = 0; i < counterArrays.eventArraySize; i++)
                restartFile >> counterArrays.events2D[i];
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            int numSpecies { 0 };
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            restartFile >> numSpecies;
            std::cout << "N species " << numSpecies << '\t';
            restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            for (unsigned itr { 0 }; itr < numSpecies; ++itr) {
                unsigned long listSize { 0 };
                restartFile >> listSize;
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                counterArrays.bindPairList.emplace_back();
                std::cout << "Specie " << itr << " N bindPairs " << listSize << '\t';
                for (unsigned itr2 { 0 }; itr2 < listSize; ++itr2) {
                    int index { 0 };
                    restartFile >> index;
                    counterArrays.bindPairList[itr].push_back(index);
                }
                restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            // for (unsigned itr { 0 }; itr < numSpecies; ++itr) {
            //     unsigned long listSize { 0 };
            //     restartFile >> listSize;
            //     restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            //     counterArrays.bindPairListIL2D.emplace_back();
            //     std::cout << "Specie " << itr << " N bindPairs " << listSize << '\t';
            //     for (unsigned itr2 { 0 }; itr2 < listSize; ++itr2) {
            //         int index { 0 };
            //         restartFile >> index;
            //         counterArrays.bindPairListIL2D[itr].push_back(index);
            //     }
            //     restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            // }
            // for (unsigned itr { 0 }; itr < numSpecies; ++itr) {
            //     unsigned long listSize { 0 };
            //     restartFile >> listSize;
            //     restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            //     counterArrays.bindPairListIL3D.emplace_back();
            //     std::cout << "Specie " << itr << " N bindPairs " << listSize << '\t';
            //     for (unsigned itr2 { 0 }; itr2 < listSize; ++itr2) {
            //         int index { 0 };
            //         restartFile >> index;
            //         counterArrays.bindPairListIL3D[itr].push_back(index);
            //     }
            //     restartFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            // }
        }
    } catch (const std::string& msg) {
        std::cerr << msg << '\n';
        exit(1);
    } catch (const std::length_error& e) {
        std::cerr << "Error in reading template vectors for " << e.what() << '\n';
        exit(1);
    }
}
