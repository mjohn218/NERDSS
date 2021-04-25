#include "io/io.hpp"
#include "tracing.hpp"
#include <chrono>
#include <ctime>
#include <iomanip>

void write_restart(long long int simItr, std::ofstream& restartFile, const Parameters& params, const SimulVolume& simulVolume,
    const std::vector<Molecule>& moleculeList, const std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, const std::vector<CreateDestructRxn>& createDestructRxns,
    const std::map<std::string, int>& observablesList, const Membrane& membraneObject, const copyCounters& counterArrays)
{
    // TRACE();
    // Write parameters
    restartFile.precision(20);
    {
        restartFile << "#Parameters--update these for restart \n";
        restartFile << "numItr = " << params.nItr << '\n';
        restartFile << "currItr = " << simItr << '\n';
        restartFile << "currSimTime (s) = " << std::scientific << (simItr - params.itrRestartFrom) * params.timeStep * 1E-6 + params.timeRestartFrom << '\n';
        restartFile << "numMolTypes = " << params.numMolTypes << '\n';
        restartFile << "numTotalSpecies = " << params.numTotalSpecies << '\n';
        restartFile << "numComplexs = " << params.numTotalComplex << '\n';
        restartFile << "numTotalUnits = " << params.numTotalUnits << '\n';
        restartFile << "numLipids = " << params.numLipids << '\n';
        restartFile << "timestep = " << std::fixed << params.timeStep << '\n';
        restartFile << "max2DRxns = " << params.max2DRxns << '\n';
        restartFile << "simulDimensions = " << membraneObject.waterBox.x << ' ' << membraneObject.waterBox.y << ' ' << membraneObject.waterBox.z
                    << '\n';
        restartFile << "membrane = " << membraneObject.implicitlipidIndex << ' ' << membraneObject.nSites << ' ' << membraneObject.nStates << ' ' << membraneObject.No_free_lipids << ' ' << membraneObject.No_protein << ' ' << membraneObject.totalSA << '\n';
        restartFile << "implicitLipidStates = ";
        for (int i = 0; i < membraneObject.nStates; i++) {
            restartFile << membraneObject.numberOfFreeLipidsEachState[i];
            if (i != membraneObject.nStates - 1) {
                restartFile << ' ';
            }
        }
        restartFile << '\n';
        restartFile << "implicitLipidsParams = " << membraneObject.implicitLipid << ' ' << membraneObject.TwoD << ' ' << membraneObject.isBox << ' ' << membraneObject.isSphere << ' ' << membraneObject.sphereR << '\n';
        restartFile << "ifaceOverlapSepLimit = " << params.overlapSepLimit << '\n';
        restartFile << "rMaxLimit = " << params.rMaxLimit << '\n';
        restartFile << "timeWrite = " << params.timeWrite << '\n';
        restartFile << "trajWrite = " << params.trajWrite << '\n';
        restartFile << "restartWrite = " << params.restartWrite << '\n';
        restartFile << "pdbWrite = " << params.pdbWrite << '\n';
        restartFile << "checkPoint = " << params.checkPoint << '\n';
        restartFile << "scaleMaxDisplace = " << params.scaleMaxDisplace << '\n';
        restartFile << "transitionWrite = " << params.transitionWrite << '\n';
        restartFile << "clusterOverlapCheck = " << params.clusterOverlapCheck << '\n';

        restartFile << Parameters::lastUpdateTransition.size();
        for (auto& index : Parameters::lastUpdateTransition)
            restartFile << ' ' << index;
        restartFile << '\n';
    }
    /*
    // write Simulation Volume
    {
        restartFile << simulVolume.numSubCells.x << ' ' << simulVolume.numSubCells.y << ' ' << simulVolume.numSubCells.z
                    << ' ' << simulVolume.numSubCells.tot << '\n';
        restartFile << simulVolume.subCellSize.x << ' ' << simulVolume.subCellSize.y << ' ' << simulVolume.subCellSize.z
                    << '\n';
        for (auto& subCell : simulVolume.subCellList) {
            restartFile << subCell.absIndex << ' ' << subCell.xIndex << ' ' << subCell.yIndex << ' ' << subCell.zIndex
                        << '\n';
            restartFile << subCell.neighborList.size();
            for (const auto& neighbor : subCell.neighborList)
                restartFile << ' ' << neighbor;
            restartFile << '\n';

            restartFile << subCell.memberMolList.size();
            for (const auto& memMol : subCell.memberMolList)
                restartFile << ' ' << memMol;
            restartFile << '\n';
        }
    }
    */
    // write MolTemplates
    {
        restartFile << "#MolTemplates \n";
        restartFile << MolTemplate::numMolTypes;
        for (auto& num : MolTemplate::numEachMolType)
            restartFile << ' ' << num;
        restartFile << '\n';
        restartFile << MolTemplate::absToRelIface.size();
        for (auto& iface : MolTemplate::absToRelIface)
            restartFile << ' ' << iface;
        restartFile << '\n';
        restartFile << Interface::State::totalNumOfStates << '\n';

        for (const auto& oneTemp : molTemplateList) {
            restartFile << oneTemp.molTypeIndex << ' ' << oneTemp.molName << '\n';
            restartFile << oneTemp.copies << ' ' << oneTemp.mass << ' ' << oneTemp.radius << '\n';
            restartFile << oneTemp.isLipid << ' ' << oneTemp.isImplicitLipid << ' ' << oneTemp.isRod << ' ' << oneTemp.isPoint << ' '
                        << oneTemp.checkOverlap << ' ' << oneTemp.countTransition << ' ' << oneTemp.transitionMatrixSize << '\n';
            restartFile << oneTemp.comCoord.x << ' ' << oneTemp.comCoord.y << ' ' << oneTemp.comCoord.z
                        << '\n';
            restartFile << oneTemp.D.x << ' ' << oneTemp.D.y << ' ' << oneTemp.D.z << '\n';
            restartFile << oneTemp.Dr.x << ' ' << oneTemp.Dr.y << ' ' << oneTemp.Dr.z << '\n';

            // reaction partners
            restartFile << oneTemp.rxnPartners.size();
            for (const auto& partner : oneTemp.rxnPartners)
                restartFile << ' ' << partner;
            restartFile << '\n';

            // optional bonds
            restartFile << oneTemp.bondList.size() << '\n';
            for (const auto& bond : oneTemp.bondList)
                restartFile << bond[0] << ' ' << bond[1] << '\n';

            // write interfaces
            restartFile << oneTemp.interfaceList.size() << '\n';
            for (const auto& oneIface : oneTemp.interfaceList) {
                restartFile << oneIface.index << ' ' << oneIface.name << '\n';
                restartFile << std::fixed << oneIface.iCoord.x << ' ' << oneIface.iCoord.y << ' ' << oneIface.iCoord.z
                            << '\n';
                restartFile << oneIface.stateList.size() << '\n';
                for (auto& oneState : oneIface.stateList) {
                    restartFile << oneState.index << ' ' << oneState.iden << '\n';

                    // partner list
                    restartFile << oneState.rxnPartners.size();
                    for (auto elem : oneState.rxnPartners)
                        restartFile << ' ' << elem;
                    restartFile << '\n';

                    // reaction lists
                    restartFile << oneState.myForwardRxns.size();
                    for (auto elem : oneState.myForwardRxns)
                        restartFile << ' ' << elem;
                    restartFile << '\n';
                    restartFile << oneState.myCreateDestructRxns.size();
                    for (auto elem : oneState.myCreateDestructRxns)
                        restartFile << ' ' << elem;
                    restartFile << '\n';
                    restartFile << oneState.stateChangeRxns.size();
                    for (auto elem : oneState.stateChangeRxns)
                        restartFile << ' ' << elem.first << ' ' << elem.second;
                    restartFile << '\n';
                }
            }

            //write ifacesWithStates
            restartFile << oneTemp.ifacesWithStates.size();
            for (auto elem : oneTemp.ifacesWithStates) {
                restartFile << ' ' << elem;
            }
            restartFile << '\n';

            //write monomerList
            restartFile << oneTemp.monomerList.size();
            for (auto elem : oneTemp.monomerList) {
                restartFile << ' ' << elem;
            }
            restartFile << '\n';

            //write lifetime
            if(oneTemp.countTransition == true){
                for(int indexOne = 0; indexOne < oneTemp.transitionMatrixSize; ++indexOne){
                    restartFile << oneTemp.lifeTime[indexOne].size();
                    for (auto elem : oneTemp.lifeTime[indexOne]) {
                        restartFile << ' ' << elem;
                    }
                    restartFile << '\n';
                }
                restartFile << '\n';
            }

            //write transition matrix
            if(oneTemp.countTransition == true){
                for(int indexOne = 0; indexOne < oneTemp.transitionMatrixSize; ++indexOne){
                    for (int indexTwo = 0; indexTwo < oneTemp.transitionMatrixSize; ++indexTwo){
                        restartFile << ' ' << oneTemp.transitionMatrix[indexOne][indexTwo];
                    }
                    restartFile << '\n';
                }
                restartFile << '\n';
            }
        }
    }

    // write Reactions
    {
        restartFile << "#Reactions \n";
        restartFile << RxnBase::numberOfRxns << ' ' << forwardRxns.size() << ' ' << backRxns.size() << ' '
                    << createDestructRxns.size() << ' ' << RxnBase::totRxnSpecies << '\n';

        // forward reactions
        for (const auto& oneRxn : forwardRxns) {
            restartFile << oneRxn.absRxnIndex << ' ' << oneRxn.relRxnIndex << ' ' << oneRxn.rxnLabel << '\n';
            restartFile << static_cast<std::underlying_type<ReactionType>::type>(oneRxn.rxnType) << ' '
                        << oneRxn.isSymmetric << ' ' << oneRxn.isOnMem << ' ' << oneRxn.hasStateChange << '\n';
            restartFile << oneRxn.isObserved << ' ' << oneRxn.observeLabel << ' ' << oneRxn.productName << '\n';
            restartFile << oneRxn.isReversible << ' ' << oneRxn.conjBackRxnIndex << ' ' << oneRxn.irrevRingClosure << ' ' << oneRxn.bindRadSameCom << ' '
                        << oneRxn.loopCoopFactor << '\n'
                        << oneRxn.length3Dto2D << '\n';
            restartFile << std::setprecision(20) << oneRxn.bindRadius << ' ' << oneRxn.assocAngles.theta1 << ' ' << oneRxn.assocAngles.theta2
                        << ' ' << oneRxn.assocAngles.phi1 << ' ' << oneRxn.assocAngles.phi2 << ' '
                        << oneRxn.assocAngles.omega << '\n';
            restartFile << oneRxn.norm1.x << ' ' << oneRxn.norm1.y << ' ' << oneRxn.norm1.z << '\n';
            restartFile << oneRxn.norm2.x << ' ' << oneRxn.norm2.y << ' ' << oneRxn.norm2.z << '\n';
            restartFile << oneRxn.excludeVolumeBound << '\n';
            restartFile << oneRxn.isCoupled;
            if (oneRxn.isCoupled)
                restartFile << ' ' << oneRxn.coupledRxn.absRxnIndex << ' ' << oneRxn.coupledRxn.relRxnIndex << ' ' << static_cast<std::underlying_type<ReactionType>::type>(oneRxn.coupledRxn.rxnType) << ' ' << oneRxn.coupledRxn.label << ' ' << std::setprecision(20) << oneRxn.coupledRxn.probCoupled;
            restartFile << '\n';

            // integer reactants
            restartFile << oneRxn.intReactantList.size();
            for (const auto& oneReact : oneRxn.intReactantList)
                restartFile << ' ' << oneReact;
            restartFile << '\n';

            // integer products
            restartFile << oneRxn.intProductList.size();
            for (const auto& oneProd : oneRxn.intProductList)
                restartFile << ' ' << oneProd;
            restartFile << '\n';

            // reactant list
            restartFile << oneRxn.reactantListNew.size() << '\n';
            for (const auto& oneReact : oneRxn.reactantListNew) {
                restartFile << oneReact.molTypeIndex << '\n';
                restartFile << oneReact.ifaceName << ' ' << oneReact.absIfaceIndex << ' ' << oneReact.relIfaceIndex
                            << '\n';
                restartFile << oneReact.requiresState << ' ' << oneReact.requiresInteraction << '\n';
            }

            // product list
            restartFile << oneRxn.productListNew.size() << '\n';
            for (const auto& oneProd : oneRxn.productListNew) {
                restartFile << oneProd.molTypeIndex << '\n';
                restartFile << oneProd.ifaceName << ' ' << oneProd.absIfaceIndex << ' ' << oneProd.relIfaceIndex
                            << '\n';
                restartFile << oneProd.requiresState << ' ' << oneProd.requiresInteraction << '\n';
            }

            // rate list
            restartFile << oneRxn.rateList.size() << '\n';
            for (auto& oneRate : oneRxn.rateList) {
                restartFile << std::setprecision(20) << oneRate.rate << '\n';
                restartFile << oneRate.otherIfaceLists.size() << '\n';
                for (const auto& otherIfaceList : oneRate.otherIfaceLists) {
                    restartFile << otherIfaceList.size() << '\n';
                    for (const auto& anccIface : otherIfaceList) {
                        restartFile << anccIface.molTypeIndex << ' ' << anccIface.ifaceName << ' '
                                    << anccIface.absIfaceIndex << ' ' << anccIface.relIfaceIndex << ' '
                                    << anccIface.requiresState << ' ' << anccIface.requiresInteraction << '\n';
                    }
                }
            }
        }

        // back reactions
        for (const auto& oneRxn : backRxns) {
            restartFile << oneRxn.absRxnIndex << ' ' << oneRxn.relRxnIndex << '\n';
            restartFile << static_cast<std::underlying_type<ReactionType>::type>(oneRxn.rxnType) << ' '
                        << oneRxn.isSymmetric << ' ' << oneRxn.isOnMem << ' ' << oneRxn.hasStateChange << '\n';
            restartFile << oneRxn.isObserved << ' ' << oneRxn.observeLabel << '\n';
            restartFile << oneRxn.conjForwardRxnIndex << '\n';
            restartFile << oneRxn.isCoupled;
            if (oneRxn.isCoupled)
                restartFile << ' ' << oneRxn.coupledRxn.absRxnIndex << ' ' << oneRxn.coupledRxn.relRxnIndex << ' ' << static_cast<std::underlying_type<ReactionType>::type>(oneRxn.coupledRxn.rxnType) << ' ' << oneRxn.coupledRxn.label << ' ' << std::setprecision(20) << oneRxn.coupledRxn.probCoupled;
            restartFile << '\n';

            // integer reactants
            restartFile << oneRxn.intReactantList.size();
            for (const auto& oneReact : oneRxn.intReactantList)
                restartFile << ' ' << oneReact;
            restartFile << '\n';

            // integer products
            restartFile << oneRxn.intProductList.size();
            for (const auto& oneProd : oneRxn.intProductList)
                restartFile << ' ' << oneProd;
            restartFile << '\n';

            // reactant list
            restartFile << oneRxn.reactantListNew.size() << '\n';
            for (const auto& oneReact : oneRxn.reactantListNew) {
                restartFile << oneReact.molTypeIndex << '\n';
                restartFile << oneReact.ifaceName << ' ' << oneReact.absIfaceIndex << ' ' << oneReact.relIfaceIndex
                            << '\n';
                restartFile << oneReact.requiresState << ' ' << oneReact.requiresInteraction << '\n';
            }

            // product list
            restartFile << oneRxn.productListNew.size() << '\n';
            for (const auto& oneProd : oneRxn.productListNew) {
                restartFile << oneProd.molTypeIndex << '\n';
                restartFile << oneProd.ifaceName << ' ' << oneProd.absIfaceIndex << ' ' << oneProd.relIfaceIndex
                            << '\n';
                restartFile << oneProd.requiresState << ' ' << oneProd.requiresInteraction << '\n';
            }

            // rate list
            restartFile << oneRxn.rateList.size() << '\n';
            for (const auto& oneRate : oneRxn.rateList) {
                restartFile << std::setprecision(20) << oneRate.rate << '\n';
                restartFile << oneRate.otherIfaceLists.size() << '\n';
                for (const auto& oneList : oneRate.otherIfaceLists) {
                    restartFile << oneList.size() << '\n';
                    for (const auto& anccIface : oneList) {
                        restartFile << anccIface.molTypeIndex << ' ' << anccIface.ifaceName << ' '
                                    << anccIface.absIfaceIndex << ' ' << anccIface.relIfaceIndex << ' '
                                    << anccIface.requiresState << ' ' << anccIface.requiresInteraction << '\n';
                    }
                }
            }
        }

        // creation and destruction reactions
        for (auto& oneRxn : createDestructRxns) {
            restartFile << oneRxn.absRxnIndex << ' ' << oneRxn.relRxnIndex << '\n';
            restartFile << static_cast<std::underlying_type<ReactionType>::type>(oneRxn.rxnType) << ' '
                        << oneRxn.isOnMem << '\n';
            restartFile << oneRxn.isObserved << ' ' << oneRxn.observeLabel << '\n';
            restartFile << oneRxn.creationRadius << '\n';

            // integer reactants
            restartFile << oneRxn.intReactantList.size();
            for (const auto& oneReact : oneRxn.intReactantList)
                restartFile << ' ' << oneReact;
            restartFile << '\n';

            // integer products
            restartFile << oneRxn.intProductList.size();
            for (const auto& oneProd : oneRxn.intProductList)
                restartFile << ' ' << oneProd;
            restartFile << std::endl;

            // reactant list
            restartFile << oneRxn.reactantMolList.size() << '\n';
            for (auto& oneReact : oneRxn.reactantMolList) {
                restartFile << oneReact.molTypeIndex << ' ' << oneReact.molName << ' ' << oneReact.interfaceList.size()
                            << '\n';
                for (auto& oneIface : oneReact.interfaceList) {
                    restartFile << oneIface.molTypeIndex << '\n';
                    restartFile << oneIface.ifaceName << ' ' << oneIface.absIfaceIndex << ' ' << oneIface.relIfaceIndex
                                << '\n';
                    restartFile << oneIface.requiresState << ' ' << oneIface.requiresInteraction << '\n';
                }
            }

            // product list
            restartFile << oneRxn.productMolList.size() << '\n';
            for (auto& oneProd : oneRxn.productMolList) {
                restartFile << oneProd.molTypeIndex << ' ' << oneProd.molName << ' ' << oneProd.interfaceList.size()
                            << '\n';
                for (auto& oneIface : oneProd.interfaceList) {
                    restartFile << oneIface.molTypeIndex << '\n';
                    restartFile << oneIface.ifaceName << ' ' << oneIface.absIfaceIndex << ' ' << oneIface.relIfaceIndex
                                << '\n';
                    restartFile << oneIface.requiresState << ' ' << oneIface.requiresInteraction << '\n';
                }
            }

            restartFile << oneRxn.rateList.size() << '\n';
            for (auto& oneRate : oneRxn.rateList) {
                restartFile << std::setprecision(20) << oneRate.rate << '\n';
                restartFile << oneRate.otherIfaceLists.size() << '\n';
                // if (oneRate.otherIfaceLists.size() != 0) {
                //     for (const auto& anccIface : oneRate.otherIfaceLists[0]) {
                //         restartFile << anccIface.molTypeIndex << ' ' << anccIface.ifaceName << ' '
                //                     << anccIface.absIfaceIndex << ' ' << anccIface.relIfaceIndex << ' '
                //                     << anccIface.requiresState << ' ' << anccIface.requiresInteraction << '\n';
                //     }
                // }
                //  restartFile << oneRate.otherIfaceLists.size() << '\n';
                for (const auto& otherIfaceList : oneRate.otherIfaceLists) {
                    restartFile << otherIfaceList.size() << '\n';
                    for (const auto& anccIface : otherIfaceList) {
                        restartFile << anccIface.molTypeIndex << ' ' << anccIface.ifaceName << ' '
                                    << anccIface.absIfaceIndex << ' ' << anccIface.relIfaceIndex << ' '
                                    << anccIface.requiresState << ' ' << anccIface.requiresInteraction << '\n';
                    }
                }
            }
        }
    }

    // write Molecules
    {
        restartFile << "#All Molecules and coordinates \n";
        restartFile << moleculeList.size() << ' ' << Molecule::numberOfMolecules << '\n';
        for (auto& oneMol : moleculeList) {
            restartFile << oneMol.index << ' ' << oneMol.isEmpty << ' ' << oneMol.myComIndex << ' '
                        << oneMol.molTypeIndex << ' ' << oneMol.mySubVolIndex << '\n';
            restartFile << oneMol.mass << ' ' << oneMol.isLipid << ' ' << oneMol.isImplicitLipid << ' ' << oneMol.linksToSurface << ' ' << oneMol.isEmpty << '\n';

            // center of mass
            restartFile << std::fixed << oneMol.comCoord.x << ' ' << oneMol.comCoord.y << ' ' << oneMol.comCoord.z
                        << '\n';

            // interface lists
            restartFile << oneMol.freelist.size();
            for (const auto& oneIface : oneMol.freelist)
                restartFile << ' ' << oneIface;
            restartFile << '\n';

            restartFile << oneMol.bndlist.size();
            for (const auto& oneIface : oneMol.bndlist)
                restartFile << ' ' << oneIface;
            restartFile << '\n';
            restartFile << oneMol.bndpartner.size();
            for (const auto& oneIface : oneMol.bndpartner)
                restartFile << ' ' << oneIface;
            restartFile << '\n';

            // interfaces
            restartFile << oneMol.interfaceList.size() << '\n';
            for (auto& oneIface : oneMol.interfaceList) {
                restartFile << oneIface.index << ' ' << oneIface.relIndex << ' ' << oneIface.molTypeIndex << ' '
                            << oneIface.stateIndex << ' ' << oneIface.stateIden << ' ' << oneIface.isBound << '\n';
                restartFile << std::fixed << oneIface.coord.x << ' ' << oneIface.coord.y << ' ' << oneIface.coord.z
                            << '\n';

                if (oneIface.isBound) {
                    restartFile << oneIface.interaction.partnerIndex << ' ' << oneIface.interaction.partnerIfaceIndex
                                << ' ' << oneIface.interaction.conjBackRxn << '\n';
                }
            }

            // reweighting lists
            restartFile << oneMol.prevlist.size();
            for (const auto& oneElem : oneMol.prevlist)
                restartFile << ' ' << oneElem;
            restartFile << '\n';
            restartFile << oneMol.prevmyface.size();
            for (const auto& oneElem : oneMol.prevmyface)
                restartFile << ' ' << oneElem;
            restartFile << '\n';
            restartFile << oneMol.prevpface.size();
            for (const auto& oneElem : oneMol.prevpface)
                restartFile << ' ' << oneElem;
            restartFile << '\n';
            restartFile << oneMol.prevnorm.size();
            for (const auto& oneElem : oneMol.prevnorm)
                restartFile << ' ' << oneElem;
            restartFile << '\n';
            restartFile << oneMol.ps_prev.size();
            for (const auto& oneElem : oneMol.ps_prev)
                restartFile << ' ' << oneElem;
            restartFile << '\n';
            restartFile << oneMol.prevsep.size();
            for (const auto& oneElem : oneMol.prevsep)
                restartFile << ' ' << oneElem;
            restartFile << '\n';
        }

        restartFile << Molecule::emptyMolList.size();
        for (auto& index : Molecule::emptyMolList)
            restartFile << ' ' << index;
        restartFile << '\n';
    }

    // write Complexes
    {
        restartFile << "#All Complexes and their components \n";
        restartFile << complexList.size() << ' ' << Complex::numberOfComplexes << '\n';
        for (const auto& oneCom : complexList) {
            restartFile << oneCom.index << ' ' << oneCom.isEmpty << ' ' << oneCom.radius << ' ' << oneCom.mass << '\n';
            restartFile << oneCom.linksToSurface << ' ' << oneCom.iLipidIndex << ' ' << oneCom.OnSurface << '\n';
            restartFile << std::fixed << oneCom.comCoord.x << ' ' << oneCom.comCoord.y << ' ' << oneCom.comCoord.z
                        << '\n';
            restartFile << std::fixed << oneCom.D.x << ' ' << oneCom.D.y << ' ' << oneCom.D.z << '\n';
            restartFile << std::fixed << oneCom.Dr.x << ' ' << oneCom.Dr.y << ' ' << oneCom.Dr.z << '\n';

            // member molecule lists
            restartFile << oneCom.memberList.size();
            for (const auto& memMol : oneCom.memberList)
                restartFile << ' ' << memMol;
            restartFile << '\n';
            restartFile << oneCom.numEachMol.size();
            for (const auto& memMol : oneCom.numEachMol)
                restartFile << ' ' << memMol;
            restartFile << '\n';

            restartFile << oneCom.lastNumberUpdateItrEachMol.size();
            for (const auto& memMol : oneCom.lastNumberUpdateItrEachMol)
                restartFile << ' ' << memMol;
            restartFile << '\n';
        }

        restartFile << Complex::emptyComList.size();
        for (auto& index : Complex::emptyComList)
            restartFile << ' ' << index;
        restartFile << '\n';
    }

    // Write observables
    {
        restartFile << "#Observables \n";
        restartFile << observablesList.size() << '\n';
        for (auto& observable : observablesList)
            restartFile << observable.first << ' ' << observable.second << '\n';
    }

    // Write counterArrays
    {
        restartFile << "#counterArrays.NLoops .nCancels\n";
        restartFile << counterArrays.nLoops << ' ' << counterArrays.nCancelOverlapPartner << ' ' << counterArrays.nCancelOverlapSystem << ' ' << counterArrays.nCancelDisplace2D << ' ' << counterArrays.nCancelDisplace3D << ' ' << counterArrays.nCancelDisplace3Dto2D << ' ' << counterArrays.nCancelSpanBox << ' ' << counterArrays.nAssocSuccess << ' ' << counterArrays.eventArraySize << '\n';
        //write events3D, 3Dto2D, and 2D
        for (auto& event : counterArrays.events3D)
            restartFile << event << ' ';
        restartFile << '\n';
        for (auto& event : counterArrays.events3Dto2D)
            restartFile << event << ' ';
        restartFile << '\n';
        for (auto& event : counterArrays.events2D)
            restartFile << event << ' ';
        restartFile << '\n';
        // Write bindPairList
        restartFile << "#counterArrays.bindPairList \n";
        restartFile << counterArrays.bindPairList.size() << '\n';
        for (auto& bindPair : counterArrays.bindPairList) {
            restartFile << bindPair.size() << '\n';
            for (auto& oneElem : bindPair)
                restartFile << ' ' << oneElem;
            restartFile << '\n';
        }
        // for (auto& bindPair : counterArrays.bindPairListIL2D) {
        //     restartFile << bindPair.size() << '\n';
        //     for (auto& oneElem : bindPair)
        //         restartFile << ' ' << oneElem;
        //     restartFile << '\n';
        // }
        // for (auto& bindPair : counterArrays.bindPairListIL3D) {
        //     restartFile << bindPair.size() << '\n';
        //     for (auto& oneElem : bindPair)
        //         restartFile << ' ' << oneElem;
        //     restartFile << '\n';
        // }
    }
}
