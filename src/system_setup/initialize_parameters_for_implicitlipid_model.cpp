#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
void initialize_paramters_for_implicitlipid_model(int& implicitlipidIndex, const Parameters& params, std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
    std::vector<Molecule>& moleculeList, std::vector<MolTemplate>& molTemplateList, std::vector<Complex>& complexList,
    Membrane& membraneObject)
{
    // check if it is a 2D simulation
    bool is2D = true;
    for (int mol { 0 }; mol < moleculeList.size(); ++mol) {
        int molTypeIndex = moleculeList[mol].molTypeIndex;
        if (molTemplateList[molTypeIndex].isImplicitLipid == true) {
            membraneObject.implicitLipid = true; // determine membraneObject.implicitLipid according to molTemp.isImplicitLipid
            membraneObject.lipidLength = molTemplateList[molTypeIndex].radius;
            continue;
        }
        if (molTemplateList[molTypeIndex].D.z > 0)
            is2D = false;
    }
    if (is2D || membraneObject.waterBox.z == 0)
        membraneObject.TwoD = true;

    // find the mol.index of the implicit-lipid, which should be implicitlipidIndex=0;
    bool systemIL = false; //Are there any IL's in this system?
    for (int mol { 0 }; mol < moleculeList.size(); ++mol) {
        int molTypeIndex = moleculeList[mol].molTypeIndex;
        if (molTemplateList[molTypeIndex].isImplicitLipid == true) {
            moleculeList[mol].isImplicitLipid = true;
            implicitlipidIndex = mol;
            membraneObject.implicitlipidIndex = mol;
            //moleculeList[mol].mass = 1; // such a large value that implicit-lipid won't move when proteins bind to it.
            systemIL = true; //yes, there are IL's in the system.
        } else {
            moleculeList[mol].isImplicitLipid = false;
        }
    }
    //    membraneObject.implicitlipidIndex = implicitlipidIndex;
    /////////////////////////////////////////////////////////////////////////////////////////////
    // set up 'bindToSurface' in molTemplateList.
    for (const auto& oneRxn : forwardRxns) {
        if (oneRxn.rxnType == ReactionType::bimolecular || oneRxn.rxnType == ReactionType::biMolStateChange) {
            int molTypeIndex = oneRxn.reactantListNew[0].molTypeIndex;
            if (molTemplateList[molTypeIndex].isImplicitLipid == true) {
                int newProteinIndex = oneRxn.reactantListNew[1].molTypeIndex;
                molTemplateList[newProteinIndex].bindToSurface = true;
            }
            molTypeIndex = oneRxn.reactantListNew[1].molTypeIndex;
            if (molTemplateList[molTypeIndex].isImplicitLipid == true) {
                int newProteinIndex = oneRxn.reactantListNew[0].molTypeIndex;
                molTemplateList[newProteinIndex].bindToSurface = true;
            }
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////////////
    if (membraneObject.isSphere) {
        membraneObject.totalSA = 4.0 * M_PI * pow(membraneObject.sphereR, 2.0);
    } else {
        membraneObject.totalSA = membraneObject.waterBox.x * membraneObject.waterBox.y;
    }

    // inititalize RS3Dvect with -1
    for (int i = 0; i < 500; i++) {
        membraneObject.RS3Dvect.push_back(-1.0);
    }

    // initialize several parameters in 'membraneObject' for implicit-lipid model
    if (systemIL == true) {
        // initial number of free lipids' interface
        membraneObject.nSites = 0;
        for (auto onetem : molTemplateList) {
            if (onetem.isImplicitLipid)
                membraneObject.nSites = onetem.copies;
            // here, we assume each implicit-lipid has only one interface.
        }

        // initial numberOfFreeLipidsEachState for non-restart sim
        if (params.fromRestart == false) {
            for (auto& iface : molTemplateList[moleculeList[membraneObject.implicitlipidIndex].molTypeIndex].interfaceList) {
                for (auto& state : iface.stateList) {
                    if (&state - &iface.stateList[0] == 0) { // if is the first state, free lipids is initial copies of IL
                        membraneObject.numberOfFreeLipidsEachState[0] = membraneObject.nSites;
                    } else { //others are zero
                    }
                }
            }
        }

        ///////////////////////////////////////////////////////////
        // initial number of proteins' interface that can bind to implicit-lipids
        const auto& implicitLipidStateList = molTemplateList[moleculeList[membraneObject.implicitlipidIndex].molTypeIndex].interfaceList[0].stateList;
        for (int mol { 0 }; mol < moleculeList.size(); ++mol) {
            if (moleculeList[mol].isImplicitLipid)
                continue;
            int molType = moleculeList[mol].molTypeIndex;
            for (int relfaceItr { 0 }; relfaceItr < moleculeList[mol].interfaceList.size(); ++relfaceItr) {
                int stateIndex = moleculeList[mol].interfaceList[relfaceItr].stateIndex;
                const Interface::State& state = molTemplateList[molType].interfaceList[relfaceItr].stateList[stateIndex];
                for (auto rxnItr : state.myForwardRxns) {
                    const ForwardRxn& oneRxn = forwardRxns[rxnItr];
                    for (int reactItr { 0 }; reactItr < oneRxn.reactantListNew.size(); ++reactItr) {
                        for (auto& implicitLipidState : implicitLipidStateList) {
                            if (implicitLipidState.index
                                == oneRxn.reactantListNew[reactItr].absIfaceIndex)
                                membraneObject.numberOfProteinEachState[static_cast<int>(&implicitLipidState - &implicitLipidStateList[0])]++;
                        }
                    }
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////
        //initialize the RS3Dvect according to sigma, ka, and D;we accord RS3D in Rs3Dvect[i+300] and the reactmol.typeindex in RS3Dvect[i+400]
        int RS3Dindex { 0 };
        //const auto& implicitLipidStateList = molTemplateList[moleculeList[membraneObject.implicitlipidIndex].molTypeIndex].interfaceList[0].stateList;
        for (const auto& oneRxn : forwardRxns) {
            for (const auto& oneMol : oneRxn.reactantListNew) {
                //for (auto& implicitLipidState : implicitLipidStateList) {
                if (oneRxn.reactantListNew.size() > 1) {
                    if (moleculeList[implicitlipidIndex].interfaceList[0].index == oneMol.absIfaceIndex) {
                        membraneObject.RS3Dvect[RS3Dindex] = oneRxn.bindRadius;
                        membraneObject.RS3Dvect[RS3Dindex + 100] = oneRxn.rateList[0].rate;
                        membraneObject.RS3Dvect[RS3Dindex + 200] = 1.0 / 3.0 * (molTemplateList[oneRxn.reactantListNew[0].molTypeIndex].D.x + molTemplateList[oneRxn.reactantListNew[1].molTypeIndex].D.x)
                            + 1.0 / 3.0 * (molTemplateList[oneRxn.reactantListNew[0].molTypeIndex].D.y + molTemplateList[oneRxn.reactantListNew[1].molTypeIndex].D.y)
                            + 1.0 / 3.0 * (molTemplateList[oneRxn.reactantListNew[0].molTypeIndex].D.z + molTemplateList[oneRxn.reactantListNew[1].molTypeIndex].D.z);
                        membraneObject.RS3Dvect[RS3Dindex + 300] = membraneObject.RS3Dvect[RS3Dindex] * membraneObject.RS3Dvect[RS3Dindex + 100] * 2.0 / (membraneObject.RS3Dvect[RS3Dindex + 100] * 2.0 + 4 * M_PI * membraneObject.RS3Dvect[RS3Dindex] * membraneObject.RS3Dvect[RS3Dindex + 200]);
                        if (molTemplateList[oneRxn.reactantListNew[0].molTypeIndex].isImplicitLipid == true)
                            membraneObject.RS3Dvect[RS3Dindex + 400] = oneRxn.reactantListNew[1].molTypeIndex;
                        else
                            membraneObject.RS3Dvect[RS3Dindex + 400] = oneRxn.reactantListNew[0].molTypeIndex;
                        RS3Dindex += 1;
                    }
                }
            }
        }
    } //only set up IL model if they exist.
    // std::cout << "membraneObject.RS3Dvect[300]: " << std::setprecision(20) << membraneObject.RS3Dvect[300] << std::endl;
    // exit(1);
}
