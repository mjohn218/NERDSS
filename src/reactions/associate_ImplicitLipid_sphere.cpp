#include "boundary_conditions/reflect_functions.hpp"
#include "classes/class_Rxns.hpp"
#include "io/io.hpp"
#include "reactions/association/association.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"
#include <cmath>
#include <iomanip>

// ifaceIndex2, reactMol2, reactCom2 are implicit-lipid's
// only works for sphere system, and reactCom1 binds to surface by 3D->2D or 2D->2D.
//
void associate_implicitlipid_sphere(
    int ifaceIndex1, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2,
    Complex& reactCom1, Complex& reactCom2, const Parameters& params,
    ForwardRxn& currRxn, std::vector<Molecule>& moleculeList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList,
    copyCounters& counterArrays, std::vector<Complex>& complexList,
    Membrane& membraneObject, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns)
{
    // TRACE();
    double RS3D { -1.0 };
    for (int RS3Di = 0; RS3Di < 100; RS3Di++) {
        if ((std::abs(membraneObject.RS3Dvect[RS3Di] - currRxn.bindRadius) < 1E-15) && (std::abs(membraneObject.RS3Dvect[RS3Di + 100] - currRxn.rateList[0].rate) < 1E-15) && std::abs(membraneObject.RS3Dvect[RS3Di + 200] - (1.0 / 3.0 * (molTemplateList[currRxn.reactantListNew[0].molTypeIndex].D.x + molTemplateList[currRxn.reactantListNew[1].molTypeIndex].D.x) + 1.0 / 3.0 * (molTemplateList[currRxn.reactantListNew[0].molTypeIndex].D.y + molTemplateList[currRxn.reactantListNew[1].molTypeIndex].D.y) + 1.0 / 3.0 * (molTemplateList[currRxn.reactantListNew[0].molTypeIndex].D.z + molTemplateList[currRxn.reactantListNew[1].molTypeIndex].D.z))) < 1E-15) {
            RS3D = membraneObject.RS3Dvect[RS3Di + 300];
            break;
        }
    }

    // mol2 is implicit-lipid, then we need to set its temporary position according to mol1.
    reactMol2.create_position_implicit_lipid(reactMol1, ifaceIndex2, currRxn.bindRadius, membraneObject);
    //reactCom2.memberList.push_back(reactMol2.index);
    reactCom2.comCoord = reactMol2.comCoord;

    // set up temporary coordinates
    for (auto& memMol : reactCom1.memberList)
        moleculeList[memMol].set_tmp_association_coords();
    for (auto& mol : reactCom2.memberList)
        moleculeList[mol].set_tmp_association_coords();
    // std::cout << " In Associate binidng to sphere surface " << std::endl;
    // write_xyz_assoc_cout(reactCom1, reactCom2, moleculeList);

    // create references to reacting interfaces
    Coord& reactIface1 = reactMol1.tmpICoords[ifaceIndex1];
    Coord& reactIface2 = reactMol2.tmpICoords[ifaceIndex2];

    // orientation corrections for membrane bound components
    bool isOnMembrane = false;
    bool transitionToSurface = false;

    /* MOVE PROTEIN TO SIGMA */
    {
        Vector sigma { reactIface1 - reactIface2 };
        Vector transVec1 { 0.0, 0.0, 0.0 };
        double displaceFrac {};
        double sigmaMag;
        // if both in 2D, ignore the z-component
        if (reactCom1.D.z < 1E-14) {
            isOnMembrane = true;
        } else { // note not on the membrane
            transitionToSurface = true;
            // std::cout << "TRANSITIONING FROM 3D->2D " << std::endl;
            sigmaMag = sqrt((sigma.x * sigma.x) + (sigma.y * sigma.y) + (sigma.z * sigma.z));
            displaceFrac = (sigmaMag - currRxn.bindRadius) / sigmaMag;
            transVec1.x = -sigma.x * displaceFrac;
            transVec1.y = -sigma.y * displaceFrac;
            transVec1.z = -sigma.z * displaceFrac;
        }
        // update the temporary coordinates
        for (auto& mp : reactCom1.memberList)
            moleculeList[mp].update_association_coords(transVec1);
        // std::cout << "Position after pushed to sigma: " << std::endl;
        // reactMol1.display_assoc_icoords("mol1");
        // reactMol2.display_assoc_icoords("mol2");
    } // Move protein to sigma
    // now rotate reactCom1
    // std::cout << "ImplicitLipid model, need to rotate the complex, and make it stick on the sphere.  \n";
    // std::cout << "Pre-MEMBRANE ROT, COORDS: " << std::endl;
    // reactMol1.display_assoc_icoords("mol1");
    if (molTemplateList[reactMol1.molTypeIndex].isPoint) {
        // If both molecules are points, no orientations to specify
        // std::cout << " Move two point particles to contact along current separation vector, NO ORIENTATION \n";
    } else {
        // THETA
        // std::cout << std::setw(8) << std::setfill('-') << ' ' << std::endl
        //           << "THETA 1" << std::endl
        //           << std::setw(8) << ' ' << std::setfill(' ') << std::endl;
        theta_rotation(reactIface1, reactIface2, reactMol1, reactMol2,
            currRxn.assocAngles.theta1, reactCom1, reactCom2, moleculeList);
        // std::cout << std::setw(30) << std::setfill('-') << ' '
        //           << std::setfill(' ') << std::endl;
        // std::cout << "THETA 2" << std::endl
        //           << std::setw(8) << std::setfill('-') << ' ' << std::setfill(' ')
        //           << std::endl;
        theta_rotation(reactIface2, reactIface1, reactMol2, reactMol1,
            currRxn.assocAngles.theta2, reactCom2, reactCom1, moleculeList);

        // OMEGA
        // if protein has theta M_PI, uses protein norm instead of com_iface vector
        // std::cout << std::setw(6) << std::setfill('-') << ' ' << std::endl
        //           << "OMEGA" << std::endl
        //           << std::setw(6) << ' ' << std::setfill(' ') << std::endl;
        if (!std::isnan(currRxn.assocAngles.omega)) {
            omega_rotation(reactIface1, reactIface2, ifaceIndex2, reactMol1,
                reactMol2, reactCom1, reactCom2,
                currRxn.assocAngles.omega, currRxn, moleculeList,
                molTemplateList);
        } //else
        // std::cout << "P1 or P2 is a rod-type protein, no dihedral for associated complex." << std::endl;

        // PHI
        // PHI 1
        // std::cout << std::setw(6) << std::setfill('-') << ' ' << std::endl
        //           << "PHI 1" << std::endl
        //           << std::setw(6) << ' ' << std::setfill(' ') << std::endl;

        if (!std::isnan(currRxn.assocAngles.phi1)) {
            phi_rotation(reactIface1, reactIface2, ifaceIndex2, reactMol1, reactMol2, reactCom1, reactCom2,
                currRxn.norm1, currRxn.assocAngles.phi1, currRxn, moleculeList, molTemplateList);
        } //else
        // std::cout << "P1 has no valid phi angle." << std::endl;

        // PHI 2
        // std::cout << std::setw(6) << std::setfill('-') << ' ' << std::endl
        //           << "PHI 2" << std::endl
        //           << std::setw(6) << ' ' << std::setfill(' ') << std::endl;

        if (!std::isnan(currRxn.assocAngles.phi2)) {
            phi_rotation(reactIface2, reactIface1, ifaceIndex1, reactMol2, reactMol1, reactCom2, reactCom1,
                currRxn.norm2, currRxn.assocAngles.phi2, currRxn, moleculeList, molTemplateList);
        } //else
        // std::cout << "P2 has no valid phi angle." << std::endl;
    } // end of if points.

    // need to adjust Lipid's orientation, make it verticle to sphere
    {
        Molecule memProtein;
        Molecule Lipid;
        set_memProtein_sphere(reactCom2, memProtein, moleculeList, membraneObject);
        find_Lipid_sphere(reactCom2, Lipid, moleculeList, membraneObject);
        Quat memRot = save_mem_orientation(memProtein, Lipid, molTemplateList[Lipid.molTypeIndex]);
        Coord pivot = Lipid.tmpComCoord;
        /*rotate the molecules and their complexes.*/
        rotate(pivot, memRot, reactCom1, moleculeList);
        rotate(pivot, memRot, reactCom2, moleculeList);
        //make sure reactIface1 is sticking on sphere.
        double lamda = (membraneObject.sphereR - reactIface1.get_magnitude()) / reactIface1.get_magnitude();
        Vector dtrans = Vector { lamda * reactIface1 };
        // std::cout << " In ImplicitLipid model, protein interface is off spherical membrane, shift up by: " << dtrans.x << " " << dtrans.y << " " << dtrans.z
        //           << std::endl;
        // update the temporary coordinates for both complexes
        for (auto& mp : reactCom1.memberList)
            moleculeList[mp].update_association_coords(dtrans);
    }

    // std::cout << "Pos-MEMBRANE ROT, COORDS: " << std::endl;
    // reactMol1.display_assoc_icoords("mol1");

    // std::cout << "OVERLAP CHECK  AND REFLECT OFF SPHERE: " << std::endl;
    // Reflect off the sphere.
    std::array<double, 3> traj;
    for (int mm = 0; mm < 3; mm++)
        traj[mm] = 0;

    /*This needs to evaluate the traj update, based on it initially being zero.
    And here, it should be called based on the tmpCoords, not the full
    coordinates. also requires updating the COM of this temporary new position
  */
    update_complex_tmp_com_crds(reactCom1, moleculeList);
    reactCom1.tmpOnSurface = true;
    reflect_traj_tmp_crds(params, moleculeList, reactCom1, traj, membraneObject, RS3D); // uses tmpCoords to calculate traj.

    if (std::abs(traj[0] + traj[1] + traj[2]) > 1E-14) {
        // update the temporary coordinates for both complexes
        Vector vtraj { traj[0], traj[1], traj[2] };
        for (auto& mp : reactCom1.memberList)
            moleculeList[mp].update_association_coords(vtraj);

        // std::cout << "CRDS after reflecting off of the SPHERE by " << traj[0] << ' '
        //           << traj[1] << ' ' << traj[2] << std::endl;
        // write_xyz_assoc_cout(reactCom1, reactCom2, moleculeList);
    }
    /* CHECKS */
    bool cancelAssoc { false };

    check_if_spans_sphere(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);
    if (cancelAssoc == true)
        counterArrays.nCancelSpanBox++;
    if (cancelAssoc == false) {
        check_for_structure_overlap_system(cancelAssoc, reactCom1, reactCom2, moleculeList, params, molTemplateList, complexList, forwardRxns, backRxns);
        if (cancelAssoc == true)
            counterArrays.nCancelOverlapSystem++;
    }
    if (cancelAssoc == false) {
        measure_complex_displacement(cancelAssoc, reactCom1, reactCom2, moleculeList, params, molTemplateList, complexList);
        if (cancelAssoc == true) {
            if (isOnMembrane)
                counterArrays.nCancelDisplace2D++;
            else if (transitionToSurface)
                counterArrays.nCancelDisplace3Dto2D++;
            else
                counterArrays.nCancelDisplace3D++;
        }
    }

    if (cancelAssoc == true) {
        // std::cout << "Canceling association, returning complexes to original state.\n";
        reactCom1.tmpOnSurface = false;
        for (auto memMol : reactCom1.memberList)
            moleculeList[memMol].clear_tmp_association_coords();
        for (auto memMol : reactCom2.memberList)
            moleculeList[memMol].clear_tmp_association_coords();
        return;
        // EXIT FROM THIS ROUTINE.
    }
    counterArrays.nAssocSuccess++; //keep track of total number of successful association moves.
    /*Keep track of the sizes of complexes that associated*/
    track_association_events(reactCom1, reactCom2, transitionToSurface, isOnMembrane, counterArrays);

    //if (reactMol2.isImplicitLipid == true)
    {
        // write temporary to real coords and clear temporary coordinates
        for (auto memMol : reactCom1.memberList) {
            moleculeList[memMol].comCoord = moleculeList[memMol].tmpComCoord;
            for (unsigned int i { 0 }; i < moleculeList[memMol].interfaceList.size(); ++i)
                moleculeList[memMol].interfaceList[i].coord = moleculeList[memMol].tmpICoords[i];

            moleculeList[memMol].clear_tmp_association_coords();
            moleculeList[memMol].trajStatus = TrajStatus::propagated;
        }
        // add link to surface for the molecule
        moleculeList[reactMol1.index].linksToSurface++;
        // update complexes
        reactCom1.iLipidIndex = reactMol2.index; // index in molecule list of the implicit lipid.
        reactCom1.linksToSurface++;
        reactCom1.update_properties(moleculeList, molTemplateList); // recalculate the properties of the first complex
        reactCom1.D.z = 0.0; // maybe redundant, since we can use linksToSurface > 0
        reactCom1.OnSurface = true; // maybe redundant, since we can use linksToSurface > 0
        /*Clear implicit lipid properties*/
        moleculeList[reactMol2.index].clear_tmp_association_coords();
        reactCom2.memberList.clear();
        reactCom2.memberList.push_back(reactMol2.index);
        // std::cout << " Size of IL's complex:" << reactCom2.memberList.size()
        //           << " Interfaces on IL: "
        //           << moleculeList[reactMol2.index].interfaceList.size()
        //           << std::endl;
        reactCom2.update_properties(moleculeList, molTemplateList); // recalculate the properties of the second complex

        // Enforce boundary conditions
        // For the sphere system, many times of reflections may need to move the complex back inside the sphere!!
        reflect_complex_rad_rot(membraneObject, reactCom1, moleculeList, RS3D);
        //------------------------START UPDATE MONOMERLIST-------------------------
        // update oneTemp.monomerList when oneTemp.canDestroy is true and mol is monomer
        // reactMol1
        {
            Molecule& oneMol { reactMol1 };
            MolTemplate& oneTemp { molTemplateList[oneMol.molTypeIndex] };
            bool isMonomer { oneMol.bndpartner.empty() };
            bool canDestroy { oneTemp.canDestroy };
            if (isMonomer && canDestroy) {
                //remove from monomerList
                std::vector<int>& oneList { oneTemp.monomerList };
                std::vector<int>::iterator result { std::find(std::begin(oneList), std::end(oneList), oneMol.index) };
                if (result != std::end(oneList)) {
                    oneList.erase(result);
                }
            }
        }
        //------------------------END UPDATE MONOMERLIST---------------------------
        // update Molecule interface statuses
        reactMol1.interfaceList[ifaceIndex1].interaction.partnerIndex = reactMol2.index;
        reactMol1.interfaceList[ifaceIndex1].interaction.partnerIfaceIndex = ifaceIndex2;
        if (currRxn.isReversible) {
            reactMol1.interfaceList[ifaceIndex1].interaction.conjBackRxn = currRxn.conjBackRxnIndex;
        }

        reactMol1.interfaceList[ifaceIndex1].isBound = true;
        reactMol1.interfaceList[ifaceIndex1].index = currRxn.productListNew[0].absIfaceIndex;

        // add to the list of bound interfaces and remove from the list of free interfaces
        reactMol1.bndlist.push_back(ifaceIndex1);
        reactMol1.bndpartner.push_back(reactMol2.index);
        {
            size_t tmpItr { 0 };
            for (unsigned i { 0 }; i < reactMol1.freelist.size(); ++i) {
                if (reactMol1.freelist[i] == ifaceIndex1)
                    tmpItr = i;
            }
            reactMol1.freelist[tmpItr] = reactMol1.freelist.back();
            reactMol1.freelist.pop_back();
        }

        // Set probability of this protein to zero in all reactions so it doesn't
        // try to react again but the partners still will avoid overlapping.
        for (unsigned crossItr { 0 }; crossItr < reactMol1.crossbase.size(); ++crossItr) {
            int skipMol { reactMol1.crossbase[crossItr] };
            if (!moleculeList[skipMol].isImplicitLipid) {
                for (unsigned crossItr2 { 0 }; crossItr2 < moleculeList[skipMol].crossbase.size(); ++crossItr2) {
                    if (moleculeList[skipMol].crossbase[crossItr2] == reactMol1.index)
                        moleculeList[skipMol].probvec[crossItr2] = 0;
                }
            }
        }

        // Update the crossed molecule lists so that the current molecules won't
        // avoid anything, but others will.
        reactCom1.ncross = -1;
        reactMol1.crossbase.clear();
        // mol2 is the implicit lipid.
    }

    // update the number of bound species
    update_Nboundpairs(reactMol1.molTypeIndex, reactMol2.molTypeIndex, 1, params,
        counterArrays);

    //Update species copy numbers
    counterArrays.copyNumSpecies[currRxn.reactantListNew[0].absIfaceIndex] -= 1; // decrement ifaceIndex1
    counterArrays.copyNumSpecies[currRxn.reactantListNew[1].absIfaceIndex] -= 1; // decrement ifaceIndex2
    //counterArrays.copyNumSpecies[reactMol1.interfaceList[ifaceIndex1].index] -= 1; // decrement ifaceIndex1
    //counterArrays.copyNumSpecies[reactMol2.interfaceList[ifaceIndex2].index] -= 1; // decrement ifaceIndex2
    counterArrays.copyNumSpecies[currRxn.productListNew[0].absIfaceIndex] += 1; // increment product state

    //-------------------------UPDATE BINDPAIRLIST----------------------------
    // auto molIndex { reactMol1.index };
    // auto comIndex { reactMol1.myComIndex };
    // if (molTemplateList[reactMol1.molTypeIndex].isImplicitLipid == true) {
    //     molIndex = reactMol2.index;
    //     comIndex = reactMol2.myComIndex;
    // }

    // if (counterArrays.canDissociate[currRxn.productListNew[0].absIfaceIndex]) {
    //     if (complexList[comIndex].linksToSurface == 1) {
    //         counterArrays.bindPairListIL3D[currRxn.productListNew[0].absIfaceIndex].emplace_back(molIndex);
    //     } else {
    //         counterArrays.bindPairListIL2D[currRxn.productListNew[0].absIfaceIndex].emplace_back(molIndex);
    //     }
    // }

    // if (complexList[comIndex].linksToSurface == 2) {
    //     // one interface in the complex bind to surface changing the linksToSurface of complex from one to two will impact the bindPairListIL of the other linking interface
    //     // need to remove the impacted interface's mol index from the bindPairListIL3D, then add it to bindPairListIL2D
    //     // try to find out the other linking interface
    //     int theOtherLinkingIfaceAbsIndex { -1 };
    //     if (moleculeList[molIndex].linksToSurface == 2) { // the other interface is in the same molecule
    //         // loop over the interfaceList to find the other interface
    //         for (auto oneIface : moleculeList[molIndex].interfaceList) {
    //             if (oneIface.interaction.partnerIndex == 0 && oneIface.index != currRxn.productListNew[0].absIfaceIndex) {
    //                 theOtherLinkingIfaceAbsIndex = oneIface.index;
    //                 if (counterArrays.canDissociate[theOtherLinkingIfaceAbsIndex]) {
    //                     counterArrays.bindPairListIL3D[theOtherLinkingIfaceAbsIndex].erase(std::find_if(counterArrays.bindPairListIL3D[theOtherLinkingIfaceAbsIndex].begin(), counterArrays.bindPairListIL3D[theOtherLinkingIfaceAbsIndex].end(), [&](const size_t& mol) { return mol == molIndex; }));
    //                     counterArrays.bindPairListIL2D[theOtherLinkingIfaceAbsIndex].emplace_back(molIndex);
    //                 }
    //             }
    //         }
    //     } else { // the other interface is in another molecule
    //         for (auto oneMember : complexList[comIndex].memberList) {
    //             if ((oneMember != molIndex) && (moleculeList[oneMember].linksToSurface == 1)) {
    //                 // loop over the interfaceList to find the other interface
    //                 for (auto oneIface : moleculeList[oneMember].interfaceList) {
    //                     if (oneIface.interaction.partnerIndex == 0) {
    //                         theOtherLinkingIfaceAbsIndex = oneIface.index;
    //                         if (counterArrays.canDissociate[theOtherLinkingIfaceAbsIndex]) {
    //                             counterArrays.bindPairListIL3D[theOtherLinkingIfaceAbsIndex].erase(std::find_if(counterArrays.bindPairListIL3D[theOtherLinkingIfaceAbsIndex].begin(), counterArrays.bindPairListIL3D[theOtherLinkingIfaceAbsIndex].end(), [&](const size_t& mol) { return mol == oneMember; }));
    //                             counterArrays.bindPairListIL2D[theOtherLinkingIfaceAbsIndex].emplace_back(oneMember);
    //                         }
    //                     }
    //                 }
    //                 break;
    //             }
    //         }
    //     }
    // }
    //-----------------------END UPDATE BINDPAIRLIST--------------------------

    //Update free number of lipids of IL for each state
    RxnIface implicitLipidState {};
    const auto& implicitLipidStateList = molTemplateList[moleculeList[membraneObject.implicitlipidIndex].molTypeIndex].interfaceList[0].stateList;
    if (molTemplateList[currRxn.reactantListNew[1].molTypeIndex].isImplicitLipid == true) {
        implicitLipidState = currRxn.reactantListNew[1];
    } else {
        implicitLipidState = currRxn.reactantListNew[0];
    }
    int relStateIndex { -1 };
    for (auto& state : implicitLipidStateList) {
        if (state.index == implicitLipidState.absIfaceIndex) {
            relStateIndex = static_cast<int>(&state - &implicitLipidStateList[0]);
            break;
        }
    }
    membraneObject.numberOfFreeLipidsEachState[relStateIndex] -= 1;

    // TODO: Insert species tracking here
    if (currRxn.isObserved) {
        auto obsItr = observablesList.find(currRxn.observeLabel);
        if (obsItr != observablesList.end())
            ++obsItr->second;
    }
    //    reactCom1.display();
}
