#include "boundary_conditions/reflect_functions.hpp"
#include "classes/class_Rxns.hpp"
#include "io/io.hpp"
#include "reactions/association/association.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "tracing.hpp"

#include <cmath>
#include <iomanip>

void perform_implicitlipid_state_change_box(int stateChangeIface, int facilitatorIface, std::array<int, 3>& rxnItr,
    Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom, Complex& facilitatorCom,
    copyCounters& counterArrays, const Parameters& params, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, Membrane& membraneObject)
{
    // TRACE();
    const auto& stateList = molTemplateList[stateChangeMol.molTypeIndex].interfaceList[stateChangeIface].stateList;
    double bindRadius { 0 };
    bool isStateChangeBackRxn { static_cast<bool>(rxnItr[2]) };
    int rxnIndex { rxnItr[0] };
    int forwardRxnIndex { rxnIndex };
    RxnIface newState {};
    ForwardRxn::Angles assocAngles {};
    if (!isStateChangeBackRxn) {
        bindRadius = forwardRxns[rxnIndex].bindRadius;
        newState = forwardRxns[rxnIndex].productListNew[1];
        assocAngles = forwardRxns[rxnIndex].assocAngles;
    } else {
        forwardRxnIndex = backRxns[rxnIndex].conjForwardRxnIndex;
        bindRadius = forwardRxns[forwardRxnIndex].bindRadius;
        newState = backRxns[rxnIndex].productListNew[1];
        assocAngles = forwardRxns[forwardRxnIndex].assocAngles;
    }
    const ForwardRxn& currRxn = forwardRxns[forwardRxnIndex];
    double RS3D { -1.0 };
    for (int RS3Di = 0; RS3Di < 100; RS3Di++) {
        if ((std::abs(membraneObject.RS3Dvect[RS3Di] - currRxn.bindRadius) < 1E-15) && (std::abs(membraneObject.RS3Dvect[RS3Di + 100] - currRxn.rateList[0].rate) < 1E-15) && std::abs(membraneObject.RS3Dvect[RS3Di + 200] - (1.0 / 3.0 * (molTemplateList[currRxn.reactantListNew[0].molTypeIndex].D.x + molTemplateList[currRxn.reactantListNew[1].molTypeIndex].D.x) + 1.0 / 3.0 * (molTemplateList[currRxn.reactantListNew[0].molTypeIndex].D.y + molTemplateList[currRxn.reactantListNew[1].molTypeIndex].D.y) + 1.0 / 3.0 * (molTemplateList[currRxn.reactantListNew[0].molTypeIndex].D.z + molTemplateList[currRxn.reactantListNew[1].molTypeIndex].D.z))) < 1E-15) {
            RS3D = membraneObject.RS3Dvect[RS3Di + 300];
            break;
        }
    }
    // stateChangeMol is implicit-lipid, then we need to set its temporary position according to facilitatorMol.
    // Coord displace = stateChangeMol.interfaceList[stateChangeIface].coord - stateChangeMol.comCoord;
    // double shift = 0.1; //do not put right underneath, so that sigma starts at non-zero.
    // stateChangeMol.comCoord.x = facilitatorMol.comCoord.x + shift;
    // stateChangeMol.comCoord.y = facilitatorMol.comCoord.y - shift;
    // stateChangeMol.comCoord.z = -membraneObject.waterBox.z / 2.0 + RS3D;

    // stateChangeMol.interfaceList[stateChangeIface].coord = stateChangeMol.comCoord + displace;

    // stateChangeCom.comCoord = stateChangeMol.comCoord;

    stateChangeMol.create_position_implicit_lipid(facilitatorMol, stateChangeIface, currRxn.bindRadius, membraneObject);
    stateChangeCom.comCoord = stateChangeMol.comCoord;

    int relStateIndex { -1 };
    for (auto& state : stateList) {
        if (state.index == newState.absIfaceIndex) {
            relStateIndex = static_cast<int>(&state - &stateList[0]);
            break;
        }
    }
    //std::cout << facilitatorCom.memberList.size() << std::endl;
    for (auto& memMol : facilitatorCom.memberList) {
        moleculeList[memMol].set_tmp_association_coords();
    }

    for (auto& mol : stateChangeCom.memberList) {
        moleculeList[mol].set_tmp_association_coords();
    }

    // create references to reacting interfaces
    Coord& reactIface1 = facilitatorMol.tmpICoords[facilitatorIface];
    Coord& reactIface2 = stateChangeMol.tmpICoords[stateChangeIface];
    Vector sigma { reactIface1 - reactIface2 };
    //    std::cout <<" ORIGINAL CRDS: "<<std::endl;
    //write_xyz_assoc_cout( stateChangeCom, facilitatorCom, moleculeList);

    /*Calculate COM of the two complexes pre-association. The COM of the new complex after should be close to this                              Here, we will force it back, as rotation can cause large displacements*/
    Coord startCOM; //=new double[3];

    com_of_two_tmp_complexes(facilitatorCom, stateChangeCom, startCOM, moleculeList); //com of c1+c2 (original coordinates).
    //    std::cout <<"INITIAL COMPLEX PAIR COM: "<<startCOM.x<<' '<<startCOM.y<<' '<<startCOM.z<<std::endl;
    //orientation corrections for membrane bound components
    // bool isOnMembrane = false;
    // bool transitionToSurface = false;
    // Molecule memProtein;
    // double tol = 1E-14;
    // int slowPro;
    // /* MOVE PROTEIN (ENZYME) TO MEMBRANE, and place implicitlipid at sigma */
    // {
    //     double DxSum { facilitatorCom.D.x + stateChangeCom.D.x };
    //     double DySum { facilitatorCom.D.y + stateChangeCom.D.y };
    //     double DzSum { facilitatorCom.D.z + stateChangeCom.D.z };
    //     Vector sigma { reactIface1 - reactIface2 };
    //     //        std::cout<<" sigma: "<<sigma.x<<' '<<sigma.y<<' '<<sigma.z<<std::endl;
    //     //std::cout <<" Dsum and components: "<<DxSum<<' '<<DySum<<' '<<DzSum<<std::endl;

    //     Vector transVec1 {};
    //     Vector transVec2 {};
    //     double displaceFrac {};
    //     // if both in 2D, ignore the z-component
    //     if (DzSum < 1E-14) {
    //         isOnMembrane = true;
    //         /*Store coordinates of one protein to recover membrane-bound orientation*/

    //         if (stateChangeCom.D.x < facilitatorCom.D.x) {
    //             slowPro = stateChangeMol.index;
    //             memProtein = stateChangeMol; //rotate relative to the slower protein.
    //         } else {
    //             slowPro = facilitatorMol.index;
    //             memProtein = facilitatorMol;
    //         }
    //         DzSum = 1; // to prevent divide by 0
    //         if (std::abs(std::abs(sigma.z) - bindRadius) < 1E-3) {
    //             // if entirety of sigma is in z-component, ignore x and y
    //             displaceFrac = 1;
    //         } else {
    //             double sigmaMag = sqrt((sigma.x * sigma.x) + (sigma.y * sigma.y));
    //             displaceFrac = (sigmaMag - bindRadius) / sigmaMag;
    //         }
    //     } else {
    //         sigma.calc_magnitude();
    //         displaceFrac = (sigma.magnitude - bindRadius) / sigma.magnitude;
    //         if (stateChangeCom.D.z < tol || facilitatorCom.D.z < tol) {
    //             transitionToSurface = true; //both can't be less than tol, or would not be in this loop.
    //             std::cout << "TRANSITIONING FROM 3D->2D " << std::endl;
    //         }
    //     }

    //     transVec1.x = -sigma.x * (facilitatorCom.D.x / DxSum) * displaceFrac;
    //     transVec1.y = -sigma.y * (facilitatorCom.D.y / DySum) * displaceFrac;
    //     transVec1.z = -sigma.z * (facilitatorCom.D.z / DzSum) * displaceFrac;

    //     transVec2.x = sigma.x * (stateChangeCom.D.x / DxSum) * displaceFrac;
    //     transVec2.y = sigma.y * (stateChangeCom.D.y / DySum) * displaceFrac;
    //     transVec2.z = sigma.z * (stateChangeCom.D.z / DzSum) * displaceFrac;
    //     //std::cout<<" translation1 to sigma: "<<transVec1.x<<' '<<transVec1.y<<' '<<transVec1.z<<std::endl;
    //     // update the temporary coordinates
    //     for (auto& mp : facilitatorCom.memberList)
    //         moleculeList[mp].update_association_coords(transVec1);
    //     for (auto& mp : stateChangeCom.memberList)
    //         moleculeList[mp].update_association_coords(transVec2);
    // } //end moving to membrane

    bool isOnMembrane = false;
    // orientation corrections for membrane bound components
    if (facilitatorCom.D.z < 1E-15 && stateChangeCom.D.z < 1E-15) {
        isOnMembrane = true;
        // std::cout << " IMPLICIT LIPID 2D ! Value of isOnMembrane: " << isOnMembrane << std::endl;
    }
    bool transitionToSurface = false;
    if (isOnMembrane == false) {
        transitionToSurface = true;
    }
    Molecule memProtein;

    int slowPro = stateChangeMol.index;

    double tol = 1E-14;

    // MOVE PROTEIN TO THE MEMBRANE, and place IL at sigma.
    {
        sigma.calc_magnitude();
        double sigmaMag = sigma.magnitude;
        Vector transVec1 = Vector(-(sigmaMag - currRxn.bindRadius) / sigmaMag * sigma);
        Vector transVec2 { 0.0, 0.0, 0.0 };
        // std::cout << " translation1 to sigma: " << transVec1.x << ' ' << transVec1.y
        //           << ' ' << transVec1.z << std::endl;

        // update the temporary coordinates
        for (auto& mp : facilitatorCom.memberList)
            moleculeList[mp].update_association_coords(transVec1);
        for (auto& mp : stateChangeCom.memberList)
            moleculeList[mp].update_association_coords(transVec2);
    }

    //    std::cout <<" CRDS TO SIGMA: "<<std::endl;
    //write_xyz_assoc_cout( stateChangeCom, facilitatorCom, moleculeList);
    // TODO: do something else if only one is a point
    if (molTemplateList[facilitatorMol.molTypeIndex].isPoint) {
        /*If both molecules are points, no orientations to specify*/
        // std::cout << " Move two point particles to contact along current separation vector, NO ORIENTATION \n";
    } else {

        /* THETA */
        // std::cout << std::setw(8) << std::setfill('-') << ' ' << std::endl
        //           << "THETA 1" << std::endl
        //           << std::setw(8) << ' ' << std::setfill(' ') << std::endl;
        theta_rotation(reactIface1, reactIface2, facilitatorMol, stateChangeMol, assocAngles.theta1, facilitatorCom,
            stateChangeCom, moleculeList);

        //        write_xyz_assoc_cout( stateChangeCom, facilitatorCom, moleculeList);

        // std::cout << std::setw(30) << std::setfill('-') << ' ' << std::setfill(' ') << std::endl;
        // std::cout << "THETA 2" << std::endl
        //           << std::setw(8) << std::setfill('-') << ' ' << std::setfill(' ') << std::endl;
        theta_rotation(reactIface2, reactIface1, stateChangeMol, facilitatorMol, assocAngles.theta2, stateChangeCom,
            facilitatorCom, moleculeList);

        //        write_xyz_assoc_cout( stateChangeCom, facilitatorCom, moleculeList);

        /* OMEGA */
        // if protein has theta M_PI, uses protein norm instead of com_iface vector
        // std::cout << std::setw(6) << std::setfill('-') << ' ' << std::endl
        //           << "OMEGA" << std::endl
        //           << std::setw(6) << ' ' << std::setfill(' ') << std::endl;
        if (!std::isnan(currRxn.assocAngles.omega)) {
            omega_rotation(reactIface1, reactIface2, stateChangeIface, facilitatorMol, stateChangeMol, facilitatorCom,
                stateChangeCom, assocAngles.omega, currRxn, moleculeList, molTemplateList);

            //            write_xyz_assoc_cout( stateChangeCom, facilitatorCom, moleculeList);

        } //else
        // std::cout << "P1 or P2 is a rod-type protein, no dihedral for associated complex." << std::endl;

        /* PHI */
        // PHI 1
        // std::cout << std::setw(6) << std::setfill('-') << ' ' << std::endl
        //           << "PHI 1" << std::endl
        //           << std::setw(6) << ' ' << std::setfill(' ') << std::endl;

        if (!std::isnan(assocAngles.phi1)) {
            phi_rotation(reactIface1, reactIface2, stateChangeIface, facilitatorMol, stateChangeMol, facilitatorCom,
                stateChangeCom, currRxn.norm1, assocAngles.phi1, currRxn, moleculeList, molTemplateList);

            //            write_xyz_assoc_cout( stateChangeCom, facilitatorCom, moleculeList);

        } //else
        // std::cout << "P1 has no valid phi angle." << std::endl;

        // PHI 2
        // std::cout << std::setw(6) << std::setfill('-') << ' ' << std::endl
        //           << "PHI 2" << std::endl
        //           << std::setw(6) << ' ' << std::setfill(' ') << std::endl;

        if (!std::isnan(assocAngles.phi2)) {
            phi_rotation(reactIface2, reactIface1, facilitatorIface, stateChangeMol, facilitatorMol, stateChangeCom,
                facilitatorCom, currRxn.norm2, assocAngles.phi2, currRxn, moleculeList, molTemplateList);

            //           write_xyz_assoc_cout( stateChangeCom, facilitatorCom, moleculeList);

        } //else
        // std::cout << "P2 has no valid phi angle." << std::endl;
    } //only rotate if they are not both points.

    /*FINISHED ROTATING, NO CONSTRAINTS APPLIED TO SURFACE REACTIONS*/
    Coord finalCOM; //=new double[3];
    com_of_two_tmp_complexes(facilitatorCom, stateChangeCom, finalCOM, moleculeList); //com of c1+c2 (final (tmp) coordinates).
    //   std::cout <<"Pre-MEMBRANE ROT: COMPLEX PAIR COM: "<<finalCOM.x<<' '<<finalCOM.y<<' '<<finalCOM.z<<std::endl;
    {

        /*return orientation of normal back to starting position*/
        // std::cout << " IS ON MEMBRANE, CORRECT ORIENTATION ! " << std::endl;
        Quat memRot;
        Coord pivot;

        if (slowPro == facilitatorMol.index) {
            //	facilitatorMol.display_assoc_icoords("CURRORIENTATION_TOROTATE1");
            memRot = save_mem_orientation(memProtein, facilitatorMol, molTemplateList[facilitatorMol.molTypeIndex]);
            pivot = facilitatorMol.tmpComCoord;

        } else {
            //stateChangeMol.display_assoc_icoords("CURRORIENTATION_TOROTATE2");
            memRot = save_mem_orientation(memProtein, stateChangeMol, molTemplateList[stateChangeMol.molTypeIndex]);
            pivot = stateChangeMol.tmpComCoord;
        }

        /*double *pivot=new double[3];
	    pivot[0]=bases[slowPro].xcom;
	    pivot[1]=bases[slowPro].ycom;
	    pivot[2]=bases[slowPro].zcom;*/
        /*rotate the molecules and their complexes.*/
        //rotate_int_quat(pivot, ind_com[c1], bases, memRot);
        //rotate_int_quat(pivot, ind_com[c2], bases, memRot);
        rotate(pivot, memRot, facilitatorCom, moleculeList);
        rotate(pivot, memRot, stateChangeCom, moleculeList);
        // std::cout <<" AFTER ROTATION1: "<<std::endl;
        //facilitatorMol.display_assoc_icoords("molecule 1");
        //stateChangeMol.display_assoc_icoords("molecule 2");

        //	    std::cout <<"CRDS after forcing back to membrane bound orientation: \n";*/
    }
    //Coord finalCOM;//=new double[3];
    com_of_two_tmp_complexes(facilitatorCom, stateChangeCom, finalCOM, moleculeList); //com of c1+c2 (final (tmp) coordinates).
    //std::cout <<"FINAL COMPLEX PAIR COM: "<<finalCOM.x<<' '<<finalCOM.y<<' '<<finalCOM.z<<std::endl;
    /*Force finalCOM to startCOM, unless transitioning from 3D->2D*/
    Vector dtrans {};
    dtrans.x = startCOM.x - finalCOM.x;
    dtrans.y = startCOM.y - finalCOM.y;
    dtrans.z = startCOM.z - finalCOM.z;

    if (transitionToSurface == true)
        dtrans.z = 0.0; //don't move in z, now they are both on membrane
    // std::cout << "TRANSLATE COMPLEX PAIR TO ORIG COM BY SHIFTING: " << dtrans.x << ' ' << dtrans.y << ' ' << dtrans.z << std::endl; // update the temporary coordinates for both complexes
    for (auto& mp : facilitatorCom.memberList)
        moleculeList[mp].update_association_coords(dtrans);
    for (auto& mp : stateChangeCom.memberList)
        moleculeList[mp].update_association_coords(dtrans);

    double zchg = 0;
    {
        /*RECHECK HERE IF ANY OF THE LIPIDS ARE SLIGHTLY BELOW THE MEMBRANE. THIS CAN HAPPEN DUE TO PRECISION ISSUES
	      always use tmpCoords in this associate routine.
	     */
        dtrans.x = 0;
        dtrans.y = 0;
        for (auto& mp : facilitatorCom.memberList) {
            if (moleculeList[mp].isLipid == true) { //this is a lipid
                if (moleculeList[mp].tmpComCoord.z < -membraneObject.waterBox.z / 2.0) {
                    double ztmp = (-membraneObject.waterBox.z / 2.0) - moleculeList[mp].tmpComCoord.z; //lipid COM is below box bottom
                    if (ztmp > zchg)
                        zchg = ztmp; //largest dip below membrane
                }
                if (moleculeList[mp].tmpComCoord.z - 0.01 > -membraneObject.waterBox.z / 2.0) {

                    double ztmp = (-membraneObject.waterBox.z / 2.0)
                        - moleculeList[mp].tmpComCoord.z; // lipid COM is ABOVE box bottom

                    // std::cout << "WARNING, during state-change associate, LIPID IS ABOVE MEMBRANE BY " << -ztmp << '\n';
                }
            }
        }
        for (auto& mp : stateChangeCom.memberList) {
            if (moleculeList[mp].isLipid == true) { //this is a lipid
                if (moleculeList[mp].tmpComCoord.z < -membraneObject.waterBox.z / 2.0) {
                    double ztmp = (-membraneObject.waterBox.z / 2.0) - moleculeList[mp].tmpComCoord.z; //lipid COM is below box bottom
                    if (ztmp > zchg)
                        zchg = ztmp; //largest dip below membrane
                }
                if (moleculeList[mp].tmpComCoord.z - 0.01 > -membraneObject.waterBox.z / 2.0) {

                    double ztmp = (-membraneObject.waterBox.z / 2.0)
                        - moleculeList[mp].tmpComCoord.z; // lipid COM is ABOVE box bottom

                    // std::cout << "WARNING, during state-change associate, LIPID IS ABOVE MEMBRANE BY " << -ztmp << '\n';
                }
            }
        }
        dtrans.z = zchg;

        // update the temporary coordinates for both complexes
        for (auto& mp : facilitatorCom.memberList)
            moleculeList[mp].update_association_coords(dtrans);
        for (auto& mp : stateChangeCom.memberList)
            moleculeList[mp].update_association_coords(dtrans);

    } //is onMembrane

    // std::cout << " FINAL COORDS PRIOR TO OVERLAP CHECK " << std::endl;
    // //write_xyz_assoc("final.xyz", facilitatorCom, stateChangeCom, moleculeList);
    // write_xyz_assoc_cout(facilitatorCom, stateChangeCom, moleculeList);

    /*Reflect off the box.*/
    //Vector traj {};
    std::array<double, 3> traj; //=new double[3];
    for (int mm = 0; mm < 3; mm++)
        traj[mm] = 0;

    /*This needs to evaluate the traj update, based on it initially being zero.
	  And here, it should be called based on the tmpCoords, not the full coordinates. 
	  also requires updating the COM of this temporary new position
	*/
    update_complex_tmp_com_crds(facilitatorCom, moleculeList);
    //update_complex_tmp_com_crds(stateChangeCom, moleculeList);

    reflect_traj_tmp_crds(params, moleculeList, facilitatorCom, traj, membraneObject, RS3D); //uses tmpCoords to calculate traj.
    //reflect_traj_tmp_crds(params, moleculeList, stateChangeCom, traj, membraneObject);

    if (std::abs(traj[0] + traj[1] + traj[2]) > 1E-50) {
        // update the temporary coordinates for both complexes
        Vector vtraj { traj[0], traj[1], traj[2] };
        for (auto& mp : facilitatorCom.memberList)
            moleculeList[mp].update_association_coords(vtraj);
        for (auto& mp : stateChangeCom.memberList)
            moleculeList[mp].update_association_coords(vtraj);

        // std::cout << "CRDS after reflecting off of the BOX by " << traj[0] << ' ' << traj[1] << ' ' << traj[2] << std::endl;
        // write_xyz_assoc_cout(facilitatorCom, stateChangeCom, moleculeList);
    }
    /* CHECKS */
    bool cancelAssoc { false };
    check_for_structure_overlap(cancelAssoc, facilitatorCom, stateChangeCom, moleculeList, params, molTemplateList);

    if (cancelAssoc == false)
        check_if_spans_box(cancelAssoc, params, facilitatorCom, stateChangeCom, moleculeList, membraneObject);
    if (cancelAssoc == false)
        check_for_structure_overlap_system(cancelAssoc, facilitatorCom, stateChangeCom, moleculeList, params, molTemplateList, complexList, forwardRxns, backRxns);

    if (cancelAssoc) {
        // std::cout << "Canceling association, returning complexes to original state.\n";
        for (auto memMol : facilitatorCom.memberList)
            moleculeList[memMol].clear_tmp_association_coords();
        for (auto memMol : stateChangeCom.memberList)
            moleculeList[memMol].clear_tmp_association_coords();
        return;
    }

    // change the state. for IL, update species copy numbers
    int oldStateIndex { stateChangeMol.interfaceList[stateChangeIface].stateIndex };
    //    stateChangeMol.interfaceList[stateChangeIface].change_state(
    //        relStateIndex, newState.absIfaceIndex, newState.requiresState);

    counterArrays.copyNumSpecies[currRxn.reactantListNew[0].absIfaceIndex] -= 1;
    counterArrays.copyNumSpecies[currRxn.reactantListNew[1].absIfaceIndex] -= 1;
    counterArrays.copyNumSpecies[currRxn.productListNew[0].absIfaceIndex] += 1;
    counterArrays.copyNumSpecies[currRxn.productListNew[1].absIfaceIndex] += 1;

    //Update free number of lipids of IL for each state
    RxnIface implicitLipidState {};
    const auto& implicitLipidStateList = molTemplateList[stateChangeMol.molTypeIndex].interfaceList[stateChangeIface].stateList;
    if (molTemplateList[currRxn.reactantListNew[1].molTypeIndex].isImplicitLipid == true) {
        implicitLipidState = currRxn.reactantListNew[1];
    } else {
        implicitLipidState = currRxn.reactantListNew[0];
    }

    for (auto& state : implicitLipidStateList) {
        if (state.index == implicitLipidState.absIfaceIndex) {
            relStateIndex = static_cast<int>(&state - &implicitLipidStateList[0]);
            break;
        }
    }
    membraneObject.numberOfFreeLipidsEachState[relStateIndex] -= 1;

    if (molTemplateList[currRxn.productListNew[1].molTypeIndex].isImplicitLipid == true) {
        implicitLipidState = currRxn.productListNew[1];
    } else {
        implicitLipidState = currRxn.productListNew[0];
    }

    for (auto& state : implicitLipidStateList) {
        if (state.index == implicitLipidState.absIfaceIndex) {
            relStateIndex = static_cast<int>(&state - &implicitLipidStateList[0]);
            break;
        }
    }
    membraneObject.numberOfFreeLipidsEachState[relStateIndex] += 1;

    // write temporary to real coords and clear temporary coordinates
    for (auto memMol : facilitatorCom.memberList) {
        moleculeList[memMol].comCoord = moleculeList[memMol].tmpComCoord;
        for (unsigned int i { 0 }; i < moleculeList[memMol].interfaceList.size(); ++i)
            moleculeList[memMol].interfaceList[i].coord = moleculeList[memMol].tmpICoords[i];
        moleculeList[memMol].clear_tmp_association_coords();
        moleculeList[memMol].trajStatus = TrajStatus::propagated;
    }
    facilitatorCom.update_properties(moleculeList, molTemplateList);

    //for (auto memMol : stateChangeCom.memberList) {
    //    moleculeList[memMol].comCoord = moleculeList[memMol].tmpComCoord;
    //    for (unsigned int i { 0 }; i < moleculeList[memMol].interfaceList.size(); ++i)
    //        moleculeList[memMol].interfaceList[i].coord = moleculeList[memMol].tmpICoords[i];
    //    moleculeList[memMol].clear_tmp_association_coords();
    //    moleculeList[memMol].trajStatus = TrajStatus::propagated;
    //}
    //stateChangeCom.update_properties(moleculeList, molTemplateList);

    // clear implicit lipid properties
    moleculeList[stateChangeMol.index].clear_tmp_association_coords();
    stateChangeCom.memberList.clear();
    stateChangeCom.memberList.push_back(stateChangeMol.index);
    stateChangeCom.update_properties(moleculeList, molTemplateList);

    // Enforce boundary conditions
    reflect_complex_rad_rot(membraneObject, facilitatorCom, moleculeList, RS3D);

    //for (unsigned crossItr { 0 }; crossItr < stateChangeMol.crossbase.size(); ++crossItr) {
    //    int skipMol { stateChangeMol.crossbase[crossItr] };
    //    for (unsigned crossItr2 { 0 }; crossItr2 < moleculeList[skipMol].crossbase.size(); ++crossItr2) {
    //        if (moleculeList[skipMol].crossbase[crossItr2] == stateChangeMol.index)
    //            moleculeList[skipMol].probvec[crossItr2] = 0;
    //    }
    //}

    // Set probability of this protein to zero in all reactions so it doesn't try to
    // react again but the partners still will avoid overlapping.
    for (unsigned crossItr { 0 }; crossItr < facilitatorMol.crossbase.size(); ++crossItr) {
        int skipMol { facilitatorMol.crossbase[crossItr] };
        for (unsigned crossItr2 { 0 }; crossItr2 < moleculeList[skipMol].crossbase.size(); ++crossItr2) {
            if (moleculeList[skipMol].crossbase[crossItr2] == facilitatorMol.index)
                moleculeList[skipMol].probvec[crossItr2] = 0;
        }
    }

    // update observables, if applicable
    // TODO: Temporarily, if backRxn, iterate down, if forwardRxn, iterate up
    if (!isStateChangeBackRxn && forwardRxns[rxnIndex].isObserved) {
        auto observeItr = observablesList.find(forwardRxns[rxnIndex].observeLabel);
        if (observeItr == observablesList.end()) {
            // std::cerr << "WARNING: Observable " << forwardRxns[rxnIndex].observeLabel << " not defined.\n";
        } else {
            ++observeItr->second;
        }
    } else if (isStateChangeBackRxn && backRxns[rxnIndex].isObserved) {
        auto observeItr = observablesList.find(backRxns[rxnIndex].observeLabel);
        if (observeItr == observablesList.end()) {
            // std::cerr << "WARNING: Observable " << backRxns[rxnIndex].observeLabel << " not defined.\n";
        } else {
            --observeItr->second;
        }
    }

    /*Update species copy numbers*/
    //counterArrays.copyNumSpecies[oldStateIndex]--; // decrement ifaceIndex1
    //counterArrays.copyNumSpecies[stateChangeMol.interfaceList[stateChangeIface].stateIndex]++; // increment product state

    // Update the crossed molecule lists so that the current molecules won't avoid anything, but others will.
    stateChangeCom.ncross = -1;
    facilitatorCom.ncross = -1;
    stateChangeMol.crossbase.clear();
    facilitatorMol.crossbase.clear();
}