#include "boundary_conditions/reflect_functions.hpp"
#include "classes/class_Rxns.hpp"
#include "io/io.hpp"
#include "reactions/association/association.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"
#include <cmath>
#include <iomanip>

void associate_sphere(long long int iter, 
    int ifaceIndex1, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2,
    Complex& reactCom1, Complex& reactCom2, const Parameters& params,
    ForwardRxn& currRxn, std::vector<Molecule>& moleculeList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList,
    copyCounters& counterArrays, std::vector<Complex>& complexList,
    Membrane& membraneObject, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns)
{
    // TRACE();
    if (reactCom1.index == reactCom2.index) {
        // skip to protein interation updates
        // std::cout << "Closing a loop, no rotations performed.\n";
        counterArrays.nLoops++;
        // update the Molecule's TrajStatus (this is done in the else, when
        // Molecules are rotated but not otherwise)
        for (auto& memMol : reactCom1.memberList)
            moleculeList[memMol].trajStatus = TrajStatus::associated;
    } else { // not in the same complex
        if (reactMol2.isImplicitLipid == true || reactMol1.isImplicitLipid == true) {
            std::cout << "WRONG: implicit-lipid binding involves wrong function file. associate_ImplicitLipid_sphere should be called instead of associate_sphere!" << std::endl;
            exit(1);
        }

        // record the previous lastNumberUpdateItrEachMol & numEachMol
        std::vector<int> numEachMolPrevious1 {};
        std::vector<long long int> lastNumberUpdateItrEachMolPrevious1 {};
        std::vector<int> numEachMolPrevious2 {};
        std::vector<long long int> lastNumberUpdateItrEachMolPrevious2 {};
        for(unsigned index=0;index<molTemplateList.size();index++){
            numEachMolPrevious1.emplace_back(complexList[reactMol1.myComIndex].numEachMol[index]);
            lastNumberUpdateItrEachMolPrevious1.emplace_back(complexList[reactMol1.myComIndex].lastNumberUpdateItrEachMol[index]);
            numEachMolPrevious2.emplace_back(complexList[reactMol2.myComIndex].numEachMol[index]);
            lastNumberUpdateItrEachMolPrevious2.emplace_back(complexList[reactMol2.myComIndex].lastNumberUpdateItrEachMol[index]);
        }

        // set up temporary coordinates
        for (auto& memMol : reactCom1.memberList)
            moleculeList[memMol].set_tmp_association_coords();
        for (auto& mol : reactCom2.memberList)
            moleculeList[mol].set_tmp_association_coords();

        // create references to reacting interfaces
        Coord& reactIface1 = reactMol1.tmpICoords[ifaceIndex1];
        Coord& reactIface2 = reactMol2.tmpICoords[ifaceIndex2];

        // orientation corrections for membrane bound components
        bool isOnMembrane = false;
        bool transitionToSurface = false;

        double tol = 1E-14;
        /*Calculate COM of the two complexes pre-association. The COM of the new
     * complex after should be close to this Here, we will force it back, as
     * rotation can cause large displacements*/
        Coord startCOM;
        com_of_two_tmp_complexes(reactCom1, reactCom2, startCOM, moleculeList); // com of c1+c2 (original coordinates).

        /* MOVE PROTEIN TO SIGMA */
        {
            double DxSum { reactCom1.D.x + reactCom2.D.x };
            double DySum { reactCom1.D.y + reactCom2.D.y };
            double DzSum { reactCom1.D.z + reactCom2.D.z };

            Vector sigma { reactIface1 - reactIface2 };

            Vector transVec1 {};
            Vector transVec2 {};
            double displaceFrac {};
            double arc1Move, arc2Move;
            double sigmaMag;
            // if both in 2D, ignore the z-component
            if (DzSum < 1E-14) {
                isOnMembrane = true;
                currRxn.bindRadius2D = calc_bindRadius2D(currRxn.bindRadius, reactIface1);
                DzSum = 1; // to prevent divide by 0
                sigmaMag = get_geodesic_distance(reactIface1, reactIface2);
                //geodesic distance (arc-length) to displace molecule 1 and molecule 2.
                arc1Move = (sigmaMag - currRxn.bindRadius2D) / DxSum * reactCom1.D.x;
                arc2Move = (sigmaMag - currRxn.bindRadius2D) / DxSum * reactCom2.D.x;
            } else { // note not on the membrane
                sigmaMag = sqrt((sigma.x * sigma.x) + (sigma.y * sigma.y) + (sigma.z * sigma.z));
                displaceFrac = (sigmaMag - currRxn.bindRadius) / sigmaMag;
                /*At least one protein is in 3D*/
                if (reactCom1.D.z < tol || reactCom2.D.z < tol) {
                    transitionToSurface = true; // both can't be less than tol, or would not be in this loop.
                    // std::cout << "TRANSITIONING FROM 3D->2D " << std::endl;
                }
            }

            // now, calculate the transVec for each complex and move these two complexes
            Coord target1Pos;
            Coord target2Pos;
            if (isOnMembrane == true) { // both complexes are on the sphere
                target1Pos = find_position_after_association(arc1Move, reactIface1, reactIface2, sigmaMag, currRxn.bindRadius2D);
                target2Pos = find_position_after_association(arc2Move, reactIface2, reactIface1, sigmaMag, currRxn.bindRadius2D);
            } else if (reactCom1.D.z < 1E-14) { //complex1 is on sphere
                transVec1 = Vector(0, 0, 0);
                sigma.calc_magnitude();
                double lamda = (sigma.magnitude - currRxn.bindRadius) / sigma.magnitude;
                transVec2 = Vector(lamda * sigma);
            } else if (reactCom2.D.z < 1E-14) { //complex2 is on sphere
                transVec2 = Vector(0, 0, 0);
                sigma.calc_magnitude();
                double lamda = (sigma.magnitude - currRxn.bindRadius) / sigma.magnitude;
                transVec1 = Vector(-lamda * sigma);
            } else { // both complexes are in solution
                transVec1.x = -sigma.x * (reactCom1.D.x / DxSum) * displaceFrac;
                transVec1.y = -sigma.y * (reactCom1.D.y / DySum) * displaceFrac;
                transVec1.z = -sigma.z * (reactCom1.D.z / DzSum) * displaceFrac;
                transVec2.x = sigma.x * (reactCom2.D.x / DxSum) * displaceFrac;
                transVec2.y = sigma.y * (reactCom2.D.y / DySum) * displaceFrac;
                transVec2.z = sigma.z * (reactCom2.D.z / DzSum) * displaceFrac;
            }
            // update the temporary coordinates
            if (isOnMembrane == true) {
                reactCom1.update_association_coords_sphere(moleculeList, reactIface1, target1Pos);
                reactCom2.update_association_coords_sphere(moleculeList, reactIface2, target2Pos);
            } else {
                for (auto& mp : reactCom1.memberList)
                    moleculeList[mp].update_association_coords(transVec1);
                for (auto& mp : reactCom2.memberList)
                    moleculeList[mp].update_association_coords(transVec2);
                // std::cout << "Position after pushed to sigma: " << std::endl;
                // reactMol1.display_assoc_icoords("mol1");
                // reactMol2.display_assoc_icoords("mol2");
            }
        } // Move protein to sigma

        // now rotate both complexes
        if (molTemplateList[reactMol1.molTypeIndex].isPoint && molTemplateList[reactMol2.molTypeIndex].isPoint) {
            /*If both molecules are points, no orientations to specify*/
            // std::cout << " Move two point particles to contact along current "
            //              "separation vector, NO ORIENTATION \n";
        } else { // both are not points
            /* THETA */
            // std::cout << std::setw(8) << std::setfill('-') << ' ' << std::endl
            //           << "THETA 1" << std::endl
            //           << std::setw(8) << ' ' << std::setfill(' ') << std::endl;
            theta_rotation(reactIface1, reactIface2, reactMol1, reactMol2,
                currRxn.assocAngles.theta1, reactCom1, reactCom2,
                moleculeList);
            // std::cout << std::setw(30) << std::setfill('-') << ' '
            //           << std::setfill(' ') << std::endl;
            // std::cout << "THETA 2" << std::endl
            //           << std::setw(8) << std::setfill('-') << ' ' << std::setfill(' ')
            //           << std::endl;
            theta_rotation(reactIface2, reactIface1, reactMol2, reactMol1,
                currRxn.assocAngles.theta2, reactCom2, reactCom1,
                moleculeList);

            /* OMEGA */
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
            // std::cout << "P1 or P2 is a rod-type protein, no dihedral for "
            //              "associated complex."
            //           << std::endl;

            /* PHI */
            // PHI 1
            // std::cout << std::setw(6) << std::setfill('-') << ' ' << std::endl
            //           << "PHI 1" << std::endl
            //           << std::setw(6) << ' ' << std::setfill(' ') << std::endl;

            if (!std::isnan(currRxn.assocAngles.phi1)) {
                phi_rotation(reactIface1, reactIface2, ifaceIndex2, reactMol1,
                    reactMol2, reactCom1, reactCom2, currRxn.norm1,
                    currRxn.assocAngles.phi1, currRxn, moleculeList,
                    molTemplateList);
            } //else
            // std::cout << "P1 has no valid phi angle." << std::endl;

            // PHI 2
            // std::cout << std::setw(6) << std::setfill('-') << ' ' << std::endl
            //           << "PHI 2" << std::endl
            //           << std::setw(6) << ' ' << std::setfill(' ') << std::endl;

            if (!std::isnan(currRxn.assocAngles.phi2)) {
                phi_rotation(reactIface2, reactIface1, ifaceIndex1, reactMol2,
                    reactMol1, reactCom2, reactCom1, currRxn.norm2,
                    currRxn.assocAngles.phi2, currRxn, moleculeList,
                    molTemplateList);
            } //else
            // std::cout << "P2 has no valid phi angle." << std::endl;
        } // end of if points.

        /*FINISHED ROTATING, NO CONSTRAINTS APPLIED TO SURFACE REACTIONS*/
        Coord finalCOM;
        com_of_two_tmp_complexes(reactCom1, reactCom2, finalCOM, moleculeList); // com of c1+c2 (final (tmp) coordinates).
        // std::cout << "Pre-MEMBRANE ROT: COMPLEX PAIR COM: " << finalCOM.x << ' '
        //           << finalCOM.y << ' ' << finalCOM.z << std::endl;
        // reactMol1.display_assoc_icoords("mol1");
        // reactMol2.display_assoc_icoords("mol2");

        if (isOnMembrane == true || transitionToSurface == true) {
            //return orientation of normal back to starting position
            // std::cout << " IS ON MEMBRANE, CORRECT ORIENTATION ! " << std::endl;
            Molecule memProtein;
            Molecule Lipid;
            if (reactCom1.D.z < reactCom2.D.z) {
                set_memProtein_sphere(reactCom1, memProtein, moleculeList, membraneObject); // rotate relative to the slower protein.
                find_Lipid_sphere(reactCom1, Lipid, moleculeList, membraneObject);
            } else {
                set_memProtein_sphere(reactCom2, memProtein, moleculeList, membraneObject); // rotate relative to the slower protein.
                find_Lipid_sphere(reactCom2, Lipid, moleculeList, membraneObject);
            }
            Quat memRot = save_mem_orientation(memProtein, Lipid, molTemplateList[Lipid.molTypeIndex]);
            Coord pivot = Lipid.tmpComCoord;
            /*rotate the molecules and their complexes.*/
            rotate(pivot, memRot, reactCom1, moleculeList);
            rotate(pivot, memRot, reactCom2, moleculeList);
        }

        // Coord finalCOM;
        com_of_two_tmp_complexes(reactCom1, reactCom2, finalCOM, moleculeList); // com of c1+c2 (final (tmp) coordinates).
        //Force finalCOM to startCOM if in solution, while on sphere or transitionToSurface make sure lipids are still on sphere
        if (isOnMembrane == true || transitionToSurface == true) {
            // For explicit-lipid model, move the lipid onto sphere surface.
            // For implicit-lipid model, move the lipid-bound interface onto surface.
            Vector dtrans { 0, 0, 0 };
            if (membraneObject.implicitLipid == false) {
                Coord targ { 0.0, 0.0, 0.0 };
                double zchg = 0.0;
                for (auto& mp : reactCom1.memberList) {
                    if (moleculeList[mp].isLipid == true) { // this is a lipid
                        double rtmp = moleculeList[mp].tmpComCoord.get_magnitude();
                        double drtmp = std::abs(rtmp - membraneObject.sphereR); // lipid COM is off spherical membrane
                        if (drtmp > zchg) {
                            zchg = drtmp;
                            targ = moleculeList[mp].tmpComCoord; // farthest point outside membrane
                        }
                    } // this is a lipid
                }
                for (auto& mp : reactCom2.memberList) {
                    if (moleculeList[mp].isLipid == true) { // this is a lipid
                        double rtmp = moleculeList[mp].tmpComCoord.get_magnitude();
                        double drtmp = std::abs(rtmp - membraneObject.sphereR); // lipid COM is off spherical membrane
                        if (drtmp > zchg) {
                            zchg = drtmp;
                            targ = moleculeList[mp].tmpComCoord; // farthest point outside membrane
                        }
                    } // this is a lipid
                }
                if (targ.get_magnitude() > 1E-8) {
                    dtrans = Vector { (membraneObject.sphereR - targ.get_magnitude()) / targ.get_magnitude() * targ };
                }
            } else {
                Coord targ { 0.0, 0.0, 0.0 };
                double dr = 0;
                for (auto& mp : reactCom1.memberList) {
                    for (int i = 0; i < moleculeList[mp].interfaceList.size(); i++) {
                        if (moleculeList[mp].interfaceList[i].isBound == true) {
                            int index = moleculeList[mp].interfaceList[i].interaction.partnerIndex;
                            if (moleculeList[index].isImplicitLipid == true) {
                                double drtmp = moleculeList[mp].tmpICoords[i].get_magnitude() - membraneObject.sphereR;
                                if (std::abs(drtmp) >= std::abs(dr)) {
                                    dr = drtmp;
                                    targ = moleculeList[mp].tmpICoords[i];
                                }
                            }
                        }
                    }
                }
                for (auto& mp : reactCom2.memberList) {
                    for (int i = 0; i < moleculeList[mp].interfaceList.size(); i++) {
                        if (moleculeList[mp].interfaceList[i].isBound == true) {
                            int index = moleculeList[mp].interfaceList[i].interaction.partnerIndex;
                            if (moleculeList[index].isImplicitLipid == true) {
                                double drtmp = moleculeList[mp].tmpICoords[i].get_magnitude() - membraneObject.sphereR;
                                if (std::abs(drtmp) >= std::abs(dr)) {
                                    dr = drtmp;
                                    targ = moleculeList[mp].tmpICoords[i];
                                }
                            }
                        }
                    }
                }
                if (targ.get_magnitude() > 1E-8) {
                    dtrans = Vector { (membraneObject.sphereR - targ.get_magnitude()) / targ.get_magnitude() * targ };
                }
            }

            // std::cout << " Lipid is off spherical membrane, shift up by: " << dtrans.x << " " << dtrans.y << " " << dtrans.z
            //           << std::endl;
            // update the temporary coordinates for both complexes
            for (auto& mp : reactCom1.memberList)
                moleculeList[mp].update_association_coords(dtrans);
            for (auto& mp : reactCom2.memberList)
                moleculeList[mp].update_association_coords(dtrans);
        } else {
            /*
	    Vector dtrans;
            dtrans.x = startCOM.x - finalCOM.x;
            dtrans.y = startCOM.y - finalCOM.y;
            dtrans.z = startCOM.z - finalCOM.z;
            std::cout << "TRANSLATE COMPLEX PAIR TO ORIG COM BY SHIFTING: "
                      << dtrans.x << ' ' << dtrans.y << ' ' << dtrans.z
                      << std::endl; // update the temporary coordinates for both complexes
            for (auto& mp : reactCom1.memberList)
                moleculeList[mp].update_association_coords(dtrans);
            for (auto& mp : reactCom2.memberList)
                moleculeList[mp].update_association_coords(dtrans);
	  */
        }

        // std::cout << " FINAL COORDS PRIOR TO OVERLAP CHECK  AND REFLECT OFF SPHERE: "
        //           << std::endl;
        // reactMol1.display_assoc_icoords("mol1");
        // reactMol2.display_assoc_icoords("mol2");

        /*Reflect off the sphere.*/
        std::array<double, 3> traj;
        for (int mm = 0; mm < 3; mm++)
            traj[mm] = 0;

        /*This needs to evaluate the traj update, based on it initially being zero.
      And here, it should be called based on the tmpCoords, not the full
      coordinates. also requires updating the COM of this temporary new position
    */
        update_complex_tmp_com_crds(reactCom1, moleculeList);
        update_complex_tmp_com_crds(reactCom2, moleculeList);

        reflect_traj_tmp_crds(params, moleculeList, reactCom1, traj, membraneObject, 0.0); // uses tmpCoords to calculate traj.
        reflect_traj_tmp_crds(params, moleculeList, reactCom2, traj, membraneObject, 0.0);

        if (std::abs(traj[0] + traj[1] + traj[2]) > 1E-15) {
            // update the temporary coordinates for both complexes
            Vector vtraj { traj[0], traj[1], traj[2] };
            for (auto& mp : reactCom1.memberList)
                moleculeList[mp].update_association_coords(vtraj);
            for (auto& mp : reactCom2.memberList)
                moleculeList[mp].update_association_coords(vtraj);

            // std::cout << "CRDS after reflecting off of the SPHERE by " << traj[0] << ' '
            //           << traj[1] << ' ' << traj[2] << std::endl;
            // reactMol1.display_assoc_icoords("mol1");
            // reactMol2.display_assoc_icoords("mol2");
        }
        /*Calculate the angles swept out by the interface to COM vectors as a result of the displacement*/
        calc_angular_displacement(ifaceIndex1, ifaceIndex2, reactMol1, reactMol2, reactCom1, reactCom2, moleculeList);

        /* CHECKS AFTER ASSOCIATION FOR STERIC COLLISIONS, FOR EXPANDING BEYOND THE
       BOX SIZE OR FOR MOVING PROTEINS A LARGE DISTANCE DUE TO SNAPPING INTO
       PLACE
    */
        bool cancelAssoc { false };
        check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2, moleculeList, params, molTemplateList);

        if (cancelAssoc == false) {
            check_if_spans_sphere(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);
            if (cancelAssoc == true)
                counterArrays.nCancelSpanBox++;
        } else
            counterArrays.nCancelOverlapPartner++; //true for structure overlap check.

        if (cancelAssoc == false) {
            check_for_structure_overlap_system(cancelAssoc, reactCom1, reactCom2,
                moleculeList, params, molTemplateList,
                complexList, forwardRxns, backRxns);
            if (cancelAssoc == true)
                counterArrays.nCancelOverlapSystem++;
        }
        if (cancelAssoc == false) {
            measure_complex_displacement(cancelAssoc, reactCom1, reactCom2,
                moleculeList, params, molTemplateList,
                complexList);
            if (cancelAssoc == true) {
                if (isOnMembrane)
                    counterArrays.nCancelDisplace2D++;
                else if (transitionToSurface)
                    counterArrays.nCancelDisplace3Dto2D++;
                else
                    counterArrays.nCancelDisplace3D++;
            }
        }
        if (cancelAssoc) {
            // std::cout << "Canceling association, returning complexes to original state.\n";
            for (auto memMol : reactCom1.memberList)
                moleculeList[memMol].clear_tmp_association_coords();
            for (auto memMol : reactCom2.memberList)
                moleculeList[memMol].clear_tmp_association_coords();
            // end routine here!
            return;
        }
        counterArrays.nAssocSuccess++; //keep track of total number of successful association moves.

        /*Keep track of the sizes of complexes that associated*/
        track_association_events(reactCom1, reactCom2, transitionToSurface, isOnMembrane, counterArrays);

        // temporary coordinates
        for (auto memMol : reactCom1.memberList) {
            moleculeList[memMol].comCoord = moleculeList[memMol].tmpComCoord;
            for (unsigned int i { 0 }; i < moleculeList[memMol].interfaceList.size(); ++i)
                moleculeList[memMol].interfaceList[i].coord = moleculeList[memMol].tmpICoords[i];
            moleculeList[memMol].clear_tmp_association_coords();
            moleculeList[memMol].trajStatus = TrajStatus::propagated;
        }

        for (auto memMol : reactCom2.memberList) {
            moleculeList[memMol].comCoord = moleculeList[memMol].tmpComCoord;
            for (unsigned int i { 0 }; i < moleculeList[memMol].interfaceList.size(); ++i)
                moleculeList[memMol].interfaceList[i].coord = moleculeList[memMol].tmpICoords[i];
            moleculeList[memMol].clear_tmp_association_coords();
            // push reactCom2 into reactCom1
            if (currRxn.rxnType != ReactionType::biMolStateChange) {
                moleculeList[memMol].myComIndex = reactCom1.index; // update their complex index
                reactCom1.memberList.push_back(memMol);
            }
            moleculeList[memMol].trajStatus = TrajStatus::propagated;
        }

        // update complexes

        // add lifetime to the previous small clusters
        for(unsigned index=0;index<molTemplateList.size();index++){
            if(molTemplateList[index].countTransition == true){
                // com2 destroty
                molTemplateList[index].lifeTime[numEachMolPrevious2[index]-1].emplace_back((iter-lastNumberUpdateItrEachMolPrevious2[index])*Parameters::dt/1E6);
            }
        }

        reactCom2.memberList.clear(); // clear the member list so the molecules don't get destroyed
        reactCom2.destroy(moleculeList, complexList); // destroy the complex
        reactCom1.update_properties(moleculeList, molTemplateList); // recalculate the properties of the first complex

        for(unsigned index=0;index<molTemplateList.size();index++){
            if(molTemplateList[index].countTransition == true && complexList[reactMol1.myComIndex].numEachMol[index]>numEachMolPrevious1[index]){
                // grow
                molTemplateList[index].lifeTime[numEachMolPrevious1[index]-1].emplace_back((iter-lastNumberUpdateItrEachMolPrevious1[index])*Parameters::dt/1E6);
                complexList[reactMol1.myComIndex].lastNumberUpdateItrEachMol[index]=iter;
            }
        }

        // update transition matrix
        // compare current cluster size with the previous ones
        for(unsigned index=0;index<molTemplateList.size();index++){
            if(molTemplateList[index].countTransition == true && complexList[reactMol1.myComIndex].numEachMol[index]>numEachMolPrevious1[index]){
                // grow
                if(numEachMolPrevious1[index]-1 >= 0)
                    molTemplateList[index].transitionMatrix[numEachMolPrevious1[index]-1][numEachMolPrevious1[index]-1] += iter-Parameters::lastUpdateTransition[index]-1;
                if(numEachMolPrevious2[index]-1 >= 0)
                    molTemplateList[index].transitionMatrix[numEachMolPrevious2[index]-1][numEachMolPrevious2[index]-1] += iter-Parameters::lastUpdateTransition[index]-1;
                if(complexList[reactMol1.myComIndex].numEachMol[index]-1 >= 0 && numEachMolPrevious1[index]-1 >= 0)
                    molTemplateList[index].transitionMatrix[complexList[reactMol1.myComIndex].numEachMol[index]-1][numEachMolPrevious1[index]-1] += 1;
                if(complexList[reactMol1.myComIndex].numEachMol[index]-1 >= 0 && numEachMolPrevious2[index]-1 >= 0)
                    molTemplateList[index].transitionMatrix[complexList[reactMol1.myComIndex].numEachMol[index]-1][numEachMolPrevious2[index]-1] += 1;

                // update diagonal elements for unchanged complexes
                for(unsigned indexCom=0;indexCom<complexList.size();indexCom++){
                    if(indexCom != reactMol1.myComIndex){
                        if(complexList[indexCom].numEachMol[index]-1 >= 0){
                            molTemplateList[index].transitionMatrix[complexList[indexCom].numEachMol[index]-1][complexList[indexCom].numEachMol[index]-1] += iter-Parameters::lastUpdateTransition[index];
                        }
                    }
                }

                // update lastUpdateTransition
                Parameters::lastUpdateTransition[index] = iter;
            }
        }

        // Enforce boundary conditions
        reflect_complex_rad_rot(membraneObject, reactCom1, moleculeList, 0.0);

    } // end of if these molecules are closing a loop or not.
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
    // reactMol2
    {
        Molecule& oneMol { reactMol2 };
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
    reactMol2.interfaceList[ifaceIndex2].interaction.partnerIndex = reactMol1.index;
    reactMol1.interfaceList[ifaceIndex1].interaction.partnerIfaceIndex = ifaceIndex2;
    reactMol2.interfaceList[ifaceIndex2].interaction.partnerIfaceIndex = ifaceIndex1;
    if (currRxn.isReversible) {
        reactMol1.interfaceList[ifaceIndex1].interaction.conjBackRxn = currRxn.conjBackRxnIndex;
        reactMol2.interfaceList[ifaceIndex2].interaction.conjBackRxn = currRxn.conjBackRxnIndex;
    }

    reactMol1.interfaceList[ifaceIndex1].isBound = true;
    reactMol1.interfaceList[ifaceIndex1].index = currRxn.productListNew[0].absIfaceIndex;
    reactMol2.interfaceList[ifaceIndex2].isBound = true;
    reactMol2.interfaceList[ifaceIndex2].index = currRxn.productListNew[0].absIfaceIndex;

    // add to the list of bound interfaces and remove from the list of free interfaces
    reactMol1.bndlist.push_back(ifaceIndex1);
    reactMol2.bndlist.push_back(ifaceIndex2);
    reactMol1.bndpartner.push_back(reactMol2.index);
    reactMol2.bndpartner.push_back(reactMol1.index);

    {
        size_t tmpItr { 0 };
        for (unsigned i { 0 }; i < reactMol1.freelist.size(); ++i) {
            if (reactMol1.freelist[i] == ifaceIndex1)
                tmpItr = i;
        }
        reactMol1.freelist[tmpItr] = reactMol1.freelist.back();
        reactMol1.freelist.pop_back();
    }
    {
        size_t tmpItr { 0 };
        for (unsigned i { 0 }; i < reactMol2.freelist.size(); ++i) {
            if (reactMol2.freelist[i] == ifaceIndex2)
                tmpItr = i;
        }
        reactMol2.freelist[tmpItr] = reactMol2.freelist.back();
        reactMol2.freelist.pop_back();
    }

    // Set probability of this protein to zero in all reactions so it doesn't try
    // to react again but the partners still will avoid overlapping.
    for (unsigned crossItr { 0 }; crossItr < reactMol1.crossbase.size(); ++crossItr) {
        int skipMol { reactMol1.crossbase[crossItr] };
        for (unsigned crossItr2 { 0 }; crossItr2 < moleculeList[skipMol].crossbase.size(); ++crossItr2) {
            if (moleculeList[skipMol].crossbase[crossItr2] == reactMol1.index)
                moleculeList[skipMol].probvec[crossItr2] = 0;
        }
    }

    // Set probability of this protein to zero in all reactions so it doesn't try
    // to react again but the partners still will avoid overlapping.
    for (unsigned crossItr { 0 }; crossItr < reactMol2.crossbase.size(); ++crossItr) {
        int skipMol { reactMol2.crossbase[crossItr] };
        for (unsigned crossItr2 { 0 }; crossItr2 < moleculeList[skipMol].crossbase.size(); ++crossItr2) {
            if (moleculeList[skipMol].crossbase[crossItr2] == reactMol2.index)
                moleculeList[skipMol].probvec[crossItr2] = 0;
        }
    }

    // Update the crossed molecule lists so that the current molecules won't avoid anything, but others will.
    reactCom1.ncross = -1;
    reactMol1.crossbase.clear();
    reactMol2.crossbase.clear();

    /*update the number of bound species*/
    update_Nboundpairs(reactMol1.molTypeIndex, reactMol2.molTypeIndex, 1, params, counterArrays);
    /*Update species copy numbers*/
    counterArrays.copyNumSpecies[currRxn.reactantListNew[0].absIfaceIndex]--; // decrement ifaceIndex1
    counterArrays.copyNumSpecies[currRxn.reactantListNew[1].absIfaceIndex]--; // decrement ifaceIndex2
    counterArrays.copyNumSpecies[currRxn.productListNew[0].absIfaceIndex]++; // increment product state

    //-------------------------UPDATE BINDPAIRLIST----------------------------
    if (counterArrays.canDissociate[currRxn.productListNew[0].absIfaceIndex]) {
        // int molIndex { reactMol1.index };
        // if (reactMol1.index > reactMol2.index)
        //     molIndex = reactMol2.index;
        int molIndex { -1 };
        if (ifaceIndex1 == currRxn.reactantListNew[0].relIfaceIndex) {
            molIndex = reactMol1.index;
        } else {
            molIndex = reactMol2.index;
        }

        counterArrays.bindPairList[currRxn.productListNew[0].absIfaceIndex].emplace_back(molIndex);
    }
    //-----------------------END UPDATE BINDPAIRLIST--------------------------

    // std::cout << " After ASSOCIATE, CHANGE COPY NUMBERS, interfaces: "
    //           << ifaceIndex1 << ' ' << ifaceIndex2
    //           << " add to product: " << currRxn.productListNew[0].absIfaceIndex
    //           << " sub from reactants: "
    //           << currRxn.reactantListNew[0].absIfaceIndex << " "
    //           << currRxn.reactantListNew[1].absIfaceIndex << std::endl;

    // TODO: Insert species tracking here
    if (currRxn.isObserved) {
        auto obsItr = observablesList.find(currRxn.observeLabel);
        if (obsItr != observablesList.end())
            ++obsItr->second;
    }
}
