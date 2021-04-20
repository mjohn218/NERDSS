#include "boundary_conditions/reflect_functions.hpp"
#include "classes/class_Rxns.hpp"
#include "io/io.hpp"
#include "reactions/association/association.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"
#include <cmath>
#include <iomanip>

/*Either 1 or 2 is implicit lipid*/
void associate_implicitlipid_box(int ifaceIndex1, int ifaceIndex2, Molecule& reactMol1, Molecule& reactMol2, Complex& reactCom1,
    Complex& reactCom2, const Parameters& params, ForwardRxn& currRxn, std::vector<Molecule>& moleculeList,
    std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays, std::vector<Complex>& complexList, Membrane& membraneObject, const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns)
{
    double RS3D { -1.0 };
    for (int RS3Di = 0; RS3Di < 100; RS3Di++) {
        if ((std::abs(membraneObject.RS3Dvect[RS3Di] - currRxn.bindRadius) < 1E-15) && (std::abs(membraneObject.RS3Dvect[RS3Di + 100] - currRxn.rateList[0].rate) < 1E-15) && std::abs(membraneObject.RS3Dvect[RS3Di + 200] - (1.0 / 3.0 * (molTemplateList[currRxn.reactantListNew[0].molTypeIndex].D.x + molTemplateList[currRxn.reactantListNew[1].molTypeIndex].D.x) + 1.0 / 3.0 * (molTemplateList[currRxn.reactantListNew[0].molTypeIndex].D.y + molTemplateList[currRxn.reactantListNew[1].molTypeIndex].D.y) + 1.0 / 3.0 * (molTemplateList[currRxn.reactantListNew[0].molTypeIndex].D.z + molTemplateList[currRxn.reactantListNew[1].molTypeIndex].D.z))) < 1E-15) {
            RS3D = membraneObject.RS3Dvect[RS3Di + 300]; //look up value of the RS3D in a table.
            break;
        }
    }

    // we need to set implicit-lipid's temporary position according to mol.
    if (reactMol2.isImplicitLipid == true) {
        reactMol2.create_position_implicit_lipid(reactMol1, ifaceIndex2, currRxn.bindRadius, membraneObject);
        reactCom2.comCoord = reactMol2.comCoord;
    } else {
        reactMol1.create_position_implicit_lipid(reactMol2, ifaceIndex1, currRxn.bindRadius, membraneObject);
        reactCom1.comCoord = reactMol1.comCoord;
    }

    // set up temporary coordinates
    for (auto& memMol : reactCom1.memberList)
        moleculeList[memMol].set_tmp_association_coords();

    for (auto& mol : reactCom2.memberList)
        moleculeList[mol].set_tmp_association_coords();

    bool isOnMembrane = false;
    if (reactCom1.D.z < 1E-15 && reactCom2.D.z < 1E-15) {
        isOnMembrane = true;
        // std::cout << " IMPLICIT LIPID 2D CASE" << std::endl;
    }
    bool transitionToSurface = false;
    if (isOnMembrane == false) {
        transitionToSurface = true;
        // std::cout << " IMPLICIT LIPID 3D->2D CASE" << std::endl;
    }
    if (isOnMembrane == false) {
        /*
      Perform reorientaion for 3D->2D case,
      and for 2D case, to correct for proteins that have dropped below their 
      need orientation corrections for implicit model
    */
        // std::cout << "Before operate: " << std::endl;
        // write_xyz_assoc_cout(reactCom1, reactCom2, moleculeList);

        // create references to reacting interfaces
        Coord& reactIface1 = reactMol1.tmpICoords[ifaceIndex1];
        Coord& reactIface2 = reactMol2.tmpICoords[ifaceIndex2];
        Vector sigma { reactIface1 - reactIface2 };
        // orientation corrections for membrane bound components
        Molecule memProtein = reactMol2;
        int slowPro = reactMol2.index;

        if (reactCom1.D.x < reactCom2.D.x) {
            slowPro = reactMol1.index;
            memProtein = reactMol1; // rotate relative to the slower protein. Clat + IL, IL is slower; Clat-IL + IL, Clat-IL is slower
        }

        double tol = 1E-14;
        /* Calculate COM of the two complexes pre-association. The COM of the new complex after should be close to this
         * Here, we will force it back, as rotation can cause large displacements
         * */
        Coord startCOM;
        com_of_two_tmp_complexes(reactCom1, reactCom2, startCOM, moleculeList);

        // MOVE PROTEIN TO THE MEMBRANE.
        {
            sigma.calc_magnitude();
            double sigmaMag = sigma.magnitude;
            Vector transVec1 = Vector(-(sigmaMag - currRxn.bindRadius - RS3D) / sigmaMag * sigma);
            Vector transVec2 { 0.0, 0.0, 0.0 };
            if (reactMol1.isImplicitLipid == true) {
                transVec2 = Vector((sigmaMag - currRxn.bindRadius - RS3D) / sigmaMag * sigma);
                transVec1.x = 0;
                transVec1.y = 0;
                transVec1.z = 0;
            }
            // if (reactMol2.isImplicitLipid == true)
            //     std::cout << " transVec1 to sigma: " << transVec1.x << ' ' << transVec1.y << ' ' << transVec1.z << std::endl;
            // if (reactMol1.isImplicitLipid == true)
            //     std::cout << " transVec2 to sigma: " << transVec2.x << ' ' << transVec2.y << ' ' << transVec2.z << std::endl;

            // update the temporary coordinates
            for (auto& mp : reactCom1.memberList)
                moleculeList[mp].update_association_coords(transVec1);
            for (auto& mp : reactCom2.memberList)
                moleculeList[mp].update_association_coords(transVec2);
        }

        // std::cout << "After move protein to the membrane: " << std::endl;
        // write_xyz_assoc_cout(reactCom1, reactCom2, moleculeList);

        if ((molTemplateList[reactMol1.molTypeIndex].isPoint && reactMol2.isImplicitLipid == true) || (molTemplateList[reactMol2.molTypeIndex].isPoint && reactMol1.isImplicitLipid == true)) {
            // binding protein is point, no orientations to specify
            // std::cout << " Move point particle to implicit lipid NO ORIENTATION." << std::endl;
        } else { // not point.
            /* THETA */
            // std::cout << std::setw(8) << std::setfill('-') << ' ' << std::endl
            //           << "THETA 1" << std::endl
            //           << std::setw(8) << ' ' << std::setfill(' ') << std::endl;
            theta_rotation(reactIface1, reactIface2, reactMol1, reactMol2, currRxn.assocAngles.theta1, reactCom1, reactCom2, moleculeList);

            // std::cout << std::setw(30) << std::setfill('-') << ' ' << std::setfill(' ') << std::endl;
            // std::cout << "THETA 2" << std::endl
            //           << std::setw(8) << std::setfill('-') << ' ' << std::setfill(' ') << std::endl;
            theta_rotation(reactIface2, reactIface1, reactMol2, reactMol1, currRxn.assocAngles.theta2, reactCom2, reactCom1, moleculeList);

            /* OMEGA */
            // if protein has theta M_PI, uses protein norm instead of com_iface vector
            // std::cout << std::setw(6) << std::setfill('-') << ' ' << std::endl
            //           << "OMEGA" << std::endl
            //           << std::setw(6) << ' ' << std::setfill(' ') << std::endl;
            if (!std::isnan(currRxn.assocAngles.omega)) {
                omega_rotation(reactIface1, reactIface2, ifaceIndex2, reactMol1, reactMol2, reactCom1, reactCom2, currRxn.assocAngles.omega, currRxn, moleculeList, molTemplateList);
            } //else
            // std::cout << "P1 or P2 is a rod-type protein, no dihedral for associated complex." << std::endl;

            /* PHI */
            // PHI 1
            // std::cout << std::setw(6) << std::setfill('-') << ' ' << std::endl
            //           << "PHI 1" << std::endl
            //           << std::setw(6) << ' ' << std::setfill(' ') << std::endl;
            if (!std::isnan(currRxn.assocAngles.phi1)) {
                phi_rotation(reactIface1, reactIface2, ifaceIndex2, reactMol1, reactMol2, reactCom1, reactCom2, currRxn.norm1, currRxn.assocAngles.phi1, currRxn, moleculeList, molTemplateList);
            } //else
            // std::cout << "P1 has no valid phi angle." << std::endl;
            // PHI 2
            // std::cout << std::setw(6) << std::setfill('-') << ' ' << std::endl
            //           << "PHI 2" << std::endl
            //           << std::setw(6) << ' ' << std::setfill(' ') << std::endl;

            if (!std::isnan(currRxn.assocAngles.phi2)) {
                phi_rotation(reactIface2, reactIface1, ifaceIndex1, reactMol2, reactMol1, reactCom2, reactCom1, currRxn.norm2, currRxn.assocAngles.phi2, currRxn, moleculeList, molTemplateList);
            } //else
            // std::cout << "P2 has no valid phi angle." << std::endl;

            /*FINISHED ROTATING, NO CONSTRAINTS APPLIED TO SURFACE REACTIONS*/
            Coord finalCOM;
            com_of_two_tmp_complexes(reactCom1, reactCom2, finalCOM, moleculeList);

            // std::cout << "After rotation: " << std::endl;
            // write_xyz_assoc_cout(reactCom1, reactCom2, moleculeList);

            if (isOnMembrane == true || transitionToSurface == true) {
                /*return orientation of normal back to starting position*/
                Quat memRot;
                Coord pivot;
                if (slowPro == reactMol1.index) {
                    memRot = save_mem_orientation(memProtein, reactMol1, molTemplateList[reactMol1.molTypeIndex]);
                    pivot = reactMol1.tmpComCoord;

                } else {
                    memRot = save_mem_orientation(memProtein, reactMol2, molTemplateList[reactMol2.molTypeIndex]);
                    pivot = reactMol2.tmpComCoord;
                }
                rotate(pivot, memRot, reactCom1, moleculeList);
                rotate(pivot, memRot, reactCom2, moleculeList);
            }

            com_of_two_tmp_complexes(reactCom1, reactCom2, finalCOM, moleculeList);
            Vector dtrans {};
            dtrans.x = startCOM.x - finalCOM.x;
            dtrans.y = startCOM.y - finalCOM.y;
            dtrans.z = startCOM.z - finalCOM.z;

            if (isOnMembrane == true || transitionToSurface == true)
                dtrans.z = 0.0;
            for (auto& mp : reactCom1.memberList)
                moleculeList[mp].update_association_coords(dtrans);
            for (auto& mp : reactCom2.memberList)
                moleculeList[mp].update_association_coords(dtrans);

            double zchg = 0;
            bool isBelowBottom = false;
            bool isAboveBottom = false;
            if (isOnMembrane == true || transitionToSurface == true) {
                /* For the impicit binding, we need to move the impicit lipid to the bottom of box */
                dtrans.x = 0;
                dtrans.y = 0;
                for (auto& mp : reactCom1.memberList) {
                    if (moleculeList[mp].isImplicitLipid == true) {
                        if (moleculeList[mp].tmpComCoord.z < -membraneObject.waterBox.z / 2.0) {
                            double ztmp = (-membraneObject.waterBox.z / 2.0)
                                - moleculeList[mp].tmpComCoord.z; // lipid COM is below box bottom
                            if (ztmp > zchg) {
                                zchg = ztmp; // largest dip below membrane
                                isBelowBottom = true;
                            }
                        }
                        if (moleculeList[mp].tmpComCoord.z - 0.01 > -membraneObject.waterBox.z / 2.0) {
                            double ztmp = (-membraneObject.waterBox.z / 2.0) - moleculeList[mp].tmpComCoord.z; // lipid COM is ABOVE box bottom, here ztmp is negtive
                            // std::cout << "WARNING, during associate, implicit LIPID IS ABOVE MEMBRANE BY " << -ztmp << '\n';
                            // move the lipid back to the bottom
                            if (-ztmp > zchg) {
                                zchg = -ztmp;
                                isAboveBottom = true;
                            }
                        }
                    }
                }
                for (auto& mp : reactCom2.memberList) {
                    if (moleculeList[mp].isImplicitLipid == true) {
                        if (moleculeList[mp].tmpComCoord.z < -membraneObject.waterBox.z / 2.0) {
                            double ztmp = (-membraneObject.waterBox.z / 2.0) - moleculeList[mp].tmpComCoord.z; // lipid COM is below box bottom
                            if (ztmp > zchg) {
                                zchg = ztmp; // largest dip below membrane
                                isBelowBottom = true;
                            }
                        }
                        if (moleculeList[mp].tmpComCoord.z - 0.01 > -membraneObject.waterBox.z / 2.0) {
                            double ztmp = (-membraneObject.waterBox.z / 2.0) - moleculeList[mp].tmpComCoord.z; // lipid COM is ABOVE box bottom, here ztmp is negtive
                            // std::cout << "WARNING, during associate, implicit LIPID IS ABOVE MEMBRANE BY " << -ztmp << '\n';
                            // move the lipid back to the bottom
                            if (-ztmp > zchg) {
                                zchg = -ztmp;
                                isAboveBottom = true;
                            }
                        }
                    }
                }
                if (isBelowBottom == true) {
                    dtrans.z = zchg;
                    // std::cout << "Implicit Lipid is below membrane, shift up by: " << zchg << std::endl;
                }
                if (isAboveBottom == true) {
                    dtrans.z = -zchg;
                    // std::cout << "Implicit Lipid is above membrane, shift down by: " << zchg << std::endl;
                }

                // update the temporary coordinates for both complexes
                for (auto& mp : reactCom1.memberList)
                    moleculeList[mp].update_association_coords(dtrans);
                for (auto& mp : reactCom2.memberList)
                    moleculeList[mp].update_association_coords(dtrans);
            }
        } //only correct orientations if both particles are not points.

        // std::cout << "Before check overlop: " << std::endl;
        // write_xyz_assoc_cout(reactCom1, reactCom2, moleculeList);

        //Reflect off the box.
        std::array<double, 3> traj; //=new double[3];
        for (int mm = 0; mm < 3; mm++)
            traj[mm] = 0;

        /* This needs to evaluate the traj update, based on it initially being zero.
         * And here, it should be called based on the tmpCoords, not the full coordinates.
         * also requires updating the COM of this temporary new position 
         * */
        update_complex_tmp_com_crds(reactCom1, moleculeList);
        reflect_traj_tmp_crds(params, moleculeList, reactCom1, traj, membraneObject, RS3D); // uses tmpCoords to calculate traj.
        if (std::abs(traj[0] + traj[1] + traj[2]) > 1E-50) {
            // update the temporary coordinates for both complexes
            Vector vtraj { traj[0], traj[1], traj[2] };
            for (auto& mp : reactCom1.memberList)
                moleculeList[mp].update_association_coords(vtraj);
            for (auto& mp : reactCom2.memberList)
                moleculeList[mp].update_association_coords(vtraj);

            // std::cout << "CRDS after reflecting off of the BOX by " << traj[0] << ' ' << traj[1] << ' ' << traj[2]
            //           << std::endl;
            // write_xyz_assoc_cout(reactCom1, reactCom2, moleculeList);
        }
        /* CHECKS */
        bool cancelAssoc { false };
        //check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2, moleculeList, params, molTemplateList);
        check_if_spans_box(cancelAssoc, params, reactCom1, reactCom2, moleculeList, membraneObject);
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
        //done deciding whether to reject move.

        if (cancelAssoc) {
            // std::cout << "Canceling association, returning complexes to original state.\n";
            for (auto memMol : reactCom1.memberList)
                moleculeList[memMol].clear_tmp_association_coords();
            for (auto memMol : reactCom2.memberList)
                moleculeList[memMol].clear_tmp_association_coords();
            return;
            //EXIT FROM THIS ROUTINE.
        }

    } //if already in 2D, no need for reorientation.

    counterArrays.nAssocSuccess++; //keep track of total number of successful association moves.
    /*Keep track of the sizes of complexes that associated*/
    track_association_events(reactCom1, reactCom2, transitionToSurface, isOnMembrane, counterArrays);

    if (reactMol2.isImplicitLipid == true) {
        // write temporary to real coords and clear temporary coordinates
        for (auto memMol : reactCom1.memberList) {
            moleculeList[memMol].comCoord = moleculeList[memMol].tmpComCoord;
            for (unsigned int i { 0 }; i < moleculeList[memMol].interfaceList.size(); ++i)
                moleculeList[memMol].interfaceList[i].coord = moleculeList[memMol].tmpICoords[i];
            moleculeList[memMol].clear_tmp_association_coords();
            moleculeList[memMol].trajStatus = TrajStatus::propagated;
        }
        //add link to surface for the molecule
        moleculeList[reactMol1.index].linksToSurface++;
        // update complexes
        reactCom1.iLipidIndex = reactMol2.index; //index in molecule list of the implicit lipid.
        reactCom1.linksToSurface++;
        reactCom1.update_properties(moleculeList, molTemplateList); // recalculate the properties of the first complex
        reactCom1.OnSurface = true; // the complex is on membrane now.
        reactCom1.D.z = 0;
        /*Clear implicit lipid properties*/
        moleculeList[reactMol2.index].clear_tmp_association_coords();
        reactCom2.memberList.clear();
        reactCom2.memberList.push_back(reactMol2.index);
        // std::cout << " Size of IL's complex:" << reactCom2.memberList.size() << " Interfaces on IL: " << moleculeList[reactMol2.index].interfaceList.size() << std::endl;
        reactCom2.update_properties(moleculeList, molTemplateList); // recalculate the properties of the second complex
        //Enforce boundary conditions
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
        //**************************************************
        //**************************************************
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

        // Set probability of this protein to zero in all reactions so it doesn't try to
        // react again but the partners still will avoid overlapping.
        for (unsigned crossItr { 0 }; crossItr < reactMol1.crossbase.size(); ++crossItr) {
            int skipMol { reactMol1.crossbase[crossItr] };
            if (!moleculeList[skipMol].isImplicitLipid) {
                for (unsigned crossItr2 { 0 }; crossItr2 < moleculeList[skipMol].crossbase.size(); ++crossItr2) {
                    if (moleculeList[skipMol].crossbase[crossItr2] == reactMol1.index)
                        moleculeList[skipMol].probvec[crossItr2] = 0;
                }
            }
        }

        // Update the crossed molecule lists so that the current molecules won't avoid anything, but others will.
        reactCom1.ncross = -1;
        reactMol1.crossbase.clear();
        //mol2 is the implicit lipid.
    } else {
        //Mol1 is implicit lipid.
        // write temporary to real coords and clear temporary coordinates
        for (auto memMol : reactCom2.memberList) {
            moleculeList[memMol].comCoord = moleculeList[memMol].tmpComCoord;
            for (unsigned int i { 0 }; i < moleculeList[memMol].interfaceList.size(); ++i)
                moleculeList[memMol].interfaceList[i].coord = moleculeList[memMol].tmpICoords[i];
            moleculeList[memMol].clear_tmp_association_coords();
            moleculeList[memMol].trajStatus = TrajStatus::propagated;
        }
        //add link to surface for the molecule
        moleculeList[reactMol2.index].linksToSurface++;
        // update complexes
        reactCom2.iLipidIndex = reactMol1.index; //index in molecule list of the implicit lipid.
        reactCom2.linksToSurface++;
        reactCom2.update_properties(moleculeList, molTemplateList); // recalculate the properties of the first complex
        reactCom2.OnSurface = true; // the complex is on membrane now.
        reactCom2.D.z = 0;
        /*Clear implicit lipid properties*/
        moleculeList[reactMol1.index].clear_tmp_association_coords();
        reactCom1.memberList.clear();
        reactCom1.memberList.push_back(reactMol1.index);
        // std::cout << " Size of IL's complex:" << reactCom1.memberList.size() << " Interfaces on IL: " << moleculeList[reactMol1.index].interfaceList.size() << std::endl;
        reactCom1.update_properties(moleculeList, molTemplateList); // recalculate the properties of the second complex

        //Enforce boundary conditions
        reflect_complex_rad_rot(membraneObject, reactCom2, moleculeList, RS3D);
        //------------------------START UPDATE MONOMERLIST-------------------------
        // update oneTemp.monomerList when oneTemp.canDestroy is true and mol is monomer
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
        reactMol2.interfaceList[ifaceIndex2].interaction.partnerIndex = reactMol1.index;
        reactMol2.interfaceList[ifaceIndex2].interaction.partnerIfaceIndex = ifaceIndex1;
        if (currRxn.isReversible) {
            reactMol2.interfaceList[ifaceIndex2].interaction.conjBackRxn = currRxn.conjBackRxnIndex;
        }
        reactMol2.interfaceList[ifaceIndex2].isBound = true;
        //**************************************************
        //**************************************************
        reactMol2.interfaceList[ifaceIndex2].index = currRxn.productListNew[0].absIfaceIndex;

        // add to the list of bound interfaces and remove from the list of free interfaces
        reactMol2.bndlist.push_back(ifaceIndex2);
        reactMol2.bndpartner.push_back(reactMol1.index);
        {
            size_t tmpItr { 0 };
            for (unsigned i { 0 }; i < reactMol2.freelist.size(); ++i) {
                if (reactMol2.freelist[i] == ifaceIndex2)
                    tmpItr = i;
            }
            reactMol2.freelist[tmpItr] = reactMol2.freelist.back();
            reactMol2.freelist.pop_back();
        }

        // Set probability of this protein to zero in all reactions so it doesn't try to
        // react again but the partners still will avoid overlapping.
        for (unsigned crossItr { 0 }; crossItr < reactMol2.crossbase.size(); ++crossItr) {
            int skipMol { reactMol2.crossbase[crossItr] };
            if (!moleculeList[skipMol].isImplicitLipid) {
                for (unsigned crossItr2 { 0 }; crossItr2 < moleculeList[skipMol].crossbase.size(); ++crossItr2) {
                    if (moleculeList[skipMol].crossbase[crossItr2] == reactMol2.index)
                        moleculeList[skipMol].probvec[crossItr2] = 0;
                }
            }
        }

        // Update the crossed molecule lists so that the current molecules won't avoid anything, but others will.
        reactCom2.ncross = -1;
        reactMol2.crossbase.clear();
    } //reactmol1 is implicit lipid

    //update the number of bound species
    update_Nboundpairs(reactMol1.molTypeIndex, reactMol2.molTypeIndex, 1, params, counterArrays);

    //Update species copy numbers
    counterArrays.copyNumSpecies[currRxn.reactantListNew[0].absIfaceIndex] -= 1; // decrement ifaceIndex1
    counterArrays.copyNumSpecies[currRxn.reactantListNew[1].absIfaceIndex] -= 1; // decrement ifaceIndex2
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

    // update free number of lipids of IL for each state
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

    // species tracking here
    if (currRxn.isObserved) {
        auto obsItr = observablesList.find(currRxn.observeLabel);
        if (obsItr != observablesList.end())
            ++obsItr->second;
    }
}
