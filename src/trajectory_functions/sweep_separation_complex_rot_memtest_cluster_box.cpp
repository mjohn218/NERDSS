#include "boundary_conditions/reflect_functions.hpp"
#include "classes/class_Cluster.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

void sweep_separation_complex_rot_memtest_cluster_box(int simItr, int pro1Index, Parameters& params, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns, const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject)
{
    // try cluster sweep first
    int comIndex1 { moleculeList[pro1Index].myComIndex };
    //determine RS3Dinput
    double RS3Dinput { 0.0 };
    Complex targCom { complexList[comIndex1] };
    for (auto& molIndex : targCom.memberList) {
        for (int RS3Dindex = 0; RS3Dindex < 100; RS3Dindex++) {
            if (std::abs(membraneObject.RS3Dvect[RS3Dindex + 400] - moleculeList[molIndex].molTypeIndex) < 1E-2) {
                RS3Dinput = membraneObject.RS3Dvect[RS3Dindex + 300];
                break;
            }
        }
    }

    std::vector<TrajStatus> movestatOrig;
    std::vector<ClusterPair> pairList;
    define_cluster_pairs(pro1Index, moleculeList, complexList, forwardRxns, pairList);

    // for the partner of the member in pairList'complex, used when nStep > 1
    std::vector<TrajStatus> movestatOrigPartnerOfPairList;
    std::vector<int> partnerOfPairList;

    int p1, p2, i1, i2, k1, k2;
    int i;
    // store original movestats
    movestatOrig.reserve(pairList.size() * 2);
    for (i = 0; i < pairList.size(); i++) {
        p1 = pairList[i].p1;
        movestatOrig.push_back(moleculeList[p1].trajStatus);
        p2 = pairList[i].p2;
        movestatOrig.push_back(moleculeList[p2].trajStatus);
    }

    // double df1, df2, df3;
    double sigma;
    int it = 0;
    int maxit = 100;
    int saveit = 0;
    int t = 0;
    int flag = 0;
    int tsave = 0;
    int flags = 0;
    double x0;
    double y0;
    double z0;

    double x02;
    double y02;
    double z02;
    int maxPairs = 30;
    int nStep = 1;
    int ns;
    double timeStepOrg = params.timeStep;

    // std::cout << "In sweep_separation: pairList.size = " << pairList.size() << std::endl;

    if (pairList.size() > maxPairs) {
        nStep = 10;
        params.timeStep = params.timeStep / (1.0 * nStep); //perform repeated smaller updates, to converge positions.

        // resample the traj for the timeStep changed
        int currStop = pairList.size(); // resample trajectories for all proteins that need keep moving
        resample_traj(currStop, pairList, moleculeList, complexList, params, membraneObject, RS3Dinput);

        // find the partner of pairList'complex'member, this is for the propagate of these molecule continuely when there is subStep
        for (i = 0; i < pairList.size(); i++) {
            p1 = pairList[i].p1;
            p2 = pairList[i].p2;
            k1 = pairList[i].k1;
            k2 = pairList[i].k2;
            for (unsigned memMolItr { 0 }; memMolItr < complexList[k1].memberList.size(); ++memMolItr) {
                int pro1Index = complexList[k1].memberList[memMolItr];
                for (int crossMemItr { 0 }; crossMemItr < moleculeList[pro1Index].crossbase.size(); ++crossMemItr) {
                    int pro2Index { moleculeList[pro1Index].crossbase[crossMemItr] };
                    if (moleculeList[pro2Index].isImplicitLipid)
                        continue;
                    partnerOfPairList.push_back(pro2Index);
                    movestatOrigPartnerOfPairList.push_back(moleculeList[pro2Index].trajStatus);
                }
            }
            for (unsigned memMolItr { 0 }; memMolItr < complexList[k2].memberList.size(); ++memMolItr) {
                int pro1Index = complexList[k2].memberList[memMolItr];
                for (int crossMemItr { 0 }; crossMemItr < moleculeList[pro1Index].crossbase.size(); ++crossMemItr) {
                    int pro2Index { moleculeList[pro1Index].crossbase[crossMemItr] };
                    if (moleculeList[pro2Index].isImplicitLipid)
                        continue;
                    partnerOfPairList.push_back(pro2Index);
                    movestatOrigPartnerOfPairList.push_back(moleculeList[pro2Index].trajStatus);
                }
            }
        }

        // resample the traj in partnerOfPairList for the timeStep changed
        {
            int i;
            int k;
            int p;
            int flag;
            std::vector<int> didMove;
            for (i = 0; i < partnerOfPairList.size(); i++) {
                p = partnerOfPairList[i];
                k = moleculeList[p].myComIndex;
                if (moleculeList[p].trajStatus == TrajStatus::none || moleculeList[p].trajStatus == TrajStatus::canBeResampled) {
                    flag = 0;
                    for (int d = 0; d < didMove.size(); d++) {
                        if (didMove[d] == k)
                            flag = 1;
                    }
                    if (flag == 0) {

                        if (membraneObject.isSphere == true && complexList[k].D.z < 1E-15) { // complex on sphere surface
                            Coord targTrans = create_complex_propagation_vectors_on_sphere(params, complexList[k1]);
                            complexList[k].trajTrans.x = targTrans.x;
                            complexList[k].trajTrans.y = targTrans.y;
                            complexList[k].trajTrans.z = targTrans.z;
                            complexList[k].trajRot.x = sqrt(2.0 * params.timeStep * complexList[k].Dr.x) * GaussV();
                            complexList[k].trajRot.y = sqrt(2.0 * params.timeStep * complexList[k].Dr.y) * GaussV();
                            complexList[k].trajRot.z = sqrt(2.0 * params.timeStep * complexList[k].Dr.z) * GaussV();
                        } else {
                            complexList[k].trajTrans.x = sqrt(2.0 * params.timeStep * complexList[k].D.x) * GaussV();
                            complexList[k].trajTrans.y = sqrt(2.0 * params.timeStep * complexList[k].D.y) * GaussV();
                            complexList[k].trajTrans.z = sqrt(2.0 * params.timeStep * complexList[k].D.z) * GaussV();
                            complexList[k].trajRot.x = sqrt(2.0 * params.timeStep * complexList[k].Dr.x) * GaussV();
                            complexList[k].trajRot.y = sqrt(2.0 * params.timeStep * complexList[k].Dr.y) * GaussV();
                            complexList[k].trajRot.z = sqrt(2.0 * params.timeStep * complexList[k].Dr.z) * GaussV();
                        }

                        reflect_traj_complex_rad_rot(params, moleculeList, complexList[k], membraneObject, RS3Dinput);

                        didMove.push_back(k);
                    }
                }
            }
        }

        // for(int index = 0;index<pairList.size();index++){
        //     resampleList[pairList[index].k1]=1;
        //     resampleList[pairList[index].k2]=1;
        // }
    }
    for (ns = 0; ns < nStep; ns++) {
        // std::cout << "In for loop: ns = " << ns << std::endl;
        // if (nStep > 1 && ns == 0) { // resample the trajectory for the timeStep changed
        //     for (auto i : resampleList) {
        //         i = 0;
        //     }
        // }

        t = 0;
        flag = 1;
        it = 0;
        i = 0;
        while (i < pairList.size() && it < maxit) {
            // std::cout << "In while loop: i = " << i << std::endl;
            p1 = pairList[i].p1;
            p2 = pairList[i].p2;
            i1 = pairList[i].i1;
            i2 = pairList[i].i2;
            k1 = pairList[i].k1;
            k2 = pairList[i].k2;

            // if (resampleList[k1] == 0) {
            //     // resample complex k1
            //     complexList[k1].trajTrans.x = sqrt(2.0 * params.timeStep * complexList[k1].D.x) * GaussV();
            //     complexList[k1].trajTrans.y = sqrt(2.0 * params.timeStep * complexList[k1].D.y) * GaussV();
            //     complexList[k1].trajTrans.z = sqrt(2.0 * params.timeStep * complexList[k1].D.z) * GaussV();
            //     complexList[k1].trajRot.x = sqrt(2.0 * params.timeStep * complexList[k1].Dr.x) * GaussV();
            //     complexList[k1].trajRot.y = sqrt(2.0 * params.timeStep * complexList[k1].Dr.y) * GaussV();
            //     complexList[k1].trajRot.z = sqrt(2.0 * params.timeStep * complexList[k1].Dr.z) * GaussV();

            //     reflect_traj_complex_rad_rot(params, moleculeList, complexList[k1], membraneObject, RS3Dinput);
            //     resampleList[k1] = 1;
            // }

            // if (resampleList[k2] == 0) {
            //     // resample complex k2
            //     complexList[k2].trajTrans.x = sqrt(2.0 * params.timeStep * complexList[k2].D.x) * GaussV();
            //     complexList[k2].trajTrans.y = sqrt(2.0 * params.timeStep * complexList[k2].D.y) * GaussV();
            //     complexList[k2].trajTrans.z = sqrt(2.0 * params.timeStep * complexList[k2].D.z) * GaussV();
            //     complexList[k2].trajRot.x = sqrt(2.0 * params.timeStep * complexList[k2].Dr.x) * GaussV();
            //     complexList[k2].trajRot.y = sqrt(2.0 * params.timeStep * complexList[k2].Dr.y) * GaussV();
            //     complexList[k2].trajRot.z = sqrt(2.0 * params.timeStep * complexList[k2].Dr.z) * GaussV();

            //     reflect_traj_complex_rad_rot(params, moleculeList, complexList[k2], membraneObject, RS3Dinput);
            //     resampleList[k2] = 1;
            // }

            // calculate separation
            if (moleculeList[p2].isImplicitLipid) {
                i++;
                // std::cout << "i++: " << i << std::endl;
                continue;
            }

            Vector iface1Vec { moleculeList[p1].interfaceList[i1].coord - complexList[k1].comCoord };
            std::array<double, 9> M = create_euler_rotation_matrix(complexList[k1].trajRot);
            iface1Vec = matrix_rotate(iface1Vec, M);

            double dx1 { complexList[k1].comCoord.x + iface1Vec.x + complexList[k1].trajTrans.x };
            double dy1 { complexList[k1].comCoord.y + iface1Vec.y + complexList[k1].trajTrans.y };
            double dz1 { complexList[k1].comCoord.z + iface1Vec.z + complexList[k1].trajTrans.z };

            Vector iface2Vec { moleculeList[p2].interfaceList[i2].coord - complexList[k2].comCoord };
            std::array<double, 9> M2 = create_euler_rotation_matrix(complexList[k2].trajRot);
            iface2Vec = matrix_rotate(iface2Vec, M2);
            double dx2 { complexList[k2].comCoord.x + iface2Vec.x + complexList[k2].trajTrans.x };
            double dy2 { complexList[k2].comCoord.y + iface2Vec.y + complexList[k2].trajTrans.y };
            double dz2 { complexList[k2].comCoord.z + iface2Vec.z + complexList[k2].trajTrans.z };

            /*separation*/
            double df1 { dx1 - dx2 };
            double df2 { dy1 - dy2 };
            double df3 { dz1 - dz2 };

            double dr2 = (df1 * df1) + (df2 * df2);
            if (pairList[i].memtest != 1)
                dr2 += (df3 * df3);

            if (dr2 < pairList[i].bindrad * pairList[i].bindrad) {
                int currStop = i + 1;
                resample_traj(currStop, pairList, moleculeList, complexList, params, membraneObject, RS3Dinput);
                it++;
                i = 0;
                // std::cout << "it++: " << it << std::endl;
                // std::cout << "reset i: " << i << std::endl;
            } else {
                i++;
                // std::cout << "i++: " << i << std::endl;
            }
        } // end looping over all pairs

        // exit(1);

        if (it == maxit) {
            // Tried cluster sweep and it failed. Now try pairs.
            // std::cout << " WARNING ***************************************************** " << '\n';
            // std::cout << "At " << (simItr - params.itrRestartFrom) * timeStepOrg * 1E-6 << " s,";
            // std::cout << "cluster sweep can't solve overlap. Now try pairs. " << '\n';
            for (i = 0; i < pairList.size(); i++) {
                p1 = pairList[i].p1;
                p2 = pairList[i].p2;
                if (moleculeList[p1].trajStatus == TrajStatus::none || moleculeList[p1].trajStatus == TrajStatus::canBeResampled) {
                    sweep_separation_complex_rot_memtest_box(simItr, p1, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);
                }
                if (moleculeList[p2].trajStatus == TrajStatus::none || moleculeList[p2].trajStatus == TrajStatus::canBeResampled) {
                    sweep_separation_complex_rot_memtest_box(simItr, p2, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);
                }
            }
        }

        // now physically move all these proteins
        for (i = 0; i < pairList.size(); i++) {
            p1 = pairList[i].p1;
            k1 = pairList[i].k1;
            p2 = pairList[i].p2;
            k2 = pairList[i].k2;
            if (moleculeList[p1].trajStatus == TrajStatus::none || moleculeList[p1].trajStatus == TrajStatus::canBeResampled) {
                complexList[k1].propagate(moleculeList, membraneObject, molTemplateList);
            }
            if (moleculeList[p2].trajStatus == TrajStatus::none || moleculeList[p2].trajStatus == TrajStatus::canBeResampled) {
                complexList[k2].propagate(moleculeList, membraneObject, molTemplateList);
            }
        }

        // if the timeStep changed, also need move physically the molecule in partnerOfPairList
        if (nStep > 1) {
            for (i = 0; i < partnerOfPairList.size(); i++) {
                int p = partnerOfPairList[i];
                int k = moleculeList[p].myComIndex;
                if (moleculeList[p].trajStatus == TrajStatus::none || moleculeList[p].trajStatus == TrajStatus::canBeResampled) {
                    complexList[k].propagate(moleculeList, membraneObject, molTemplateList);
                }
            }
        }

        if (ns < nStep - 1) {
            // we have broken up the full time_step into multiple smaller steps. need these proteins to continue moving in the next substep
            for (i = 0; i < pairList.size(); i++) {
                p1 = pairList[i].p1;
                p2 = pairList[i].p2;
                if (movestatOrig[i * 2] == TrajStatus::none || movestatOrig[i * 2] == TrajStatus::canBeResampled) {
                    moleculeList[p1].trajStatus = movestatOrig[i * 2]; // keep p1 moving
                }
                if (movestatOrig[i * 2 + 1] == TrajStatus::none || movestatOrig[i * 2 + 1] == TrajStatus::canBeResampled) {
                    moleculeList[p2].trajStatus = movestatOrig[i * 2 + 1]; // keep p2 moving
                }
            }
            int currStop = pairList.size(); // resample trajectories for all proteins that need keep moving
            resample_traj(currStop, pairList, moleculeList, complexList, params, membraneObject, RS3Dinput);

            for (i = 0; i < partnerOfPairList.size(); i++) {
                int p = partnerOfPairList[i];
                if (movestatOrigPartnerOfPairList[i] == TrajStatus::none || movestatOrigPartnerOfPairList[i] == TrajStatus::canBeResampled) {
                    moleculeList[p].trajStatus = movestatOrigPartnerOfPairList[i]; // keep moving
                }
            }
            // resample
            {
                int i;
                int k;
                int p;
                int flag;
                std::vector<int> didMove;
                for (i = 0; i < partnerOfPairList.size(); i++) {
                    p = partnerOfPairList[i];
                    k = moleculeList[p].myComIndex;
                    if (moleculeList[p].trajStatus == TrajStatus::none || moleculeList[p].trajStatus == TrajStatus::canBeResampled) {
                        flag = 0;
                        for (int d = 0; d < didMove.size(); d++) {
                            if (didMove[d] == k)
                                flag = 1;
                        }
                        if (flag == 0) {

                            if (membraneObject.isSphere == true && complexList[k].D.z < 1E-15) { // complex on sphere surface
                                Coord targTrans = create_complex_propagation_vectors_on_sphere(params, complexList[k1]);
                                complexList[k].trajTrans.x = targTrans.x;
                                complexList[k].trajTrans.y = targTrans.y;
                                complexList[k].trajTrans.z = targTrans.z;
                                complexList[k].trajRot.x = sqrt(2.0 * params.timeStep * complexList[k].Dr.x) * GaussV();
                                complexList[k].trajRot.y = sqrt(2.0 * params.timeStep * complexList[k].Dr.y) * GaussV();
                                complexList[k].trajRot.z = sqrt(2.0 * params.timeStep * complexList[k].Dr.z) * GaussV();
                            } else {
                                complexList[k].trajTrans.x = sqrt(2.0 * params.timeStep * complexList[k].D.x) * GaussV();
                                complexList[k].trajTrans.y = sqrt(2.0 * params.timeStep * complexList[k].D.y) * GaussV();
                                complexList[k].trajTrans.z = sqrt(2.0 * params.timeStep * complexList[k].D.z) * GaussV();
                                complexList[k].trajRot.x = sqrt(2.0 * params.timeStep * complexList[k].Dr.x) * GaussV();
                                complexList[k].trajRot.y = sqrt(2.0 * params.timeStep * complexList[k].Dr.y) * GaussV();
                                complexList[k].trajRot.z = sqrt(2.0 * params.timeStep * complexList[k].Dr.z) * GaussV();
                            }

                            reflect_traj_complex_rad_rot(params, moleculeList, complexList[k], membraneObject, RS3Dinput);

                            didMove.push_back(k);
                        }
                    }
                }
            }
        }
    } // end looping over Nsteps

    // Done moving, set all trajStatus to propagated
    for (i = 0; i < pairList.size(); i++) {
        p1 = pairList[i].p1;
        p2 = pairList[i].p2;

        moleculeList[p1].trajStatus = TrajStatus::propagated;
        moleculeList[p2].trajStatus = TrajStatus::propagated;
        complexList[moleculeList[p1].myComIndex].trajTrans.zero_crds();
        complexList[moleculeList[p1].myComIndex].trajRot.zero_crds();
        complexList[moleculeList[p2].myComIndex].trajTrans.zero_crds();
        complexList[moleculeList[p2].myComIndex].trajRot.zero_crds();
    }
    for (i = 0; i < partnerOfPairList.size(); i++) {
        int p = partnerOfPairList[i];
        moleculeList[p].trajStatus = TrajStatus::propagated;
        complexList[moleculeList[p].myComIndex].trajTrans.zero_crds();
        complexList[moleculeList[p].myComIndex].trajRot.zero_crds();
    }

    params.timeStep = timeStepOrg; // recover timeStep
}