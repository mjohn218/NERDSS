#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "tracing.hpp"
#include <sstream>

void determine_2D_bimolecular_reaction_probability(int simItr, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    unsigned& DDTableIndex, double* tableIDs, BiMolData& biMolData, const Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, Membrane& membraneObject, std::vector<gsl_matrix*>& normMatrices,
    std::vector<gsl_matrix*>& survMatrices, std::vector<gsl_matrix*>& pirMatrices)
{
    // TRACE();
    double Dr1 {};
    {
        double cf { cos(sqrt(2.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep)) };
        Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
    }

    double Dr2 {};
    {
        double cf = cos(sqrt(2.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
        Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
    }

    biMolData.Dtot += (Dr1 + Dr2) / (4.0 * params.timeStep); // add in contributions from rotation

    {
        // Only allow 2D. diffusion at certain intervals, to avoid generating too many 2D. Tables
        // Keep only one sig fig for <0.1, 2 for 0.1<d<10, 3 for 10<d<100, etc
        double dtmp;
        if (biMolData.Dtot < 0.0001)
            dtmp = biMolData.Dtot * 100000;
        else if (biMolData.Dtot < 0.001)
            dtmp = biMolData.Dtot * 10000;
        else if (biMolData.Dtot < 0.01)
            dtmp = biMolData.Dtot * 1000;
        else if (biMolData.Dtot < 0.1)
            dtmp = biMolData.Dtot * 100;
        else
            dtmp = biMolData.Dtot * 100;

        int d_ones = int(round(dtmp));

        if (biMolData.Dtot < 0.0001)
            biMolData.Dtot = d_ones * 0.00001;
        else if (biMolData.Dtot < 0.001)
            biMolData.Dtot = d_ones * 0.0001;
        else if (biMolData.Dtot < 0.01)
            biMolData.Dtot = d_ones * 0.001;
        else if (biMolData.Dtot < 0.1)
            biMolData.Dtot = d_ones * 0.01;
        else
            biMolData.Dtot = d_ones * 0.01;

        if (biMolData.Dtot < 1E-50)
            biMolData.Dtot = 0;
    }

    double RMax { 3.5 * sqrt(4.0 * biMolData.Dtot * params.timeStep) + forwardRxns[rxnIndex].bindRadius };
    double sep {};
    double R1 {};
    bool withinRmax { get_distance(biMolData.pro1Index, biMolData.pro2Index, biMolData.relIface1, biMolData.relIface2,
        rxnIndex, rateIndex, isStateChangeBackRxn, sep, R1, RMax, complexList, forwardRxns[rxnIndex], moleculeList, membraneObject) };
    if (withinRmax) {
        // in case they dissociated
        moleculeList[biMolData.pro1Index].probvec.push_back(0);
        moleculeList[biMolData.pro2Index].probvec.push_back(0);
    }

    if (moleculeList[biMolData.pro1Index].trajStatus != TrajStatus::propagated
        && moleculeList[biMolData.pro2Index].trajStatus != TrajStatus::propagated) {
        /*This movestat check is if you allow just dissociated proteins to avoid
         * overlap*/
        if (withinRmax && forwardRxns[rxnIndex].rateList[rateIndex].rate > 0) {
            bool probValExists { false };
            int probMatrixIndex { 0 };
            // get_distance, because ncross is > 0 iff R1 < RMax
            /*Evaluate probability of reaction, with reweighting*/
            // Generate 2D tables unless they were not before
            double ktemp { forwardRxns[rxnIndex].rateList[rateIndex].rate / forwardRxns[rxnIndex].length3Dto2D };
            if (forwardRxns[rxnIndex].isSymmetric == true)
                ktemp *= 2.0; // for A(a)+A(a)->A(a!).A(a!) case

            for (int l = 0; l < DDTableIndex; ++l) {
                if (std::abs(tableIDs[l] - ktemp) < 1e-8 && std::abs(tableIDs[params.max2DRxns + l] - biMolData.Dtot) < 1E-4) {
                    probValExists = true;
                    probMatrixIndex = l;
                    break;
                }
            }

            if (!probValExists) {
                // first dimension (+0*params.max2DRxns)
                tableIDs[DDTableIndex] = ktemp;
                // second dimension (+1*params.max2DRxns)
                tableIDs[DDTableIndex + params.max2DRxns] = biMolData.Dtot;
                size_t veclen { size_lookup(forwardRxns[rxnIndex].bindRadius, biMolData.Dtot, params, RMax) };
                // std::cout << "Create new 2D table: " << ktemp << ", Dtot: " << biMolData.Dtot << " size: " << veclen
                //           << '\n';
                survMatrices.resize(DDTableIndex + 1);
                normMatrices.resize(DDTableIndex + 1);
                pirMatrices.resize(DDTableIndex + 1);
                survMatrices[DDTableIndex] = gsl_matrix_alloc(2, veclen);
                normMatrices[DDTableIndex] = gsl_matrix_alloc(2, veclen);
                pirMatrices[DDTableIndex] = gsl_matrix_alloc(veclen, veclen);

                create_DDMatrices(survMatrices[DDTableIndex], normMatrices[DDTableIndex], pirMatrices[DDTableIndex],
                    forwardRxns[rxnIndex].bindRadius, biMolData.Dtot, RMax, ktemp, params);
                probMatrixIndex = DDTableIndex;
                DDTableIndex += 1;
                if (DDTableIndex == params.max2DRxns) {
                    std::cout << "You have hit the maximum number of unique 2D reactions "
                                 "allowed: "
                              << params.max2DRxns << '\n';
                    std::cout << "terminating...." << '\n';
                    exit(1);
                }
            }
            probValExists = false; // reset

            double ratio = forwardRxns[rxnIndex].bindRadius / R1;
            if (sep < 0) {
                if (biMolData.com1Index != biMolData.com2Index) {
                    // && i1!=0 && i2!=0) {
                    // std::cout << "*****************************************************\n"
                    //           << " WARNING AT ITERATION " << simItr << "\n";
                    // std::cout << "SEPARATION BETWEEN INTERFACE " << biMolData.relIface1 << " ON MOLECULE "
                    //           << biMolData.pro1Index << " AND INTERFACE " << biMolData.relIface2 << " ON MOLECULE "
                    //           << biMolData.pro2Index << " IS LESS THAN 0\n";
                    // std::cout << "separation: " << sep << " r1: " << R1 << " p1: " << biMolData.pro1Index
                    //           << " p2: " << biMolData.pro2Index << " it " << simItr << " i1: " << biMolData.absIface1
                    //           << " i2: " << biMolData.absIface2 << '\n';
                    // std::cout << "MOL1 COM: " << moleculeList[biMolData.pro1Index].comCoord
                    //           << " freelist.size(): " << moleculeList[biMolData.pro1Index].freelist.size() << '\n';
                    // std::cout << "IFACE1: "
                    //           << moleculeList[biMolData.pro1Index].interfaceList[biMolData.relIface1].coord << '\n';
                    // std::cout << "MOL2 COM: " << moleculeList[biMolData.pro2Index].comCoord
                    //           << " freelist.size(): " << moleculeList[biMolData.pro2Index].freelist.size() << '\n';
                    // std::cout << "IFACE2: "
                    //           << moleculeList[biMolData.pro2Index].interfaceList[biMolData.relIface2].coord << '\n';
                    // std::cout << "*****************************************************\n";
                } else {
                    // std::cout << " WARNING "
                    //              "***************************************************** "
                    //           << '\n';
                    // std::cout << "Protein interfaces within a complex, trying to bind. "
                    //              "Separation <0: "
                    //           << sep << " r1 " << R1 << " p1: " << biMolData.pro1Index << " p2: " << biMolData.pro2Index
                    //           << " it " << simItr << " i1: " << biMolData.absIface1 << " i2: " << biMolData.absIface2
                    //           << '\n';
                    // std::cout << "P1 COORDS: " << moleculeList[biMolData.pro1Index].comCoord
                    //           << " freelist.size(): " << moleculeList[biMolData.pro1Index].freelist.size() << '\n';
                    // std::cout << "P2 COORDS: " << moleculeList[biMolData.pro2Index].comCoord
                    //           << " freelist.size(): " << moleculeList[biMolData.pro2Index].freelist.size() << '\n';
                }
                sep = 0;
                ratio = 1;
                R1 = forwardRxns[rxnIndex].bindRadius;
            }
            double currnorm { 1.0 };
            double p0_ratio { 1.0 };

            /*protein pro1 is molTypeIndex wprot and pro2 is wprot2*/
            int proA = biMolData.pro1Index;
            int ifaceA = biMolData.relIface1;
            int proB = biMolData.pro2Index;
            int ifaceB = biMolData.relIface2;

            if (biMolData.pro1Index > biMolData.pro2Index) {
                proA = biMolData.pro2Index;
                ifaceA = biMolData.relIface2;
                proB = biMolData.pro1Index;
                ifaceB = biMolData.relIface1;
            }

            double rxnProb {};
            for (int s { 0 }; s < moleculeList[proA].prevlist.size(); ++s) {
                if (moleculeList[proA].prevlist[s] == proB && moleculeList[proA].prevmyface[s] == ifaceA
                    && moleculeList[proA].prevpface[s] == ifaceB) {
                    if (moleculeList[proA].prevsep[s] >= RMax) {
                        // BEcause previous reweighting was for 3D, now
                        p0_ratio = 1.0;
                        // restart reweighting in 2D.
                        currnorm = 1.0;
                    } else {
                        p0_ratio = DDpirr_pfree_ratio_ps(pirMatrices[probMatrixIndex], survMatrices[probMatrixIndex],
                            normMatrices[probMatrixIndex], R1, biMolData.Dtot, params.timeStep,
                            moleculeList[proA].prevsep[s], moleculeList[proA].ps_prev[s], 1E-10,
                            forwardRxns[rxnIndex].bindRadius);
                        currnorm = moleculeList[proA].prevnorm[s] * p0_ratio;
                    }
                    break;
                }
            }
            rxnProb = get_prevSurv(
                survMatrices[probMatrixIndex], biMolData.Dtot, params.timeStep, R1, forwardRxns[rxnIndex].bindRadius);
            moleculeList[biMolData.pro1Index].probvec.back() = rxnProb * currnorm;
            moleculeList[biMolData.pro2Index].probvec.back() = rxnProb * currnorm;
            if (rxnProb > 1.000001) {
                std::cerr << "Error: prob of reaction is: " << rxnProb << " > 1. Avoid this using a smaller time step." << std::endl;
                exit(1);
            }
            if (rxnProb > 0.5) {
                // std::cout << "WARNING: prob of reaction > 0.5. If this is a reaction for a bimolecular binding with multiple binding sites, please use a smaller time step." << std::endl;
            }

            //------------------Calculate 2D bimolecular prob---------------------
            // if (std::abs(biMolData.Dtot - 0.32) < 1e-10) {
            //     std::cout << "RMax: " << RMax << std::endl;
            //     std::cout << "Dtot: " << biMolData.Dtot << std::endl;
            //     std::cout << "Dr1: " << Dr1 << std::endl;
            //     std::cout << "Dr2: " << Dr2 << std::endl;
            //     std::cout << "Com1: " << std::endl;
            //     for (auto& tmp : complexList[biMolData.com1Index].memberList) {
            //         std::cout << tmp << "\t";
            //     }
            //     std::cout << std::endl;
            //     std::cout << "Com2: " << std::endl;
            //     for (auto& tmp : complexList[biMolData.com2Index].memberList) {
            //         std::cout << tmp << "\t";
            //     }
            //     std::cout << std::endl;
            //     double z = 0.0;
            //     int total = 10000;
            //     double zmin = forwardRxns[rxnIndex].bindRadius;
            //     double zmax = RMax;
            //     double prob = 0.0;
            //     std::ofstream probout("prob.dat");
            //     for (int i = 1; i < total; i++) {
            //         z = zmin + (zmax - zmin) * i / total;
            //         prob = get_prevSurv(
            //             survMatrices[probMatrixIndex], biMolData.Dtot, params.timeStep, z, forwardRxns[rxnIndex].bindRadius);
            //         probout << z << "\t" << prob << std::endl;
            //     }
            //     probout.close();
            //     exit(1);
            // }
            //------------------------------------------------------------------------

            /*Store all the reweighting numbers for next step.*/
            // below used to be moleculeList[proA].vector[s] = value
            //                                            s =
            //                                            moleculeList[proA].currlist.size();
            moleculeList[proA].currprevsep.push_back(R1);
            moleculeList[proA].currlist.push_back(proB);
            moleculeList[proA].currmyface.push_back(ifaceA);
            moleculeList[proA].currpface.push_back(ifaceB);
            moleculeList[proA].currprevnorm.push_back(currnorm);
            moleculeList[proA].currps_prev.push_back(1.0 - rxnProb * currnorm);
        } // Within reaction zone
    }
}
