#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "tracing.hpp"

void determine_3D_implicitlipid_reaction_probability(int simItr, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    BiMolData& biMolData, const Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, Membrane& membraneObject, const int& relStateIndex)
{
    // TRACE();
    /*3D reaction*/
    double Dr1 {};
    if (std::abs(complexList[biMolData.com1Index].D.z - 0) < 1E-15) {
        double cf = cos(sqrt(2.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep));
        Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
        biMolData.Dtot += Dr1 / (4.0 * params.timeStep);
    } else {
        double cf = cos(sqrt(4.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep));
        Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
        biMolData.Dtot += Dr1 / (6.0 * params.timeStep);
    }

    double Dr2;
    if (std::abs(complexList[biMolData.com2Index].D.z - 0) < 1E-15) {
        double cf = cos(sqrt(2.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
        Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
        biMolData.Dtot += Dr2 / (4.0 * params.timeStep);
    } else {
        double cf = cos(sqrt(4.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
        Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
        biMolData.Dtot += Dr2 / (6.0 * params.timeStep);
    }

    double Rmax { 3.0 * sqrt(6.0 * biMolData.Dtot * params.timeStep) + forwardRxns[rxnIndex].bindRadius };
    // std::cout << "Rmax: " << std::setprecision(20) << Rmax << std::endl;
    // exit(1);

    double sep {};
    double R1 {};
    bool withinRmax = get_distance_to_surface(biMolData.pro1Index, biMolData.pro2Index, biMolData.relIface1, biMolData.relIface2,
        rxnIndex, rateIndex, isStateChangeBackRxn, sep, R1, Rmax, complexList,
        forwardRxns[rxnIndex], moleculeList, membraneObject);

    if (withinRmax) {
        // in case the molecule dissociated
        moleculeList[biMolData.pro1Index].probvec.push_back(0);
        //moleculeList[biMolData.pro2Index].probvec.push_back(0);
    }

    // if (moleculeList[biMolData.pro1Index].trajStatus != TrajStatus::propagated && moleculeList[biMolData.pro2Index].trajStatus != TrajStatus::propagated)
    if (moleculeList[biMolData.pro1Index].trajStatus != TrajStatus::propagated) {
        // This movestat check is if you allow just dissociated proteins to avoid overlap
        if (withinRmax && forwardRxns[rxnIndex].rateList[rateIndex].rate > 0) {
            // declare intrinsic binding rate of 3D->2D case.
            double ktemp { 2.0 * forwardRxns[rxnIndex].rateList[rateIndex].rate };

            // Evaluate probability of reaction, implicit-lipid method doesn't need reweighting
            if (sep < 0) {
                // std::cout << "*****************************************************\n"
                //           << " WARNING AT ITERATION " << simItr << "\n";
                // std::cout << "SEPARATION BETWEEN INTERFACE " << biMolData.relIface1 << " ON MOLECULE "
                //           << biMolData.pro1Index << " AND THE MEMBRANE SURFACE IS LESS THAN 0\n";
                // std::cout << "separation: " << sep << " r1: " << R1 << " p1: " << biMolData.pro1Index
                //           << " it " << simItr << " i1: " << biMolData.absIface1 << '\n';
                // std::cout << "MOL1 COM: " << moleculeList[biMolData.pro1Index].comCoord
                //           << " freelist.size(): " << moleculeList[biMolData.pro1Index].freelist.size() << '\n';
                // std::cout << "IFACE1: "
                //           << moleculeList[biMolData.pro1Index].interfaceList[biMolData.relIface1].coord << '\n';
                // std::cout << "*****************************************************\n";
                sep = 0;
                R1 = forwardRxns[rxnIndex].bindRadius;
            }
            int proA = biMolData.pro1Index;
            int ifaceA = biMolData.relIface1;
            int proB = biMolData.pro2Index;
            int ifaceB = biMolData.relIface2;
            double currnorm { 1.0 };
            double rxnProb {};

            paramsIL params3D {};
            params3D.R2D = 0.0;
            params3D.sigma = forwardRxns[rxnIndex].bindRadius;
            params3D.Dtot = biMolData.Dtot;
            params3D.ka = ktemp;
            //params3D.kb = backRxns[rxnIndex].rateList[rateIndex].rate;
            params3D.area = membraneObject.totalSA;
            params3D.dt = params.timeStep;
            params3D.Nlipid = membraneObject.numberOfFreeLipidsEachState[relStateIndex];
            double rho = 1.0 * membraneObject.numberOfFreeLipidsEachState[relStateIndex] / params3D.area;
            double z = R1;
            rxnProb = rho * pimplicitlipid_3D(z, params3D);
            if (rxnProb > 1.000001) {
                std::cerr << "Error: prob of reaction is: " << rxnProb << " > 1. Avoid this using a smaller time step." << std::endl;
                exit(1);
            }
            if (rxnProb > 0.5) {
                // std::cout << "WARNING: prob of reaction > 0.5. If this is a reaction for a bimolecular binding with multiple binding sites, please use a smaller time step." << std::endl;
            }
            if (sep < 0)
                rxnProb = 1.0;

            moleculeList[biMolData.pro1Index].probvec.back() = rxnProb * currnorm;

            moleculeList[proA].currprevsep.push_back(R1);
            moleculeList[proA].currlist.push_back(proB);
            moleculeList[proA].currmyface.push_back(ifaceA);
            moleculeList[proA].currpface.push_back(ifaceB);
            moleculeList[proA].currprevnorm.push_back(currnorm);
            moleculeList[proA].currps_prev.push_back(1.0 - rxnProb * currnorm);

            // std::cout << "bind prob: " << std::setprecision(20) << rxnProb << std::endl;
            // std::cout << "rho: " << std::setprecision(20) << rho << std::endl;
            // std::cout << "z: " << std::setprecision(20) << z << std::endl;
            // std::cout << "params3D.sigma: " << std::setprecision(20) << params3D.sigma << std::endl;
            // std::cout << "params3D.Dtot: " << std::setprecision(20) << params3D.Dtot << std::endl;
            // std::cout << "params3D.ka: " << std::setprecision(20) << params3D.ka << std::endl;
            // std::cout << "params3D.area: " << std::setprecision(20) << params3D.area << std::endl;
            // std::cout << "params3D.dt: " << std::setprecision(20) << params3D.dt << std::endl;
            // exit(1);
        } // Within reaction zone
    }
}
