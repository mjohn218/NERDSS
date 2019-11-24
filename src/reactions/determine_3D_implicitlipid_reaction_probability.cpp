#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"

void determine_3D_implicitlipid_reaction_probability(int simItr, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    BiMolData& biMolData, const Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, Membrane &membraneObject)
{
    /*3D reaction*/
    double Dr1 {};
    if (complexList[biMolData.com1Index].D.z == 0) {
        double cf = cos(sqrt(2.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep));
        Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
        biMolData.Dtot += Dr1 / (4.0 * params.timeStep);
    } else {
        double cf = cos(sqrt(4.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep));
        Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
        biMolData.Dtot += Dr1 / (6.0 * params.timeStep);
    }
    
    double Dr2;
    if (complexList[biMolData.com2Index].D.z == 0) {
        double cf = cos(sqrt(2.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
        Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
        biMolData.Dtot += Dr2 / (4.0 * params.timeStep);
    } else {
        double cf = cos(sqrt(4.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
        Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
        biMolData.Dtot += Dr2 / (6.0 * params.timeStep);
    }
    
    double Rmax { 3.0 * sqrt(6.0 * biMolData.Dtot * params.timeStep) + forwardRxns[rxnIndex].bindRadius };

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
    if (moleculeList[biMolData.pro1Index].trajStatus != TrajStatus::propagated) 
    {
        // This movestat check is if you allow just dissociated proteins to avoid overlap
        if (withinRmax && forwardRxns[rxnIndex].rateList[rateIndex].rate > 0) {
        	   // declare intrinsic binding rate of 3D->2D case.
        	   double ktemp { 2.0 * forwardRxns[rxnIndex].rateList[rateIndex].rate };
        	   
            // Evaluate probability of reaction, implicit-lipid method doesn't need reweighting
            if (sep < 0) {
                std::cout << "*****************************************************\n"
                          << " WARNING AT ITERATION " << simItr << "\n";
                std::cout << "SEPARATION BETWEEN INTERFACE " << biMolData.relIface1 << " ON MOLECULE "
                          << biMolData.pro1Index << " AND THE MEMBRANE SURFACE IS LESS THAN 0\n";
                std::cout << "separation: " << sep << " r1: " << R1 << " p1: " << biMolData.pro1Index
                          << " it " << simItr << " i1: " << biMolData.absIface1 << '\n';
                std::cout << "MOL1 COM: " << moleculeList[biMolData.pro1Index].comCoord
                          << " freelist.size(): " << moleculeList[biMolData.pro1Index].freelist.size() << '\n';
                std::cout << "IFACE1: "
                          << moleculeList[biMolData.pro1Index].interfaceList[biMolData.relIface1].coord << '\n';
                std::cout << "*****************************************************\n";
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
            params3D.Nlipid = membraneObject.No_free_lipids; // the initial number of lipids' interfaces on the surface
            double rho = params3D.Nlipid/params3D.area; 
            double z = R1;
            rxnProb = rho * pimplicitlipid_3D(z,params3D);
            if (sep < 0)
               rxnProb = 1.0;
                  
            moleculeList[biMolData.pro1Index].probvec.back() = rxnProb * currnorm;

            moleculeList[proA].currprevsep.push_back(R1);
            moleculeList[proA].currlist.push_back(proB);
            moleculeList[proA].currmyface.push_back(ifaceA);
            moleculeList[proA].currpface.push_back(ifaceB);
            moleculeList[proA].currprevnorm.push_back(currnorm);
            moleculeList[proA].currps_prev.push_back(1.0 - rxnProb * currnorm);
        } // Within reaction zone
    } 
}
