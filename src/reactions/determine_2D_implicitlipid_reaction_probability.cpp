#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "tracing.hpp"

void determine_2D_implicitlipid_reaction_probability(int simItr, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    std::vector<double>& ILTableIDs, BiMolData& biMolData, const Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, std::vector<double>& IL2DbindingVec, std::vector<double>& IL2DUnbindingVec, Membrane& membraneObject, const int& relStateIndex)
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
    double sep = 0.0;
    double R1 = 0.0;

    // in case they dissociated
    moleculeList[biMolData.pro1Index].probvec.push_back(0);
    //moleculeList[biMolData.pro2Index].probvec.push_back(0);

    if (moleculeList[biMolData.pro1Index].trajStatus != TrajStatus::propagated) {
        // This movestat check is if you allow just dissociated proteins to avoid overlap
        if (forwardRxns[rxnIndex].rateList[rateIndex].rate > 0) {
            bool probValExists { false };
            int probMatrixIndex { 0 };

            // declare intrinsic binding rate of 2D->2D case.
            double ktemp { forwardRxns[rxnIndex].rateList[rateIndex].rate / forwardRxns[rxnIndex].length3Dto2D };
            //int backIndex = forwardRxns[rxnIndex].conjBackRxnIndex;
            double kb { 0 };

            if (forwardRxns[rxnIndex].isReversible == true) {
                kb = backRxns[forwardRxns[rxnIndex].conjBackRxnIndex].rateList[rateIndex].rate;
            }

            for (int l = 0; l < IL2DbindingVec.size(); ++l) {
                if (std::abs(ILTableIDs[l * 3] - ktemp) < 1e-8 && std::abs(ILTableIDs[l * 3 + 1] - biMolData.Dtot) < 1E-4 && std::abs(ILTableIDs[l * 3 + 2] - kb) < 1e-8) {
                    probValExists = true;
                    probMatrixIndex = l;
                    break;
                }
            }

            if (!probValExists) {
                // first dimension out of i elements (2*i)
                ILTableIDs.push_back(ktemp);
                // second dimension (2*i+1)
                ILTableIDs.push_back(biMolData.Dtot);

                ILTableIDs.push_back(kb);
                paramsIL params2D {};
                params2D.kb = kb;
                params2D.R2D = 0.0;
                params2D.sigma = forwardRxns[rxnIndex].bindRadius;
                params2D.Dtot = biMolData.Dtot;
                params2D.ka = ktemp;

                params2D.area = membraneObject.totalSA;
                params2D.dt = params.timeStep;
                params2D.Nlipid = membraneObject.numberOfFreeLipidsEachState[relStateIndex];
                params2D.Na = membraneObject.numberOfProteinEachState[relStateIndex]; // the initial number of protein's interfaces that can bind to surface

                probMatrixIndex = IL2DbindingVec.size();
                IL2DbindingVec.push_back(pimplicitlipid_2D(params2D));
            }
            probValExists = false; // reset

            int proA = biMolData.pro1Index;
            int ifaceA = biMolData.relIface1;
            int proB = biMolData.pro2Index;
            int ifaceB = biMolData.relIface2;
            double currnorm { 1.0 };
            double rho = 1.0 * membraneObject.numberOfFreeLipidsEachState[relStateIndex] / membraneObject.totalSA;

            double rxnProb = rho * IL2DbindingVec[probMatrixIndex];
            if (rxnProb > 1.000001) {
                std::cerr << "Error: prob of reaction is: " << rxnProb << " > 1. Avoid this using a smaller time step." << std::endl;
                exit(1);
            }
            if (rxnProb > 0.5) {
                // std::cout << "WARNING: prob of reaction > 0.5. If this is a reaction for a bimolecular binding with multiple binding sites, please use a smaller time step." << std::endl;
            }
            // std::cout <<" IN DETERMIN 2D BINDING, RHO: "<<rho<<" N free lipids: "<<membraneObject.No_free_lipids<<" Binding PROB: "<<rxnProb<<" ";
            moleculeList[biMolData.pro1Index].probvec.back() = rxnProb * currnorm;
            //std::cout <<" SIZE OF PROBVEC FOR MOL: "<<biMolData.pro1Index<<" size:" <<moleculeList[biMolData.pro1Index].probvec.size()<<std::endl;
            moleculeList[proA].crossbase.push_back(proB);
            moleculeList[proA].mycrossint.push_back(ifaceA);
            moleculeList[proA].crossrxn.push_back(std::array<int, 3> { rxnIndex, rateIndex, isStateChangeBackRxn });
            ++complexList[moleculeList[proA].myComIndex].ncross;
            //moleculeList[biMolData.pro2Index].probvec.back() = rxnProb * currnorm;

            moleculeList[proA].currprevsep.push_back(R1);
            moleculeList[proA].currlist.push_back(proB);
            moleculeList[proA].currmyface.push_back(ifaceA);
            moleculeList[proA].currpface.push_back(ifaceB);
            moleculeList[proA].currprevnorm.push_back(currnorm);
            moleculeList[proA].currps_prev.push_back(1.0 - rxnProb * currnorm);
        } // Within reaction zone
    }
}
