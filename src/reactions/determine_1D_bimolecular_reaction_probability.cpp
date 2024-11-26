#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "tracing.hpp"
#include <iostream>
#include <sstream>

void determine_1D_bimolecular_reaction_probability(
    int simItr, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    BiMolData &biMolData, const Parameters &params,
    std::vector<Molecule> &moleculeList, std::vector<Complex> &complexList,
    const std::vector<ForwardRxn> &forwardRxns,
    const std::vector<BackRxn> &backRxns) {

  // rotation
  // double Dr1{};
  // double cf{cos(sqrt(2.0 * complexList[biMolData.com1Index].Dr.z * params.timeStep))};
  // Dr1 = 2.0 * biMolData.magMol1 * (1.0 - cf);
  // double Dr2{};
  // cf = cos(sqrt(2.0 * complexList[biMolData.com2Index].Dr.z * params.timeStep));
  // Dr2 = 2.0 * biMolData.magMol2 * (1.0 - cf);
  // biMolData.Dtot += (Dr1 + Dr2) / (4.0 * params.timeStep); // add in contributions from rotation

  double RMax{4.0 * sqrt(2.0 * biMolData.Dtot * params.timeStep) +
              forwardRxns[rxnIndex].bindRadius};
  double sep{};
  double R1{};
  bool withinRmax{get_distance(
      biMolData.pro1Index, biMolData.pro2Index, biMolData.relIface1,
      biMolData.relIface2, rxnIndex, rateIndex, isStateChangeBackRxn, sep, R1,
      RMax, complexList, forwardRxns[rxnIndex], moleculeList, false)};
  if (withinRmax) {
    // in case they dissociated
    // moleculeList[biMolData.pro1Index].display_all();
    // moleculeList[biMolData.pro2Index].display_all();
    moleculeList[biMolData.pro1Index].probvec.push_back(0);
    moleculeList[biMolData.pro2Index].probvec.push_back(0);
  }

  if (moleculeList[biMolData.pro1Index].isDissociated != true &&
      moleculeList[biMolData.pro2Index].isDissociated != true) {
    /*This movestat check is if you allow just dissociated proteins to avoid
     * overlap*/
    if (withinRmax && forwardRxns[rxnIndex].rateList[rateIndex].rate > 0) {
      // int probMatrixIndex{0};
      /*Evaluate probability of reaction, with reweighting*/
      // double ratio = forwardRxns[rxnIndex].bindRadius / R1;
      if (sep < 0) {
        // (biMolData.com1Index != biMolData.com2Index) 
        sep = 0;
        // ratio = 1;
        R1 = forwardRxns[rxnIndex].bindRadius;
      }

      // declare intrinsic binding rate of 1D->1D case.
      double kact{forwardRxns[rxnIndex].rateList[rateIndex].rate / forwardRxns[rxnIndex].area3Dto1D};
      if (forwardRxns[rxnIndex].isSymmetric == false)
        kact /= 2.0; // for A(a)+B(b)->A(a!).B(b!) case

      double currnorm{1.0};
      double p0_ratio{1.0};

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

      double rxnProb{};
      for (int s{0}; s < moleculeList[proA].prevlist.size(); ++s) {
        if (moleculeList[proA].prevlist[s] == proB &&
            moleculeList[proA].prevmyface[s] == ifaceA &&
            moleculeList[proA].prevpface[s] == ifaceB) {
          if (moleculeList[proA].prevsep[s] >= RMax) {
            // BEcause previous reweighting was for 3D, now
            p0_ratio = 1.0;
            // restart reweighting in 2D.
            currnorm = 1.0;
          } else {
            p0_ratio = pirr_pfree_ratio_psF_1D(
                R1, moleculeList[proA].prevsep[s], params.timeStep,
                biMolData.Dtot, forwardRxns[rxnIndex].bindRadius, kact,
                moleculeList[proA].ps_prev[s]);
            currnorm = moleculeList[proA].prevnorm[s] * p0_ratio;
          }
          break;
        }
      }
      rxnProb = passocF_1D(R1, params.timeStep, biMolData.Dtot,
                           forwardRxns[rxnIndex].bindRadius, kact);
      // std::cout << "In determine_1D_prob proB " << proB << " proA " << proA
      //           << " R1 " << R1 << " kact " << kact << " rxnProb " << rxnProb
      //           << " p0_ratio " << p0_ratio << " currnorm " << currnorm << std::endl;
      moleculeList[biMolData.pro1Index].probvec.back() = rxnProb * currnorm;
      moleculeList[biMolData.pro2Index].probvec.back() = rxnProb * currnorm;
      if (rxnProb > 1.000001) {
        std::cerr << "Error: prob of reaction is: " << rxnProb
                  << " > 1. Avoid this using a smaller time step." << std::endl;
        // exit(1);
      }
      if (rxnProb > 0.5) {
        std::cout << "WARNING: prob of reaction > 0.5. If this is a reaction "
                     "for a bimolecular binding with multiple binding sites, "
                     "please use a smaller time step."
                  << std::endl;
      }
      moleculeList[proA].currprevsep.push_back(R1);
      moleculeList[proA].currlist.push_back(proB);
      moleculeList[proA].currmyface.push_back(ifaceA);
      moleculeList[proA].currpface.push_back(ifaceB);
      moleculeList[proA].currprevnorm.push_back(currnorm);
      moleculeList[proA].currps_prev.push_back(1.0 - rxnProb * currnorm);
    } // Within reaction zone
  }
}