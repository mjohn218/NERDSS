#include "classes/class_Coord.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Vector.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include <iostream>

double calc_distance_1d(std::vector<Molecule> &moleculeList, int pro, int iface, const ForwardRxn& currRxn) {
  double coordx{moleculeList[pro].interfaceList[iface].coord.x};
  // Only consider proteins that only binds to a promoter
  if (moleculeList[pro].bndpartner.size() == 1) {
    int idx_bndpartner{moleculeList[pro].bndpartner[0]};
    if (moleculeList[idx_bndpartner].isPromoter) {
      // the relative coordinate of this interface to its COM
      double rel_self_x{moleculeList[pro].interfaceList[iface].coord.x - moleculeList[pro].comCoord.x};
      // find the correct position of "pro" COM
      ///// find the interface bound with promoter
      Molecule::Iface iface_self{moleculeList[pro].interfaceList[0]};
      if (!(iface_self.interaction.partnerIndex == idx_bndpartner)) {
        std::cout << "WRONG INTERFACE!" << std::endl;
      }
      int idx_iface_partner{iface_self.interaction.partnerIfaceIndex};
      Coord promoterIfaceCoord{moleculeList[idx_bndpartner].interfaceList[idx_iface_partner].coord};
      ///// find the correct position of self interface
      Vector iface2iface;
      iface2iface.x = iface_self.coord.x - promoterIfaceCoord.x;
      iface2iface.y = iface_self.coord.y - promoterIfaceCoord.y;
      iface2iface.z = iface_self.coord.z - promoterIfaceCoord.z;
      double iface2ifaceMag{sqrt((iface2iface.x * iface2iface.x) +
                                 (iface2iface.y * iface2iface.y) +
                                 (iface2iface.z * iface2iface.z))};
      double iface_x{currRxn.bindRadius / iface2ifaceMag * iface2iface.x + promoterIfaceCoord.x};
      double com_x{iface_x + moleculeList[pro].comCoord.x - iface_self.coord.x};
      coordx = com_x + rel_self_x;
    }
  }
  return coordx;
}

bool get_distance(int pro1, int pro2, int iface1, int iface2, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    double& sep, double& R1, double Rmax, std::vector<Complex>& complexList, const ForwardRxn& currRxn,
    std::vector<Molecule>& moleculeList, bool isSphere)
{
    bool is2D = false;
    bool is1D = false;
    if (complexList[moleculeList[pro1].myComIndex].onFiber && complexList[moleculeList[pro2].myComIndex].onFiber) {
        is1D = true;
    // } else if (complexList[moleculeList[pro1].myComIndex].D.z < 1E-8 && complexList[moleculeList[pro2].myComIndex].D.z < 1E-8) {
    } else if (complexList[moleculeList[pro1].myComIndex].OnSurface &&
               complexList[moleculeList[pro2].myComIndex].OnSurface) {
      is2D = true;
    }

    if (isSphere == true && is2D == true) {
        Coord iface11 = moleculeList[pro1].interfaceList[iface1].coord;
        Coord iface22 = moleculeList[pro2].interfaceList[iface2].coord;
        double r1 = iface11.get_magnitude();
        double r2 = iface22.get_magnitude();
        double r = (r1 + r2) / 2.0; //membraneObject.sphereR; //
        double theta = acos((iface11.x * iface22.x + iface11.y * iface22.y + iface11.z * iface22.z) / r1 / r2);
        R1 = r * theta;
        sep = R1 - currRxn.bindRadius;
    } else if (is1D) {
        double coordx1{moleculeList[pro1].interfaceList[iface1].coord.x};
        double coordx2{moleculeList[pro2].interfaceList[iface2].coord.x};
        R1 = abs(coordx1 - coordx2);
        sep = R1 - currRxn.bindRadius;
        } else {
          double dx = moleculeList[pro1].interfaceList[iface1].coord.x -
                      moleculeList[pro2].interfaceList[iface2].coord.x;
          double dy = moleculeList[pro1].interfaceList[iface1].coord.y -
                      moleculeList[pro2].interfaceList[iface2].coord.y;
          // double dz {
          // (std::abs(complexList[moleculeList[pro1].myComIndex].D.z - 0) <
          // 1E-10
          //                 &&
          //                 std::abs(complexList[moleculeList[pro2].myComIndex].D.z
          //                 - 0) < 1E-10)
          //         ? 0
          //         : moleculeList[pro1].interfaceList[iface1].coord.z -
          //         moleculeList[pro2].interfaceList[iface2].coord.z };
          double dz{};
          if (is2D == true) {
            dz = 0;
        } else {
            dz = moleculeList[pro1].interfaceList[iface1].coord.z - moleculeList[pro2].interfaceList[iface2].coord.z;
        }
        R1 = sqrt((dx * dx) + (dy * dy) + (dz * dz));
        sep = R1 - currRxn.bindRadius;
        }
    /*Rmax should be the binding radius plus ~max diffusion distance, using
   * 3*sqrt(6*Dtot*deltat)*/
    if (R1 < Rmax) {
        /*in this case we evaluate the probability of this reaction*/
        moleculeList[pro1].crossbase.push_back(pro2);
        moleculeList[pro2].crossbase.push_back(pro1);
        moleculeList[pro1].mycrossint.push_back(iface1);
        moleculeList[pro2].mycrossint.push_back(iface2);
        moleculeList[pro1].crossrxn.push_back(
            std::array<int, 3> { rxnIndex, rateIndex, isStateChangeBackRxn });
        moleculeList[pro2].crossrxn.push_back(
            std::array<int, 3> { rxnIndex, rateIndex, isStateChangeBackRxn });
        ++complexList[moleculeList[pro1].myComIndex].ncross;
        ++complexList[moleculeList[pro2].myComIndex].ncross;
        return true;
    }
    return false;
}
