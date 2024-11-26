
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"
#include <iostream>
#include <sstream>
#include <vector>

void correct_structure(const std::vector<Molecule> &moleculeList, Complex &complex,
                       const std::vector<ForwardRxn> &forwardRxns) {
  int numPromoter{0};
  int numProtein{0};
  int promoterIndex;
  int proteinIndex;
  for (int molIndex : complex.memberList) {
    if (moleculeList[molIndex].isPromoter) {
      numPromoter += 1;
      promoterIndex = molIndex;
    } else if (moleculeList[molIndex].isLipid == false && moleculeList[molIndex].isImplicitLipid == false) {
      // a simple protein
      numProtein += 1;
      proteinIndex = molIndex;
    }
  }
  if (numPromoter == 1 && numProtein == 1) {
    Molecule::Iface ifaceProtein{moleculeList[proteinIndex].interfaceList[0]};
    if (ifaceProtein.interaction.partnerIndex == promoterIndex) {
      // std::cout << "Try to correct a sturcture!" << std::endl;
      int ifacePromoterIndex{ifaceProtein.interaction.partnerIfaceIndex};
      // std::cout << "Found interface " << ifacePromoterIndex << std::endl;
      Coord ifacePromoterCoord{moleculeList[promoterIndex].interfaceList[ifacePromoterIndex].coord};
      Vector iface2iface{};
      iface2iface.x = ifacePromoterCoord.x - ifaceProtein.coord.x;
      iface2iface.y = ifacePromoterCoord.y - ifaceProtein.coord.y;
      iface2iface.z = ifacePromoterCoord.z - ifaceProtein.coord.z;
      // std::cout << "Got the displacement between two interfaces." << std::endl;
      double iface2ifaceMag{sqrt((iface2iface.x * iface2iface.x) +
                                 (iface2iface.y * iface2iface.y) +
                                 (iface2iface.z * iface2iface.z))};
      Vector dispVec{};
      // std::cout << "Try to find the binding radius by bndRxnList" << std::endl;
      double bindRadius{forwardRxns[moleculeList[proteinIndex].bndRxnList[0]].bindRadius};
      // std::cout << "Got the binding radius of this bond." << std::endl;
      dispVec.x = iface2iface.x * (1 - bindRadius / iface2ifaceMag);
      dispVec.y = iface2iface.y * (1 - bindRadius / iface2ifaceMag);
      dispVec.z = iface2iface.z * (1 - bindRadius / iface2ifaceMag);
      // std::cout << "Got the displacement vector." << std::endl;
      Molecule movingMol{moleculeList[proteinIndex]};
      // std::cout << "Moving the protein." << std::endl;
      movingMol.comCoord = dispVec + movingMol.comCoord;
      for (unsigned int j{0}; j < movingMol.interfaceList.size(); ++j)
        movingMol.interfaceList[j].coord = dispVec + movingMol.interfaceList[j].coord;
      // std::cout << "Structure correction finished +1!" << std::endl;
    }
  }
}
