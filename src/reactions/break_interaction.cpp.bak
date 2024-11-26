#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "math/rand_gsl.hpp"
#include "tracing.hpp"

bool break_interaction(long long int iter, size_t relIface1, size_t relIface2, Molecule& reactMol1, Molecule& reactMol2,
    const BackRxn& currRxn, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, std::vector<MolTemplate>& molTemplateList, int ILindexMol, 
    const ForwardRxn& conjForwardRxn, bool& breakLinkComplex, double timeStep) {
  bool cancelDissociation = false; // this can be true during loop breaking due to correction term.
  breakLinkComplex = false; // default is full dissociation

  reactMol1.bndpartner.erase(remove(reactMol1.bndpartner.begin(), reactMol1.bndpartner.end(), reactMol2.index), reactMol1.bndpartner.end());
  if (reactMol2.isImplicitLipid == false)
      reactMol2.bndpartner.erase(remove(reactMol2.bndpartner.begin(), reactMol2.bndpartner.end(), reactMol1.index), reactMol2.bndpartner.end());
  
  bool keepSameComplex;

  unsigned newComIndex = complexList.size();

  if (Complex::emptyComList.size() != 0 && complexList[Complex::emptyComList.back()].isEmpty) {
      // if there is an empty complex slot, make the new Complex in it
      newComIndex = Complex::emptyComList.back();
      Complex::emptyComList.pop_back(); // remove the index from the list
  } else {
      // if we're making a new Complex, create an empty one at the end (index complexList.size())
      complexList.emplace_back();
  }

  if (reactMol2.isImplicitLipid == false)
      keepSameComplex = determine_parent_complex_IL(reactMol1.index, reactMol2.index, newComIndex, moleculeList, complexList, ILindexMol);
  else
      keepSameComplex = false;

  /*if they are within the same complex, use the loop correction rate to decide whether to go forward with the dissociation*/
  if (keepSameComplex) {
      double rateForward = conjForwardRxn.rateList[0].rate;
      double largestRate = -1;
      for (auto oneForwardRate : conjForwardRxn.rateList) {
          if (oneForwardRate.rate > largestRate)
              largestRate = oneForwardRate.rate;
      }
      double c0_nm3 = 0.602; // standard state 1M in units of /nm^3
      double coop = conjForwardRxn.loopCoopFactor;
      double rateClose = largestRate * c0_nm3 * coop;
      // rateClose will be in units of /us (due to ka units), so no need for 1E-6 factor, since timeStep has
      // units of us

      double poisson = timeStep * rateClose;
      double correctionRatio{(1 - exp(-poisson)) / poisson};
  
    //std::cout <<"Correction Ratio: "<<correctionRatio<<std::endl;

      if (1.0 * rand_gsl() > correctionRatio) {
          /*Cancel the dissociation!!*/
          cancelDissociation = true;
      }
  }

  if (cancelDissociation) {
        /*cancel the dissociation, move both interfaces back into the bndpartner list, which is the only change we've made so far.*/
        reactMol1.bndpartner.push_back(reactMol2.index);
        reactMol2.bndpartner.push_back(reactMol1.index);

        // reset empty complexList
        if (newComIndex + 1 == complexList.size())
            complexList.pop_back();
        else
            Complex::emptyComList.push_back(newComIndex);

        return true;
    } else {
      std::vector<int> numEachMolPrevious {};
      std::vector<long long int> lastNumberUpdateItrEachMolPrevious {};
      for (unsigned index = 0; index < molTemplateList.size(); index++) {
          numEachMolPrevious.emplace_back(complexList[reactMol1.myComIndex].numEachMol[index]);
          lastNumberUpdateItrEachMolPrevious.emplace_back(complexList[reactMol1.myComIndex].lastNumberUpdateItrEachMol[index]);
      }

      int absIface1 { -1 };
      int absIface2 { -1 };
      if (currRxn.isSymmetric) {
          absIface1 = currRxn.productListNew[0].absIfaceIndex;
          absIface2 = currRxn.productListNew[1].absIfaceIndex;
      } else {
          /*Here, if both proteins are the same protein, same interface, but distinct states, need to correct for that.*/
          // std::cout << "State of mol1: " << reactMol1.interfaceList[relIface1].stateIden << " of mol2: " << reactMol2.interfaceList[relIface2].stateIden << "\n";
          // std::cout << " Product requires state: " << currRxn.productListNew[0].requiresState << " " << currRxn.productListNew[1].requiresState << std::endl;
          if (reactMol1.molTypeIndex == currRxn.productListNew[0].molTypeIndex
              && relIface1 == currRxn.productListNew[0].relIfaceIndex && currRxn.productListNew[0].requiresState == reactMol1.interfaceList[relIface1].stateIden) {
              // matched protein, interface, and state of product[0] to mol1
              absIface1 = currRxn.productListNew[0].absIfaceIndex;
              absIface2 = currRxn.productListNew[1].absIfaceIndex;

          } else if (reactMol1.molTypeIndex == currRxn.productListNew[1].molTypeIndex
              && relIface1 == currRxn.productListNew[1].relIfaceIndex && currRxn.productListNew[1].requiresState == reactMol1.interfaceList[relIface1].stateIden) {
              // matched protein, interface and state of product[1] to mol1
              absIface2 = currRxn.productListNew[0].absIfaceIndex;
              absIface1 = currRxn.productListNew[1].absIfaceIndex;
          } else {
              std::cout << " IN BREAK INTERACTION, DID NOT MATCH protein, interface, and state of a product to Mol1 " << reactMol1.index << std::endl;
              std::cout << "reactMol1.molTypeIndex: " << reactMol1.molTypeIndex << std::endl;
              std::cout << "currRxn.productListNew[0].molTypeIndex: " << currRxn.productListNew[0].molTypeIndex << std::endl;
              std::cout << "currRxn.productListNew[1].molTypeIndex: " << currRxn.productListNew[1].molTypeIndex << std::endl;
              std::cout << "reactMol1.molTypeIndex == currRxn.productListNew[0].molTypeIndex: " << (reactMol1.molTypeIndex == currRxn.productListNew[0].molTypeIndex) << std::endl;
              std::cout << "reactMol1.molTypeIndex == currRxn.productListNew[1].molTypeIndex: " << (reactMol1.molTypeIndex == currRxn.productListNew[1].molTypeIndex) << std::endl;
              std::cout << "relIface1: " << relIface1 << std::endl;
              std::cout << "currRxn.productListNew[0].relIfaceIndex: " << currRxn.productListNew[0].relIfaceIndex << std::endl;
              std::cout << "currRxn.productListNew[1].relIfaceIndex: " << currRxn.productListNew[1].relIfaceIndex << std::endl;
              std::cout << "relIface1 == currRxn.productListNew[0].relIfaceIndex: " << (relIface1 == currRxn.productListNew[0].relIfaceIndex) << std::endl;
              std::cout << "relIface1 == currRxn.productListNew[1].relIfaceIndex: " << (relIface1 == currRxn.productListNew[1].relIfaceIndex) << std::endl;
              std::cout << "currRxn.productListNew[0].requiresState: " << currRxn.productListNew[0].requiresState << std::endl;
              std::cout << "currRxn.productListNew[1].requiresState: " << currRxn.productListNew[1].requiresState << std::endl;
              std::cout << "reactMol1.interfaceList[relIface1].stateIden: " << reactMol1.interfaceList[relIface1].stateIden << std::endl;
              std::cout << "currRxn.productListNew[0].requiresState == reactMol1.interfaceList[relIface1].stateIden: " << (currRxn.productListNew[0].requiresState == reactMol1.interfaceList[relIface1].stateIden) << std::endl;
              std::cout << "currRxn.productListNew[1].requiresState == reactMol1.interfaceList[relIface1].stateIden: " << (currRxn.productListNew[1].requiresState == reactMol1.interfaceList[relIface1].stateIden) << std::endl;
              exit(1);
          }
      }

      reactMol1.interfaceList[relIface1].index = absIface1;
      reactMol1.interfaceList[relIface1].interaction.clear();
      reactMol1.interfaceList[relIface1].isBound = false;
      if (reactMol2.isImplicitLipid == false) {
          reactMol2.interfaceList[relIface2].index = absIface2;
          reactMol2.interfaceList[relIface2].interaction.clear();
          reactMol2.interfaceList[relIface2].isBound = false;
      }
      // Add these protein into the bimolecular association list
      reactMol1.freelist.push_back(relIface1);
      reactMol1.bndlist.erase(remove(reactMol1.bndlist.begin(), reactMol1.bndlist.end(), relIface1), reactMol1.bndlist.end());
      // reactMol1.bndlist.erase(std::find_if(reactMol1.bndlist.begin(), reactMol1.bndlist.end(), [&](const size_t& iface) { return iface == relIface1; }));
      // reactMol1.bndpartner.erase(std::find_if(reactMol1.bndpartner.begin(), reactMol1.bndpartner.end(), [&](const size_t& mol) { return mol == reactMol2.index; }));
      if (reactMol2.isImplicitLipid == false) {
          reactMol2.freelist.push_back(relIface2);
          reactMol2.bndlist.erase(remove(reactMol2.bndlist.begin(), reactMol2.bndlist.end(), relIface2), reactMol2.bndlist.end());
          // reactMol2.bndlist.erase(std::find_if(reactMol2.bndlist.begin(), reactMol2.bndlist.end(), [&](const size_t& iface) { return iface == relIface2; }));
          //	  reactMol2.bndpartner.erase(std::find_if(reactMol2.bndpartner.begin(), reactMol2.bndpartner.end(), [&](const size_t& mol) { return mol == reactMol1.index; }));
      }

      if (!keepSameComplex) {
          /*continue on with the dissociation that creates two complexes*/
          complexList[reactMol1.myComIndex].update_properties(moleculeList, molTemplateList);

          // add a lifetime to the previous larger cluster
          // compare current cluster size with the previous ones
          for (unsigned index = 0; index < molTemplateList.size(); index++) {
              if (molTemplateList[index].countTransition == true && complexList[reactMol1.myComIndex].numEachMol[index] < numEachMolPrevious[index]) {
                  // shrink
                  molTemplateList[index].lifeTime[numEachMolPrevious[index] - 1].emplace_back((iter - lastNumberUpdateItrEachMolPrevious[index]) * Parameters::dt / 1E6);
                  complexList[reactMol1.myComIndex].lastNumberUpdateItrEachMol[index] = iter;
              }
          }

          if (reactMol2.isImplicitLipid == false) {
              complexList[reactMol2.myComIndex].update_properties(moleculeList, molTemplateList);

              // update lastNumberUpdateItrEachMol for new complex
              complexList[reactMol2.myComIndex].lastNumberUpdateItrEachMol.resize(molTemplateList.size());
              for (unsigned index = 0; index < molTemplateList.size(); index++) {
                  complexList[reactMol2.myComIndex].lastNumberUpdateItrEachMol[index] = iter;
              }

              complexList[reactMol2.myComIndex].isEmpty = false;
          }

          // update transition matrix
          // compare current cluster size with the previous ones
          for (unsigned index = 0; index < molTemplateList.size(); index++) {
              if (molTemplateList[index].countTransition == true && complexList[reactMol1.myComIndex].numEachMol[index] < numEachMolPrevious[index]) {
                  // shrink
                  if (numEachMolPrevious[index] - 1 >= 0)
                      molTemplateList[index].transitionMatrix[numEachMolPrevious[index] - 1][numEachMolPrevious[index] - 1] += iter - Parameters::lastUpdateTransition[index] - 1;
                  if (complexList[reactMol1.myComIndex].numEachMol[index] - 1 >= 0 && numEachMolPrevious[index] - 1 >= 0)
                      molTemplateList[index].transitionMatrix[complexList[reactMol1.myComIndex].numEachMol[index] - 1][numEachMolPrevious[index] - 1] += 1;
                  if (complexList[reactMol2.myComIndex].numEachMol[index] - 1 >= 0 && numEachMolPrevious[index] - 1 >= 0)
                      molTemplateList[index].transitionMatrix[complexList[reactMol2.myComIndex].numEachMol[index] - 1][numEachMolPrevious[index] - 1] += 1;

                  // update diagonal elements for unchanged complexes
                  for (unsigned indexCom = 0; indexCom < complexList.size(); indexCom++) {
                      if ((indexCom != reactMol1.myComIndex) && (indexCom != reactMol2.myComIndex)) {
                          if (complexList[indexCom].numEachMol[index] - 1 >= 0) {
                              molTemplateList[index].transitionMatrix[complexList[indexCom].numEachMol[index] - 1][complexList[indexCom].numEachMol[index] - 1] += iter - Parameters::lastUpdateTransition[index];
                          }
                      }
                  }

                  // update lastUpdateTransition
                  Parameters::lastUpdateTransition[index] = iter;
              }
          }

          ++Complex::numberOfComplexes;

          complexList[newComIndex].id = Complex::maxID++;

          for (auto i : complexList[newComIndex].memberList) {
            moleculeList[i].complexId = complexList[newComIndex].id;
          }
      } else {
          /*All proteins remain in complex c1, dissociation
                will break the product state of the two proteins that dissociated but here they
                are linked in a closed loop so it will not create a new complex.
                positions don't change
          */
          // reset empty complexList
          breakLinkComplex = true;
          if (newComIndex + 1 == complexList.size())
              complexList.pop_back();
          else
              Complex::emptyComList.push_back(newComIndex);
      }
      //------------------------START UPDATE MONOMERLIST-------------------------
      // update oneTemp.monomerList when oneTemp.canDestroy is true and mol is monomer, add to monomerList if new monomer produced
      // reactMol1
      // {
      //     Molecule& oneMol { reactMol1 };
      //     MolTemplate& oneTemp { molTemplateList[oneMol.molTypeIndex] };
      //     bool isMonomer { oneMol.bndpartner.empty() };
      //     bool canDestroy { oneTemp.canDestroy };
      //     if (isMonomer && canDestroy) {
      //         // add to monomerList
      //         oneTemp.monomerList.emplace_back(oneMol.index);
      //     }

      //     // std::cout << "For mol " << oneMol.index << ": "
      //     //           << "canDestory is " << oneTemp.canDestroy << "\t"
      //     //           << "isMonomer is " << isMonomer << std::endl;
      //     // std::cout << "Now the monomerList is: ";
      //     // for (auto one : oneTemp.monomerList) {
      //     //     std::cout << one << "\t";
      //     // }
      //     // std::cout << std::endl;
      // }
      // reactMol2
      // {
      //     Molecule& oneMol { reactMol2 };
      //     MolTemplate& oneTemp { molTemplateList[oneMol.molTypeIndex] };
      //     bool isMonomer { oneMol.bndpartner.empty() };
      //     bool canDestroy { oneTemp.canDestroy };
      //     if (isMonomer && canDestroy) {
      //         // add to monomerList
      //         oneTemp.monomerList.emplace_back(oneMol.index);
      //     }
      //     // std::cout << "For mol " << oneMol.index << ": "
      //     //           << "canDestory is " << oneTemp.canDestroy << "\t"
      //     //           << "isMonomer is " << isMonomer << std::endl;
      //     // std::cout << "Now the monomerList is: ";
      //     // for (auto one : oneTemp.monomerList) {
      //     //     std::cout << one << "\t";
      //     // }
      //     // std::cout << std::endl;
      // }
      //------------------------END UPDATE MONOMERLIST---------------------------
  } // end if to actually perform dissociation, so not cancelled
  return cancelDissociation;
}
