#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"
#include "tracing.hpp"

using namespace std;

// write the header for file tracking all species
int init_speciesFile(ofstream& speciesFile, copyCounters& counterArrays,
                     std::vector<MolTemplate>& molTemplateList,
                     std::vector<ForwardRxn>& forwardRxns, Parameters& params) {
  // TRACE();
  int nSpecies = 0;
  speciesFile << "Time (s)";
  for (const auto& oneTemp : molTemplateList) {
    for (auto& iface : oneTemp.interfaceList) {
      if (iface.stateList.size() == 1) {
        speciesFile << ',' << oneTemp.molName << '(' << iface.name << ')';
        counterArrays.singleDouble.push_back(1);
        counterArrays.implicitDouble.push_back(false);
        counterArrays.canDissociate.push_back(false);
        nSpecies++;
      } else {
        for (auto& state : iface.stateList) {
          speciesFile << ',' << oneTemp.molName << '(' << iface.name << '~'
                      << state.iden << ')';
          counterArrays.singleDouble.push_back(1);
          counterArrays.implicitDouble.push_back(false);
          counterArrays.canDissociate.push_back(false);
          nSpecies++;
        }
      }
    }
  }
  for (const auto& oneRxn : forwardRxns) {
    if (oneRxn.rxnType == ReactionType::bimolecular) {
      speciesFile << ',' << oneRxn.productName;
      counterArrays.singleDouble.push_back(2);  // contains two species.
      if (molTemplateList[oneRxn.reactantListNew[0].molTypeIndex]
              .isImplicitLipid ||
          molTemplateList[oneRxn.reactantListNew[1].molTypeIndex]
              .isImplicitLipid) {
        counterArrays.implicitDouble.push_back(true);
      } else {
        counterArrays.implicitDouble.push_back(false);
      }
      if (oneRxn.isReversible == true) {
        counterArrays.canDissociate.push_back(true);
      } else {
        counterArrays.canDissociate.push_back(false);
      }
      nSpecies++;
    }
  }
  speciesFile << std::endl;

  cout << " Total species calculated from init_speciesFile.cpp: " << nSpecies
       << endl;

  return nSpecies;
}
