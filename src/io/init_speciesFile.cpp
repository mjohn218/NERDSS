#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_copyCounters.hpp"
#include "classes/class_Rxns.hpp"

using namespace std;


// write the header for file tracking all species
void init_speciesFile(ofstream &speciesFile, copyCounters& counterArrays, std::vector<MolTemplate>& molTemplateList, std::vector<ForwardRxn>&forwardRxns)
{
  int nSpecies=0;
  speciesFile << "Itr";
  for (const auto& oneTemp : molTemplateList) {
    for (auto& iface : oneTemp.interfaceList) {
                if (iface.stateList.size() == 1) {
                    speciesFile << ',' << oneTemp.molName << '(' << iface.name << ')';
		    counterArrays.singleDouble.push_back(1);
		    nSpecies++;
                } else {
                    for (auto& state : iface.stateList) {
                        speciesFile << ',' << oneTemp.molName << '(' << iface.name << '~' << state.iden << ')';
			counterArrays.singleDouble.push_back(1);
			nSpecies++;
                    }
                }
            }
        }
        for (const auto& oneRxn : forwardRxns) {
	  if (oneRxn.rxnType == ReactionType::bimolecular){
                speciesFile << ',' << oneRxn.productName;
		counterArrays.singleDouble.push_back(2);//contains two species.
		nSpecies++;
	  }
        }
        speciesFile <<std::endl;
	
	cout <<" Total species calculated from init_speciesFile.cpp: "<<nSpecies<<endl;

}

