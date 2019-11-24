#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_Membrane.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
void initialize_paramters_for_implicitlipid_model(int& implicitlipidIndex, const Parameters& params, std::vector<ForwardRxn>& forwardRxns, std::vector<BackRxn>& backRxns,
                                                  std::vector<Molecule>& moleculeList, std::vector<MolTemplate>& molTemplateList, std::vector<Complex>& complexList, 
                                                  Membrane& membraneObject)
{
// check if it is a 2D simulation
   bool is2D = true;
   for (int mol {0}; mol < moleculeList.size(); ++mol){
    	  int molTypeIndex = moleculeList[mol].molTypeIndex;
    	  if (molTemplateList[molTypeIndex].isImplicitLipid == true)
             continue;
        if (molTemplateList[molTypeIndex].D.z > 0)
            is2D = false;         
    }
    if (is2D || membraneObject.waterBox.z == 0)
        membraneObject.TwoD = true;
  
// find the mol.index of the implicit-lipid, which should be implicitlipidIndex=0;
    bool systemIL = false;//Are there any IL's in this system?
   for (int mol {0}; mol < moleculeList.size(); ++mol){
    	  int molTypeIndex = moleculeList[mol].molTypeIndex;
    	  if (molTemplateList[molTypeIndex].isImplicitLipid == true){
            moleculeList[mol].isImplicitLipid = true;
            implicitlipidIndex = mol;
	    std::cout<<" INDEX OF IMPLICIT LIPID MOLECULE: "<<implicitlipidIndex<<std::endl;
	    membraneObject.implicitlipidIndex=mol;
            moleculeList[mol].mass = 1;  // such a large value that implicit-lipid won't move when proteins bind to it.
	    systemIL = true;//yes, there are IL's in the system.
        }else{
            moleculeList[mol].isImplicitLipid = false;
        }
    }
   //    membraneObject.implicitlipidIndex = implicitlipidIndex;
/////////////////////////////////////////////////////////////////////////////////////////////
    // set up 'bindToSurface' in molTemplateList. 
    for (const auto& oneRxn : forwardRxns) {
       if (oneRxn.rxnType == ReactionType::bimolecular) {
           int molTypeIndex = oneRxn.reactantListNew[0].molTypeIndex;
           if (molTemplateList[molTypeIndex].isImplicitLipid == true) {
                int newProteinIndex = oneRxn.reactantListNew[1].molTypeIndex;
                molTemplateList[newProteinIndex].bindToSurface = true;
           }
           molTypeIndex = oneRxn.reactantListNew[1].molTypeIndex;
          if (molTemplateList[molTypeIndex].isImplicitLipid == true) {
               int newProteinIndex = oneRxn.reactantListNew[0].molTypeIndex;
               molTemplateList[newProteinIndex].bindToSurface = true;
          } 
       }
    }
/////////////////////////////////////////////////////////////////////////////////////////////
    membraneObject.totalSA = membraneObject.waterBox.x * membraneObject.waterBox.y; 
    // here, we assume that the membrane surface is on the Z-axis bottom of the cubic box.
      
// initialize several parameters in 'membraneObject' for implicit-lipid model  
    if(systemIL==true){
      // initial number of free lipids' interface
      membraneObject.nSites = 0;
      for (auto onetem:molTemplateList){
        if(onetem.isImplicitLipid)
	  membraneObject.nSites = onetem.copies;
	// here, we assume each implicit-lipid has only one interface. 
      }
      //   membraneObject.No_free_lipids = membraneObject.nSites;
      ///////////////////////////////////////////////////////////   
      // initial number of proteins' interface that can bind to implicit-lipids
      membraneObject.No_protein = 0;
      for(int mol {0}; mol < moleculeList.size(); ++mol){
	if(moleculeList[mol].isImplicitLipid)
	  continue;   	     
        int molType = moleculeList[mol].molTypeIndex;
        for(int relfaceItr {0}; relfaceItr<moleculeList[mol].interfaceList.size(); ++relfaceItr){
	  int stateIndex = moleculeList[mol].interfaceList[relfaceItr].stateIndex;
	  const Interface::State& state = molTemplateList[molType].interfaceList[relfaceItr].stateList[stateIndex];
	  for(auto rxnItr : state.myForwardRxns) {
	    const ForwardRxn& oneRxn = forwardRxns[rxnItr];
	    for(int reactItr { 0 }; reactItr < oneRxn.reactantListNew.size(); ++reactItr) {
	      if (moleculeList[implicitlipidIndex].interfaceList[0].index 
		  == oneRxn.reactantListNew[reactItr].absIfaceIndex)
		membraneObject.No_protein ++ ;
	    }
	  }
        }      
      }
      ///////////////////////////////////////////////////////////
      // area of the membrane surface where the implicit-lipids are on.
      ///////////////////////////////////////////////////////////////////////////////////////////////////    
      // calculate the reflecting-surface RS3D
      //////////////////// to calculate Dtot;   
      int com1Index = moleculeList[implicitlipidIndex].myComIndex;
      int com2Index;
      for(auto& mol : moleculeList){
	bool thisone = false;
	if(mol.isImplicitLipid)
	  continue; 	 
	int moltype = mol.molTypeIndex;
	for (auto partner : molTemplateList[moltype].rxnPartners) {
          if (partner == moleculeList[implicitlipidIndex].molTypeIndex) {
	    thisone = true;
	    break;
          }
	}
	if(thisone){
	  com2Index = mol.myComIndex;
	  break;
	}   	     
      }
      double Dtot = 1.0/3.0 * (complexList[com1Index].D.x + complexList[com2Index].D.x)
	+ 1.0/3.0 * (complexList[com1Index].D.y + complexList[com2Index].D.y)
	+ 1.0/3.0 * (complexList[com1Index].D.z + complexList[com2Index].D.z); 
      ///////////////// to find sigma and ka;
      double RS3D = 0; 
      int rxnIndex { -1 };
      int rateIndex { -1 };
      double ka;
      double sigma;
      if (membraneObject.implicitLipid && !membraneObject.TwoD){
	for(auto& mol : moleculeList){
	  if (mol.isImplicitLipid)
      	    continue;    
    	  int proMolType = mol.molTypeIndex;
	  for (int relIface1 = 0; relIface1 < mol.interfaceList.size(); ++relIface1) {
            int absIface1 { mol.interfaceList[relIface1].index };
            int stateIndex1 { mol.interfaceList[relIface1].stateIndex };
            const Interface::State& state {molTemplateList[proMolType].interfaceList[relIface1].stateList[stateIndex1]};
            for (auto& statePartner : state.rxnPartners) {
	      int relIface2 { 0 };
	      int absIface2 { moleculeList[implicitlipidIndex].interfaceList[relIface2].index };
	      if (absIface2 == statePartner) {
		bool isStateChangeBackRxn { false };
		find_which_reaction(relIface1, relIface2, rxnIndex, rateIndex, isStateChangeBackRxn, state,
				    mol, moleculeList[implicitlipidIndex], forwardRxns, backRxns, molTemplateList); 
		if (rxnIndex != -1 && rateIndex != -1) {
		  ka = 2.0 * forwardRxns[rxnIndex].rateList[rateIndex].rate;
		  sigma = forwardRxns[rxnIndex].bindRadius;
		  RS3D = sigma*ka/(ka + 4*M_PI*sigma*Dtot);  
		  break;
		}
	      }
            }
            if (RS3D>0){break;}
	  }
	  if (RS3D>0){break;} 
	}
	if (RS3D <= 0){
	  std::cout<<"Error: reflecting surface RS3D is not successfully calculated. Exiting..."<<std::endl;
	  exit(1);
	}
      }
      membraneObject.RS3D = RS3D;
      std::cout <<" RS3D value: "<<RS3D<<std::endl;
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //calculate unbinding probabilities of 2D->3D, and 2D->2D
      double kb = 0;
      int rxnIndexback, rateIndexback;
      if (membraneObject.implicitLipid){   	
	for (auto& mol : moleculeList){
	  if (mol.isImplicitLipid)
    	       continue;
	  for (int relIface1 = 0; relIface1 < mol.interfaceList.size(); ++relIface1){
	    rxnIndexback  = rxnIndex;
	    rateIndexback = rateIndex;
	    kb = backRxns[rxnIndexback].rateList[rateIndexback].rate; // 
	    //std::cout<<"kb ="<<kb<<std::endl; 
	    if (kb > 0)
	      break;
          }
          if (kb > 0) break;
	}
	double h = params.timeStep;  
	membraneObject.Punbinding3D = dissociate3D(h,Dtot,sigma,ka,kb); 
	if (membraneObject.Punbinding3D <= 0) {
	  std::cout<<"Error: unbinding probability (2D->3D) is not successfully calculated !"<<std::endl;
	  exit(1);
        }
        std::cout<<"membraneObject.Punbinding3D"<<std::endl;
	///////////////////////////////////////////////////////////////////////////
	double Dtot2D = 1;//1.0/3.0 * (complexList[com1Index].D.x + complexList[com1Index].D.y + complexList[com1Index].D.z);
	paramsIL params2D {};
	params2D.R2D = 0.0;
	params2D.sigma = forwardRxns[rxnIndex].bindRadius;
	params2D.Dtot = Dtot2D;
	params2D.ka = forwardRxns[rxnIndex].rateList[rateIndex].rate / 2.0 / forwardRxns[rxnIndex].bindRadius;
	params2D.kb = backRxns[rxnIndexback].rateList[rateIndexback].rate;       
	params2D.area = membraneObject.totalSA;
	params2D.dt = params.timeStep;
	params2D.Nlipid = membraneObject.nSites; // the initial number of lipids' interfaces on the surface
	params2D.Na = membraneObject.No_protein; // the initial number of protein's interfaces that can bind to surface
	membraneObject.Pbinding2D = pimplicitlipid_2D(params2D);// binding probability, must multiply the lipid density
	membraneObject.Punbinding2D = dissociate2D(params2D); // unbinding probability.
	if (membraneObject.Punbinding2D <= 0) {
	  std::cout<<"Error: unbinding probability (2D->2D) is not successfully calculated !"<<std::endl;
	  exit(1);
	}
      }
      std::cout<<"unbinding3D="<<membraneObject.Punbinding3D
	       <<", unbinding2D="<<membraneObject.Punbinding2D<<", binding2D="<<membraneObject.Pbinding2D<<std::endl;
    }//only set up IL model if they exist. 

}
 
