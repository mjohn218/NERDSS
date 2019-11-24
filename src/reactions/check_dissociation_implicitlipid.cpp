#include "boundary_conditions/reflect_functions.hpp"
#include "io/io.hpp"
#include "math/constants.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"

void check_dissociation_implicitlipid(unsigned int simItr, const Parameters& params, SimulVolume& simulVolume,
    const std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, unsigned int molItr,
    std::vector<int>& emptyMolList, std::vector<int>& emptyComList, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, const std::vector<BackRxn>& backRxns, const std::vector<ForwardRxn>& forwardRxns,
				      const std::vector<CreateDestructRxn>& createDestructRxns, copyCounters& counterArrays, const Membrane &membraneObject, std::vector<double> & IL2DbindingVec, std::vector<double> & IL2DUnbindingVec, std::vector<double> &ILTableIDs)
{
    double Dtot;
    for (int relIface1Itr { 0 }; relIface1Itr < moleculeList[molItr].bndlist.size(); ++relIface1Itr) {
        int relIface1 { moleculeList[molItr].bndlist[relIface1Itr] };
        if (moleculeList[molItr].interfaceList[relIface1].isBound) { // make sure it's actually bound
            int pro2Index { moleculeList[molItr].interfaceList[relIface1].interaction.partnerIndex };
            int relIface2 { moleculeList[molItr].interfaceList[relIface1].interaction.partnerIfaceIndex };
            
            int com1Index { moleculeList[molItr].myComIndex };
            int com2Index { moleculeList[pro2Index].myComIndex };
            // only consider when pro2 is implicit-lipid 
            if (complexList[com1Index].OnSurface ==false || moleculeList[pro2Index].isImplicitLipid == false)
                continue;
                
            int mu { moleculeList[molItr].interfaceList[relIface1].interaction.conjBackRxn };    
            // if the reaction is irreversible, the conjBackRxn index will be -1, so continue
            if (mu == -1)
                  continue;

            int rateItr { find_reaction_rate_state(simItr, relIface1, relIface2, moleculeList[molItr],
                    moleculeList[pro2Index], backRxns[mu], molTemplateList) };
            if (rateItr == -1)
                    continue;

            double kb { backRxns[mu].rateList[rateItr].rate }; // <- kr[mu]
            int kfIndex= backRxns[mu].conjForwardRxnIndex ; //!< the index of this reaction's ForwardRxn counterpart
	    double ka=forwardRxns[kfIndex].rateList[rateItr].rate;
	    double prob;         
	    /*Calculate Dtot for the binding reaction that gave rise to this complex. Assume Dcomplex=(1/D1+1/D2)-1 (two components)
	      D1=molTemplate[implicitLipidIndex].D, and Dcomplex is current D, solve for D2, and Dtot=D1+D2.
	     */
	    	    //std::cout <<" In check implicit lipid Unbinding Complex D: "<<Dcomplex_x<<' '<<Dcomplex_y<<' '<<Dcomplex_z<<" And DIL:" <<DIL_x<<' '<<DIL_y<<' '<<DIL_z<<std::endl;
	    //std::cout <<" In check implicit lipid Unbinding, calculated Dtot of protein bound to IL: "<<D2_x<<' '<<D2_y<<' '<<D2_z<<" And Dtot:" <<Dtot<<std::endl;
            // the expression of unbinding probability varies as implicit-lipid model (2D/3D) and explicit-lipid model
            // if the complex has only one protein-implicit.lipid bond, then the unbinding will make the complex into solution, 2D->3D, 
            // then the reflect-surface RS3D needs to be introduced; otherwise, it is just 2D->2D, we donnot need use RS3D.
            if (complexList[com1Index].linksToSurface == 1 && membraneObject.TwoD == false) { 
                // 2D->3D
		double Dcomplex_x=complexList[com1Index].D.x;
		double DIL_x=molTemplateList[moleculeList[pro2Index].molTypeIndex].D.x;
		double D2_x=1.0/(1.0/Dcomplex_x-1.0/DIL_x);
		//If it is a protein bound to a lipid only, then DIL might be== Dcomplex due to forced rounding
		if(Dcomplex_x==DIL_x)D2_x=molTemplateList[moleculeList[molItr].molTypeIndex].D.x;
		
		double Dcomplex_y=complexList[com1Index].D.y;
		double DIL_y=molTemplateList[moleculeList[pro2Index].molTypeIndex].D.y;
		double D2_y=1.0/(1.0/Dcomplex_y-1.0/DIL_y);
		if(Dcomplex_y==DIL_y)D2_y=molTemplateList[moleculeList[molItr].molTypeIndex].D.y;
		//For z, both are currently on the membrane, so Dcomplex_z and DIL_z==0
		double Dcomplex_z=complexList[com1Index].D.z;
		double DIL_z=molTemplateList[moleculeList[pro2Index].molTypeIndex].D.z;
		double D2_z=0;
		//If Links>1, then it is stuck in 2D even if this link breaks, so Dz=0.//otherwise, use the bound protein's Dz
		
		if(complexList[com1Index].linksToSurface == 1) D2_z=molTemplateList[moleculeList[molItr].molTypeIndex].D.z;
		
		Dtot=1.0/3.0 * (D2_x+DIL_x) +1.0/3.0 *(D2_y+DIL_y) +1.0/3.0* (D2_z+DIL_z);
		
		double sigma=forwardRxns[kfIndex].bindRadius;
		double h = params.timeStep;  
                prob = dissociate3D(h,Dtot,sigma,ka,kb); //membraneObject.Punbinding3D; 
            }else{ 
                // 2D->2D
		bool probValExists { false };
		int probMatrixIndex { 0 };
		double Dcomplex_x=complexList[com1Index].D.x;
		double DIL_x=molTemplateList[moleculeList[pro2Index].molTypeIndex].D.x;
		double D2_x=1.0/(1.0/Dcomplex_x-1.0/DIL_x);
		//If both are on membrane, D2 should be close to DIL
		if(Dcomplex_x==DIL_x)D2_x=DIL_x;//molTemplateList[moleculeList[molItr].molTypeIndex].D.x;
	    
		double Dcomplex_y=complexList[com1Index].D.y;
		double DIL_y=molTemplateList[moleculeList[pro2Index].molTypeIndex].D.y;
		double D2_y=1.0/(1.0/Dcomplex_y-1.0/DIL_y);
		if(Dcomplex_y==DIL_y)D2_y=DIL_y;//molTemplateList[moleculeList[molItr].molTypeIndex].D.y;
		//For z, both are currently on the membrane, so Dcomplex_z and DIL_z==0
		double Dcomplex_z=complexList[com1Index].D.z;
		double DIL_z=molTemplateList[moleculeList[pro2Index].molTypeIndex].D.z;
		double D2_z=0;
		//If Links>1, then it is stuck in 2D even if this link breaks, so Dz=0.//otherwise, use the bound protein's Dz
		
		Dtot=1.0/3.0 * (D2_x+DIL_x) +1.0/3.0 *(D2_y+DIL_y) +1.0/3.0* (D2_z+DIL_z);
		
		// declare intrinsic binding rate of 2D->2D case.
		double ktemp { ka /forwardRxns[kfIndex].length3Dto2D };
		for (int l = 0; l < IL2DbindingVec.size(); ++l) {
		    if (std::abs(ILTableIDs[l*3] - ktemp) < 1e-8 && std::abs(ILTableIDs[l*3+1] - Dtot) < 1E-4 && std::abs(ILTableIDs[l*3+2]- kb) < 1e-8) {
			probValExists = true;
			probMatrixIndex = l;
			break;
		    }
		}
		
		if (!probValExists) {
		    // first dimension out of i elements (2*i)
		    ILTableIDs.push_back(ktemp);
		    // second dimension (2*i+1)
		    ILTableIDs.push_back(Dtot);
		    ILTableIDs.push_back(kb);
		    paramsIL params2D {};
		    params2D.R2D = 0.0;
		    params2D.sigma = forwardRxns[kfIndex].bindRadius;
		    params2D.Dtot = Dtot;
		    params2D.ka = ktemp;
		    int backIndex = forwardRxns[kfIndex].conjBackRxnIndex;
		    params2D.kb = kb;//backRxns[backIndex].rateList[rateIndex].rate;       
		    params2D.area = membraneObject.totalSA;
		    params2D.dt = params.timeStep;
		    params2D.Nlipid = membraneObject.nSites; // the initial number of lipids' interfaces on the surface
		    params2D.Na = membraneObject.No_protein; // the initial number of protein's interfaces that can bind to surface
		    //		membraneObject.Pbinding2D = pimplicitlipid_2D(params2D);// binding probability, must multiply the lipid density
		    probMatrixIndex=IL2DbindingVec.size();
		    IL2DbindingVec.push_back(pimplicitlipid_2D(params2D));
		    IL2DUnbindingVec.push_back(dissociate2D(params2D)); // unbinding probability.
		    
		    std::cout.precision(12);
		    std::cout << "Create new IL 2D Un/binding prob (diss check): " << ktemp << ", Dtot: " << Dtot <<" kb: "<<kb<<" binding prob: "<<IL2DbindingVec[probMatrixIndex]<< " Unbinding prob: "<<IL2DUnbindingVec[probMatrixIndex]<<'\n';
		    

		}
		probValExists = false; // reset
		
                prob = IL2DUnbindingVec[probMatrixIndex]; //2D-2D unbinding          
            } 
            if (params.debugParams.forceDissoc)
                prob = 1.0;
	    
            double rnum { rand_gsl() };
            if (prob > rnum) {
            double rnum2 { rnum + rand_gsl() * Constants::iRandMax }; // to get higher resolution                   
            if (prob > rnum2){                          
                        std::cout << "Dissociation at iteration: " << simItr << " protein: " << molItr
                                  << " from the membrane, with probability: " << prob<<'\n';
                        std::cout << "Complex " << moleculeList[molItr].myComIndex << ", composed of "
                                  << complexList[moleculeList[molItr].myComIndex].memberList.size() << " molecules\n";
			std::cout <<" protein coords for : "<<molItr<<std::endl;
			moleculeList[molItr].display_my_coords("proteinonsurface");
                        /*std::cout << "Complex members:";
                        for (const auto& memMol : complexList[moleculeList[molItr].myComIndex].memberList)
                            std::cout << ' ' << memMol;
                        std::cout << '\n';
			*/

                        // 'break_interaction' frees certain protein and its interface 
                        break_interaction_implicitlipid(relIface1, relIface2, moleculeList[molItr], moleculeList[pro2Index],
                               backRxns[mu], emptyComList, moleculeList, complexList, molTemplateList); 

                        // Change the number of bound pairs in the system.
                        update_Nboundpairs(moleculeList[molItr].molTypeIndex, moleculeList[pro2Index].molTypeIndex, -1,
                            params, counterArrays);
                        //Update species copy numbers
			//std::cout <<"Add new to index:" <<backRxns[mu].productListNew[0].absIfaceIndex<<' '<<backRxns[mu].productListNew[1].absIfaceIndex;
			//std::cout <<"Subtract from index: "<<backRxns[mu].reactantListNew[0].absIfaceIndex<<std::endl;
                        counterArrays.copyNumSpecies[backRxns[mu].reactantListNew[0].absIfaceIndex] -= 1;
                        counterArrays.copyNumSpecies[backRxns[mu].productListNew[0].absIfaceIndex]  += 1;
                        counterArrays.copyNumSpecies[backRxns[mu].productListNew[1].absIfaceIndex]  += 1;
                        //counterArrays.copyNumSpecies[moleculeList[molItr].interfaceList[relIface1].index]  += 1;
                        //counterArrays.copyNumSpecies[moleculeList[pro2Index].interfaceList[relIface2].index]  += 1;
                        
                        // update the number of bonds that this complex has connected to the membrane surface.
			//this also needs to be done for the individual proteins.
                        complexList[com1Index].linksToSurface -= 1;
			moleculeList[molItr].linksToSurface--;
                        // consider the reflecting-surface movement   
                        Vector transVec1 {};
                        if (complexList[com1Index].linksToSurface == 0 && membraneObject.TwoD == false){ 
                            transVec1.x = 0;
                            transVec1.y = 0;
                            //transVec1.z = 0;
			    transVec1.z = (-membraneObject.waterBox.z/2.0+membraneObject.RS3D)-moleculeList[molItr].interfaceList[relIface1].coord.z;//move this interface to the RS3D position in z.
			    std::cout <<" unbinding POSITION UPDATE IN Z: "<<transVec1.z<<std::endl;
                            complexList[com1Index].OnSurface = false; // the complex is not bound to the membrane anymore.
                        }else {
			  transVec1.x = 0;
                            transVec1.y = 0;
                            transVec1.z = 0;
                            complexList[com1Index].OnSurface = true; 
                        }
			//BINDING CURRENT PUTS PROTEINS AT RS3D, SO NO DISPLACEMENT AFTER DISSOCIATIONS
                        // update the temporary coordinates. If binding places it at RS, do not move here!!
                        for (auto& mp : complexList[com1Index].memberList)
                             moleculeList[mp].update_association_coords(transVec1);  
                            
                        for (auto memMol : complexList[com1Index].memberList) {
                             moleculeList[memMol].comCoord = moleculeList[memMol].tmpComCoord;
                             for (unsigned int i { 0 }; i < moleculeList[memMol].interfaceList.size(); ++i)
                                  moleculeList[memMol].interfaceList[i].coord = moleculeList[memMol].tmpICoords[i];
                             moleculeList[memMol].clear_tmp_association_coords();
                         }
			//*/
                         complexList[com1Index].update_properties(moleculeList, molTemplateList); // recalculate the properties of the first complex
                        //reflect_complex_rad_rot(membraneObject, complexList[moleculeList[molItr].myComIndex], moleculeList);                         
                             
			 
			 std::cout << "Coords of p1 (COM) after dissociation: \n";// << moleculeList[molItr].comCoord << '\n';
			 moleculeList[molItr].display_my_coords("proteinatRS");
                        std::cout << "The partner was an implicit lipid" << '\n';
			complexList[com1Index].display();
                        // change the traj Status
                        //complexList[moleculeList[molItr].myComIndex].trajStatus = TrajStatus::propagated;
                        for (auto memMol : complexList[com1Index].memberList)
                            moleculeList[memMol].trajStatus = TrajStatus::propagated;
 
                        // TODO: Temporary implementation for destruction coupled to dissociation
                        /*Must be a C->A+B
                          Find out if A or B is being destroyed.
                         */
                        if (backRxns[mu].isCoupled) {
                        	 int destroyProIndex { -1 };
                            const CreateDestructRxn& coupledRxn
                                = createDestructRxns[backRxns[mu].coupledRxn.relRxnIndex]; // which reaction is being
                                                                                           // performed.
                            if (moleculeList[molItr].molTypeIndex == coupledRxn.reactantMolList[0].molTypeIndex) {
                                destroyProIndex = molItr; // Is it A or B?
                            } else{
                                destroyProIndex = pro2Index;
                            }

                            std::cout << "Performing coupled reaction.\n";
                            // decrement the copy number array for everything in complex
                            for (auto& memMol : complexList[moleculeList[destroyProIndex].myComIndex].memberList) {
                                for (auto& iface : moleculeList[memMol].interfaceList) {
                                    --counterArrays.copyNumSpecies[iface.index];
                                }
                            }
                            complexList[moleculeList[destroyProIndex].myComIndex].destroy(moleculeList, emptyMolList,
                                emptyComList); // destroying the entire complex that this molecule is a part of.

                            // remove the molecule from the SimulVolume subsCellList
                            // have this here to avoid circular header calls with SimulVolume and
                            // Molecule_Complex
                            for (auto itr = simulVolume.subCellList[moleculeList[destroyProIndex].mySubVolIndex].memberMolList.begin();
                                      itr!= simulVolume.subCellList[moleculeList[destroyProIndex].mySubVolIndex].memberMolList.end();
                                      ++itr) {
                                if (*itr == moleculeList[destroyProIndex].index) {
                                    simulVolume.subCellList[moleculeList[destroyProIndex].mySubVolIndex].memberMolList.erase(itr);
                                    break;
                                }
                            }
                            moleculeList[destroyProIndex].mySubVolIndex = -1; // reinitialize index

                            if (coupledRxn.isObserved) {
                                auto observeItr = observablesList.find(coupledRxn.observeLabel);
                                if (observeItr == observablesList.end()) {
                                    std::cerr << "WARNING: Observable " << coupledRxn.observeLabel << " not defined.\n";
                                } else {
                                    --observeItr->second;
                                }
                            }
                        } // finished with IsCoupled?

                        if (backRxns[mu].isObserved) {
                        	 auto obsItr = observablesList.find(backRxns[mu].observeLabel);
                            if (obsItr != observablesList.end())
                                --obsItr->second;
                        }
                  }
             }
        }
    }
}
