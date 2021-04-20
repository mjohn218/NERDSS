#include "classes/class_bngl_parser.hpp"
#include "reactions/bimolecular/2D_reaction_table_functions.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"

#include <algorithm>

void write_debug_info(int pro1Index, int pro2Index, int iface1Index, int iface2Index, const ForwardRxn& oneRxn,
    const copyCounters& counterArrays, double R1, double Keq, double probvec1, int Nc, double pdiss, double self,
    double Volume)
{
    // TRACE();
    // just writes debug info to cout, function call is easier to comment out
    // std::cout << " WITHIN COMPLEX: Proteins are R1=: " << R1 << " < bindRadSameCom: " << oneRxn.bindRadSameCom
    //           << " Proteins: " << pro1Index << ' ' << pro2Index << " interfaces: " << iface1Index << ' ' << iface2Index
    //           << " set probvec: " << probvec1 << " pdiss: " << pdiss << " Keq/V: " << Keq / Volume
    //           << " copy A: " << counterArrays.copyNumSpecies[iface2Index] << " copy C: " << Nc << " self: " << self
    //           << '\n';
}

void evaluate_binding_within_complex(int pro1Index, int pro2Index, int iface1Index, int iface2Index, int rxnIndex,
    int rateIndex, bool isBiMolStateChange, const Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, const std::vector<MolTemplate>& molTemplateList, const ForwardRxn& oneRxn,
    const std::vector<BackRxn>& backRxns, Membrane& membraneObject, copyCounters& counterArrays)
{
    if (moleculeList[pro1Index].trajStatus != TrajStatus::propagated
        && moleculeList[pro2Index].trajStatus != TrajStatus::propagated) {
        if (oneRxn.rateList[rateIndex].rate > 0) {

            /*Don't use Rmax, because it is based on diffusion*/

            /*Separate 3D and 2D?*/
            bool withinRmax {};
            double sep { 0 };
            double R1 { 0 };
            double bindSep { oneRxn.bindRadSameCom * oneRxn.bindRadius };
            withinRmax = get_distance(pro1Index, pro2Index, iface1Index, iface2Index, rxnIndex, rateIndex,
                isBiMolStateChange, sep, R1, bindSep, complexList, oneRxn, moleculeList, membraneObject);
            /*Proteins are in the same complex!*/
            /*If these proteins are at contact (sep=0) or close to at contact (<bindrad) then
            they will associate.
            Otherwise, they are in the same complex, but the binding interfaces are not close enought together, and they
            are not diffusing, so probvec should be zero!
            */
            /*Try using a detailed balance expression, where p_u->b=p_b->u *p_bound/p_unbound. Ratio of states
            could be N_C/(N_A*N_B) for equilibrium. p_b->u = 1-exp(-kb*deltat)
            p_b=Keq/V/ (1+Keq/V). Thus p_b/p_u=Keq/V/
            For SELF, Keq=ka/2/kb! FOR DISTINCT, Keq=ka/kb.
            */
            // cout <<" Proteins are in the same complex evaluate binding: sep is: "<<sep<<endl;

            /*Because this is closing a loop, introduce cooperativity in that the Keq is not just due to forming the one
             * interface, but could be stronger.
             */
            double coop = oneRxn.loopCoopFactor;
            if (withinRmax) {
                int backIndex = oneRxn.conjBackRxnIndex;

                double ka = oneRxn.rateList[rateIndex].rate; // units of nm^3/us
                // this uses the rate corresponding to the same forward reaction rate.
                double kb = backRxns[backIndex].rateList[rateIndex].rate; // units of per s
                // copy numbers of the product state
                // int Nc = counterArrays.copyNumSpecies[prod];

                // double pdiss = 1 - exp(-params.timeStep * Nc * kb * 1E-6);
                // If SELF, KD=kb*2/ka. Else, KD=kb/ka.
                double self = 1.0;
                // double Volume, Ka;
                if (iface1Index == iface2Index)
                    self = 2.0;
                // if (complexList[moleculeList[pro1Index].myComIndex].D.z == 0) { // This is a reaction on the membrane
                //                     // here again, assuming the forward reaction has only 1 rate
                //                     Keq = ka / (2.0 * oneRxn.bindRadius) / kb * 1E6 / self * coop;
                //                     Volume = params.waterBox.area; // nm2
                //                 } else {
                //                     Keq = ka / kb / self * 1E6 * coop;
                //                     Volume = params.waterBox.volume; // nm3
                //                 }
                /*Passoc needs to be
                  pdiss=1-exp(-kb*dt*Nc)
                  passoc=pdiss*ndiss*B(i)/C(i)*KaV
                  KaV=ka/kb/V;
                  ndiss max(1, Nevents_last_step). So either 1 or >1.
                  If you do this for each A molecule, then multiple by Nb. If you do it for each B molecule, multiply by
                  Na. Here we need copy numbers of each reactant.
                */
                double c0_nm3 = 0.602; // standard state 1M in units of /nm^3
                //double m = kb / (ka * 1E6 * c0_nm3) * 2.0; // 1E6 is to convert ka to units of /s
                double rateClose = ka * c0_nm3 * coop / self; //(1.0 + m); // ka*c0/(1+m) where m is defined above

                // rateClose will be in units of /us (due to ka units), so no need for 1E-6 factor, since timeStep has
                // units of us!!!!!
                double probvec1 = 1 - exp(-params.timeStep * rateClose);
                if (oneRxn.irrevRingClosure) {
                    // std::cout << "Ring closure reaction probability for species " << pro1Index << " and " << pro2Index
                    //           << " in complex " << moleculeList[pro1Index].myComIndex << " is set to 1.0.\n";
                    probvec1 = 1.0;
                }
                // pdiss * Keq * counterArrays.copyNumSpecies[iface2Index] / Nc / Volume;
                // if (Nc == 0) {
                //                     probvec1 = Keq / Volume * counterArrays.copyNumSpecies[iface2Index]; //
                //                     pdiss/Nc=1
                //                 }

                // std::cout<<"LOOP CLOSURE RATE, /us: "<<rateClose<<" m value: "<<m<<" c0 value: "<<c0_nm3<<" reaction
                // prob: "<<probvec1<<std::endl; write_debug_info(pro1Index, pro2Index, iface1Index, iface2Index,
                // oneRxn, counterArrays, R1, Keq, probvec1, Nc, pdiss, self, Volume);

                int proA = pro1Index;
                int ifaceA = iface1Index;
                int proB = pro2Index;
                int ifaceB = iface2Index;
                if (pro1Index > pro2Index) {
                    // only store for one of the two proteins.
                    proA = pro2Index;
                    ifaceA = iface2Index;
                    proB = pro1Index;
                    ifaceB = iface1Index;
                }

                moleculeList[proA].probvec.push_back(probvec1 * 1.0);
                moleculeList[proB].probvec.push_back(probvec1 * 1.0);

                /*Store all the reweighting numbers for next step.*/
                moleculeList[proA].currprevsep.push_back(R1);
                moleculeList[proA].currlist.push_back(proB);
                moleculeList[proA].currmyface.push_back(ifaceA);
                moleculeList[proA].currpface.push_back(ifaceB);
                moleculeList[proA].currprevnorm.push_back(1.0);
                moleculeList[proA].currps_prev.push_back(1.0 - probvec1 * 1.0);
            }
        } // rate is >0
    } // did not just dissociate
}
