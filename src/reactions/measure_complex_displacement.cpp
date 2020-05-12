#include "classes/class_Rxns.hpp"
#include "reactions/association/association.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"

/*
  tmpComCoord hold temporary crds. 
  ComCoord hold pre-association crds.

  Measure how much each protein has moved from previous to new positions, and the complex_com as well.
  includes all proteins, not just checkOverlap proteins.
  CURRENTLY IS NOT RETURNING THE CANCELATION FLAG.
*/
void measure_complex_displacement(bool& flag, Complex& reactCom1, Complex& reactCom2,
    std::vector<Molecule>& moleculeList, const Parameters& params,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<Complex>& complexList)
{
    // TRACE();
    double dx, dy, dz;
    double R2;
    double LARGE_DISP = 1000.0; //nm
    double LDISPSQ = LARGE_DISP * LARGE_DISP;

    int mp;

    /*Distance between ind_com[c1] and all_com[c1Main], and ind_com[c2] and all_com[c2Main]
    protein and the complex COMs.
  */
    update_complex_tmp_com_crds(reactCom1, moleculeList);
    update_complex_tmp_com_crds(reactCom2, moleculeList); //update COM of tmp complex coords

    //measure distance between initial and final.
    dx = reactCom1.tmpComCoord.x - reactCom1.comCoord.x;
    dy = reactCom1.tmpComCoord.y - reactCom1.comCoord.y;
    dz = reactCom1.tmpComCoord.z - reactCom1.comCoord.z;

    R2 = dx * dx + dy * dy + dz * dz;
    if (R2 > LDISPSQ) {
        std::cout << "Cancel Association: NEW COMPLEX COM FOR ASSOCIATING " << reactCom1.index << " DISPLACED BY: " << sqrt(R2) << std::endl;
        flag = true;
    }
    dx = reactCom2.tmpComCoord.x - reactCom2.comCoord.x;
    dy = reactCom2.tmpComCoord.y - reactCom2.comCoord.y;
    dz = reactCom2.tmpComCoord.z - reactCom2.comCoord.z;

    R2 = dx * dx + dy * dy + dz * dz;
    if (R2 > LDISPSQ) {
        std::cout << "Cancel Association: NEW COMPLEX COM FOR ASSOCIATING " << reactCom2.index << " DISPLACED BY: " << sqrt(R2) << std::endl;
        flag = true;
    }

    //NOW LOOP OVER ALL PROTEINS in Com1
    for (int i = 0; i < reactCom1.memberList.size(); i++) {
        //proteins are in the same order in both lists.
        mp = reactCom1.memberList[i];
        dx = moleculeList[mp].tmpComCoord.x - moleculeList[mp].comCoord.x;
        dy = moleculeList[mp].tmpComCoord.y - moleculeList[mp].comCoord.y;
        dz = moleculeList[mp].tmpComCoord.z - moleculeList[mp].comCoord.z;

        R2 = dx * dx + dy * dy + dz * dz;
        // double moleculeRad=molTemplateList[moleculeList[mp].molTypeIndex].radius*2.0;//they are the exact same protein at different positions
        // double molRadSq=moleculeRad*moleculeRad;//don't use this-for lipids it will be close to zero.

        if (R2 > LDISPSQ) {
            std::cout << "Cancel Association: NEW PROTEIN COM FOR ASSOCIATING PROTEIN " << mp << " DISPLACED BY: " << sqrt(R2) << std::endl;
            flag = true;
        }

    } //all proteins in Com1
    //NOW LOOP OVER ALL PROTEINS in C2

    for (int i = 0; i < reactCom2.memberList.size(); i++) {
        //proteins are in the same order in both lists.
        mp = reactCom2.memberList[i];
        dx = moleculeList[mp].tmpComCoord.x - moleculeList[mp].comCoord.x;
        dy = moleculeList[mp].tmpComCoord.y - moleculeList[mp].comCoord.y;
        dz = moleculeList[mp].tmpComCoord.z - moleculeList[mp].comCoord.z;

        R2 = dx * dx + dy * dy + dz * dz;
        // double moleculeRad=molTemplateList[moleculeList[mp].molTypeIndex].radius*2.0;//they are the exact same protein at different positions
        // double molRadSq=moleculeRad*moleculeRad;//don't use this-for lipids it will be close to zero.

        if (R2 > LDISPSQ) {
            std::cout << "Cancel Association: NEW PROTEIN COM FOR ASSOCIATING PROTEIN " << mp << " DISPLACED BY: " << sqrt(R2) << std::endl;
            flag = true;
        }

    } //all proteins in Com2
}
