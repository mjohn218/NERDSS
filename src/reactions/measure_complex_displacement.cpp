#include "classes/class_Rxns.hpp"
#include "reactions/association/association.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"

/*
  tmpComCoord hold temporary crds. 
  ComCoord hold pre-association crds.

  Measure how much each protein has moved from previous to new positions, and the complex_com as well.
  includes all proteins, not just checkOverlap proteins.

*/
void measure_complex_displacement(bool& flag, Complex& reactCom1, Complex& reactCom2,
    std::vector<Molecule>& moleculeList, const Parameters& params,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<Complex>& complexList)
{
    // TRACE();
    double dx, dy, dz;
    double R2;
    double Dtot1, Dtot2, cf, Dr1, Dr2;
    bool in2D = false;
    double dim1 = 3; //dimensionality
    double dim2 = 3; //dimensionality

    /*Assume diffusion is isotropic!*/
    Dtot1 = reactCom1.D.x; //Dtot1 = 1.0 / 3.0 * (reactCom1.D.x) + 1.0 / 3.0 * (reactCom1.D.y) + 1.0 / 3.0 * (reactCom1.D.z);
    Dtot2 = reactCom2.D.x; //Dtot2 = 1.0 / 3.0 * ( reactCom2.D.x) + 1.0 / 3.0 * ( reactCom2.D.y) + 1.0 / 3.0 * ( reactCom2.D.z);

    /*Calculate a mean square displacement for each complex given D and dimensionality.*/
    if (reactCom1.D.z < 1e-16) {
        //complex 1 is in 2D
        dim1 = 2;
    }
    if (reactCom2.D.z < 1e-16) {
        //complex 2 is in 2D
        dim2 = 2;
    }
    if (dim1 == 2 && dim2 == 2)
        in2D = true;

    /*rotational displacement*/
    cf = cos(sqrt(2.0 * (dim1 - 1) * reactCom1.Dr.z * params.timeStep));
    Dr1 = 2.0 * reactCom1.radius * reactCom1.radius * (1.0 - cf);
    cf = cos(sqrt(2.0 * (dim2 - 1) * reactCom2.Dr.z * params.timeStep));
    Dr2 = 2.0 * reactCom2.radius * reactCom2.radius * (1.0 - cf);

    Dtot1 += Dr1 / (2.0 * dim1 * params.timeStep); //in 2D, use 4
    Dtot2 += Dr2 / (2.0 * dim2 * params.timeStep); //in 2D, use 4, 3D, use 6

    double avgDisp1 = sqrt(2.0 * dim1 * Dtot1 * params.timeStep); //from Einstein relation.
    double avgDisp2 = sqrt(2.0 * dim2 * Dtot2 * params.timeStep); //from Einstein relation.

    /*calculated average displacement of both complexes in the step, due to both translational and rotational diffusion. Use this distance multiplied by params.scaleMaxDisplace to decide if motion is too large */
    double LARGE_DISP1 = params.scaleMaxDisplace * avgDisp1; //nm
    double LDISP1SQ = LARGE_DISP1 * LARGE_DISP1;
    double LARGE_DISP2 = params.scaleMaxDisplace * avgDisp2; //nm
    double LDISP2SQ = LARGE_DISP2 * LARGE_DISP2;

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
    if (R2 > LDISP1SQ) {
        // std::cout << "Cancel Association: NEW COMPLEX COM1 FOR ASSOCIATING " << reactCom1.index << " DISPLACED BY: " << sqrt(R2) << " max displace is : " << LARGE_DISP1 << std::endl;
        flag = true;
        return;
    } else {
        // std::cout << "ASSOCIATION COMPLEX1 DISPLACED BY: " << sqrt(R2) << " max displace is: " << LARGE_DISP1 << "Dtot is : "<<Dtot1<<" Dr1: "<<Dr1<<" complex_rad: "<<reactCom1.radius<<std::endl;
    }
    dx = reactCom2.tmpComCoord.x - reactCom2.comCoord.x;
    dy = reactCom2.tmpComCoord.y - reactCom2.comCoord.y;
    dz = reactCom2.tmpComCoord.z - reactCom2.comCoord.z;

    R2 = dx * dx + dy * dy + dz * dz;
    if (R2 > LDISP2SQ) {
        // std::cout << "Cancel Association: NEW COMPLEX COM FOR ASSOCIATING " << reactCom2.index << " DISPLACED BY: " << sqrt(R2) << " max displace is : " << LARGE_DISP2 << std::endl;
        flag = true;
        return;
    } else {
        // std::cout << "ASSOCIATION COMPLEX DISPLACED BY: " << sqrt(R2) << " max displace is " << LARGE_DISP2 <<" Dtot is : "<<Dtot2<<" Dr2: "<<Dr2<<" complex_rad: "<<reactCom2.radius<<std::endl;
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

        if (R2 > LDISP1SQ) {
            // std::cout << "Cancel Association: NEW PROTEIN COM1 FOR ASSOCIATING PROTEIN " << mp << " DISPLACED BY: " << sqrt(R2) << " max displace is : " << LARGE_DISP1 << std::endl;
            flag = true;
            return;
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

        if (R2 > LDISP2SQ) {
            // std::cout << "Cancel Association: NEW PROTEIN COM2 FOR ASSOCIATING PROTEIN " << mp << " DISPLACED BY: " << sqrt(R2) << " max displace is " << LARGE_DISP2 << std::endl;
            flag = true;
            return;
        }

    } //all proteins in Com2
}
