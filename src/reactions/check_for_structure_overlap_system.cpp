#include "classes/class_Rxns.hpp"
#include "reactions/association/association.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"

/*  Checks to see if the centers of masses of any of the molecules that are undergoing physical association
 * overlap with any of the other molecules in the system! Only checks molecules that are flagged with checkOverlap=1
 *
 * If so, cancels association.
 *
 *For reactCom1 and reactCom2, their current coordinates are in tmpComCoord. For all other complexes, use their actual
 *coordinates
 */
void check_for_structure_overlap_system(bool& flag, const Complex& reactCom1, const Complex& reactCom2,
    std::vector<Molecule>& moleculeList, const Parameters& params,
    const std::vector<MolTemplate>& molTemplateList, const std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns)
{
    // TRACE();

    // int measure_overlap_system(int c1, int c2, Complex *ind_com, Fullmol *bases, Parms plist, Protein *wholep,
    // Complex *all_com, Fullmol *allbases, int c1Main, int c2Main) {

    int t1 = 0;
    int t2 = 0;
    int k;
    int j, mp, mp2;
    //    bool flag = false;
    int s1 = reactCom1.memberList.size(); // ind_com[c1].mysize;
    int s2 = reactCom2.memberList.size();
    double tol2 = params.overlapSepLimit * params.overlapSepLimit;
    double dx, dy, dz;
    double r2;
    /*Evaluate whether the two complexes that are binding to one another will, once they are oriented properly,
     generate structures with overlapping interfaces, against all other complexes in the system!
     It is possible that for a large complex, a large rotation will cause it to overlap with other complexes, not just
     the one it binds to.
    */

    /*No overlap between the two complexes found. But, now evaluate whether the new structure overlaps significantly
      with other structures that are in the simulation, as a result of large orientational changes*/
    for (int c = 0; c < complexList.size(); c++) {
        /*c is looping over all complexes, which index matches allbases
         */
        if (!complexList[c].isEmpty) {
            if (c != reactCom1.index && c != reactCom2.index) {
                // no self, and c1 vs c2 was done in check_for_structure_overlap()

                int sAll = complexList[c].memberList.size();
                for (int i = 0; i < sAll; i++) {
                    // loop over all proteins in complex c
                    int pp = complexList[c].memberList[i];
                    // only check overlap bewteen distinct types of proteins.
                    if (molTemplateList[moleculeList[pp].molTypeIndex].checkOverlap) {
                        double xm = moleculeList[pp].comCoord.x; // for proteins that were not just associated, use
                            // actual cooordinates
                        double ym = moleculeList[pp].comCoord.y;
                        double zm = moleculeList[pp].comCoord.z;

                        /*measure distance between the proteins in complex c and proteins in reactCom1 */

                        for (j = 0; j < s1; j++) {
                            mp = reactCom1.memberList[j]; // these are the newly rotated proteins.
                            if (molTemplateList[moleculeList[mp].molTypeIndex].checkOverlap) {

                                dx = moleculeList[mp].tmpComCoord.x
                                    - xm; // for proteins that just associated, still use tmp coords
                                dy = moleculeList[mp].tmpComCoord.y
                                    - ym; // for proteins that just associated, still use tmp coords
                                dz = moleculeList[mp].tmpComCoord.z
                                    - zm; // for proteins that just associated, still use tmp coords
                                r2 = dx * dx + dy * dy + dz * dz;
                                if (r2 < tol2) {
                                    flag = true;
                                    // i = sAll;//break i loop
                                    // j = s1;//break j loop
                                    // c=plist.ntotalcomplex;//break full loop
                                    // std::cout << " WARNING: Cancel association, overlap with other proteins COM in SYSTEM! "
                                    //           << mp << ' ' << pp << " SEPARATION: " << sqrt(r2) << std::endl;
                                    // std::cout << " complex moved size: " << s1 << '\n';
                                    // moleculeList[mp].display_assoc_icoords("moved_protein");
                                    // std::cout << "Complex in system, size:" << sAll << '\n';
                                    // moleculeList[pp].display_my_coords("pro_in_system");
                                    return;
                                }

                                double moleculeRad = molTemplateList[moleculeList[mp].molTypeIndex].radius
                                    + molTemplateList[moleculeList[pp].molTypeIndex].radius;
                                double molRadSq = moleculeRad * moleculeRad; //this is only used to cut down the number of times evaluating interface overlap.

                                if (r2 < molRadSq) {
                                    // The COMs are not close, but the binding interfaces might be. Check if these two
                                    // proteins have overlapping interfaces, not just COMs.
                                    //measure_overlap_free_protein_interfaces(moleculeList[pp], moleculeList[mp], flag, molTemplateList, forwardRxns, backRxns);
                                    measure_overlap_protein_interfaces(moleculeList[pp], moleculeList[mp], flag); // first one is actual coords, second one is tempCoords.
                                    if (flag == true) {
                                        // std::cout << " WARNING, CANCEL ASSOC: Protein iface in association overlaps protein in SYSTEM! " << mp
                                        //           << ' ' << pp << std::endl;
                                        // std::cout << " complex moved size: " << s1 << '\n';
                                        // moleculeList[mp].display_assoc_icoords("moved_protein");
                                        // std::cout << "Complex in system, size:" << sAll << '\n';
                                        // moleculeList[pp].display_my_coords("pro_in_system");

                                        return;
                                    }
                                }

                            } // if overlapcheck
                        } // all proteins in s1

                        /*measure distance between the proteins in complex c and proteins in reactCom2 */

                        for (j = 0; j < s2; j++) {
                            mp = reactCom2.memberList[j]; // these are the newly rotated proteins.
                            if (molTemplateList[moleculeList[mp].molTypeIndex].checkOverlap == 1) {

                                dx = moleculeList[mp].tmpComCoord.x
                                    - xm; // for proteins that just associated, still use tmp coords
                                dy = moleculeList[mp].tmpComCoord.y
                                    - ym; // for proteins that just associated, still use tmp coords
                                dz = moleculeList[mp].tmpComCoord.z
                                    - zm; // for proteins that just associated, still use tmp coords
                                r2 = dx * dx + dy * dy + dz * dz;
                                if (r2 < tol2) {
                                    flag = true;
                                    // i = sAll;//break i loop
                                    // j = s1;//break j loop
                                    // c=plist.ntotalcomplex;//break full loop
                                    // std::cout << " WARNING: Cancel association, overlap with other proteins COM in SYSTEM! "
                                    //           << mp << ' ' << pp << " SEPARATION: " << sqrt(r2) << std::endl;
                                    // std::cout << " complex moved size: " << s2 << '\n';
                                    // moleculeList[mp].display_assoc_icoords("moved_protein");
                                    // std::cout << "Complex in system, size:" << sAll << '\n';
                                    // moleculeList[pp].display_my_coords("pro_in_system");

                                    return;
                                }

                                double moleculeRad = molTemplateList[moleculeList[mp].molTypeIndex].radius
                                    + molTemplateList[moleculeList[pp].molTypeIndex].radius;
                                double molRadSq = moleculeRad * moleculeRad;

                                if (r2 < molRadSq) {
                                    // The COMs are not close, but the binding interfaces might be. Check if these two
                                    // proteins have overlapping interfaces, not just COMs.
                                    measure_overlap_free_protein_interfaces(moleculeList[pp], moleculeList[mp], flag, molTemplateList, forwardRxns, backRxns);
                                    //measure_overlap_protein_interfaces(moleculeList[pp], moleculeList[mp],flag); // first one is actual coords, second one is tempCoords.
                                    if (flag == true) {
                                        // std::cout << " WARNING, CANCEL ASSOC: Protein iface in association overlaps protein in SYSTEM! " << mp
                                        //           << ' ' << pp << std::endl;
                                        // std::cout << " complex moved size: " << s2 << '\n';
                                        // moleculeList[mp].display_assoc_icoords("moved_protein");
                                        // std::cout << "Complex in system, size:" << sAll << '\n';
                                        // moleculeList[pp].display_my_coords("pro_in_system");

                                        return;
                                    }
                                }

                            } // if overlapcheck
                        } // all proteins in s2

                    } // if overlapcheck, proteins in c
                } // all proteins in c
            } // if not c2 and not c1
        } // if complex is empty
    } // all complexes in the system

    return;
}
