#include "reactions/association/association.hpp"
#include "tracing.hpp"

/*Evaluates overlap between newly bound complex, based on distance only of COMs,
  If distance between them is < systems' params.overlapSepLimit (default 10),
  then cancel the move.
*/
void check_for_structure_overlap(bool& cancelAssoc, const Complex& reactCom1, const Complex& reactCom2,
    const std::vector<Molecule>& moleculeList, const Parameters& params,
    const std::vector<MolTemplate>& molTemplateList)
{
    // TRACE();
    double overlapTolerance { params.overlapSepLimit * params.overlapSepLimit };
    for (int memMol : reactCom1.memberList) {
        if (molTemplateList[moleculeList[memMol].molTypeIndex].checkOverlap) {

            for (int memMol2 : reactCom2.memberList) {
                if (molTemplateList[moleculeList[memMol2].molTypeIndex].checkOverlap) {
                    double dx = moleculeList[memMol].tmpComCoord.x - moleculeList[memMol2].tmpComCoord.x;
                    double dy = moleculeList[memMol].tmpComCoord.y - moleculeList[memMol2].tmpComCoord.y;
                    double dz = moleculeList[memMol].tmpComCoord.z - moleculeList[memMol2].tmpComCoord.z;
                    double r2 = dx * dx + dy * dy + dz * dz;
                    if (r2 < (overlapTolerance)) {
                        //   std::cout << "WARNING: Canceling association, complexes overlap. Distance COM-COM: " <<sqrt(r2)<<" mol1: "<<memMol<<" mol2 "<<memMol2<<"\n";
                        cancelAssoc = true;
                        return;
                    }
                }
            }
        }
    }
}
