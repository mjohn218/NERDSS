#include "classes/class_Rxns.hpp"
#include "reactions/association/association.hpp"


void update_complex_tmp_com_crds(Complex &reactCom,  std::vector<Molecule>& moleculeList)
{

    // update center of mass
    double totMass { 0 };
    reactCom.tmpComCoord.x=0.0;
    reactCom.tmpComCoord.y=0.0;
    reactCom.tmpComCoord.z=0.0;
    
    for (auto& memMol : reactCom.memberList) {
        totMass += moleculeList[memMol].mass;
        reactCom.tmpComCoord.x += moleculeList[memMol].tmpComCoord.x * moleculeList[memMol].mass;
        reactCom.tmpComCoord.y += moleculeList[memMol].tmpComCoord.y * moleculeList[memMol].mass;
        reactCom.tmpComCoord.z += moleculeList[memMol].tmpComCoord.z * moleculeList[memMol].mass;
    }
    reactCom.tmpComCoord /= totMass;
    

}
