#include "boundary_conditions/reflect_functions.hpp"
#include "classes/class_Rxns.hpp"
#include "io/io.hpp"
#include "reactions/association/association.hpp"
#include "reactions/shared_reaction_functions.hpp"

/*Calculate COM of two complexes that have not yet associated
  applies to the temporrary protein coordinates, not final!
 */
void com_of_two_tmp_complexes(
    Complex& reactCom1, Complex& reactCom2, Coord& vectorCOM, std::vector<Molecule>& moleculeList)
{

    /*Calculate the COM of two separate complexes*/
    double totalmass = 0;
    vectorCOM.x = 0.0;
    vectorCOM.y = 0.0;
    vectorCOM.z = 0.0;
    int mp;
    for (int i = 0; i < reactCom1.memberList.size(); i++) {
        mp = reactCom1.memberList[i];
        totalmass += moleculeList[mp].mass; // mass of a protein is currently set based on its radius!

        vectorCOM.x += moleculeList[mp].tmpComCoord.x * moleculeList[mp].mass;
        vectorCOM.y += moleculeList[mp].tmpComCoord.y * moleculeList[mp].mass;
        vectorCOM.z += moleculeList[mp].tmpComCoord.z * moleculeList[mp].mass;
    }
    for (int i = 0; i < reactCom2.memberList.size(); i++) {
        mp = reactCom2.memberList[i];
        totalmass += moleculeList[mp].mass;
        vectorCOM.x += moleculeList[mp].tmpComCoord.x * moleculeList[mp].mass;
        vectorCOM.y += moleculeList[mp].tmpComCoord.y * moleculeList[mp].mass;
        vectorCOM.z += moleculeList[mp].tmpComCoord.z * moleculeList[mp].mass;
    }
    vectorCOM /= totalmass;
}
