#include "classes/class_Rxns.hpp"
#include "reactions/association/association.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"

/*
  Loop over all interfaces of base1 and then base2.
  If they are within bindrad of each other (just set here to 1!!!) then set flagCancel=1, so the association will be
  Cancelled. base1 is not associating in this step, it is part of the system, base2 is performing association this step.
  For baseTmp, access its Tmp coords, as it is testing its new orientation!
 */
void measure_overlap_protein_interfaces(Molecule base1, Molecule baseTmp, bool& flagCancel)
{
    // TRACE();
    for (unsigned i { 0 }; i < base1.interfaceList.size(); i++) {
        // just check all interfaces.
        for (int m = 0; m < baseTmp.interfaceList.size(); m++) {
            /*distance between these two interfaces*/
            double d2;
            double bindrad2
                = 1.0; // THIS SHOULD BE AN ACTUAL BINDING RADIUS squared BETWEEN INTERFACES WITHIN THIS PROTEIN
            double dx = base1.interfaceList[i].coord.x - baseTmp.tmpICoords[m].x;
            double dy = base1.interfaceList[i].coord.y - baseTmp.tmpICoords[m].y;
            double dz = base1.interfaceList[i].coord.z - baseTmp.tmpICoords[m].z;

            d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < bindrad2) {
                // std::cout << " WARNING: Cancel Association, overlap with INTERFACES in SYSTEM ! " << i << ' ' << m
                //           << " separation : " << sqrt(d2) << " from proteins: " << '\t';
                // cout <<" UNBOUND PROTEINSs "<<mp<<' '<<mp2 <<" Interfaces: "<<myIfaceIndex<<' '<<m<<  " WITHIN A
                // COMPLEX ARE CLOSE TOGETHER: "<<d<<endl;
                flagCancel = true;
            }

        } // check other protein interfaces

    } // loop over all interfaces pro1
}
