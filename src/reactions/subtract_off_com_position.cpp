#include "reactions/association/association.hpp"

//for one molecule, set COM to zero, and all interfaces relative to this.
void subtract_off_com_position(Molecule &base1)
{
  
    for(int i=0;i<base1.interfaceList.size();i++){
      	base1.interfaceList[i].coord=base1.interfaceList[i].coord-base1.comCoord;
	      
    }
    base1.comCoord=Coord {0.0, 0.0, 0.0};//set to zero
    
    
}
