#include "system_setup/system_setup.hpp"

Complex initialize_complex(const Molecule& mol, const MolTemplate& molTemp)
{
  /*Creates a new complex (end of complexList array) to hold the newly generated molecule, assigns coords*/
  Complex tmp;
    tmp.comCoord = mol.comCoord;
    tmp.D = molTemp.D;
    tmp.Dr = molTemp.Dr;
    tmp.radius = molTemp.radius;
    tmp.index = Complex::numberOfComplexes;
    tmp.mass = molTemp.mass;
    tmp.memberList.push_back(mol.index);
    tmp.isEmpty = false;
    tmp.numEachMol = std::vector<int>(MolTemplate::numMolTypes);
    ++tmp.numEachMol[molTemp.molTypeIndex];
    tmp.lastNumberUpdateItrEachMol.resize(MolTemplate::numMolTypes);
    
    ++Complex::numberOfComplexes;
    
    return tmp;
}
