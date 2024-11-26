#include <numeric>

#include "system_setup/system_setup.hpp"

Molecule initialize_molecule(int comIndex, const Parameters& params,
                             const MolTemplate& molTemplate,
                             const Membrane& membraneObject) {
  Molecule tmp{};
  tmp.molTypeIndex = molTemplate.molTypeIndex;
  tmp.mass = molTemplate.mass;
  tmp.isLipid = molTemplate.isLipid;
  tmp.myComIndex = comIndex;
  tmp.id = Molecule::maxID++;

  // Set up interface state vectors
  tmp.freelist = std::vector<int>(molTemplate.interfaceList.size());
  std::iota(tmp.freelist.begin(), tmp.freelist.end(), 0);

  // Create center of mass
  tmp.interfaceList =
      std::vector<Molecule::Iface>(molTemplate.interfaceList.size());
  tmp.create_random_coords(molTemplate, membraneObject);

  // iterate number of molecules in the system and set index
  tmp.index = Molecule::numberOfMolecules;
  ++Molecule::numberOfMolecules;

  // keep track of molecule types
  ++MolTemplate::numEachMolType[molTemplate.molTypeIndex];

  return tmp;
}
