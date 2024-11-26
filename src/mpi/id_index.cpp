#include <iostream>
#include <unordered_map>

#include "macro.hpp"
#include "mpi/mpi_function.hpp"

using namespace std;

/*
// Hash maps to help translating IDs to indices faster:
unordered_map<size_t, size_t> Molecule::mapIdToIndex {};
unordered_map<size_t, size_t> Complex::mapIdToIndex {};

int find_complex(vector<Complex> &complexList, int id){
    size_t iCom = 0;
    // Fist check if the given id can already be found in the mapIdToIndex:
    if(Complex::mapIdToIndex.find(id) != Complex::mapIdToIndex.end()){
        iCom = Complex::mapIdToIndex[id];
        if( (iCom < complexList.size()) && (complexList[iCom].id == id) )
            return iCom;
    }
    // Search complexList for a complex with the same ID:
    iCom = 0;
    while( (iCom < complexList.size()) && (complexList[iCom].id != id) )
        iCom++;
    if( iCom == complexList.size() )
        return -1;
    // Store for future reference:
    Complex::mapIdToIndex[id] = iCom;
    //TODO: Empty old values once in e.g. 100000 iterations.
    return iCom;
}


int find_molecule(vector<Molecule> &moleculeList, int id){
    size_t iMol = 0;
    // Fist check if the given id can already be found in the mapIdToIndex:
    if(Molecule::mapIdToIndex.find(id) != Molecule::mapIdToIndex.end()){
        iMol = Molecule::mapIdToIndex[id];
        if( (iMol < moleculeList.size()) && (moleculeList[iMol].id == id) )
            return iMol;
    }
    // Search moleculeList for a molecule with the same ID:
    iMol = 0;
    while( (iMol < moleculeList.size()) && (moleculeList[iMol].id != id) )
        iMol++;
    if( iMol == moleculeList.size() )
        return -1;
    // Store for future reference:
    Molecule::mapIdToIndex[id] = iMol;
    return iMol;
}
*/

int find_complex(vector<Complex> &complexList, int id) {
  size_t iCom = 0;
  // Search complexList for a complex with the same ID:
  while ((iCom < complexList.size()) && (complexList[iCom].id != id)) iCom++;
  if (iCom == complexList.size()) return -1;
  return iCom;
}

int find_molecule(vector<Molecule> &moleculeList, int id) {
  size_t iMol = 0;
  // Search moleculeList for a molecule with the same ID:
  while ((iMol < moleculeList.size()) && (moleculeList[iMol].id != id)) iMol++;
  if (iMol == moleculeList.size()) return -1;
  return iMol;
}

void IDs_to_indices(MpiContext &mpiContext, vector<Molecule> &moleculeList,
                    vector<int> &indices) {
  if (VERBOSE) cout << "IDs_to_indices begins" << endl;
  for (auto &index : indices) {
    Molecule &mol = moleculeList[index];
    //        DEBUG_MOL("4A");
    // Update interfaces to store also partner-molecule ID and
    // re-populate bndpartner and bndlist:
    mol.bndpartner.clear();
    mol.bndlist.clear();
    for (size_t i = 0; i < mol.interfaceList.size(); i++) {
      Molecule::Interaction &inter = mol.interfaceList[i].interaction;
      // Update only interfaces whose partner is visible on this rank:
      if (inter.partnerId != -1) {  // if it can be seen
        // cout << "4A0 + " << inter.partnerId << endl;
        if (inter.partnerId == -3) {
          // partner is implicit lipid
          inter.partnerIndex = 0;
          inter.partnerId = 0;
          mol.bndpartner.push_back(inter.partnerIndex);
          mol.bndlist.push_back(i);
          continue;
        }
        inter.partnerIndex = find_molecule(moleculeList, inter.partnerId);
        if (inter.partnerIndex != -1) {  // if found local partner
                                         // cout << "4AA" << endl;
          //  Both have to be updated:
          mol.bndpartner.push_back(inter.partnerIndex);
          mol.bndlist.push_back(i);

          // Update interfaces of a partner as well:
          Molecule &partnerMol = moleculeList[inter.partnerIndex];
          for (size_t iPartner = 0; iPartner < partnerMol.interfaceList.size();
               iPartner++) {
            Molecule::Interaction &interPartner =
                partnerMol.interfaceList[iPartner].interaction;
            if (interPartner.partnerId ==
                mol.id) {  // if it is connected to mol
              interPartner.partnerIndex = mol.index;
              //  Add this mol as a partner only if it doesn't exist:
              size_t iBndpartner = 0;
              while ((iBndpartner < partnerMol.bndlist.size()) &&
                     (partnerMol.bndlist[iBndpartner] != iPartner))
                iBndpartner++;
              if (iBndpartner == partnerMol.bndpartner.size()) {
                partnerMol.bndpartner.push_back(mol.index);
                partnerMol.bndlist.push_back(iPartner);
              } else {
                partnerMol.bndpartner[iBndpartner] = mol.index;
              }
            }
          }
        }
      }
    }

    // map the prevlist of the molecule back to index
    // for (auto &prevId : mol.prevlist) {
    //   prevId = find_molecule(moleculeList, prevId);
    // }
  }
  if (VERBOSE) cout << "IDs_to_indices ends" << endl;
}

void indices_to_IDs(Molecule &mol, vector<Molecule> &moleculeList,
                    vector<Complex> &complexList) {
  if (VERBOSE) cout << "indices_to_IDs begins" << endl;
  // Update interfaces to store also partner-molecule ID:
  for (size_t i = 0; i < mol.interfaceList.size(); i++) {
    Molecule::Interaction &inter = mol.interfaceList[i].interaction;
    if (inter.partnerIndex != -1) {  // if it can be seen
      if (moleculeList[inter.partnerIndex].isImplicitLipid == false) {
        // Update only interfaces whose partner is visible on this rank:
        inter.partnerId = moleculeList[inter.partnerIndex].id;
        if (VERBOSE) cout << "Assigned: " << inter.partnerId << ". ";
      } else {
        inter.partnerId = -3;  // its partner is implicit lipid
      }
    }
  }
  // Set myComIndex to the ID of the complex:
  mol.myComIndex = complexList[mol.myComIndex].id;
  // Map the prevlist of the molecule to IDs:
  // for (auto &prevIndex : mol.prevlist) {
  //   prevIndex = moleculeList[prevIndex].id;
  // }
  if (VERBOSE) cout << "indices_to_IDs ends" << endl;
}
