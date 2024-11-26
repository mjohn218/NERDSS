#include <cstring>
#include <vector>

#include "debug/debug.hpp"
#include "error/error.hpp"
#include "macro.hpp"
#include <stdexcept>

using namespace std;

// Function debug_molecule_complex_missmatch should find inconsistencies.
// It should be called only in the debugging mode, because of the amount of
// processing it includes.
void debug_molecule_complex_missmatch(MpiContext &mpiContext,
                                      vector<Molecule> &moleculeList,
                                      vector<Complex> &complexList,
                                      string errorMessage) {
  if (DEBUG) {
    // Check each molecule:
    for (auto &mol : moleculeList) {
      if (mol.isImplicitLipid) continue;  // no checking for potential implicit lipid
      if (mol.isEmpty) continue;
      // myComIndex has to be within complexList:
      if (mol.myComIndex > complexList.size()) {
        cout << "!!!!!!!!!!!!!!!!! debug_molecule_complex_missmatch "
                "!!!!!!!!!!!!!!!!!\n"
             << errorMessage << endl;
        cout << "mol.myComIndex = " << mol.myComIndex
             << ", complexList.size() == " << complexList.size() << endl;
        error(mpiContext, "3.1: mol(id=" + to_string(mol.id) +
                              ").myComIndex > complexList.size(), " +
                              to_string(mol.myComIndex) + ">" +
                              to_string(complexList.size()) + errorMessage);
      }
      // Find corresponding complex of the molecule
      // and check whether it has mol.index within memberList:
      Complex &com = complexList[mol.myComIndex];
      bool found = false;
      for (auto &it : com.memberList) {
        if (it == mol.index) found = true;
      }
      if (!found) {
        cout << "!!!!!!!!!!!!!!!!! debug_molecule_complex_missmatch "
                "!!!!!!!!!!!!!!!!!\n"
             << errorMessage << endl;
        cout << "mol.myComIndex = " << mol.myComIndex << endl;
        cout << "com.memberList = ";
        for (auto &it : com.memberList) {
          cout << it << ", ";
        }
        cout << endl;
        error(mpiContext, "3.5: mol(id=" + to_string(mol.id) +
                              ").myComIndex not found in complex.memberList com.id= " +
                              to_string(com.id) + errorMessage);
      }
    }
    // Check each complex:
    for (auto &com : complexList) {
      if (com.isEmpty) continue;  // no checking for deleted complexes
      // Check each memberList element:
      for (auto &molIndex : com.memberList) {
        if (moleculeList[molIndex].isImplicitLipid) continue;  // no checking for potential implicit lipid
        // Molecule index has to be within moleculeList:
        if (molIndex > moleculeList.size()) {
          cout << "!!!!!!!!!!!!!!!!! debug_molecule_complex_missmatch "
                  "!!!!!!!!!!!!!!!!!\n"
               << errorMessage << endl;
          cout << "com.memberList element = " << molIndex
               << ", moleculeList.size() == " << moleculeList.size() << endl;
          error(mpiContext,
                "3.2: complex memberList element > moleculeList.size(), " +
                    to_string(molIndex) + ">" + to_string(moleculeList.size()) +
                    errorMessage);
        }
      }
    }
    // Check each complex:
    for (auto &com : complexList) {
      if (com.isEmpty) continue;  // no checking for deleted complexes
      // Check each memberList element:
      for (auto &molIndex : com.memberList) {
        if (!molIndex) continue;  // no checking for potential implicit lipid
        // Ensure myComIndex refers to com:
        Molecule &mol = moleculeList[molIndex];
        if (mol.isEmpty) continue;  // no checking for deleted molecules
        if (mol.myComIndex != com.index) {
          cout << "!!!!!!!!!!!!!!!!! debug_molecule_complex_missmatch "
                  "!!!!!!!!!!!!!!!!!\n"
               << errorMessage << endl;

          cout << "EEE: com.index = " << com.index << ", com.id = " << com.id
               << ", member mol index = " << molIndex << ", mol.id = " << mol.id
               << ", moleculeList[" << molIndex
               << "].myComIndex = " << moleculeList[molIndex].myComIndex
               << endl;

          error(mpiContext,
                "3.6: mol.myComIndex != com.index, " +
                    to_string(mol.myComIndex) + "!=" + to_string(com.index) +
                    errorMessage);
        }
      }
    }
  }
}

void check_complex_index(MpiContext &mpiContext,
                                      vector<Molecule> &moleculeList,
                                      vector<Complex> &complexList) {
  // check each molecule has the correct complex index
  for (auto &mol : moleculeList) {
    int complexNum = complexList.size();
    if (mol.isEmpty) continue;
    if (mol.myComIndex < 0 || mol.myComIndex >= complexNum) {
      mol.print(mpiContext);
      std::cout.flush();
      // throw error
      throw std::runtime_error("Molecule has an incorrect complex index");
    }
  }
}