#include "error/error.hpp"

#include <cstring>
#include <iostream>

using namespace std;

void error(string errorString) {
  cerr << "!!! ################################## !!! "
          "################################## !!!"
       << endl;
  cerr << "Error: " << errorString << "!!!" << endl;
  cerr << "!!! ################################## !!! "
          "################################## !!!"
       << endl;
  exit(1);
}

void error(MpiContext &mpiContext, string errorString) {
  error(errorString + ": Rank=" + to_string(mpiContext.rank) + ":\n");
}

void error(MpiContext &mpiContext, Molecule &mol, string errorString) {
  cerr << "!!! ################################## !!! "
          "################################## !!!"
       << endl;
  cerr << "mol.id = " << mol.id << endl;
  cerr << "mol.index = " << mol.index << endl;
  cerr << "moleculeList.size() = " << (*(mpiContext.moleculeList)).size()
       << endl;
  cerr << "mol.myComIndex = " << mol.myComIndex << endl;
  cerr << "mol.my complex id = "
       << (*(mpiContext.complexList))[mol.myComIndex].id << endl;
  error(mpiContext, errorString);
}