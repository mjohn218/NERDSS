/*! \file debug.hpp

 * ###
 * ### Purpose
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */

#pragma once

#include <cstring>
#include <vector>

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Observable.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_SimulVolume.hpp"
#include "classes/class_copyCounters.hpp"
#include "split.cpp"

using namespace std;

/*! \defgroup debug
 * \brief Functions for debugging.
 */

/*! \ingroup debug
 * \brief write the debug information to debugFileName
 */
void write_debug_information(MpiContext &mpiContext, long long int simItr,
                             ofstream &debugFile,
                             vector<Molecule> &moleculeList,
                             vector<Complex> &complexList,
                             vector<MolTemplate> &molTemplateList,
                             copyCounters &counterArrays, double simTime,
                             string s);

/*! \ingroup debug
 * \brief get x_bin of one molecule
 */
inline int get_x_bin(MpiContext &mpiContext, Molecule &mol) {
  return int((mol.comCoord.x + (*(mpiContext.membraneObject)).waterBox.x / 2) /
             (*(mpiContext.simulVolume)).subCellSize.x) -
         mpiContext.xOffset;
}

/*! \ingroup debug
 * \brief print the information about the target molecule
 */
void debug_print(MpiContext &mpiContext, Molecule &mol, string s);

/*! \ingroup debug
 * \brief print the information about the target complex
 */
void debug_print_complex(MpiContext &mpiContext, Complex &com, string s);

/*! \ingroup debug
 * \brief check if index of Molecule::emptyMolList is larger than
 * moleculeList.size()
 */
void debug_firstEmptyIndex(MpiContext &mpiContext, string errorMessage);

/*! \ingroup debug
 * \brief check consistence between Iface.interaction.partnerIndex and
 * mol.bndpartner
 */
void debug_bndpartner_interface(MpiContext &mpiContext, string errorMessage);

/*! \ingroup debug
 * \brief print a vector to standard output
 */
template <typename T>
void print_vector(vector<T> vec) {
  if (!vec.size()) {
    cout << "/";
    return;
  }
  cout << "(";
  bool first = true;
  for (auto &it : vec) {
    if (!first) cout << ",";
    first = false;
    cout << it;
  }
  cout << ")";
}

/*! \ingroup debug
 * \brief Check consistence of mol.myComIndex and com.memberList
 */
void debug_molecule_complex_missmatch(MpiContext &mpiContext,
                                      vector<Molecule> &moleculeList,
                                      vector<Complex> &complexList,
                                      string errorMessage);

void check_complex_index(MpiContext &mpiContext,
                                      vector<Molecule> &moleculeList,
                                      vector<Complex> &complexList);

/*! \ingroup debug
 * \brief Testing serialization and deserialization is done just for debugging
 * purposes, solely on rank 0:
 */
void test_serialization(MpiContext &mpiContext, vector<Molecule> &moleculeList,
                        SimulVolume &simulVolume, Membrane &membraneObject,
                        vector<MolTemplate> &molTemplateList,
                        Parameters &params, vector<ForwardRxn> &forwardRxns,
                        vector<BackRxn> &backRxns,
                        vector<CreateDestructRxn> &createDestructRxns,
                        copyCounters &counterArrays,
                        vector<Complex> &complexList, unsigned char *arrayRank);

/*! \ingroup debug
 * \brief print Code section (number, text), rank and simulation iteration
 */
void LOC(float no, char *text, int proc, int simItr = -1);
