#include "debug/debug.hpp"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "error/error.hpp"
#include "macro.hpp"

using namespace std;

// Function called by macros: DEBUG_MOL and DEBUG_FIND_MOL
void debug_print(MpiContext &mpiContext, Molecule &mol, string s) {
  // Limit printing to iterations of interest:
  // if (mpiContext.simItr < 113000) return;

  // Do not print for deleted molecules:
  // if (mol.myComIndex == -1) return;

  vector<int> targetMolIds = {};  // 29166, 31989};

  const char *TrajStatusTypes[] = {"none",       "dissociate",     "associated",
                                   "propagated", "canBeResampled", "empty"};

  // Condition for printing data related molecule/molecules:
  if (find(targetMolIds.begin(), targetMolIds.end(), mol.id) !=
      targetMolIds.end()) {
    mol.print(mpiContext);
    // By uncommenting lines, one can choose what to print out:
    cout << "OOOOOOOO: " << s << "simItr= " << mpiContext.simItr
         << ", rank=" << mpiContext.rank << ", mol.id="
         << mol.id
         //<< ", mol.index = " << mol.index
         //<< ", mol.receivedFromNeighborRank=" << mol.receivedFromNeighborRank
        //  << ", partnerId = " << mol.interfaceList[0].interaction.partnerId
         << ", trajStatus= "
         << TrajStatusTypes[static_cast<std::underlying_type<TrajStatus>::type>(
                mol.trajStatus)]
         << ", get_x_bin=" << get_x_bin(mpiContext, mol) << "/"
         << mpiContext.endGhosted - mpiContext.xOffset
         << ", mol.comCoord.x=" << mol.comCoord.x
         << ", mol.comCoord.y=" << mol.comCoord.y
         << ", mol.comCoord.z=" << mol.comCoord.z << ", complex.comCoord.x="
         << (*mpiContext.complexList)[mol.myComIndex].comCoord.x
         << ", complex.comCoord.y="
         << (*mpiContext.complexList)[mol.myComIndex].comCoord.y
         << ", complex.comCoord.z="
         << (*mpiContext.complexList)[mol.myComIndex].comCoord.z
         //<< ", mol.bndpartner.size()=" << mol.bndpartner.size()
         << ", mol.isGhosted=" << mol.isGhosted
         << ", mol.isEmpty=" << mol.isEmpty
         << ", mol.myComIndex=" << mol.myComIndex << ", mol.my complex id="
         << (*mpiContext.complexList)[mol.myComIndex].id
         << ", mol.my complex receivedFromNeighborRank="
         << (*mpiContext.complexList)[mol.myComIndex].receivedFromNeighborRank
         << ", myComIndex's memberList size="
         << (*(mpiContext.complexList))[mol.myComIndex].memberList.size()
         << endl;

    // loop the interfaceList to print all the partnerIds
    cout << "partnerIds: (";
    for (auto &tmpIface : mol.interfaceList) {
      cout << tmpIface.interaction.partnerId << ", ";
    }
    cout << ")" << endl;

    if (mol.isEmpty) return;

    //*     // Printing complex members:
    // if(1){
    if (s[0] != '1') {  // if not just deserialized, before updating indices
                        // based on IDs
      std::cout << " # Complex members: ";
      if (s[0] != '1' && mol.myComIndex > (*(mpiContext.complexList)).size())
        error(
            mpiContext, mol,
            s + ": 27: mol.myComIndex > (*(mpiContext.complexList)).size(), " +
                to_string(mol.myComIndex) + ">" +
                to_string((*(mpiContext.complexList)).size()));
      for (auto &it : (*(mpiContext.complexList))[mol.myComIndex].memberList) {
        if (s[0] != '1' && it >= (*(mpiContext.moleculeList)).size())
          error(s + ": 28: it >= (*(mpiContext.moleculeList)).size(), " +
                to_string(it) +
                ">=" + to_string((*(mpiContext.moleculeList)).size()));
        std::cout << (*(mpiContext.moleculeList))[it].id << "; ";
      }
      std::cout << endl;
    }
    //*/
  }
}

// Function called by macros: DEBUG_COMPLEX and DEBUG_FIND_COMPLEX
void debug_print_complex(MpiContext &mpiContext, Complex &com, string s) {
  // Limit printing to iterations of interest:
  // if (mpiContext.simItr < 3500) return;

  // target complex id:
  // initialize a vector targetComplexIds storing the target complex ids
  vector<int> targetComplexIds = {};  // 29166};

  const char *TrajStatusTypes[] = {"none",       "dissociate",     "associated",
                                   "propagated", "canBeResampled", "empty"};

  // Condition for printing data related to complex/complexes:
  // check if com.id is in the targetComplexIds
  if (find(targetComplexIds.begin(), targetComplexIds.end(), com.id) !=
      targetComplexIds.end()) {
    // By uncommenting lines, one can choose what to print out:
    cout << "OOOOOOOO: " << s << "simItr= " << mpiContext.simItr
         << ", rank=" << mpiContext.rank << ", com.id=" << com.id
         << ", com.index = " << com.index
         << ", com.receivedFromNeighborRank=" << com.receivedFromNeighborRank
         << ", com.deleteIfNotReceivedBack=" << com.deleteIfNotReceivedBack
         << ", trajStatus= "
         << TrajStatusTypes[static_cast<std::underlying_type<TrajStatus>::type>(
                com.trajStatus)]
         << ", com.comCoord.x=" << com.comCoord.x
         << ", com.comCoord.y=" << com.comCoord.y
         << ", com.comCoord.z=" << com.comCoord.z
         << ", com.isEmpty=" << com.isEmpty << endl;

    if (com.isEmpty) return;

    if (1) {
      std::cout << " # Complex members id (index): ";
      for (auto &it : com.memberList) {
        std::cout << (*(mpiContext.moleculeList))[it].id << "(" << it << "); ";
      }
      std::cout << endl;
    }
  }
}

void debug_bndpartner_interface(MpiContext &mpiContext, string errorMessage) {
  if (!DEBUG) return;
  auto &moleculeList = *mpiContext.moleculeList;
  for (auto &mol : moleculeList) {
    if (mol.isImplicitLipid == true) {
      continue;
    }
    vector<int> vec;
    for (auto &tmpIface : mol.interfaceList) {
      if ((tmpIface.interaction.partnerIndex >= moleculeList.size()) &&
          (tmpIface.interaction.partnerIndex != -1))
        error("35: tmpIface.interaction.partnerIndex >= moleculeList.size(), " +
              to_string(tmpIface.interaction.partnerIndex) +
              ">=" + to_string(moleculeList.size()));
      vec.push_back(tmpIface.interaction.partnerIndex);

      if (tmpIface.interaction.partnerIndex != -1) {
        bool found = false;
        for (auto &it : mol.bndpartner) {
          if (it == tmpIface.interaction.partnerIndex) found = true;
        }
        if (!found) {
          error(
              "29: mol[id=" + to_string(mol.id) +
              "].tmpIface.interaction.partnerIndex not found mol.bndpartner!" +
              errorMessage);
        }
      }
    }
    for (auto &partner : mol.bndpartner) {
      if (partner >= moleculeList.size() && partner != -1)
        error("36: mol.bndpartner element >= moleculeList.size(), " +
              to_string(partner) + ">=" + to_string(moleculeList.size()));
      bool found = false;
      for (auto &it : vec)
        if (it == partner) found = true;
      if (!found) {
        error("30: mol[id=" + to_string(mol.id) +
              "].bndpartner not found in interfaceList "
              "interaction.partnerIndex!" +
              errorMessage);
      }
    }
  }
}

void debug_firstEmptyIndex(MpiContext &mpiContext, string errorMessage) {
  if (!DEBUG) return;
  auto &moleculeList = *mpiContext.moleculeList;
  for (auto &firstEmptyIndex : Molecule::emptyMolList) {
    if (firstEmptyIndex >= moleculeList.size())
      error(errorMessage + ": firstEmptyIndex >= moleculeList.size() ");
  }
}

void test_serialization(MpiContext &mpiContext, vector<Molecule> &moleculeList,
                        SimulVolume &simulVolume, Membrane &membraneObject,
                        vector<MolTemplate> &molTemplateList,
                        Parameters &params, vector<ForwardRxn> &forwardRxns,
                        vector<BackRxn> &backRxns,
                        vector<CreateDestructRxn> &createDestructRxns,
                        copyCounters &counterArrays,
                        vector<Complex> &complexList,
                        unsigned char *arrayRank) {
  if (TEST_SERIALIZATION_AND_DESERIALIZATION && (!mpiContext.rank)) {
    test_abstract_vector_serialization<Molecule>(moleculeList, arrayRank, true);
    test_abstract_vector_serialization<Complex>(complexList, arrayRank, true);
    test_object_serialization<SimulVolume>(simulVolume, arrayRank, true);
    test_object_serialization<Membrane>(membraneObject, arrayRank, true);
    test_abstract_vector_serialization<MolTemplate>(molTemplateList, arrayRank,
                                                    true);
    test_object_serialization<Parameters>(params, arrayRank, true);
    test_abstract_vector_serialization<ForwardRxn>(forwardRxns, arrayRank,
                                                   true);
    test_abstract_vector_serialization<BackRxn>(backRxns, arrayRank, true);
    test_abstract_vector_serialization<CreateDestructRxn>(createDestructRxns,
                                                          arrayRank, true);
    test_object_serialization<copyCounters>(counterArrays, arrayRank, true);
    // sleep(5); // so that one could have a glance on results
  }
}

// void LOC(float no, char *text, int proc, int simItr=-1){
void LOC(float no, char *text, int proc, int simItr) {
  using namespace std;
  /*
        This routine was created to define sections of the code with a location
     number and string. When an error is encountered, turn on printing to
     determine the location, iteration and rank for the error.

        Syntax: LOC(const location_number, text_about_location, rank [,simItr]);
        Syntax: LOC(const location_number, text_about_location, PROC_ALL
     [,simItr]); Syntax: LOC(const location_number, text_about_location, PROC_0
     [,simItr]);

        PROC_ALL macro just substitutes mpiContext.rank for rank.
        PROC_0   macro for printing only for rank == 0

        Best to use macros: Allows to turn of all output by changing
        #define RANK rank   to #define RANK -2 (at header space  of main)

         Prints location number/text, rank and simItr,
             e.g LOC(2.1, "Deserialization", mpiContext.rank,simItr);

             -> RNK: 0000  ITR: [10000] LOC: 2.1 @ DESERIALIZTION
             -> ...
             -> RNK: 0004  ITR: [10000] LOC: 2.1 @ DESERIALIZTION



         Prints location number/text, rank and NO simItr  (omit simItr argument)
             e.g LOC(2.1, "Deserialization", PROC_ALL);

             -> RNK: 0000 LOC: 2.1 @ DESERIALIZTION
                ...
             -> RNK: 0003 LOC: 2.1 @ DESERIALIZTION


         Prints location number/text,
                for ONLY RANK0 and     (use PROC_0 macro)
                NO simItr
             e.g LOC(2.1, "Deserialization", PROC_0);

             ->            LOC: 2.1 @ DESERIALIZTION
  */

  stringstream ss2;
  if (proc >= 0)
    ss2 << "RNK: " << setw(4) << setfill('0') << proc << " ";
  else
    ss2 << "          ";      // rank -1 by RANK_0 , means rank is 0,
  string rnk_no = ss2.str();  // but print without rank number

  if (proc == -2)
    return;  // rank = -2 by RANK_0, Don't print non-zero ranks
             // or by RANK macro turning of all printing (rank argument
             // given fake rank of -2 by PROC_0 or RANK macro).

  // Format Iteration Number for loop & nonloop (-1)
  stringstream ss1;  // If LOC() has no simItr arg, default value is -1
  if (simItr >= 0)
    ss1 << " ITR:"
        << " [" << setw(6) << simItr << "] ";
  else
    ss1 << " ";               // No simItr arg, hence, not in sim loop,
  string itr_no = ss1.str();  // print without simItr number

  std::cout << " ->" << rnk_no << itr_no << "LOC:" << fixed << setw(4)
            << std::setprecision(1) << setfill(' ') << no << text << std::endl;
}
