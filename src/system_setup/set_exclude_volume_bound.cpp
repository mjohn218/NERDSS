#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>

#include "boundary_conditions/reflect_functions.hpp"
#include "debug/debug.hpp"
#include "error/error.hpp"
#include "io/io.hpp"
#include "macro.hpp"
#include "math/constants.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "mpi.h"
#include "mpi/mpi_function.hpp"
#include "parser/parser_functions.hpp"
#include "reactions/association/association.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "split.cpp"
#include "system_setup/system_setup.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

using namespace std;

void set_exclude_volume_bound(std::vector<ForwardRxn>& forwardRxns,
                              std::vector<MolTemplate>& molTemplateList) {
  // Set the excludeVolumeBoundList for each molTemplate according to the
  // reaction list
  for (auto& oneReaction : forwardRxns) {
    if (oneReaction.excludeVolumeBound == true) {
      // this reaction declare that excludeVolumeBound, now add the partner's
      // index to the excludeVolumeBoundList
      int react1MolIndex = oneReaction.reactantListNew[0].molTypeIndex;
      int react1InfIndex = oneReaction.reactantListNew[0].relIfaceIndex;
      int react1InfAbsIndex = oneReaction.reactantListNew[0].absIfaceIndex;
      int react2MolIndex = oneReaction.reactantListNew[1].molTypeIndex;
      int react2InfIndex = oneReaction.reactantListNew[1].relIfaceIndex;
      int react2InfAbsIndex = oneReaction.reactantListNew[1].absIfaceIndex;
      molTemplateList[react1MolIndex]
          .interfaceList[react1InfIndex]
          .excludeVolumeBoundList.push_back(react2MolIndex);
      molTemplateList[react2MolIndex]
          .interfaceList[react2InfIndex]
          .excludeVolumeBoundList.push_back(react1MolIndex);
      molTemplateList[react1MolIndex]
          .interfaceList[react1InfIndex]
          .excludeVolumeBoundIfaceList.push_back(react2InfIndex);
      molTemplateList[react2MolIndex]
          .interfaceList[react2InfIndex]
          .excludeVolumeBoundIfaceList.push_back(react1InfIndex);
      molTemplateList[react1MolIndex]
          .interfaceList[react1InfIndex]
          .excludeRadiusList.push_back(oneReaction.bindRadius);
      molTemplateList[react2MolIndex]
          .interfaceList[react2InfIndex]
          .excludeRadiusList.push_back(oneReaction.bindRadius);
      molTemplateList[react1MolIndex]
          .interfaceList[react1InfIndex]
          .excludeVolumeBoundReactList.push_back(oneReaction.relRxnIndex);
      molTemplateList[react2MolIndex]
          .interfaceList[react2InfIndex]
          .excludeVolumeBoundReactList.push_back(oneReaction.relRxnIndex);
      molTemplateList[react1MolIndex].excludeVolumeBound = true;
      molTemplateList[react2MolIndex].excludeVolumeBound = true;
    }
  }
}