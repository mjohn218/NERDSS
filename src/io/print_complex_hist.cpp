#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "classes/class_Molecule_Complex.hpp"
#include "debug/debug.hpp"
#include "io/io.hpp"
#include "mpi/mpi_function.hpp"
#include "tracing.hpp"

double print_complex_hist(std::vector<Complex>& complexList,
                          std::ofstream& outfile, int it, Parameters params,
                          std::vector<MolTemplate>& molTemplateList,
                          std::vector<Molecule>& moleculeList,
                          int nImplicitLipids, MpiContext& mpiContext,
                          SimulVolume& simulVolume) {
  int n_mol_types = molTemplateList.size();
  int n_complexes = complexList.size();
  std::unordered_map<std::string, int> complex_count_dict;
  std::vector<std::string> mol_type_names;

  for (auto i = 0; i < n_mol_types; i++) {
    mol_type_names.push_back(molTemplateList[i].molName);
  }

  for (auto i = 0; i < n_complexes; i++) {
    auto& com = complexList[i];
    if (com.isEmpty) {
      continue;
    }

    std::vector<int> numEachMol(n_mol_types, 0);

    int n_members = com.memberList.size();

    int largest_id = -1;
    int largest_id_index = -1;

    for (auto j = 0; j < n_members; j++) {
      auto& mol = moleculeList[com.memberList[j]];
      numEachMol[mol.molTypeIndex] += 1;

      if (mol.id > largest_id && mol.isImplicitLipid == false) {
        largest_id = mol.id;
        largest_id_index = mol.index;
      }
    }

    if (moleculeList[largest_id_index].isLeftGhost ||
        moleculeList[largest_id_index].isRightGhost) {
      continue;
    }

    if (molTemplateList[0].isImplicitLipid) {
      numEachMol[0] = com.linksToSurface;
    }

    std::string components_string = "";

    for (auto j = 0; j < n_mol_types; j++) {
      if (numEachMol[j] != 0) {
        components_string +=
            mol_type_names[j] + ": " + std::to_string(numEachMol[j]) + ". ";
      }
    }

    // update the count of this complex
    if(components_string!=""){
      if (complex_count_dict.find(components_string) ==
          complex_count_dict.end()) {
        complex_count_dict[components_string] = 1;
      } else {
        complex_count_dict[components_string] += 1;
      }
    }
  }

  // write the complex counts to the output file
  outfile << "Time (s): "
          << (it - params.itrRestartFrom) * params.timeStep * 1E-6 +
                 params.timeRestartFrom
          << "\n";

  if(molTemplateList[0].isImplicitLipid){
    outfile << nImplicitLipids << '\t' << mol_type_names[0] << ": 1." << "\n";
  }

  for (auto& pair : complex_count_dict) {
    auto& components_string = pair.first;
    auto& count = pair.second;
    outfile << count << '\t' << components_string << "\n";
  }

  outfile << std::flush;

  return 0.0;
}