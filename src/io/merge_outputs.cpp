#include <dirent.h>  // For directory iteration in POSIX systems
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

struct ComplexKey {
  std::string components;
  bool operator==(const ComplexKey& other) const {
    return components == other.components;
  }
};

namespace std {
template <>
struct hash<ComplexKey> {
  std::size_t operator()(const ComplexKey& k) const {
    return std::hash<std::string>()(k.components);
  }
};
}  // namespace std

void processFile(const std::string& filename,
                 std::map<double, std::unordered_map<ComplexKey, int>>& data) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Failed to open " << filename << std::endl;
    return;
  }

  std::string line;
  double currentTime = -1;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    if (line.find("Time (s):") != std::string::npos) {
      // New time block
      std::string dummy;
      iss >> dummy >> dummy >> currentTime;
    } else {
      // Complex data
      int count;
      iss >> count;
      std::string components, component;
      while (iss >> component) {
        components += component + " ";
      }
      data[currentTime][{components}] += count;
    }
  }
  file.close();
}

void merge_outputs(int totalrank, int molTemplateNum) {
  // MERGE THE PDB FILES IN THE PDB FOLDER AND STORE THE MERGED FILES IN
  // MERGEPDB FOLDER
  int n = totalrank;
  int skiplines = molTemplateNum + 1;
  std::string pdb_folder = "PDB";
  std::string merge_folder = "mergePDB";

  // COLLECT ALL .PDB FILENAMES
  std::vector<std::string> pdb_filenames;

  DIR* dir;
  struct dirent* ent;
  if ((dir = opendir(pdb_folder.c_str())) != NULL) {
    while ((ent = readdir(dir)) != NULL) {
      std::string filename = ent->d_name;
      if (filename.find(".pdb") != std::string::npos) {
        pdb_filenames.push_back(filename);
      }
    }
    closedir(dir);
  } else {
    std::cerr << "Could not open directory: " << pdb_folder << std::endl;
  }

  // EXTRACT UNIQUE TIMEFRAMES
  std::set<int> timeframes;
  for (const auto& filename : pdb_filenames) {
    size_t underscorePos = filename.find("_");
    size_t dotPos = filename.find(".");
    if (underscorePos != std::string::npos && dotPos != std::string::npos) {
      int timeframe = std::stoi(
          filename.substr(underscorePos + 1, dotPos - underscorePos - 1));
      timeframes.insert(timeframe);
    }
  }

  // MERGE
  for (int timeframe : timeframes) {
    std::string merged_filename = std::to_string(timeframe) + ".pdb";
    std::string merged_path = merge_folder + "/" + merged_filename;

    std::ofstream merged_file(merged_path);
    if (!merged_file.is_open()) {
      std::cerr << "Failed to open " << merged_path << std::endl;
      continue;
    }

    for (int rank = 0; rank < n; ++rank) {
      std::string input_filename =
          std::to_string(rank) + "_" + std::to_string(timeframe) + ".pdb";
      std::string input_path = pdb_folder + "/" + input_filename;

      std::ifstream input_file(input_path);
      if (!input_file.is_open()) {
        std::cerr << "Failed to open " << input_path << std::endl;
        continue;
      }

      if (rank == 0) {
        merged_file << input_file.rdbuf();
      } else {
        std::string line;
        for (int i = 0; i < skiplines; ++i) {
          std::getline(input_file, line);  // Skip lines
        }
        merged_file << input_file.rdbuf();
      }
      input_file.close();
    }
    merged_file.close();
  }

  // MERGE THE COPY NUMBER FILES FROM ALL THE PROCESSORS AND STORE IN THE
  // MERGEOUT FOLDER
  std::string input_prefix = "copy_numbers_time_";
  std::string output_folder = "mergeOUT";
  std::string output_filename = output_folder + "/copy_numbers_time.dat";
  // Open output file
  std::ofstream output_file(output_filename);
  if (!output_file.is_open()) {
    std::cerr << "Failed to open output file for writing." << std::endl;
  }

  // PROCESS EACH FILE
  std::string line;
  std::vector<std::ifstream> input_files(n);
  bool headerProcessed = false;
  int numColumns = 0;
  for (int i = 0; i < n; ++i) {
    input_files[i].open(input_prefix + std::to_string(i) + ".dat");
    if (!input_files[i].is_open()) {
      std::cerr << "Failed to open input file " << i << std::endl;
    }

    if (!headerProcessed) {
      if (std::getline(input_files[i], line)) {
        output_file << line << std::endl;
        numColumns = std::count(line.begin(), line.end(), ',');
        headerProcessed = true;
      }
    } else {
      std::getline(input_files[i], line);
    }
  }

  while (true) {
    std::map<double, std::vector<double>> summed_values;
    bool end_of_files = false;
    for (int i = 0; i < n; ++i) {
      if (std::getline(input_files[i], line)) {
        std::istringstream iss(line);
        std::string value;
        double time;
        getline(iss, value, ',');
        time = std::stod(value);
        if (summed_values.find(time) == summed_values.end()) {
          summed_values[time] = std::vector<double>(numColumns, 0.0);
        }
        for (int j = 0; j < numColumns; ++j) {
          if (getline(iss, value, ',')) {
            summed_values[time][j] += std::stod(value);
          }
        }
      } else {
        end_of_files = true;
        break;
      }
    }
    if (end_of_files) break;
    for (const auto& entry : summed_values) {
      output_file << entry.first;
      for (const auto& sum : entry.second) {
        output_file << "," << sum;
      }
      output_file << std::endl;
    }
  }
  for (auto& file : input_files) {
    file.close();
  }
  output_file.close();

  // MERGE HISTOGRAM FILES FROM ALL PROCESSORS
  input_prefix = "histogram_complexes_time_";
  output_filename = output_folder + "/histogram_complexes_time.dat";

  std::map<double, std::unordered_map<ComplexKey, int>> mergedData;

  for (int i = 0; i < n; ++i) {
    processFile(input_prefix + std::to_string(i) + ".dat", mergedData);
  }

  std::ofstream output_file_hist(output_filename);
  if (!output_file_hist.is_open()) {
    std::cerr << "Failed to open output file for writing." << std::endl;
  }

  for (const auto& timeData : mergedData) {
    output_file_hist << "Time (s): " << timeData.first << std::endl;
    for (const auto& complexCount : timeData.second) {
      output_file_hist << complexCount.second << "\t"
                       << complexCount.first.components << std::endl;
    }
  }
  output_file_hist.close();
}