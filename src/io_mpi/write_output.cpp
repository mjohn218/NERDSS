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
#ifdef mpi_
#include "mpi.h"
#endif
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

using MDTimer = std::chrono::system_clock;

void write_output(long long int simItr, Parameters& params,
                  std::string trajFileName, std::vector<Molecule>& moleculeList,
                  std::vector<MolTemplate>& molTemplateList,
                  Membrane& membraneObject, MpiContext& mpiContext,
                  std::string transitionFileName, copyCounters& counterArrays,
                  std::ofstream& speciesFile1, std::ofstream& debugFile,
                  std::string debugFileName, std::string restartFileName,
                  std::vector<Complex>& complexList, SimulVolume& simulVolume,
                  std::vector<ForwardRxn>& forwardRxns,
                  std::vector<BackRxn>& backRxns,
                  std::vector<CreateDestructRxn>& createDestructRxns,
                  std::vector<std::chrono::duration<double>>& durationList,
                  MDTimer::time_point& startStep,
                  MDTimer::time_point& totalTimeStart,
                  std::map<std::string, int>& observablesList,
                  std::ofstream& assemblyfile) {
  // if (simItr % params.trajWrite == 0) {
  //   std::ofstream trajFile{trajFileName, std::ios::app};  // for append
  //   write_traj(simItr, trajFile, params, moleculeList, molTemplateList,
  //              membraneObject);
  //   trajFile.close();
  // }

  if (params.pdbWrite != -1) {
    if (simItr % params.pdbWrite == 0) {
      // std::cout << "Writing PDB file for current frame...\n";
      write_pdb(simItr, simItr, params, moleculeList, molTemplateList,
                membraneObject, mpiContext.rank, mpiContext, simulVolume);
    }
  }

  // if (simItr % params.restartWrite) { // turn off the restart now
  if (0) {
    std::ofstream restartFile{restartFileName,
                              std::ios::out};  // to show different from append
    write_rng_state(mpiContext.rank);          // write the current RNG state
    std::vector<TransmissionRxn> transmissionRxns{};
    write_restart(simItr, restartFile, params, simulVolume, moleculeList,
                  complexList, molTemplateList, forwardRxns, backRxns,
                  createDestructRxns,transmissionRxns, observablesList, membraneObject,
                  counterArrays);

    restartFile.close();
  }

  if (simItr % params.timeWrite == 0) {
    write_all_species((simItr - params.itrRestartFrom) * params.timeStep *
                              Constants::usToSeconds +
                          params.timeRestartFrom,
                      moleculeList, speciesFile1, counterArrays, membraneObject,
                      mpiContext, simulVolume);

    int number_of_lipids = 0;  // sum of all states of IL
    for (int i = 0; i < membraneObject.numberOfFreeLipidsEachState.size();
         i++) {
      number_of_lipids += membraneObject.numberOfFreeLipidsEachState[i];
    }
    int meanComplexSize = print_complex_hist(
        complexList, assemblyfile, simItr, params, molTemplateList,
        moleculeList, number_of_lipids, mpiContext, simulVolume);

    char charTime[24];

    using duration = std::chrono::duration<double>;
    durationList.erase(durationList.begin());
    durationList.emplace_back(MDTimer::now() - startStep);
    if (simItr % params.timeWrite == 0) {
      double timeSimulated{(simItr - params.itrRestartFrom) * params.timeStep *
                               Constants::usToSeconds +
                           params.timeRestartFrom};
      std::cout << linebreak;
      std::cout << "End iteration: " << simItr << ", simulation time: ";
      std::cout << std::scientific << timeSimulated << " seconds.\n";
      auto endTime = MDTimer::now();
      auto endTimeFormat = MDTimer::to_time_t(endTime);
      std::cout << "System time: ";
      if (0 < strftime(charTime, sizeof(charTime), "%F %T",
                       std::localtime(&endTimeFormat)))
        std::cout << charTime << '\n';
      // elapsed time in minutes
      auto elapsedTime = std::chrono::duration_cast<std::chrono::minutes>(
                             MDTimer::now() - totalTimeStart)
                             .count();
      std::cout << "Elapsed time: " << elapsedTime << " minutes\n";

      // std::cout << "Number of molecules: " << Molecule::numberOfMolecules
      //           << '\n';
      // std::cout << "Number of complexes: " << Complex::numberOfComplexes
      //           << '\n';
      // std::cout << "Total reaction matches: " << totMatches << '\n';
      // if (params.debugParams.printSystemInfo) {
      //   std::cout << "Printing full system information...\n";
      //   std::ofstream systemInfoFile{"system_information.dat",
      //   std::ios::app}; print_system_information(simItr, systemInfoFile,
      //   moleculeList,
      //                            complexList, molTemplateList);
      //   systemInfoFile.close();
      // }

      int numSavedDurations{1000};

      // Estimate time remaining
      duration avgTimeStepDuration =
          std::accumulate(durationList.begin(), durationList.end(),
                          duration{0}) /
          numSavedDurations;
      duration timeLeft = (params.nItr - simItr) * avgTimeStepDuration;
      auto timeLeftMinutes = 1.0 * elapsedTime /
                             (simItr - params.itrRestartFrom) *
                             (params.nItr - simItr);
      // round to integer
      timeLeftMinutes = std::round(timeLeftMinutes);
      std::cout
          // << "Avg timestep duration: " << avgTimeStepDuration.count()
          << "Iterations remaining: " << params.nItr - simItr
          << ", Time left: "
          // <<
          // std::chrono::duration_cast<std::chrono::minutes>(timeLeft).count()
          // show in float, not scientific notation
          << std::fixed << std::setprecision(0) << timeLeftMinutes
          << " minutes\n";
      auto estTimeLeft =
          MDTimer::now() + timeLeftMinutes * std::chrono::minutes{1};
      auto estTimeEnd = std::chrono::system_clock::to_time_t(
          std::chrono::system_clock::from_time_t(
              std::time_t(std::chrono::duration_cast<std::chrono::seconds>(
                              estTimeLeft.time_since_epoch())
                              .count())));
      std::cout << "Estimated end time: ";
      //<< std::put_time(std::localtime(&estTimeEnd), "%F %T") << '\n';
      if (0 < strftime(charTime, sizeof(charTime), "%F %T",
                       std::localtime(&estTimeEnd)))
        std::cout << charTime << '\n';
      std::cout << llinebreak;
    }
  }
}