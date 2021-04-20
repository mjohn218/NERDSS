#include "io/io.hpp"
#include "tracing.hpp"
#include <chrono>
#include <ctime>

void write_observables(
    double simTime, std::ofstream& observablesFile, const std::map<std::string, int>& observablesList)
{
    // TRACE();
    if (observablesList.size() == 1)
        observablesFile << simTime << ',' << observablesList.begin()->second << '\n';
    else {
        observablesFile << simTime;
        for (auto obsItr = observablesList.begin(); obsItr != observablesList.end(); ++obsItr)
            observablesFile << ',' << obsItr->second;
        observablesFile << std::endl;
    }
}
