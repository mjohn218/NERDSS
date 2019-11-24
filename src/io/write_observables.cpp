#include "io/io.hpp"
#include <chrono>
#include <ctime>

void write_observables(
    unsigned simItr, std::ofstream& observablesFile, const std::map<std::string, int>& observablesList)
{
    if (observablesList.size() == 1)
        observablesFile << simItr << ',' << observablesList.begin()->second << '\n';
    else {
        observablesFile << simItr;
        for (auto obsItr = observablesList.begin(); obsItr != observablesList.end(); ++obsItr)
            observablesFile << ',' << obsItr->second;
        observablesFile << '\n';
    }
}
