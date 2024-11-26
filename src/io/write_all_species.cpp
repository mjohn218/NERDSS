/*! \file write_all_species.cpp
 *
 * \brief
 *
 * ### Created on 2019-06-05 by Matthew Varga
 */
#include "io/io.hpp"
#include "tracing.hpp"
#include "debug/debug.hpp"
#include "mpi/mpi_function.hpp"

void write_all_species(double simTime, std::ofstream& speciesFile, const copyCounters& counterArray)
{
    // TRACE();
    speciesFile << simTime;
    for (auto elem : counterArray.copyNumSpecies)
        speciesFile << ',' << elem;
    speciesFile << std::endl;
}
