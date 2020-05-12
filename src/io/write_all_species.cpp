/*! \file write_all_species.cpp
 *
 * \brief
 *
 * ### Created on 2019-06-05 by Matthew Varga
 */
#include "io/io.hpp"
#include "tracing.hpp"

void write_all_species(unsigned simItr, std::ofstream& speciesFile, const copyCounters& counterArray)
{
    // TRACE();
    speciesFile << simItr;
    for (auto elem : counterArray.copyNumSpecies)
        speciesFile << ',' << elem;
    speciesFile << '\n';
}