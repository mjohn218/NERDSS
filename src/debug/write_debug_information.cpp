#include <cstring>
#include <iostream>
#include <unordered_set>

#include "debug/debug.hpp"
#include "error/error.hpp"

using namespace std;

void write_debug_information(MpiContext& mpiContext, long long int simItr,
                             ofstream& debugFile,
                             vector<Molecule>& moleculeList,
                             vector<Complex>& complexList,
                             vector<MolTemplate>& molTemplateList,
                             copyCounters& counterArrays, double simTime,
                             string s) {}