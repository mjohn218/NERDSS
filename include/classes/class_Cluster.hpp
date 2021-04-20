#pragma once

// #include "classes/class_Membrane.hpp"
// #include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"

class ClusterPair {
public:
    int p1;
    int p2;
    int i1;
    int i2;
    int k1;
    int k2;
    int priority;
    int memtest;
    double bindrad;
    //int flagOverlap;

    ClusterPair(void); //Constructor
    ClusterPair(int setp1, int setp2); //Constructor
};

void cluster_one_complex(int k1, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns, std::vector<ClusterPair>& pairList, std::vector<int>& finished);
void define_cluster_pairs(int p1, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns, std::vector<ClusterPair>& pairList);
void resample_traj(int currStop, std::vector<ClusterPair>& pairList, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const Parameters& params, const Membrane& membraneObject, double RS3Dinput);