#include "boundary_conditions/reflect_functions.hpp"
#include "classes/class_Cluster.hpp"
#include "math/rand_gsl.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

//Definition of constructor
ClusterPair::ClusterPair()
{
    p1 = -1;
}

ClusterPair::ClusterPair(int setp1, int setp2)
{
    p1 = setp1;
    p2 = setp2;
}

void cluster_one_complex(int k1, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns, std::vector<ClusterPair>& pairList, std::vector<int>& finished)
{
    // Figure out all pairs for one complex, k1
    // This is for checking overlap, so do not add any pairs to the list that within th same complex
    int csize = complexList[k1].memberList.size();
    int c, i, t;
    t = 0;
    int p1, p2, i1, i2;
    int k2, rxn;
    int flag;
    for (c = 0; c < csize; c++) {
        p1 = complexList[k1].memberList[c];
        int p1basesize = moleculeList[p1].crossbase.size();
        for (i = 0; i < p1basesize; i++) {
            t++;
            p2 = moleculeList[p1].crossbase[i];
            ClusterPair newPair(p1, p2);
            k2 = moleculeList[p2].myComIndex;
            if (k2 != k1) {
                // no pairs within the same complex
                i1 = moleculeList[p1].mycrossint[i];
                std::array<int, 3> rxnItr = moleculeList[p1].crossrxn[i];
                rxn = rxnItr[0];
                newPair.i1 = i1;
                newPair.bindrad = forwardRxns[rxn].bindRadius;
                // get the partner interface
                newPair.i2 = (forwardRxns[rxn].reactantListNew[0].relIfaceIndex == i1)
                    ? forwardRxns[rxn].reactantListNew[1].relIfaceIndex
                    : forwardRxns[rxn].reactantListNew[0].relIfaceIndex;

                flag = 0;
                for (int j = 0; j < pairList.size(); j++) {
                    if (pairList[j].p1 == p1 && pairList[j].p2 == p2) {
                        if (pairList[j].i1 == i1 && pairList[j].i2 == i2) {
                            flag = 1;
                        }
                    } else if (pairList[j].p1 == p2 && pairList[j].p2 == p1) { //check for p2 (i2), p1 (i1)
                        if (pairList[j].i1 == i2 && pairList[j].i2 == i1) {
                            flag = 1;
                        }
                    }
                }
                if (flag == 0) {
                    newPair.priority = 1;
                    if (complexList[k1].D.z < 1E-15 || complexList[k2].D.z < 1E-15)
                        newPair.priority = 0; //low prioirity to solve versus lipid, since no structure

                    if (complexList[k2].D.z < 1E-15 && complexList[k1].D.z < 1E-15) {
                        newPair.memtest = 1;
                    } else
                        newPair.memtest = 0;
                    newPair.k1 = k1;
                    newPair.k2 = k2;
                    pairList.push_back(newPair);
                }
            }
        }
    }
    finished.push_back(k1);
}

void define_cluster_pairs(int p1, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns, std::vector<ClusterPair>& pairList)
{
    // start off with one protein, and descend through all its' partners and their partners
    int c, i;
    int k1 = moleculeList[p1].myComIndex;
    std::vector<int> finished;
    int k2;
    int flag1, flag2;

    // Add all pairs involving complex k1 and its cross partners
    cluster_one_complex(k1, moleculeList, complexList, forwardRxns, pairList, finished);

    // Below, loop over all current pairs.
    for (i = 0; i < pairList.size(); i++) {
        k1 = pairList[i].k1;
        k2 = pairList[i].k2;
        flag1 = 0;
        flag2 = 0;
        for (int f = 0; f < finished.size(); f++) {
            if (k1 == finished[f])
                flag1 = 1;
            if (k2 == finished[f])
                flag2 = 1;
        }
        if (flag1 == 0)
            cluster_one_complex(k1, moleculeList, complexList, forwardRxns, pairList, finished);
        if (flag2 == 0)
            cluster_one_complex(k2, moleculeList, complexList, forwardRxns, pairList, finished);
    }
}

void resample_traj(int currStop, std::vector<ClusterPair>& pairList, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const Parameters& params, const Membrane& membraneObject, double RS3Dinput)
{
    int i;
    int k1, k2;
    int p1, p2;
    int flag;
    std::vector<int> didMove;
    for (i = 0; i < currStop; i++) {
        p1 = pairList[i].p1;
        p2 = pairList[i].p2;
        k1 = pairList[i].k1;
        k2 = pairList[i].k2;

        if (moleculeList[p1].trajStatus == TrajStatus::none || moleculeList[p1].trajStatus == TrajStatus::canBeResampled) {
            flag = 0;
            for (int d = 0; d < didMove.size(); d++) {
                if (didMove[d] == k1)
                    flag = 1;
            }
            if (flag == 0) {

                if (membraneObject.isSphere == true && complexList[k1].D.z < 1E-15) { // complex on sphere surface
                    Coord targTrans = create_complex_propagation_vectors_on_sphere(params, complexList[k1]);
                    complexList[k1].trajTrans.x = targTrans.x;
                    complexList[k1].trajTrans.y = targTrans.y;
                    complexList[k1].trajTrans.z = targTrans.z;
                    complexList[k1].trajRot.x = sqrt(2.0 * params.timeStep * complexList[k1].Dr.x) * GaussV();
                    complexList[k1].trajRot.y = sqrt(2.0 * params.timeStep * complexList[k1].Dr.y) * GaussV();
                    complexList[k1].trajRot.z = sqrt(2.0 * params.timeStep * complexList[k1].Dr.z) * GaussV();
                } else {
                    complexList[k1].trajTrans.x = sqrt(2.0 * params.timeStep * complexList[k1].D.x) * GaussV();
                    complexList[k1].trajTrans.y = sqrt(2.0 * params.timeStep * complexList[k1].D.y) * GaussV();
                    complexList[k1].trajTrans.z = sqrt(2.0 * params.timeStep * complexList[k1].D.z) * GaussV();
                    complexList[k1].trajRot.x = sqrt(2.0 * params.timeStep * complexList[k1].Dr.x) * GaussV();
                    complexList[k1].trajRot.y = sqrt(2.0 * params.timeStep * complexList[k1].Dr.y) * GaussV();
                    complexList[k1].trajRot.z = sqrt(2.0 * params.timeStep * complexList[k1].Dr.z) * GaussV();
                }

                reflect_traj_complex_rad_rot(params, moleculeList, complexList[k1], membraneObject, RS3Dinput);

                didMove.push_back(k1);
            }
        }
        if (moleculeList[p2].trajStatus == TrajStatus::none || moleculeList[p2].trajStatus == TrajStatus::canBeResampled) {
            flag = 0;
            for (int d = 0; d < didMove.size(); d++) {
                if (didMove[d] == k2)
                    flag = 1;
            }
            if (flag == 0) {

                if (membraneObject.isSphere == true && complexList[k2].D.z < 1E-15) { // complex on sphere surface
                    Coord targTrans = create_complex_propagation_vectors_on_sphere(params, complexList[k2]);
                    complexList[k2].trajTrans.x = targTrans.x;
                    complexList[k2].trajTrans.y = targTrans.y;
                    complexList[k2].trajTrans.z = targTrans.z;
                    complexList[k2].trajRot.x = sqrt(2.0 * params.timeStep * complexList[k2].Dr.x) * GaussV();
                    complexList[k2].trajRot.y = sqrt(2.0 * params.timeStep * complexList[k2].Dr.y) * GaussV();
                    complexList[k2].trajRot.z = sqrt(2.0 * params.timeStep * complexList[k2].Dr.z) * GaussV();
                } else {
                    complexList[k2].trajTrans.x = sqrt(2.0 * params.timeStep * complexList[k2].D.x) * GaussV();
                    complexList[k2].trajTrans.y = sqrt(2.0 * params.timeStep * complexList[k2].D.y) * GaussV();
                    complexList[k2].trajTrans.z = sqrt(2.0 * params.timeStep * complexList[k2].D.z) * GaussV();
                    complexList[k2].trajRot.x = sqrt(2.0 * params.timeStep * complexList[k2].Dr.x) * GaussV();
                    complexList[k2].trajRot.y = sqrt(2.0 * params.timeStep * complexList[k2].Dr.y) * GaussV();
                    complexList[k2].trajRot.z = sqrt(2.0 * params.timeStep * complexList[k2].Dr.z) * GaussV();
                }

                reflect_traj_complex_rad_rot(params, moleculeList, complexList[k2], membraneObject, RS3Dinput);

                didMove.push_back(k2);
            }
        }
    }
}
