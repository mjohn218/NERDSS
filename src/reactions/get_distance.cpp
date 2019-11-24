#include "reactions/bimolecular/bimolecular_reactions.hpp"

bool get_distance(int pro1, int pro2, int iface1, int iface2, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    double& sep, double& R1, double Rmax, std::vector<Complex>& complexList, const ForwardRxn& currRxn,
    std::vector<Molecule>& moleculeList)
{
    double dx = moleculeList[pro1].interfaceList[iface1].coord.x - moleculeList[pro2].interfaceList[iface2].coord.x;
    double dy = moleculeList[pro1].interfaceList[iface1].coord.y - moleculeList[pro2].interfaceList[iface2].coord.y;

    double dz { (complexList[moleculeList[pro1].myComIndex].D.z == 0
                    && complexList[moleculeList[pro2].myComIndex].D.z == 0)
            ? 0
            : moleculeList[pro1].interfaceList[iface1].coord.z - moleculeList[pro2].interfaceList[iface2].coord.z };

    R1 = sqrt((dx * dx) + (dy * dy) + (dz * dz));
    sep = R1 - currRxn.bindRadius;

    /*Rmax should be the binding radius plus ~max diffusion distance, using 3*sqrt(6*Dtot*deltat)*/
    if (R1 < Rmax) {
        /*in this case we evaluate the probability of this reaction*/
        moleculeList[pro1].crossbase.push_back(pro2);
        moleculeList[pro2].crossbase.push_back(pro1);
        moleculeList[pro1].mycrossint.push_back(iface1);
        moleculeList[pro2].mycrossint.push_back(iface2);
        moleculeList[pro1].crossrxn.push_back(std::array<int, 3> { rxnIndex, rateIndex, isStateChangeBackRxn });
        moleculeList[pro2].crossrxn.push_back(std::array<int, 3> { rxnIndex, rateIndex, isStateChangeBackRxn });
        ++complexList[moleculeList[pro1].myComIndex].ncross;
        ++complexList[moleculeList[pro2].myComIndex].ncross;
        return true;
    }
    return false;
}
