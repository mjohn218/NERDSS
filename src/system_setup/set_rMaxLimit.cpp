#include "system_setup/system_setup.hpp"

void set_rMaxLimit( Parameters& params, const std::vector<MolTemplate>& molTemplateList, const std::vector<ForwardRxn>& forwardRxns)
{
    /*For each reaction, need distance from the interface to the COM for both partners, plus the
     * bindrad+sqrt(6*Dtot*deltat)
     */
    double rMaxTot { 0 };
    for (const auto& oneRxn : forwardRxns) {
        if (oneRxn.rxnType == ReactionType::bimolecular) {
            const RxnIface& rxnIface1 = oneRxn.reactantListNew.at(0);
            const RxnIface& rxnIface2 = oneRxn.reactantListNew.at(1);

            // MolTemplate Interfaces which are involved in the reaction
            const Interface& iface1 = molTemplateList.at(rxnIface1.molTypeIndex)
                                          .interfaceList.at(MolTemplate::absToRelIface.at(rxnIface1.absIfaceIndex));
            const Interface& iface2 = molTemplateList.at(rxnIface2.molTypeIndex)
                                          .interfaceList.at(MolTemplate::absToRelIface.at(rxnIface2.absIfaceIndex));

            // MolTemplates of the interfaces involved in the reaction
            const MolTemplate& pro1Temp = molTemplateList.at(rxnIface1.molTypeIndex);
            const MolTemplate& pro2Temp = molTemplateList.at(rxnIface2.molTypeIndex);

            //            double Dtot = pro1.D.x + pro2.D.x;
            double scal { 1.0 / 3.0 };
            double Dtot { (scal * (pro1Temp.D.x + pro2Temp.D.x)) + (scal * (pro1Temp.D.y + pro2Temp.D.y))
                + (scal * (pro1Temp.D.z + pro2Temp.D.z)) };
            // TODO: add rotational diffusion contribution

            /*Now calculate distance from the interface to the protein COM.*/
            double pro1R1 { 0 };
            double pro2R1 { 0 };
            // pro1
            if (pro1Temp.Dr.z == 0) {
                double R2 = (iface1.iCoord.x * iface1.iCoord.x) + (iface1.iCoord.y * iface1.iCoord.y);
                pro1R1 = sqrt(R2);
                double einsStks = cos(sqrt(2.0 * pro1Temp.Dr.z * params.timeStep));
                double Dr = 2.0 * R2 * (1.0 - einsStks);
                Dtot += Dr / (4.0 * params.timeStep);
            } else {
                double R2 = (iface1.iCoord.x * iface1.iCoord.x) + (iface1.iCoord.y * iface1.iCoord.y)
                    + (iface1.iCoord.z * iface1.iCoord.z);
                pro1R1 = sqrt(R2);
                double einsStks = cos(sqrt(4.0 * pro1Temp.Dr.z * params.timeStep));
                double Dr = 2.0 * R2 * (1.0 - einsStks);
                Dtot += Dr / (6.0 * params.timeStep);
            }

            // pro2
            if (pro2Temp.Dr.z == 0) {
                double R2 = (iface2.iCoord.x * iface2.iCoord.x) + (iface2.iCoord.y * iface2.iCoord.y);
                pro2R1 = sqrt(R2);
                double einsStks = cos(sqrt(2.0 * pro2Temp.Dr.z * params.timeStep));
                double Dr = 2.0 * R2 * (1.0 - einsStks);
                Dtot += Dr / (4.0 * params.timeStep);
            } else {
                double R2 = (iface2.iCoord.x * iface2.iCoord.x) + (iface2.iCoord.y * iface2.iCoord.y)
                    + (iface2.iCoord.z * iface2.iCoord.z);
                pro2R1 = sqrt(R2);
                double einsStks = cos(sqrt(4.0 * pro2Temp.Dr.z * params.timeStep));
                double Dr = 2.0 * R2 * (1.0 - einsStks);
                Dtot += Dr / (6.0 * params.timeStep);
            }

            double RmaxDiff = 3.0 * sqrt(6.0 * Dtot * params.timeStep) + oneRxn.bindRadius;
            rMaxTot = RmaxDiff + pro1R1 + pro2R1;
            if (rMaxTot > params.rMaxLimit) {
                params.rMaxLimit = rMaxTot;
                params.rMaxRadius = pro1R1 + pro2R1;
            }
        }
    }
    std::cout << "Rmaxlimit: " << params.rMaxLimit << std::endl;
}
