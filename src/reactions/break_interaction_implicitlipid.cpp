#include "reactions/implicitlipid/implicitlipid_reactions.hpp"

void break_interaction_implicitlipid(size_t relIface1, size_t relIface2, Molecule& reactMol1, Molecule& reactMol2,
                       const BackRxn& currRxn, std::vector<int>& emptyComList, std::vector<Molecule>& moleculeList,
                       std::vector<Complex>& complexList, const std::vector<MolTemplate>& molTemplateList)
{
    /*Find the new absIface for mol1*/
    int absIface1 { -1 };
    int absIface2 { -1 };
    if (currRxn.isSymmetric) {
        absIface1 = currRxn.productListNew[0].absIfaceIndex;
        absIface2 = currRxn.productListNew[1].absIfaceIndex;
    } else {
        if (reactMol1.molTypeIndex == currRxn.productListNew[0].molTypeIndex
            && relIface1 == currRxn.productListNew[0].relIfaceIndex) {
            absIface1 = currRxn.productListNew[0].absIfaceIndex;
            absIface2 = currRxn.productListNew[1].absIfaceIndex;
        } else if (reactMol1.molTypeIndex == currRxn.productListNew[1].molTypeIndex
            && relIface1 == currRxn.productListNew[1].relIfaceIndex) {
            absIface2 = currRxn.productListNew[0].absIfaceIndex;
            absIface1 = currRxn.productListNew[1].absIfaceIndex;
        }
    }

    reactMol1.interfaceList[relIface1].interaction.clear();
    reactMol1.interfaceList[relIface1].isBound = false;
    reactMol1.interfaceList[relIface1].index = absIface1;

    //Add these protein into the bimolecular association list
    reactMol1.freelist.push_back(relIface1);
    reactMol1.bndlist.erase(std::find_if(reactMol1.bndlist.begin(), reactMol1.bndlist.end(), [&](const size_t& iface) { return iface == relIface1; }));  
    reactMol1.bndpartner.erase(std::find_if(reactMol1.bndpartner.begin(), reactMol1.bndpartner.end(),[&](const size_t& mol) { return mol == reactMol2.index; }));
}
