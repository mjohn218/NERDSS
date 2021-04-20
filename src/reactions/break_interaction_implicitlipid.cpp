#include "reactions/implicitlipid/implicitlipid_reactions.hpp"
#include "tracing.hpp"

void break_interaction_implicitlipid(size_t relIface1, size_t relIface2, Molecule& reactMol1, Molecule& reactMol2,
    const BackRxn& currRxn, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, std::vector<MolTemplate>& molTemplateList)
{
    // TRACE();
    /*Find the new absIface for mol1*/
    int absIface1 { -1 };
    int absIface2 { -1 };
    if (currRxn.isSymmetric) {
        absIface1 = currRxn.productListNew[0].absIfaceIndex;
        absIface2 = currRxn.productListNew[1].absIfaceIndex;
    } else {
        /*Here, if both proteins are the same protein, same interface, but distinct states, need to correct for that.*/
        // std::cout << "State of mol1: " << reactMol1.interfaceList[relIface1].stateIden << " of mol2: " << reactMol2.interfaceList[relIface2].stateIden << "\n";
        // std::cout << " Product requires state: " << currRxn.productListNew[0].requiresState << " " << currRxn.productListNew[1].requiresState << std::endl;
        if (reactMol1.molTypeIndex == currRxn.productListNew[0].molTypeIndex
            && relIface1 == currRxn.productListNew[0].relIfaceIndex && currRxn.productListNew[0].requiresState == reactMol1.interfaceList[relIface1].stateIden) {
            //matched protein, interface, and state of product[0] to mol1
            absIface1 = currRxn.productListNew[0].absIfaceIndex;
            absIface2 = currRxn.productListNew[1].absIfaceIndex;

        } else if (reactMol1.molTypeIndex == currRxn.productListNew[1].molTypeIndex
            && relIface1 == currRxn.productListNew[1].relIfaceIndex && currRxn.productListNew[1].requiresState == reactMol1.interfaceList[relIface1].stateIden) {
            //matched protein, interface and state of product[1] to mol1
            absIface2 = currRxn.productListNew[0].absIfaceIndex;
            absIface1 = currRxn.productListNew[1].absIfaceIndex;
        } else {
            std::cout << " IN BREAK INTERACTION, DID NOT MATCH protein, interface, and state of a product to Mol1 " << reactMol1.index << std::endl;
            exit(1);
        }
    }

    reactMol1.interfaceList[relIface1].interaction.clear();
    reactMol1.interfaceList[relIface1].isBound = false;
    reactMol1.interfaceList[relIface1].index = absIface1;

    //Add these protein into the bimolecular association list
    reactMol1.freelist.push_back(relIface1);
    reactMol1.bndlist.erase(std::find_if(reactMol1.bndlist.begin(), reactMol1.bndlist.end(), [&](const size_t& iface) { return iface == relIface1; }));
    reactMol1.bndpartner.erase(std::find_if(reactMol1.bndpartner.begin(), reactMol1.bndpartner.end(), [&](const size_t& mol) { return mol == reactMol2.index; }));

    //------------------------START UPDATE MONOMERLIST-------------------------
    // update oneTemp.monomerList when oneTemp.canDestroy is true and mol is monomer, add to monomerList if new monomer produced
    // reactMol1
    {
        Molecule& oneMol { reactMol1 };
        MolTemplate& oneTemp { molTemplateList[oneMol.molTypeIndex] };
        bool isMonomer { oneMol.bndpartner.empty() };
        bool canDestroy { oneTemp.canDestroy };
        if (isMonomer && canDestroy) {
            //add to monomerList
            oneTemp.monomerList.emplace_back(oneMol.index);
        }
    }
    //------------------------END UPDATE MONOMERLIST---------------------------
}
