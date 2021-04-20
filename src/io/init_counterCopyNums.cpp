#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_copyCounters.hpp"
#include "tracing.hpp"

using namespace std;

/*Initialize the array of copyNumSpecies. Based off the copy numbers looped over all Molecules, thus should work
  for either default start or for restart.
 */
void init_counterCopyNums(copyCounters& counterArrays, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject, int totalSpeciesNum, Parameters& params)
{
    // TRACE();
    int i, j;

    int index;
    int p1, p2;
    // initialize the copyNums array to zero, for all possible interfaces and states (includes product species)
    /*For products containing two components, it will double count them!
     */
    //int totalSpecies = RxnBase::totRxnSpecies; //Interface::State::totalNumOfStates;
    int totalSpecies = totalSpeciesNum;
    counterArrays.copyNumSpecies.reserve(totalSpecies);
    for (i = 0; i < totalSpecies; i++)
        counterArrays.copyNumSpecies.push_back(0);

    for (p1 = 0; p1 < moleculeList.size(); p1++) {
        int numIfaces = moleculeList[p1].interfaceList.size();
        for (j = 0; j < numIfaces; j++) {
            // find out which state each interface on the molecule is in, and increment copyNumSpecies array.
            index = moleculeList[p1].interfaceList[j].index;
            if (moleculeList[p1].isImplicitLipid == false) {
                counterArrays.copyNumSpecies[index]++;
                if (params.fromRestart == false) {
                    if (counterArrays.canDissociate[index] == true) {
                        // add the mol index to the bindPairList
                        // if (complexList[moleculeList[p1].myComIndex].linksToSurface > 1)
                        //     counterArrays.bindPairListIL2D[index].emplace_back(p1);
                        // if (complexList[moleculeList[p1].myComIndex].linksToSurface == 1)
                        //     counterArrays.bindPairListIL3D[index].emplace_back(p1);
                        if (complexList[moleculeList[p1].myComIndex].linksToSurface == 0) {
                            // for the explicit case, only add when the partner does not exsit
                            int partnerMolIndex { moleculeList[p1].interfaceList[j].interaction.partnerIndex };
                            std::vector<int>::iterator itr = counterArrays.bindPairList[index].begin();
                            for (itr = counterArrays.bindPairList[index].begin(); itr != counterArrays.bindPairList[index].end(); ++itr) {
                                if (*itr == partnerMolIndex)
                                    break;
                            }
                            if (itr == std::end(counterArrays.bindPairList[index])) {
                                counterArrays.bindPairList[index].emplace_back(p1);
                            }
                        }
                    }
                }
            } else {
                //For implicit lipid, set copy numbers based on read in template.
                //int molTypeIndex = moleculeList[p1].molTypeIndex;
                for (int tmpStateIndex = 0; tmpStateIndex < membraneObject.nStates; tmpStateIndex++) {
                    counterArrays.copyNumSpecies[index + tmpStateIndex] = membraneObject.numberOfFreeLipidsEachState[tmpStateIndex]; //IL mol must be the first place
                }
            }
        } // end all interfaces
    } // end all current molecules
    cout << "Initialized copyNumSpecies. Printing all values. total species: " << totalSpecies << endl;
    for (i = 0; i < totalSpecies; i++) {
        if (counterArrays.singleDouble[i] == 2 && counterArrays.implicitDouble[i] == false) {
            //product state, contains two proteins, so will be double counted above.
            counterArrays.copyNumSpecies[i] *= 0.5;
        }
        cout << " absolute index: " << i << " Copy nums: " << counterArrays.copyNumSpecies[i] << endl;
    }
}
