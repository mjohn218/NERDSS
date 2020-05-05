#include "boundary_conditions/reflect_functions.hpp"
#include "io/io.hpp"
#include "reactions/association/association.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"

#include <iomanip>

void perform_bimolecular_state_change(int stateChangeIface, int facilitatorIface, std::array<int, 3>& rxnItr,
    Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom, Complex& facilitatorCom,
    copyCounters& counterArrays, const Parameters& params, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, Membrane& membraneObject)
{
    if (membraneObject.isSphere == true) {
        perform_bimolecular_state_change_sphere(stateChangeIface, facilitatorIface, rxnItr,
            stateChangeMol, facilitatorMol, stateChangeCom, facilitatorCom,
            counterArrays, params, forwardRxns,
            backRxns, moleculeList, complexList,
            molTemplateList, observablesList, membraneObject);
    } else {
        perform_bimolecular_state_change_box(stateChangeIface, facilitatorIface, rxnItr,
            stateChangeMol, facilitatorMol, stateChangeCom, facilitatorCom,
            counterArrays, params, forwardRxns,
            backRxns, moleculeList, complexList,
            molTemplateList, observablesList, membraneObject);
    }
}
