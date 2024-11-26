#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

void sweep_separation_complex_rot_fiber(int simItr, int pro1Index, Parameters& params, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns, const std::vector<MolTemplate>& molTemplateList, const Membrane& membraneObject)
{
    // TRACE();
    // TODO: if (membraneObject.isSphere)
        // sweep_separation_complex_rot_memtest_sphere(simItr, pro1Index, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);
    // else
        sweep_separation_complex_rot_fiber_box(simItr, pro1Index, params, moleculeList, complexList, forwardRxns, molTemplateList, membraneObject);
}
