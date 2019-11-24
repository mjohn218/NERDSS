#include "trajectory_functions/trajectory_functions.hpp"
#include "math/rand_gsl.hpp"
#include "math/matrix.hpp"
#include "boundary_conditions/reflect_functions.hpp"

void create_complex_propagation_vectors(const Parameters& params, Complex& targCom, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, 
                                        const std::vector<MolTemplate>& molTemplateList, const Membrane &membraneObject)
{ 
    // Create Gaussian distributed random translational motion
    targCom.trajTrans.x = sqrt(2.0 * params.timeStep * targCom.D.x) * GaussV();
    targCom.trajTrans.y = sqrt(2.0 * params.timeStep * targCom.D.y) * GaussV();
    targCom.trajTrans.z = sqrt(2.0 * params.timeStep * targCom.D.z) * GaussV();
    // create Gaussian distributed random rotational motion
    targCom.trajRot.x = sqrt(2.0 * params.timeStep * targCom.Dr.x) * GaussV();
    targCom.trajRot.y = sqrt(2.0 * params.timeStep * targCom.Dr.y) * GaussV();
    targCom.trajRot.z = sqrt(2.0 * params.timeStep * targCom.Dr.z) * GaussV();
    std::array<double, 9> M = create_euler_rotation_matrix(targCom.trajRot);
    //reflect_traj_complex_rad_rot(params, moleculeList, targCom, M, membraneObject);
    reflect_traj_complex_rad_rot_new(params, moleculeList, targCom, M, membraneObject);
}
