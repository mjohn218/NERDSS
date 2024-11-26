#include <cstdio>

#include "boundary_conditions/reflect_functions.hpp"
#include "error/error.hpp"
#include "macro.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

void create_complex_propagation_vectors(
    const Parameters& params, Complex& targCom,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    const std::vector<MolTemplate>& molTemplateList,
    const Membrane& membraneObject) {
  if (targCom.D.z < 1E-14 &&
      membraneObject.isSphere == true) {  // on sphere surface
    Coord targTrans =
        create_complex_propagation_vectors_on_sphere(params, targCom);
    // now setup translational motion vector
    targCom.trajTrans.x = targTrans.x;
    targCom.trajTrans.y = targTrans.y;
    targCom.trajTrans.z = targTrans.z;

    // now setup the rotation
    targCom.trajRot.x = sqrt(2.0 * params.timeStep * targCom.Dr.x) * GaussV();
    targCom.trajRot.y = sqrt(2.0 * params.timeStep * targCom.Dr.y) * GaussV();
    targCom.trajRot.z = sqrt(2.0 * params.timeStep * targCom.Dr.z) * GaussV();

  } else {  // box system or inside the sphere (not on the sphere)
    // Create Gaussian distributed random translational motion
    targCom.trajTrans.x = sqrt(2.0 * params.timeStep * targCom.D.x) * GaussV();
    targCom.trajTrans.y = sqrt(2.0 * params.timeStep * targCom.D.y) * GaussV();
    targCom.trajTrans.z = sqrt(2.0 * params.timeStep * targCom.D.z) * GaussV();
    // create Gaussian distributed random rotational motion
    targCom.trajRot.x = sqrt(2.0 * params.timeStep * targCom.Dr.x) * GaussV();
    targCom.trajRot.y = sqrt(2.0 * params.timeStep * targCom.Dr.y) * GaussV();
    targCom.trajRot.z = sqrt(2.0 * params.timeStep * targCom.Dr.z) * GaussV();
  }

  if (VERBOSE) {
    printf("Prop Vec has been created for Complex (%d) \n", targCom.id);
    printf("COM of Complex: x (%.3f), y (%.3f), z (%.3f)\n", targCom.comCoord.x,
           targCom.comCoord.y, targCom.comCoord.z);
    printf("Members of Complex (%d): ", targCom.id);
    for (auto i : targCom.memberList) {
      printf("mol (%d), ", moleculeList[i].id);
      printf("COM of mol: x (%.3f), y (%.3f), z (%.3f)\n",
             moleculeList[i].comCoord.x, moleculeList[i].comCoord.y,
             moleculeList[i].comCoord.z);
    }
    printf(
        "\nProp Vec: Trans.x (%.3f), Trans.y (%.3f), Trans.z (%.3f),\n\t Rot.x "
        "(%.3f), Rot.y (%.3f), Rot.z (%.3f), \n",
        targCom.trajTrans.x, targCom.trajTrans.y, targCom.trajTrans.z,
        targCom.trajRot.x, targCom.trajRot.y, targCom.trajRot.z);
  }

  // determine RS3Dinput
  double RS3Dinput{0.0};

  for (auto& molIndex : targCom.memberList) {
    for (int RS3Dindex = 0; RS3Dindex < 100; RS3Dindex++) {
      if (std::abs(membraneObject.RS3Dvect[RS3Dindex + 400] -
                   moleculeList[molIndex].molTypeIndex) < 1E-2) {
        RS3Dinput = membraneObject.RS3Dvect[RS3Dindex + 300];
        break;
      }
    }
  }

  if (VERBOSE) {
    printf("RS3D value (%.3f)\n", RS3Dinput);
  }

  //  reflect the boundary and also check_if_span
  reflect_traj_complex_rad_rot(params, moleculeList, targCom, membraneObject,
                               RS3Dinput);
}
