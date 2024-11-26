#include <algorithm>

#include "boundary_conditions/reflect_functions.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "tracing.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

void sweep_separation_complex_rot_fiber_box(
    int simItr, int selfIndex, Parameters &params,
    std::vector<Molecule> &moleculeList, std::vector<Complex> &complexList,
    const std::vector<ForwardRxn> &forwardRxns,
    const std::vector<MolTemplate> &molTemplateList,
    const Membrane &membraneObject) {
  /*
    In this version, complex comIndex1 is on the membrane.
    If both proteins are on the membrane (Dz==0), evaluate only xy displacement,
   not z. In this _complex_ version, it tests overlap not just for each protein,
   but for each complex, so all the proteins in a complex, before performing
   position updates. NEW 2018: IN THIS VERSION, IT DOES NOT ATTEMPT TO SOLVE
   OVERLAP FOR PROTEINS WITHIN THE SAME COMPLEX, SINCE THEY CANNOT DIFFUSE
   RELATIVE TO ONE ANOTHER! This routine checks whether protein p1 is
   overlapping any partners in its reaction zone at its new position that is
   given by its current position +traj. If it does overlap, the displacement by
   traj is rejected and a new position for itself and any overlapping partners
   are selected. Once it no longer overlaps anyone, this protein and its complex
   are moved and the partners retain their stored new displacements. If a
   protein has already updated its position (done in sequential order) then it
   cannot resample a new position, the current protein must still continue to
   avoid overlapping, however.
  */
  // TRACE();
  int comIndex1{moleculeList[selfIndex].myComIndex};
  int startProIndex = selfIndex;
  int com1Size = complexList[comIndex1].memberList.size();

  int maxRows{1};
  for (auto memMol : complexList[comIndex1].memberList) {
    if (moleculeList[memMol].crossbase.size() > maxRows)
      maxRows = moleculeList[memMol].crossbase.size();
  }

  int ifaceList[maxRows * com1Size];
  int overlapList[maxRows * com1Size];
  int memCheckList[maxRows * com1Size];

  // int reflectList[complexList.size()]; // if this is 0, we need call
  // reflect_traj; if this is 1, we do not need call reflect_traj for (auto i :
  // reflectList) { // initialize
  //     i = 0;
  // }

  /*The sampled displacement for p1 is stored in traj. the position from the
   previous step is still stored in bases[p1].xcom, etc, and will be updated
   at the end of this routine*/

  /*figure out i2*/
  for (int c{0}; c < com1Size; ++c) {
    int proIndex = complexList[comIndex1].memberList[c];
    for (int i{0}; i < moleculeList[proIndex].crossbase.size(); ++i) {
      int p2{moleculeList[proIndex].crossbase[i]};
      int k2{moleculeList[p2].myComIndex};
      // if (complexList[k2].D.z < 1E-15) {
      if (complexList[k2].onFiber) {
        memCheckList[maxRows * c + i] = 1;
      } else
        memCheckList[maxRows * c + i] = 0;
      int i1{moleculeList[proIndex].mycrossint[i]};
      std::array<int, 3> rxnItr = moleculeList[proIndex].crossrxn[i];

      // get the partner interface
      ifaceList[maxRows * c + i] =
          (forwardRxns[rxnItr[0]].reactantListNew[0].relIfaceIndex == i1)
              ? forwardRxns[rxnItr[0]].reactantListNew[1].relIfaceIndex
              : forwardRxns[rxnItr[0]].reactantListNew[0].relIfaceIndex;
    }
  }

  // determine RS3Dinput
  double RS3Dinput{0.0};
  Complex targCom{complexList[comIndex1]};
  if (membraneObject.implicitLipid) {
    for (auto &molIndex : targCom.memberList) {
      for (int RS3Dindex = 0; RS3Dindex < 100; RS3Dindex++) {
        if (std::abs(membraneObject.RS3Dvect[RS3Dindex + 400] -
                     moleculeList[molIndex].molTypeIndex) < 1E-2) {
          RS3Dinput = membraneObject.RS3Dvect[RS3Dindex + 300];
          break;
        }
      }
    }
  }

  int tsave = 0;
  // if (reflectList[comIndex1] == 0) {
  //     reflect_traj_complex_rad_rot(params, moleculeList,
  //     complexList[comIndex1], membraneObject, RS3Dinput);
  //     reflectList[comIndex1] = 1;
  // }

  int itr{0};
  int maxItr{20};
  int saveit{0};
  while (itr < maxItr) {
    int numOverlap{0};
    bool hasOverlap{false};
    int comIndex2{};
    double dr2{};
    for (unsigned memMolItr{0};
         memMolItr < complexList[comIndex1].memberList.size(); ++memMolItr) {
      int pro1Index = complexList[comIndex1].memberList[memMolItr];
      for (int crossMemItr{0};
           crossMemItr < moleculeList[pro1Index].crossbase.size();
           ++crossMemItr) {
        int pro2Index{moleculeList[pro1Index].crossbase[crossMemItr]};
        if (moleculeList[pro2Index].isImplicitLipid) continue;

        comIndex2 = moleculeList[pro2Index].myComIndex;
        /*Do not sweep for overlap if proteins are in the same complex, they
         * cannot move relative to one another!
         */
        if (comIndex1 != comIndex2) {
          int relIface1{moleculeList[pro1Index].mycrossint[crossMemItr]};
          int rxnItr{moleculeList[pro1Index].crossrxn[crossMemItr][0]};
          int relIface2{ifaceList[maxRows * memMolItr + crossMemItr]};

          Vector iface1Vec{
              moleculeList[pro1Index].interfaceList[relIface1].coord -
              complexList[comIndex1].comCoord};
          std::array<double, 9> M =
              create_euler_rotation_matrix(complexList[comIndex1].trajRot);
          iface1Vec = matrix_rotate(iface1Vec, M);

          double dx1{complexList[comIndex1].comCoord.x + iface1Vec.x +
                     complexList[comIndex1].trajTrans.x};
          double dy1{complexList[comIndex1].comCoord.y + iface1Vec.y +
                     complexList[comIndex1].trajTrans.y};
          double dz1{complexList[comIndex1].comCoord.z + iface1Vec.z +
                     complexList[comIndex1].trajTrans.z};

          /*Now complex 2*/
          // if (reflectList[comIndex2] == 0) {
          //     reflect_traj_complex_rad_rot(params, moleculeList,
          //     complexList[comIndex2], membraneObject, RS3Dinput);
          //     reflectList[comIndex2] = 1;
          // }

          Vector iface2Vec{
              moleculeList[pro2Index].interfaceList[relIface2].coord -
              complexList[comIndex2].comCoord};
          std::array<double, 9> M2 =
              create_euler_rotation_matrix(complexList[comIndex2].trajRot);
          iface2Vec = matrix_rotate(iface2Vec, M2);
          double dx2{complexList[comIndex2].comCoord.x + iface2Vec.x +
                     complexList[comIndex2].trajTrans.x};
          double dy2{complexList[comIndex2].comCoord.y + iface2Vec.y +
                     complexList[comIndex2].trajTrans.y};
          double dz2{complexList[comIndex2].comCoord.z + iface2Vec.z +
                     complexList[comIndex2].trajTrans.z};

          /*separation*/
          double dfx{dx1 - dx2};
          double dfy{dy1 - dy2};
          double dfz{dz1 - dz2};

          double sepLimit_sq =
              forwardRxns[rxnItr].bindRadius * forwardRxns[rxnItr].bindRadius;

          dr2 = (dfx * dfx);
          if (memCheckList[maxRows * memMolItr + crossMemItr] != 1) {
            dr2 += (dfy * dfy) + (dfz * dfz);
            if (dr2 < sepLimit_sq) {
              /*reselect positions for protein p2*/
              overlapList[numOverlap] = pro2Index;
              numOverlap++;
              hasOverlap = true;
            }
          } else {
            double dfx_origin = complexList[comIndex1].comCoord.x + iface1Vec.x -
                                complexList[comIndex2].comCoord.x - + iface2Vec.x;
            // no hop and no overlap
            if (dfx_origin * dfx < 0 || dr2 < sepLimit_sq) {
              // std::cout << "Iteration:" << simItr << " resample:" << std::endl;
              // std::cout << "dr^2:" << dr2 << std::endl;
              // std::cout << "minimum:" << sepLimit_sq << std::endl;
              // std::cout << "dx0:" << dfx_origin << std::endl;
              // std::cout << "dx1:" << dfx << std::endl;
              /*reselect positions for protein p2*/
              overlapList[numOverlap] = pro2Index;
              numOverlap++;
              hasOverlap = true;
            }
          }
          //   std::cout << "dr^2:" << dr2 << std::endl;
          //   std::cout << "sigma:" << forwardRxns[rxnItr].bindRadius <<
          //   std::endl; std::cout << "overlapSep:" << params.overlapSepLimit
          //   << std::endl;
        }  // ignore proteins within the same complex.
      }
    }
    /*Now resample positions of p1 and overlapList, if t>0, otherwise no
     overlap, so break from loop*/
    if (hasOverlap) {
      // std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" <<
      // std::endl; std::cout << complexList[comIndex1].comCoord.x +
      // complexList[comIndex1].trajTrans.x << ","
      //           << complexList[comIndex2].comCoord.x +
      //           complexList[comIndex2].trajTrans.x << std::endl;
      ++itr;
      complexList[comIndex1].trajTrans.x =
          sqrt(2.0 * params.timeStep * complexList[comIndex1].D.x) * GaussV();
      complexList[comIndex1].trajTrans.y =
          sqrt(2.0 * params.timeStep * complexList[comIndex1].D.y) * GaussV();
      complexList[comIndex1].trajTrans.z =
          sqrt(2.0 * params.timeStep * complexList[comIndex1].D.z) * GaussV();
      complexList[comIndex1].trajRot.x =
          sqrt(2.0 * params.timeStep * complexList[comIndex1].Dr.x) * GaussV();
      complexList[comIndex1].trajRot.y =
          sqrt(2.0 * params.timeStep * complexList[comIndex1].Dr.y) * GaussV();
      complexList[comIndex1].trajRot.z =
          sqrt(2.0 * params.timeStep * complexList[comIndex1].Dr.z) * GaussV();

      // reflectList[comIndex1] = 0;
      reflect_traj_complex_rad_rot(params, moleculeList, complexList[comIndex1],
                                   membraneObject, RS3Dinput,
                                   false);  // set isInsideCompartment = false.
      // reflectList[comIndex1] = 1;

      int resampleList[complexList.size()];  // if this is 0, we need resample
      for (auto &i : resampleList) {         // initialize
        i = 0;
      }

      for (int checkMolItr{0}; checkMolItr < numOverlap; checkMolItr++) {
        int p2{overlapList[checkMolItr]};
        comIndex2 = moleculeList[p2].myComIndex;

        // avoid repeated resample
        if (resampleList[comIndex2] == 0) {
          // if (p2 > startProIndex && (moleculeList[p2].trajStatus ==
          // TrajStatus::none || moleculeList[p2].trajStatus ==
          // TrajStatus::canBeResampled)) {
          if ((moleculeList[p2].trajStatus == TrajStatus::none ||
               moleculeList[p2].trajStatus == TrajStatus::canBeResampled)) {
            /*
         We loop over proteins sequentially, so earlier proteins have already
         moved and avoided their neighbors and should not be moved again. These
         new positions selected for proteins not yet moved will be stored and
         then used when they test for overlap themselves.
         */

            /*If p2 just dissociated, also don't try to move again*/
            complexList[comIndex2].trajTrans.x =
                sqrt(2.0 * params.timeStep * complexList[comIndex2].D.x) *
                GaussV();
            complexList[comIndex2].trajTrans.y =
                sqrt(2.0 * params.timeStep * complexList[comIndex2].D.y) *
                GaussV();
            complexList[comIndex2].trajTrans.z =
                sqrt(2.0 * params.timeStep * complexList[comIndex2].D.z) *
                GaussV();
            complexList[comIndex2].trajRot.x =
                sqrt(2.0 * params.timeStep * complexList[comIndex2].Dr.x) *
                GaussV();
            complexList[comIndex2].trajRot.y =
                sqrt(2.0 * params.timeStep * complexList[comIndex2].Dr.y) *
                GaussV();
            complexList[comIndex2].trajRot.z =
                sqrt(2.0 * params.timeStep * complexList[comIndex2].Dr.z) *
                GaussV();

            // reflectList[comIndex2] = 0;
            reflect_traj_complex_rad_rot(
                params, moleculeList, complexList[comIndex2], membraneObject,
                RS3Dinput, false);  // set isInsideCompartment = false.
            // reflectList[comIndex2] = 1;
            resampleList[comIndex2] = 1;
          }
        }
      }
      tsave = numOverlap;
    } else {
      saveit = itr + 1;
      itr = maxItr;  // break from loop
    }
  }  // end maximum iterations
     // cancel propagation if no available place found
  // if (saveit == 0) {
  //   // std::cout << "Overlap Check Reached Iteration Limitation" << std::endl;
  //   complexList[comIndex1].trajTrans.zero_crds();
  //   complexList[comIndex1].trajRot.zero_crds();
  // }

  // accept any propagation
  complexList[comIndex1].propagate(moleculeList, membraneObject,
                                   molTemplateList);

  // Reset displacements to zero so distance is measured to your current
  // updated position that won't change again this turn
  complexList[comIndex1].trajTrans.zero_crds();
  complexList[comIndex1].trajRot.zero_crds();
}