#include "boundary_conditions/reflect_functions.hpp"
#include "classes/class_Rxns.hpp"
#include "io/io.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/association/association.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>

void associate_box(long long int iter, int ifaceIndex1, int ifaceIndex2,
                   Molecule &reactMol1, Molecule &reactMol2, Complex &reactCom1,
                   Complex &reactCom2, const Parameters &params,
                   ForwardRxn &currRxn, std::vector<Molecule> &moleculeList,
                   std::vector<MolTemplate> &molTemplateList,
                   std::map<std::string, int> &observablesList,
                   copyCounters &counterArrays,
                   std::vector<Complex> &complexList, Membrane &membraneObject,
                   const std::vector<ForwardRxn> &forwardRxns,
                   const std::vector<BackRxn> &backRxns, std::ofstream& assocDissocFile) {
  
  // bool monitor{false};
  // if (iter == 200000295){
  //   if (reactMol1.index == 32 || reactMol2.index == 32) {
  //     monitor = true;
  //   }
  // }

  bool isInsideCompartment {molTemplateList[reactMol1.molTypeIndex].insideCompartment};
  if (reactCom1.index == reactCom2.index) {
    // skip to protein interation updates
    // std::cout << "Closing a loop, no rotations performed.\n";
    counterArrays.nLoops++;
    // update the Molecule's TrajStatus (this is done in the else, when
    // Molecules are rotated but not otherwise)
    for (auto &memMol : reactCom1.memberList)
      moleculeList[memMol].trajStatus = TrajStatus::propagated;
  } else { // not in the same complex

    // record the previous lastNumberUpdateItrEachMol & numEachMol
    std::vector<int> numEachMolPrevious1{};
    std::vector<long long int> lastNumberUpdateItrEachMolPrevious1{};
    std::vector<int> numEachMolPrevious2{};
    std::vector<long long int> lastNumberUpdateItrEachMolPrevious2{};
    for (unsigned index = 0; index < molTemplateList.size(); index++) {
      numEachMolPrevious1.emplace_back(
          complexList[reactMol1.myComIndex].numEachMol[index]);
      lastNumberUpdateItrEachMolPrevious1.emplace_back(
          complexList[reactMol1.myComIndex].lastNumberUpdateItrEachMol[index]);
      numEachMolPrevious2.emplace_back(
          complexList[reactMol2.myComIndex].numEachMol[index]);
      lastNumberUpdateItrEachMolPrevious2.emplace_back(
          complexList[reactMol2.myComIndex].lastNumberUpdateItrEachMol[index]);
    }

    // set up temporary coordinates (= current coordinates)
    for (auto &memMol : reactCom1.memberList)
      moleculeList[memMol].set_tmp_association_coords();
    for (auto &mol : reactCom2.memberList)
      moleculeList[mol].set_tmp_association_coords();

    // create references to reacting interfaces
    Coord &reactIface1 = reactMol1.tmpICoords[ifaceIndex1];
    Coord &reactIface2 = reactMol2.tmpICoords[ifaceIndex2];

    // orientation corrections for membrane bound components
    bool isOnMembrane = false;
    bool transitionToSurface = false;
    bool isOnFiber = false;
    bool transitionToFiber = false;
    Molecule memProtein;

    // reactMol1.display_all();
    // reactMol2.display_all();

    int slowPro = reactMol2.index;
    if (reactCom1.D.x + reactCom1.D.y + reactCom1.D.z <
        reactCom2.D.x + reactCom2.D.y + reactCom2.D.z) {
      slowPro = reactMol1.index;
      memProtein = reactMol1; // rotate relative to the slower protein.
    } else
      memProtein = reactMol2;

    double tol = 1E-14;
    /* Calculate COM of the two complexes pre-association. The COM of the
     * new complex after should be close to this. Here, we will force it
     * back, as rotation can cause large displacements*/
    Coord startCOM; // COM of two reactants
    com_of_two_tmp_complexes(reactCom1, reactCom2, startCOM, moleculeList);
    Coord startCOM1; // COM of protein 1
    com_of_two_tmp_complexes(reactCom1, reactCom1, startCOM1, moleculeList);
    Coord startCOM2; // COM of protein 2
    com_of_two_tmp_complexes(reactCom2, reactCom2, startCOM2, moleculeList);

    double zCom1Temp{};
    double yCom1Temp{};
    double zCom2Temp{};
    double yCom2Temp{};
    // Assume membrane's normal on z direction, and fiber on x direction
    if (reactCom1.OnSurface) {
      // com1 is the 2D complex
      zCom1Temp = startCOM1.z;
    } else if (reactCom1.onFiber) {
      // com2 is the 1D complex
      zCom1Temp = startCOM1.z;
      yCom1Temp = startCOM1.y;
    }
    if (reactCom2.OnSurface) {
      // com2 is the 2D complex
      zCom2Temp = startCOM2.z;
    } else if (reactCom2.onFiber) {
      // com2 is the 1D complex
      zCom2Temp = startCOM2.z;
      yCom2Temp = startCOM2.y;
    }

    /*
    // calibrate for some fixed proteins
    double xCom1Temp{};
    double xCom2Temp{};
    if (reactCom1.D.x < tol) {
      xCom1Temp = startCOM1.x;
    }
    if (reactCom2.D.x < tol) {
        xCom2Temp = startCOM2.x;
    }
    double yCom1Temp{};
    double yCom2Temp{};
    if (reactCom1.D.y < tol) {
      yCom1Temp = startCOM1.y;
    }
    if (reactCom2.D.y < tol) {
      yCom2Temp = startCOM2.y;
    }
    double zCom1Temp{};
    double zCom2Temp{};
    if (reactCom1.D.z < tol) {
      zCom1Temp = startCOM1.z;
    }
    if (reactCom2.D.z < tol) {
      zCom2Temp = startCOM2.z;
    }
    */

    /* MOVE PROTEIN TO SIGMA */

    double DxSum{reactCom1.D.x + reactCom2.D.x};
    double DySum{reactCom1.D.y + reactCom2.D.y};
    double DzSum{reactCom1.D.z + reactCom2.D.z};

    Vector sigma{reactIface1 - reactIface2};
    Vector transVec1{};
    Vector transVec2{};
    double displaceFrac{};

    double sigmaMag3D =
        sqrt((sigma.x * sigma.x) + (sigma.y * sigma.y) + (sigma.z * sigma.z));
    // sigma.calc_magnitude();
    double displaceFrac3D = (sigmaMag3D - currRxn.bindRadius) / sigmaMag3D;

    // if both in 2D, ignore the z-component
    if (reactCom1.OnSurface == true && reactCom2.OnSurface == true) {
      isOnMembrane = true;
      /*Store coordinates of one protein to recover membrane-bound orientation*/

      DzSum = 1; // to prevent divide by 0
      if (std::abs(std::abs(sigma.z) - currRxn.bindRadius) < 1E-6) {
        // if entirety of sigma is in z-component, ignore x and y
        displaceFrac = 1;
      } else {
        double sigmaMag = sqrt((sigma.x * sigma.x) + (sigma.y * sigma.y));
        displaceFrac = (sigmaMag - currRxn.bindRadius) / sigmaMag;
      }
      // if both in 1D, ignore the y,z-component
    } 
    // if both in 1D
    else if (reactCom1.onFiber == true && reactCom2.onFiber == true) {
      isOnFiber = true;
      DzSum = 1; // to prevent divide by 0
      DySum = 1; // to prevent divide by 0
      if (std::abs(sigma.x) < 1E-6) {
        // if there is no seperation of sigma is in x-component
        displaceFrac = 1;
      } else {
        double sigmaMag = sqrt((sigma.x * sigma.x));
        // randomly move to the left or the right side
        double randDouble = rand_gsl() < 0.5 ? 1.0 : -1.0;
        displaceFrac = (sigmaMag - randDouble * currRxn.bindRadius) / sigmaMag;
        // std::cout << "frac: " << displaceFrac << " rand: " << randDouble << std::endl; 
      }
    }
    // not both in 2D nor both in 1D
    else {
      // Supposed to have at least one protein in 3D
      // TODO: Validate if at least one protein in 3D?
      double sigmaMag = sqrt((sigma.x * sigma.x) + (sigma.y * sigma.y) + (sigma.z * sigma.z));
      // sigma.calc_magnitude();
      displaceFrac = (sigmaMag - currRxn.bindRadius) / sigmaMag;
      if (reactCom1.OnSurface || reactCom2.OnSurface) {
        transitionToSurface = true; // both can't be less than tol, or would not be in this loop.
        // std::cout << "TRANSITIONING FROM 3D->2D " << std::endl;
      }
      if (reactCom1.onFiber || reactCom2.onFiber) {
        transitionToFiber = true;
      }
    }

    // if (monitor){
    //   std::cout << "displacement frac" << std::endl;
    //   std::cout << displaceFrac << std::endl;
    //   std::cout << "On fiber?" << isOnFiber << std::endl;
    // }

    if (isOnFiber == true) {
      // Complex conformation has first priority. Fiber binding is relaxed.
      // keep the slow protein on fiber (only move along x direction)
      // move the other reactant (fast protein) freely in 3 dimension
      // The fast protein will be pushed back to fiber in orientation correction.
      //TODO: For explicit model, if the Complex conformation requires one protein to detach from the fiber,
      //      the binding site (which should be always on fiber) will also be moved away
      if (slowPro == reactMol1.index){
        if (abs(reactCom1.D.x)<1e-6){
          transVec1.x = 0.0;
          // std::cout << "found Dx=0" << std::endl;
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// This is commented out when merging
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// void associate_box(long long int iter, int ifaceIndex1, int ifaceIndex2, 
//     Molecule& reactMol1, Molecule& reactMol2, Complex& reactCom1,
//     Complex& reactCom2, const Parameters& params, 
//     ForwardRxn& currRxn, std::vector<Molecule>& moleculeList,
//     std::vector<MolTemplate>& molTemplateList,
//     std::map<std::string, int>& observablesList, 
//     copyCounters& counterArrays, 
//     std::vector<Complex>& complexList, Membrane& membraneObject, 
//     const std::vector<ForwardRxn>& forwardRxns, 
//     const std::vector<BackRxn>& backRxns)
// {
//     bool isInsideCompartment {molTemplateList[reactMol1.molTypeIndex].insideCompartment};
//     if (reactCom1.index == reactCom2.index) {
//         // skip to protein interation updates
//         //std::cout << "Closing a loop, no rotations performed.\n";
//         counterArrays.nLoops++;
//         // update the Molecule's TrajStatus (this is done in the else, when Molecules are rotated but not otherwise)
//         for (auto& memMol : reactCom1.memberList)
//             moleculeList[memMol].trajStatus = TrajStatus::associated;
//     } else { //not in the same complex

//         // record the previous lastNumberUpdateItrEachMol & numEachMol
//         std::vector<int> numEachMolPrevious1 {};
//         std::vector<long long int> lastNumberUpdateItrEachMolPrevious1 {};
//         std::vector<int> numEachMolPrevious2 {};
//         std::vector<long long int> lastNumberUpdateItrEachMolPrevious2 {};


//         for (unsigned index = 0; index < molTemplateList.size(); index++) {
//             numEachMolPrevious1.emplace_back(complexList[reactMol1.myComIndex].numEachMol[index]);
//             lastNumberUpdateItrEachMolPrevious1.emplace_back(complexList[reactMol1.myComIndex].lastNumberUpdateItrEachMol[index]);
//             numEachMolPrevious2.emplace_back(complexList[reactMol2.myComIndex].numEachMol[index]);
//             lastNumberUpdateItrEachMolPrevious2.emplace_back(complexList[reactMol2.myComIndex].lastNumberUpdateItrEachMol[index]);
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        }
        else{
          transVec1.x = -sigma.x * (reactCom1.D.x / DxSum) * displaceFrac;
        }
        transVec1.y = 0.0;
        transVec1.z = 0.0;

        if (abs(reactCom2.D.x)<1e-6){
          transVec2.x = 0.0;
          // std::cout << "found Dx=0" << std::endl;
        }
        else{
          transVec2.x = sigma.x * (reactCom2.D.x / DxSum) * displaceFrac;
        }
        transVec2.y = sigma.y * displaceFrac3D;
        transVec2.z = sigma.z * displaceFrac3D;
      } else {
        if (abs(reactCom1.D.x) < 1e-6) {
          transVec1.x = 0.0;
          // std::cout << "found Dx=0" << std::endl;
        } else {
          transVec1.x = -sigma.x * (reactCom1.D.x / DxSum) * displaceFrac;
        }
        transVec1.y = -sigma.y * displaceFrac3D;
        transVec1.z = -sigma.z * displaceFrac3D;

        if (abs(reactCom2.D.x) < 1e-6) {
          transVec2.x = 0.0;
          // std::cout << "found Dx=0" << std::endl;
        } else {
          transVec2.x = sigma.x * (reactCom2.D.x / DxSum) * displaceFrac;
        }
        transVec2.y = 0.0;
        transVec2.z = 0.0;
      }
    } else {
      transVec1.x = -sigma.x * (reactCom1.D.x / DxSum) * displaceFrac;
      transVec1.y = -sigma.y * (reactCom1.D.y / DySum) * displaceFrac;
      transVec1.z = -sigma.z * (reactCom1.D.z / DzSum) * displaceFrac;

      transVec2.x = sigma.x * (reactCom2.D.x / DxSum) * displaceFrac;
      transVec2.y = sigma.y * (reactCom2.D.y / DySum) * displaceFrac;
      transVec2.z = sigma.z * (reactCom2.D.z / DzSum) * displaceFrac;
    }

    // if (monitor){
    //   std::cout << "translation vectors" << std::endl;
    //   std::cout << reactMol1.index << std::endl;
    //   std::cout << transVec1.x << "," << transVec1.y << "," << transVec1.z << std::endl;
    //   std::cout << "&" << std::endl; /////////////////////////
    //   std::cout << reactMol2.index << std::endl;
    //   std::cout << transVec2.x << "," << transVec2.y << "," << transVec2.z << std::endl;
    // }

    // update the temporary coordinates
    for (auto &mp : reactCom1.memberList)
      moleculeList[mp].update_association_coords(transVec1);
    for (auto &mp : reactCom2.memberList)
      moleculeList[mp].update_association_coords(transVec2);

    // matches move protein to sigma
    Coord afterSigmaCOM1;
    com_of_two_tmp_complexes(reactCom1, reactCom1, afterSigmaCOM1, moleculeList);
    Coord afterSigmaCOM2;
    com_of_two_tmp_complexes(reactCom2, reactCom2, afterSigmaCOM2, moleculeList);

    // std::cout << "Diplace during pushing to sigma: " << std::endl;
    // std::cout << "COM1: " << std::endl;
    // std::cout << sqrt((afterSigmaCOM1.x - startCOM1.x) *
    // (afterSigmaCOM1.x - startCOM1.x) + (afterSigmaCOM1.y - startCOM1.y)
    // * (afterSigmaCOM1.y - startCOM1.y) + (afterSigmaCOM1.z -
    // startCOM1.z) * (afterSigmaCOM1.z - startCOM1.z)) << std::endl;
    // std::cout << "COM2: " << std::endl;
    // std::cout << sqrt((afterSigmaCOM2.x - startCOM2.x) *
    // (afterSigmaCOM2.x - startCOM2.x) + (afterSigmaCOM2.y - startCOM2.y)
    // * (afterSigmaCOM2.y - startCOM2.y) + (afterSigmaCOM2.z -
    // startCOM2.z) * (afterSigmaCOM2.z - startCOM2.z)) << std::endl;

    if (molTemplateList[reactMol1.molTypeIndex].isPoint && molTemplateList[reactMol2.molTypeIndex].isPoint) {
    //   if (monitor){
    //   std::cout << "Correct: both are points" << std::endl;
    // }
      /*If both molecules are points, no orientations to specify*/
      //   std::cout << " Move two point particles to contact along current
      //   separation vector, NO ORIENTATION \n";
    } else { 
      // Not both are points
      // rotate the proteins ignoring all geometry restrictions
      /* THETA */
      // std::cout << std::setw(8) << std::setfill('-') << ' ' << std::endl
      //           << "THETA 1" << std::endl
      //           << std::setw(8) << ' ' << std::setfill(' ') << std::endl;
      // if (iter == 818823) {
        // std::cout << "Beofre ROT: " << std::endl;
        // reactMol1.display_assoc_icoords("mol1");
        // reactMol2.display_assoc_icoords("mol2");
        // reactCom1.display();
        // std::cout << reactMol1.index << std::endl;
        // std::cout << "Corresponding Iface" << reactIface1 << std::endl;
        // reactCom2.display();
        // std::cout << reactMol2.index << std::endl;
        // std::cout << "Corresponding Iface" << reactIface2 << std::endl;
      // }

      if (!std::isnan(currRxn.assocAngles.theta1))
        theta_rotation(reactIface1, reactIface2, reactMol1, reactMol2,
                       currRxn.assocAngles.theta1, reactCom1, reactCom2,
                       moleculeList);

      // if (iter == 818823) {
      //   std::cout << "after theta1 rotation: " << std::endl;
      //   reactMol1.display_assoc_icoords("mol1");
      //   reactMol2.display_assoc_icoords("mol2");
      // }
      // else
      //   std::cout <<"No THETA1 !"<<std::endl;
      //     std::cout << std::setw(30) << std::setfill('-') << ' ' <<
      //     std::setfill(' ') << std::endl; std::cout << "THETA 2" << std::endl
      //               << std::setw(8) << std::setfill('-') << ' ' <<
      //               std::setfill(' ') << std::endl;
      if (!std::isnan(currRxn.assocAngles.theta2))
        theta_rotation(reactIface2, reactIface1, reactMol2, reactMol1,
                       currRxn.assocAngles.theta2, reactCom2, reactCom1,
                       moleculeList);

      // if (iter == 818823) {
      //   std::cout << "after theta2 rotation: " << std::endl;
      //   reactMol1.display_assoc_icoords("mol1");
      //   reactMol2.display_assoc_icoords("mol2");
      // }
      // else
      //   std::cout <<" NO THETA 2 "<<std::endl;
      /* OMEGA */
      // if protein has theta M_PI, uses protein norm instead of com_iface
      // vector std::cout << std::setw(6) << std::setfill('-') << ' ' <<
      // std::endl
      //           << "OMEGA" << std::endl
      //           << std::setw(6) << ' ' << std::setfill(' ') << std::endl;
      if (!std::isnan(currRxn.assocAngles.omega)) {
        omega_rotation(reactIface1, reactIface2, ifaceIndex2, reactMol1,
                       reactMol2, reactCom1, reactCom2,
                       currRxn.assocAngles.omega, currRxn, moleculeList,
                       molTemplateList);
      } // else
      // if (iter == 818823) {
      //   std::cout << "after omega rotation: " << std::endl;
      //   reactMol1.display_assoc_icoords("mol1");
      //   reactMol2.display_assoc_icoords("mol2");
      // }
      // std::cout << "P1 or P2 is a rod-type protein, no dihedral for
      // associated complex." << std::endl;

      /* PHI */
      // PHI 1
      // std::cout << std::setw(6) << std::setfill('-') << ' ' << std::endl
      //           << "PHI 1" << std::endl
      //           << std::setw(6) << ' ' << std::setfill(' ') << std::endl;
      if (!std::isnan(currRxn.assocAngles.phi1)) {
        phi_rotation(reactIface1, reactIface2, ifaceIndex2, reactMol1,
                     reactMol2, reactCom1, reactCom2, currRxn.norm1,
                     currRxn.assocAngles.phi1, currRxn, moleculeList,
                     molTemplateList);
      } // else
      // if (iter == 818823) {
      //   std::cout << "after phi1 rotation: " << std::endl;
      //   reactMol1.display_assoc_icoords("mol1");
      //   reactMol2.display_assoc_icoords("mol2");
      // }
      // std::cout << "P1 has no valid phi angle." << std::endl;

      // PHI 2
      // std::cout << std::setw(6) << std::setfill('-') << ' ' << std::endl
      //           << "PHI 2" << std::endl
      //           << std::setw(6) << ' ' << std::setfill(' ') << std::endl;
      if (!std::isnan(currRxn.assocAngles.phi2)) {
        phi_rotation(reactIface2, reactIface1, ifaceIndex1, reactMol2,
                     reactMol1, reactCom2, reactCom1, currRxn.norm2,
                     currRxn.assocAngles.phi2, currRxn, moleculeList,
                     molTemplateList);
      } // else
      // if (iter == 818823) {
      //   std::cout << "after phi2 rotation: " << std::endl;
      //   reactMol1.display_assoc_icoords("mol1");
      //   reactMol2.display_assoc_icoords("mol2");
      // }
      // std::cout << "P2 has no valid phi angle." << std::endl;
    } // end of if points.

    // Now the orientation of interfaces is correct. However, geometry
    // restrictions have been lifted and need correction.

    Coord afterRotateCOM1;
    com_of_two_tmp_complexes(reactCom1, reactCom1, afterRotateCOM1, moleculeList);
    Coord afterRotateCOM2;
    com_of_two_tmp_complexes(reactCom2, reactCom2, afterRotateCOM2, moleculeList);
    
    // if (monitor){
    //   std::cout << "After rotation" << std::endl;
    //   std::cout << reactMol1.index << std::endl;
    //   std::cout << afterRotateCOM1.x << "," << afterRotateCOM1.y << ","
    //             << afterRotateCOM1.z << std::endl;
    //   std::cout << "&" << std::endl; /////////////////////////
    //   std::cout << reactMol2.index << std::endl;
    //   std::cout << afterRotateCOM2.x << "," << afterRotateCOM2.y << ","
    //             << afterRotateCOM2.z << std::endl;
    // }

    // std::cout << "Diplace during rotating: " << std::endl;
    // std::cout << "COM1: " << std::endl;
    // std::cout << sqrt((afterSigmaCOM1.x - afterRotateCOM1.x) *
    // (afterSigmaCOM1.x - afterRotateCOM1.x) + (afterSigmaCOM1.y -
    // afterRotateCOM1.y) * (afterSigmaCOM1.y - afterRotateCOM1.y) +
    // (afterSigmaCOM1.z - afterRotateCOM1.z) * (afterSigmaCOM1.z -
    // afterRotateCOM1.z)) << std::endl; std::cout << "COM2: " << std::endl;
    // std::cout << sqrt((afterSigmaCOM2.x - afterRotateCOM2.x) *
    // (afterSigmaCOM2.x - afterRotateCOM2.x) + (afterSigmaCOM2.y -
    // afterRotateCOM2.y) * (afterSigmaCOM2.y - afterRotateCOM2.y) +
    // (afterSigmaCOM2.z - afterRotateCOM2.z) * (afterSigmaCOM2.z -
    // afterRotateCOM2.z)) << std::endl;

    /*FINISHED ROTATING, NO CONSTRAINTS APPLIED TO SURFACE NOR FIBER REACTIONS*/
    Coord finalCOM;
    com_of_two_tmp_complexes(reactCom1, reactCom2, finalCOM, moleculeList); 
    // com of c1+c2 (final (tmp) coordinates).

    Coord preCOM;
    if (isOnMembrane == true || transitionToSurface == true || isOnFiber || transitionToFiber) {
      // bool allpoints{true}; // check if all molecules are points
      // for (auto &mp : reactCom1.memberList){
      //   if (!molTemplateList[moleculeList[mp].molTypeIndex].isPoint) {
      //     allpoints = false;
      //   }
      // }
      // for (auto &mp : reactCom2.memberList) {
      //   if (!molTemplateList[moleculeList[mp].molTypeIndex].isPoint) {
      //     allpoints = false;
      //   }
      // }
      // if (allpoints) {
        // Both complexes are treated as points
      // } else {
        /*return orientation of normal back to starting position*/
        /*This part is equivalent for fiber and membrane*/
        // std::cout << " IS ON MEMBRANE, CORRECT ORIENTATION ! " << std::endl;
        Quat memRot;
        Coord pivot;
        // also translate the slowPro back to its same COM.

        if (slowPro == reactMol1.index) {
          memRot = save_mem_orientation(
              memProtein, reactMol1, molTemplateList[reactMol1.molTypeIndex]);
          pivot = reactMol1.tmpComCoord;
          preCOM = afterSigmaCOM1;
        } else {
          memRot = save_mem_orientation(
              memProtein, reactMol2, molTemplateList[reactMol2.molTypeIndex]);
          pivot = reactMol2.tmpComCoord;
          preCOM = afterSigmaCOM2;
        }

        rotate(pivot, memRot, reactCom1, moleculeList);
        rotate(pivot, memRot, reactCom2, moleculeList);
      // }
    }

    com_of_two_tmp_complexes(reactCom1, reactCom2, finalCOM, moleculeList);
    Coord afterOrientationCOM1;
    com_of_two_tmp_complexes(reactCom1, reactCom1, afterOrientationCOM1, moleculeList);
    Coord afterOrientationCOM2;
    com_of_two_tmp_complexes(reactCom2, reactCom2, afterOrientationCOM2, moleculeList);

    // if (monitor){
    //   std::cout << "After orientation correction" << std::endl;
    //   std::cout << reactMol1.index << std::endl;
    //   std::cout << afterOrientationCOM1.x << "," << afterOrientationCOM1.y
    //             << "," << afterOrientationCOM1.z << std::endl;
    //   std::cout << "&" << std::endl; /////////////////////////
    //   std::cout << reactMol2.index << std::endl;
    //   std::cout << afterOrientationCOM2.x << "," << afterOrientationCOM2.y
    //             << "," << afterOrientationCOM2.z << std::endl;
    // }

    Coord postCOM;
    if (slowPro == reactMol1.index) {
      postCOM = afterOrientationCOM1;
    } else {
      postCOM = afterOrientationCOM2;
    }
    // std::cout << "Diplace during orientation to the membrane: " << std::endl;
    // std::cout << "COM1: " << std::endl;
    // std::cout << sqrt((afterOrientationCOM1.x - afterRotateCOM1.x) *
    // (afterOrientationCOM1.x - afterRotateCOM1.x) + (afterOrientationCOM1.y -
    // afterRotateCOM1.y) * (afterOrientationCOM1.y - afterRotateCOM1.y) +
    // (afterOrientationCOM1.z - afterRotateCOM1.z) * (afterOrientationCOM1.z -
    // afterRotateCOM1.z)) << std::endl; std::cout << "COM2: " << std::endl;
    // std::cout << sqrt((afterOrientationCOM2.x - afterRotateCOM2.x) *
    // (afterOrientationCOM2.x - afterRotateCOM2.x) + (afterOrientationCOM2.y -
    // afterRotateCOM2.y) * (afterOrientationCOM2.y - afterRotateCOM2.y) +
    // (afterOrientationCOM2.z - afterRotateCOM2.z) * (afterOrientationCOM2.z -
    // afterRotateCOM2.z)) << std::endl;
    Vector dtrans{};
    bool correctPosition{true};
    if (abs(DxSum)<1e-6 && isOnFiber){
      correctPosition = false;
    }
    // std::cout << "correctPosition: " << correctPosition << std::endl;
    if ((isOnMembrane || transitionToSurface || isOnFiber || transitionToFiber)&&correctPosition) {
      /*Force finalCOM to afterSigmaCOM, just for the slowPro.*/
      /*This part is equivalent to fiber and membrane*/

      dtrans.x = preCOM.x - postCOM.x;
      dtrans.y = preCOM.y - postCOM.y;
      // dtrans.z = (float)preCOM.z - (float)postCOM.z;

      dtrans.z = 0.0; // don't move in z, now they are both on membrane

      //   std::cout << "TRANSLATE SLOWPRO TO ORIG SIGMA COM BY SHIFTING: " <<
      //   dtrans.x << ' ' << dtrans.y << ' ' << dtrans.z << std::endl; //
      //   update the temporary coordinates for both complexes
      for (auto &mp : reactCom1.memberList)
        moleculeList[mp].update_association_coords(dtrans);
      for (auto &mp : reactCom2.memberList)
        moleculeList[mp].update_association_coords(dtrans);

      // if (iter == 818823) {
        // std::cout << "NEW AFTER TRANSLATIONAL SHIFT: " << std::endl;
        // reactMol1.display_assoc_icoords("mol1");
        // reactMol2.display_assoc_icoords("mol2");
      // }
      Coord afterTranShiftCOM1;
      com_of_two_tmp_complexes(reactCom1, reactCom1, afterTranShiftCOM1, moleculeList);
      Coord afterTranShiftCOM2;
      com_of_two_tmp_complexes(reactCom2, reactCom2, afterTranShiftCOM2, moleculeList);

      /*
      FOR 2D AND 3D->2D, the SLOWPRO IS NOW BACK TO ITS ORIGINAL POSITION,
      MEANING THE OTHER PROTEIN DID ALL THE DISPLACING. BASED ON THE AMOUNT THAT
      PRO2 HAD TO REORIENT, USE FRACTIONS OF DR TO ROTATE IT A BIT BACK. SINCE
      THEY ARE BOTH NOW IN 2D, ONLY WANT TO EVALUATE THE DISPLACEMENT THAT HAS
      OCCURED IN 2D.
      */
      double dispAngle, rotAng, skip;
      if (slowPro == reactMol1.index) {
        dispAngle = calc_one_angular_displacement(ifaceIndex2, reactMol2, reactCom2); // pro2 is the one that moved.
        determine_rotation_angles(dispAngle, 0, rotAng, skip, reactCom1, reactCom2); // positive angle

      } else {
        dispAngle = calc_one_angular_displacement(
            ifaceIndex1, reactMol1, reactCom1); // pro1 is the one that moved
        determine_rotation_angles(0, dispAngle, skip, rotAng, reactCom1, reactCom2); // positive angle
      }
      Coord origin{0.5 * (reactIface1 + reactIface2)}; // halfway along the sigma vector, or
                                         // midway between the interfaces.
      /*axis of rotation is the normal to the plane, at the origin position*/
      Vector rotAxis{0, 0, 1}; // For a Box, we can just use the z-axis.
      Quat rotQuatPos(cos(rotAng / 2), sin(rotAng / 2) * rotAxis.x, sin(rotAng / 2) * rotAxis.y, sin(rotAng / 2) * rotAxis.z);
      rotQuatPos = rotQuatPos.unit();

      // rotate the two complexes, using same angle here, so no orientations should change!
      rotate(origin, rotQuatPos, reactCom1, moleculeList);
      rotate(origin, rotQuatPos, reactCom2, moleculeList);
      Coord afterRotShiftCOM1;
      com_of_two_tmp_complexes(reactCom1, reactCom1, afterRotShiftCOM1, moleculeList);
      Coord afterRotShiftCOM2;
      com_of_two_tmp_complexes(reactCom2, reactCom2, afterRotShiftCOM2, moleculeList);

      //   std::cout << "Displace during rotation shift to move COM1: " <<
      //   std::endl; std::cout << "COM1: " << std::endl;
      Vector d1{afterRotShiftCOM1 - afterTranShiftCOM1};
      d1.calc_magnitude();
      //   std::cout <<d1.magnitude<<std::endl;
      //   std::cout << "COM2: " << std::endl;
      Vector d2{afterRotShiftCOM2 - afterTranShiftCOM2};
      d2.calc_magnitude();
      //   std::cout <<d2.magnitude<<std::endl;
      // if (iter == 818823) {
        // std::cout << "NEW AFTER ROTATIONAL ALIGN: " << std::endl;
        // reactMol1.display_assoc_icoords("mol1");
        // reactMol2.display_assoc_icoords("mol2");
      // }
    }
    Coord currCOM1;
    com_of_two_tmp_complexes(reactCom1, reactCom1, currCOM1, moleculeList);
    Coord currCOM2;
    com_of_two_tmp_complexes(reactCom2, reactCom2, currCOM2, moleculeList);

    /*
      CORRECT FIBER PROTEIN POSITION
    */
    // Find the slow protein. All the molecules in this protein should be at
    // correct position.
    // Loop over all molecules. If the molecule is a intrinsic fiber protein,
    // then move it to the x axis ( y=y_waterbox/2 and z=z_waterbox/2 ).
    for (auto &mp : reactCom2.memberList) {
      if (moleculeList[mp].isPromoter == true) {
        Vector dtrans{0, -moleculeList[mp].tmpComCoord.y, -moleculeList[mp].tmpComCoord.z};
        moleculeList[mp].update_association_coords(dtrans);
      }
    }
    for (auto &mp : reactCom1.memberList) {
      if (moleculeList[mp].isPromoter == true) {
        Vector dtrans{0, -moleculeList[mp].tmpComCoord.y, -moleculeList[mp].tmpComCoord.z};
        moleculeList[mp].update_association_coords(dtrans);
      }
    }

    /*
      CHECK IF BELOW BOX
    */
    double zchg = 0;
    bool isBelowBottom = false;
    bool isAboveBottom = false;
    if (isOnMembrane == true || transitionToSurface == true) {
      /* RECHECK HERE IF ANY OF THE LIPIDS ARE SLIGHTLY BELOW THE MEMBRANE. THIS
         CAN HAPPEN DUE TO PRECISION ISSUES always use tmpCoords in this
         associate routine. */
      dtrans.x = 0;
      dtrans.y = 0;
      for (auto &mp : reactCom1.memberList) {
        if (moleculeList[mp].isLipid == true) {
          if (moleculeList[mp].tmpComCoord.z <
              -membraneObject.waterBox.z / 2.0) {
            double ztmp =
                (-membraneObject.waterBox.z / 2.0) -
                moleculeList[mp].tmpComCoord.z; // lipid COM is below box
                                                // bottom, here ztmp is positive
            if (ztmp > zchg) {
              zchg = ztmp; // largest dip below membrane
              isBelowBottom = true;
            }
          }
          if (moleculeList[mp].tmpComCoord.z - 0.01 >
              -membraneObject.waterBox.z / 2.0) {
            double ztmp =
                (-membraneObject.waterBox.z / 2.0) -
                moleculeList[mp].tmpComCoord.z; // lipid COM is ABOVE box
                                                // bottom, here ztmp is negtive
            // std::cout << "WARNING, during associate, LIPID IS ABOVE MEMBRANE
            // BY " << -ztmp << '\n'; move the lipid back to the bottom
            if (-ztmp > zchg) {
              zchg = -ztmp;
              isAboveBottom = true;
            }
          }
        } // this is a lipid
      }
      // if (reactCom1.D.z < tol) { // for the implicit case, we need to compare
      // the z-value of the 2D complex to keep it unchanged
      if (reactCom1.OnSurface) { // Use the molecule crds, since the tmpComCoord
                                 // is not updated
        if (currCOM1.z - zCom1Temp > 0.0) {
          if (currCOM1.z - zCom1Temp > zchg) {
            zchg = currCOM1.z - zCom1Temp;
            isAboveBottom = true;
          }
        }
        if (currCOM1.z - zCom1Temp < 0.0) {
          if (-currCOM1.z + zCom1Temp > zchg) {
            zchg = -currCOM1.z + zCom1Temp;
            isBelowBottom = true;
          }
        }
      }
      for (auto &mp : reactCom2.memberList) {
        if (moleculeList[mp].isLipid == true) {
          if (moleculeList[mp].tmpComCoord.z <
              -membraneObject.waterBox.z / 2.0) {
            double ztmp = (-membraneObject.waterBox.z / 2.0) - moleculeList[mp].tmpComCoord.z; 
            // lipid COM is below box bottom
            if (ztmp > zchg) {
              zchg = ztmp; // largest dip below membrane
              isBelowBottom = true;
            }
          }
          if (moleculeList[mp].tmpComCoord.z - 0.01 >
              -membraneObject.waterBox.z / 2.0) {
            double ztmp = (-membraneObject.waterBox.z / 2.0) - moleculeList[mp].tmpComCoord.z; 
            // lipid COM is ABOVE box bottom, here ztmp is negtive
            // std::cout << "WARNING, during associate, LIPID IS ABOVE MEMBRANE
            // BY " << -ztmp << '\n'; move the lipid back to the bottom
            if (-ztmp > zchg) {
              zchg = -ztmp;
              isAboveBottom = true;
            }
          }
        } // this is a lipid
      }
      // if (reactCom2.D.z < tol) { // for the implicit case, we need to compare
      // the z-value of the 2D complex to keep it unchanged
      if (reactCom2.OnSurface) { // for the implicit case, we need to compare the
                                 // z-value of the 2D complex to keep it unchanged
        if (currCOM2.z - zCom2Temp > 0.0) {
          if (currCOM2.z - zCom2Temp > zchg) {
            zchg = currCOM2.z - zCom2Temp;
            isAboveBottom = true;
          }
        }
        if (currCOM2.z - zCom2Temp < 0.0) {
          if (-currCOM2.z + zCom2Temp > zchg) {
            zchg = -currCOM2.z + zCom2Temp;
            isBelowBottom = true;
          }
        }
      }
      if (isBelowBottom == true) {
        dtrans.z = zchg;
        // std::cout << " Lipid is below membrane, shift up by: " << zchg <<
        // std::endl;
      }
      if (isAboveBottom == true) {
        dtrans.z = -zchg;
        // std::cout << " Lipid is above membrane, shift down by: " << zchg <<
        // std::endl;
      }
      // update the temporary coordinates for both complexes
      for (auto &mp : reactCom1.memberList)
        moleculeList[mp].update_association_coords(dtrans);
      for (auto &mp : reactCom2.memberList)
        moleculeList[mp].update_association_coords(dtrans);
    } // is on membrane

    Coord afterBackCOM1;
    com_of_two_tmp_complexes(reactCom1, reactCom1, afterBackCOM1, moleculeList);
    Coord afterBackCOM2;
    com_of_two_tmp_complexes(reactCom2, reactCom2, afterBackCOM2, moleculeList);


    // std::cout << "Diplace during back to membrane, ignoring trans and rot
    // shift: " << std::endl; std::cout << "COM1: " << std::endl; std::cout <<
    // sqrt((afterOrientationCOM1.x - afterBackCOM1.x) * (afterOrientationCOM1.x
    // - afterBackCOM1.x) + (afterOrientationCOM1.y - afterBackCOM1.y) *
    // (afterOrientationCOM1.y - afterBackCOM1.y) + (afterOrientationCOM1.z -
    // afterBackCOM1.z) * (afterOrientationCOM1.z - afterBackCOM1.z)) <<
    // std::endl; std::cout << "COM2: " << std::endl; std::cout <<
    // sqrt((afterOrientationCOM2.x - afterBackCOM2.x) * (afterOrientationCOM2.x
    // - afterBackCOM2.x) + (afterOrientationCOM2.y - afterBackCOM2.y) *
    // (afterOrientationCOM2.y - afterBackCOM2.y) + (afterOrientationCOM2.z -
    // afterBackCOM2.z) * (afterOrientationCOM2.z - afterBackCOM2.z)) <<
    // std::endl;

    // if (iter == 818823) {
      // std::cout << " FINAL COORDS PRIOR TO OVERLAP CHECK  AND REFLECT OFF BOX: "
      //           << std::endl;
      // reactMol1.display_assoc_icoords("mol1");
      // reactMol2.display_assoc_icoords("mol2");
    // }

    std::array<double, 3> traj; //=new double[3];
    for (int mm = 0; mm < 3; mm++)
      traj[mm] = 0;

    /* This needs to evaluate the traj update, based on it initially being zero.
       And here, it should be called based on the tmpCoords, not the full
       coordinates. also requires updating the COM of this temporary new
       position */
    update_complex_tmp_com_crds(reactCom1, moleculeList);
    update_complex_tmp_com_crds(reactCom2, moleculeList);

    reflect_traj_tmp_crds(params, moleculeList, reactCom1, traj, membraneObject,
                          0.0, isInsideCompartment); // uses tmpCoords to calculate traj.
    reflect_traj_tmp_crds(params, moleculeList, reactCom2, traj, membraneObject,
                          0.0, isInsideCompartment);

    if (std::abs(traj[0] + traj[1] + traj[2]) > 1E-15) {
      // update the temporary coordinates for both complexes
      Vector vtraj{traj[0], traj[1], traj[2]};
      for (auto &mp : reactCom1.memberList)
        moleculeList[mp].update_association_coords(vtraj);
      for (auto &mp : reactCom2.memberList)
        moleculeList[mp].update_association_coords(vtraj);

      // if (iter == 818823) {
      //   std::cout << "CRDS after reflecting off of the BOX by " << traj[0]
      //             << ' ' << traj[1] << ' ' << traj[2] << std::endl;
      //   reactMol1.display_assoc_icoords("mol1");
      //   reactMol2.display_assoc_icoords("mol2");
      // }
    }
    /*Calculate the angles swept out by the interface to COM vectors as a result
     * of the displacement*/
    calc_angular_displacement(ifaceIndex1, ifaceIndex2, reactMol1, reactMol2,
                              reactCom1, reactCom2, moleculeList);

    /* CHECKS AFTER ASSOCIATION FOR STERIC COLLISIONS, FOR EXPANDING BEYOND THE
       BOX SIZE OR FOR MOVING PROTEINS A LARGE DISTANCE DUE TO SNAPPING INTO
       PLACE */
    bool cancelAssoc{false};
    check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2, moleculeList,
                                params, molTemplateList);
    if (cancelAssoc == false) {
      check_if_spans_box(cancelAssoc, params, reactCom1, reactCom2,
                         moleculeList, membraneObject);
      if (cancelAssoc == true)
        counterArrays.nCancelSpanBox++;
    } else
      counterArrays.nCancelOverlapPartner++; // true for structure overlap
                                             // check.
    if (cancelAssoc == false) {
      check_for_structure_overlap_system(cancelAssoc, reactCom1, reactCom2,
                                         moleculeList, params, molTemplateList,
                                         complexList, forwardRxns, backRxns);
      if (cancelAssoc == true)
        counterArrays.nCancelOverlapSystem++;
    }
    if (cancelAssoc == false) {
      measure_complex_displacement(cancelAssoc, reactCom1, reactCom2,
                                   moleculeList, params, molTemplateList,
                                   complexList);
      if (cancelAssoc == true) {
        if (isOnMembrane)
          counterArrays.nCancelDisplace2D++;
        else if (transitionToSurface)
          counterArrays.nCancelDisplace3Dto2D++;
        else
          counterArrays.nCancelDisplace3D++;
      }
    }

    if (cancelAssoc) {
      // std::cout << "Canceling association, returning complexes to original
      // state.\n";
      for (auto memMol : reactCom1.memberList)
        moleculeList[memMol].clear_tmp_association_coords();
      for (auto memMol : reactCom2.memberList)
        moleculeList[memMol].clear_tmp_association_coords();
      // end routine here!
      return;
    }
    counterArrays.nAssocSuccess++; // keep track of total number of successful
                                   // association moves.
    /*Keep track of the sizes of complexes that associated*/
    track_association_events(reactCom1, reactCom2, transitionToSurface, isOnMembrane, counterArrays);

    // else if cancelAssoc==false, write temporary to real coords and clear
    // temporary coordinates
    for (auto memMol : reactCom1.memberList) {
      moleculeList[memMol].comCoord = moleculeList[memMol].tmpComCoord;
      for (unsigned int i{0}; i < moleculeList[memMol].interfaceList.size(); ++i)
        moleculeList[memMol].interfaceList[i].coord = moleculeList[memMol].tmpICoords[i];
      moleculeList[memMol].clear_tmp_association_coords();
      moleculeList[memMol].trajStatus = TrajStatus::propagated;
      moleculeList[memMol].need_to_send = true;
    }

    for (auto memMol : reactCom2.memberList) {
      moleculeList[memMol].comCoord = moleculeList[memMol].tmpComCoord;
      for (unsigned int i{0}; i < moleculeList[memMol].interfaceList.size();
           ++i)
        moleculeList[memMol].interfaceList[i].coord = moleculeList[memMol].tmpICoords[i];
      moleculeList[memMol].clear_tmp_association_coords();
      if (currRxn.rxnType != ReactionType::biMolStateChange) {
        moleculeList[memMol].myComIndex = reactCom1.index; // update their complex index
        reactCom1.memberList.push_back(memMol);
      }
      moleculeList[memMol].trajStatus = TrajStatus::propagated;
      moleculeList[memMol].need_to_send = true;
    }

    // update complexes

    // add lifetime to the previous small clusters
    for (unsigned index = 0; index < molTemplateList.size(); index++) {
      if (molTemplateList[index].countTransition == true) {
        // com2 destroty
        if (numEachMolPrevious2[index] - 1 >= 0)
          molTemplateList[index].lifeTime[numEachMolPrevious2[index] - 1].emplace_back(
                  (iter - lastNumberUpdateItrEachMolPrevious2[index]) * Parameters::dt / 1E6);
      }
    }

    reactCom2.memberList.clear(); // clear the member list so the molecules don't get destroyed
    reactCom2.destroy(moleculeList, complexList); // destroy the complex
    reactCom1.update_properties(moleculeList, molTemplateList); // recalculate the properties of the first complex

    for (unsigned index = 0; index < molTemplateList.size(); index++) {
      if (molTemplateList[index].countTransition == true &&
          complexList[reactMol1.myComIndex].numEachMol[index] > numEachMolPrevious1[index]) {
        // grow
        if (numEachMolPrevious1[index] - 1 >= 0) {
          molTemplateList[index].lifeTime[numEachMolPrevious1[index] - 1].emplace_back(
                  (iter - lastNumberUpdateItrEachMolPrevious1[index]) * Parameters::dt / 1E6);
          complexList[reactMol1.myComIndex].lastNumberUpdateItrEachMol[index] = iter;
        }
      }
    }

    // update transition matrix
    // compare current cluster size with the previous ones
    for (unsigned index = 0; index < molTemplateList.size(); index++) {
      if (molTemplateList[index].countTransition == true &&
          complexList[reactMol1.myComIndex].numEachMol[index] >
              numEachMolPrevious1[index]) {
        // grow
        if (numEachMolPrevious1[index] - 1 >= 0)
          molTemplateList[index].transitionMatrix[numEachMolPrevious1[index] - 1][numEachMolPrevious1[index] - 1] +=
              iter - Parameters::lastUpdateTransition[index] - 1;
        if (numEachMolPrevious2[index] - 1 >= 0)
          molTemplateList[index].transitionMatrix[numEachMolPrevious2[index] - 1][numEachMolPrevious2[index] - 1] +=
              iter - Parameters::lastUpdateTransition[index] - 1;
        if (complexList[reactMol1.myComIndex].numEachMol[index] - 1 >= 0 &&
            numEachMolPrevious1[index] - 1 >= 0)
          molTemplateList[index].transitionMatrix
              [complexList[reactMol1.myComIndex].numEachMol[index] - 1]
              [numEachMolPrevious1[index] - 1] += 1;
        if (complexList[reactMol1.myComIndex].numEachMol[index] - 1 >= 0 &&
            numEachMolPrevious2[index] - 1 >= 0)
          molTemplateList[index].transitionMatrix
              [complexList[reactMol1.myComIndex].numEachMol[index] - 1]
              [numEachMolPrevious2[index] - 1] += 1;

        // update diagonal elements for unchanged complexes
        for (unsigned indexCom = 0; indexCom < complexList.size(); indexCom++) {
          if (indexCom != reactMol1.myComIndex) {
            if (complexList[indexCom].numEachMol[index] - 1 >= 0) {
              molTemplateList[index].transitionMatrix[complexList[indexCom].numEachMol[index] - 1]
                                   [complexList[indexCom].numEachMol[index] - 1] +=
                  iter - Parameters::lastUpdateTransition[index];
            }
          }
        }

        // update lastUpdateTransition
        Parameters::lastUpdateTransition[index] = iter;
      }
    }

    //-----------------------END UPDATE BINDPAIRLIST--------------------------

    // std::cout << " After ASSOCIATE, CHANGE COPY NUMBERS, interfaces: " << ifaceIndex1 << ' ' << ifaceIndex2
    //           << " add to product: " << currRxn.productListNew[0].absIfaceIndex
    //           << " sub from reactants: " << currRxn.reactantListNew[0].absIfaceIndex << " "
    //           << currRxn.reactantListNew[1].absIfaceIndex << std::endl;

    // Enforce boundary conditions
    reflect_complex_rad_rot(membraneObject, reactCom1, moleculeList, 0.0, isInsideCompartment);
  } // end of if these molecules are closing a loop or not.

  //------------------------START UPDATE MONOMERLIST-------------------------
  //TODO: why is it possible that the bndpartner is empty right after binding?
  // update oneTemp.monomerList when oneTemp.canDestroy is true and mol is
  // monomer reactMol1
  {
    Molecule &oneMol{reactMol1};
    MolTemplate &oneTemp{molTemplateList[oneMol.molTypeIndex]};
    bool isMonomer{oneMol.bndpartner.empty()};
    bool canDestroy{oneTemp.canDestroy};
    if (isMonomer && canDestroy) {
      // remove from monomerList
      std::vector<int> &oneList{oneTemp.monomerList};
      std::vector<int>::iterator result{
          std::find(std::begin(oneList), std::end(oneList), oneMol.index)};
      if (result != std::end(oneList)) {
        oneList.erase(result);
      }
    }
  }
  // reactMol2
  {
    Molecule &oneMol{reactMol2};
    MolTemplate &oneTemp{molTemplateList[oneMol.molTypeIndex]};
    bool isMonomer{oneMol.bndpartner.empty()};
    bool canDestroy{oneTemp.canDestroy};
    if (isMonomer && canDestroy) {
      // remove from monomerList
      std::vector<int> &oneList{oneTemp.monomerList};
      std::vector<int>::iterator result{
          std::find(std::begin(oneList), std::end(oneList), oneMol.index)};
      if (result != std::end(oneList)) {
        oneList.erase(result);
      }
    }
  }
  //------------------------END UPDATE MONOMERLIST---------------------------

  // update Molecule interface statuses
  reactMol1.interfaceList[ifaceIndex1].interaction.partnerIndex = reactMol2.index;
  reactMol2.interfaceList[ifaceIndex2].interaction.partnerIndex = reactMol1.index;
  reactMol1.interfaceList[ifaceIndex1].interaction.partnerIfaceIndex = ifaceIndex2;
  reactMol2.interfaceList[ifaceIndex2].interaction.partnerIfaceIndex = ifaceIndex1;
  if (currRxn.isReversible) {
    reactMol1.interfaceList[ifaceIndex1].interaction.conjBackRxn = currRxn.conjBackRxnIndex;
    reactMol2.interfaceList[ifaceIndex2].interaction.conjBackRxn = currRxn.conjBackRxnIndex;
  }
  if (assocDissocFile.is_open()) {
      assocDissocFile << "ITR:" << iter << "," << "BOND," 
      << molTemplateList[reactMol1.molTypeIndex].molName << "," << reactMol1.index << "," << ifaceIndex1 << "," 
      << molTemplateList[reactMol2.molTypeIndex].molName << "," << reactMol2.index << "," << ifaceIndex2 << std::endl;
  }
  reactMol1.interfaceList[ifaceIndex1].isBound = true;
  reactMol1.interfaceList[ifaceIndex1].index = currRxn.productListNew[0].absIfaceIndex;
  reactMol2.interfaceList[ifaceIndex2].isBound = true;
  reactMol2.interfaceList[ifaceIndex2].index = currRxn.productListNew[0].absIfaceIndex;

  // add to the list of bound interfaces and remove from the list of free
  // interfaces
  reactMol1.bndlist.push_back(ifaceIndex1);
  reactMol2.bndlist.push_back(ifaceIndex2);
  reactMol1.bndpartner.push_back(reactMol2.index);
  reactMol2.bndpartner.push_back(reactMol1.index);
  reactMol1.bndRxnList.push_back(currRxn.relRxnIndex);
  reactMol2.bndRxnList.push_back(currRxn.relRxnIndex);

  {
    size_t tmpItr{0};
    for (unsigned i{0}; i < reactMol1.freelist.size(); ++i) {
      if (reactMol1.freelist[i] == ifaceIndex1)
        tmpItr = i;
    }
    reactMol1.freelist[tmpItr] = reactMol1.freelist.back();
    reactMol1.freelist.pop_back();
  }
  {
    size_t tmpItr{0};
    for (unsigned i{0}; i < reactMol2.freelist.size(); ++i) {
      if (reactMol2.freelist[i] == ifaceIndex2)
        tmpItr = i;
    }
    reactMol2.freelist[tmpItr] = reactMol2.freelist.back();
    reactMol2.freelist.pop_back();
  }

  // Set probability of this protein to zero in all reactions so it doesn't try
  // to react again but the partners still will avoid overlapping.
  for (unsigned crossItr{0}; crossItr < reactMol1.crossbase.size();
       ++crossItr) {
    int skipMol{reactMol1.crossbase[crossItr]};
    for (unsigned crossItr2{0};
         crossItr2 < moleculeList[skipMol].crossbase.size(); ++crossItr2) {
      if (moleculeList[skipMol].crossbase[crossItr2] == reactMol1.index)
        moleculeList[skipMol].probvec[crossItr2] = 0;
    }
  }

  // Set probability of this protein to zero in all reactions so it doesn't try
  // to react again but the partners still will avoid overlapping.
  for (unsigned crossItr{0}; crossItr < reactMol2.crossbase.size(); ++crossItr) {
    int skipMol{reactMol2.crossbase[crossItr]};
    for (unsigned crossItr2{0};
         crossItr2 < moleculeList[skipMol].crossbase.size(); ++crossItr2) {
      if (moleculeList[skipMol].crossbase[crossItr2] == reactMol2.index)
        moleculeList[skipMol].probvec[crossItr2] = 0;
    }
  }

  // Update the crossed molecule lists so that the current molecules won't avoid
  // anything, but others will.
  reactCom1.ncross = -1;
  reactMol1.crossbase.clear();
  reactMol2.crossbase.clear();

  /*update the number of bound species*/
  update_Nboundpairs(reactMol1.molTypeIndex, reactMol2.molTypeIndex, 1, params, counterArrays);
  /*Update species copy numbers*/
  counterArrays.copyNumSpecies[currRxn.reactantListNew[0].absIfaceIndex]--; // decrement ifaceIndex1
  counterArrays.copyNumSpecies[currRxn.reactantListNew[1].absIfaceIndex]--; // decrement ifaceIndex2
  counterArrays.copyNumSpecies[currRxn.productListNew[0].absIfaceIndex]++; // increment product state

  //-------------------------UPDATE BINDPAIRLIST----------------------------
  if (counterArrays.canDissociate[currRxn.productListNew[0].absIfaceIndex]) {
    // int molIndex { reactMol1.index };
    // if (reactMol1.index > reactMol2.index)
    //     molIndex = reactMol2.index;
    int molIndex{-1};
    if (ifaceIndex1 == currRxn.reactantListNew[0].relIfaceIndex) {
      molIndex = reactMol1.index;
    } else {
      molIndex = reactMol2.index;
    }

    counterArrays.bindPairList[currRxn.productListNew[0].absIfaceIndex]
        .emplace_back(molIndex);
  }
  //-----------------------END UPDATE BINDPAIRLIST--------------------------

  // std::cout << " After ASSOCIATE, CHANGE COPY NUMBERS, interfaces: " <<
  // ifaceIndex1 << ' ' << ifaceIndex2
  //           << " add to product: " << currRxn.productListNew[0].absIfaceIndex
  //           << " sub from reactants: " <<
  //           currRxn.reactantListNew[0].absIfaceIndex << " "
  //           << currRxn.reactantListNew[1].absIfaceIndex << std::endl;

  // species tracking here
  if (currRxn.isObserved) {
    auto obsItr = observablesList.find(currRxn.observeLabel);
    if (obsItr != observablesList.end())
      ++obsItr->second;
  }
}
