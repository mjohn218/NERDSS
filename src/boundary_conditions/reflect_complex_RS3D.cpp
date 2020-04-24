/*! \file reflect_complex_targCom.radius_rot.cpp

 * ### Created on 11/8/18 by Matthew Varga
 * ### Purpose
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */
#include "boundary_conditions/reflect_functions.hpp"

void reflect_complex_RS3D(double RS3Dinput, Complex& targCom, std::vector<Molecule>& moleculeList)
{
    Vector vec {};
    vec.x = 0;
    vec.y = 0;
    vec.z = RS3Dinput;
    targCom.comCoord += vec;
    for (auto& mol : targCom.memberList) {
        moleculeList[mol].comCoord += vec;
        for (auto& iface : moleculeList[mol].interfaceList) {
            iface.coord += vec;
        }
    }
    /*
    double currx{ targCom.comCoord.x + targCom.trajTrans.x };
    double curry{ targCom.comCoord.y + targCom.trajTrans.y };
    double currz{ targCom.comCoord.z + targCom.trajTrans.z }; 
    bool propagate = true;
    if (currx > membraneObject.waterBox.x/2.0 || currx < -membraneObject.waterBox.x/2.0) 
        propagate = false; 
    if (curry > membraneObject.waterBox.y/2.0 || curry < -membraneObject.waterBox.y/2.0) 
        propagate = false; 
    if (currz > membraneObject.waterBox.z/2.0 || currz < -membraneObject.waterBox.z/2.0) 
        propagate = false; 
    if (propagate = false)
    {   targCom.trajTrans.x =0;
        targCom.trajTrans.y =0;
        targCom.trajTrans.z =0;
    }
    */
}
