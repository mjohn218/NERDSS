/*! \file calculate_update_position_interface.cpp
 * ### Created on 02/28/2020 by Yiben Fu
 * ### Purpose: to calculate the new position after the translation and rotation on the sphere surface
 * ***
 *
 * ### Notes
 * ***
 *
 * ### TODO List
 * ***
 */
#include "reactions/association/functions_for_spherical_system.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

Coord calculate_update_position_interface(const Complex targCom, const Coord ifacecrds) // iface is cardesian coords
{
    Coord finalcrds; // for output
    Coord trajTrans = targCom.trajTrans;
    double dangle = targCom.trajRot.x;
    if (trajTrans.get_magnitude() < 1E-14) {
        finalcrds = ifacecrds;
        return finalcrds;
    }
    /*
    // define reference points
        Coord COM = targCom.comCoord;
        Coord COMnew = targCom.comCoord + targCom.trajTrans;
	     Coord REF = COMnew;
	     Coord dist = COMnew - COM;
	     double l = dist.get_magnitude();
	     double R = COM.get_magnitude();
	     double lnew = l + l*R*R/(R*R - l*l);
	     Coord dREFnew = (lnew/l) * dist;
	     Coord REFnew = COM + dREFnew;
	     REFnew = (R/REFnew.get_magnitude()) * REFnew;
	
	     // trans into spherical coords
	     COM = find_spherical_coords(COM);
	     COMnew = find_spherical_coords(COMnew);
	     REF = find_spherical_coords(REF);
	     REFnew = find_spherical_coords(REFnew);
    // get the rotation angle: dangle
    double dangle = targCom.trajRot.x; // we select the rotation angle x as the angle that the complex rotate on the sphere
  
    Coord temp = find_spherical_coords( ifacecrds );
    Coord targ = translate_on_sphere(temp, COM, REF, COMnew, REFnew);
    Coord finalcrds1 = rotate_on_sphere(targ, COMnew, dangle);
    finalcrds = find_cardesian_coords(finalcrds1);
    */
    Coord COM = targCom.comCoord;
    Coord COMnew = targCom.comCoord + targCom.trajTrans;
    std::array<double, 9> Crdset = inner_coord_set(COM, COMnew);
    std::array<double, 9> Crdsetnew = inner_coord_set_new(COM, COMnew);
    finalcrds = translate_on_sphere(ifacecrds, COM, COMnew, Crdset, Crdsetnew);
    finalcrds = rotate_on_sphere(finalcrds, COMnew, Crdsetnew, dangle);
    return finalcrds;
}