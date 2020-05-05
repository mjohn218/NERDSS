/*
 * ### Created on 2/05/2020 by Yiben Fu
 * ### Purpose
 * ***
 * all functions that are used in spherical systems
 */

#include "classes/class_Membrane.hpp"
#include "classes/class_MolTemplate.hpp"
#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Vector.hpp"

#include <array>
#include <cmath>

double radius(Coord mol);

Coord find_spherical_coords(Coord mol); // mol: cardesian coords, output spherical coords

Coord find_cardesian_coords(Coord mol); // mol: spherical coords, output cardesian coords

double theta_plus(double theta1, double theta2); // sum of two theta

double phi_plus(double phi1, double phi2); // sum of two phi

Coord angle_plus(Coord angle1, Coord angle2); // two spherical coords sum

Coord find_position_after_association(double alpha1, Coord Iface1, Coord Iface2, double alpha_total, double bindRadius); // on sphere, when association, the new position of Iface1. alpha1 is the geodesic angle that Iface1 moves.

//Coord translate_on_sphere(Coord targ, Coord COM, Coord REF, Coord COMnew, Coord REFnew); // translation on sphere, output cardesian coords.
std::array<double, 9> inner_coord_set(Coord com, Coord comnew);
std::array<double, 9> inner_coord_set_new(Coord com, Coord comnew);
std::array<double, 3> calculate_inner_coord_coefficients(Coord TARG, Coord COM, std::array<double, 9> crdset);
Coord translate_on_sphere(Coord targ, Coord COM, Coord COMnew, std::array<double, 9> crdset, std::array<double, 9> crdsetnew);
//Coord rotate_on_sphere(Coord targ, Coord COM, double dangle);
Coord rotate_on_sphere(Coord Targ, Coord COM, std::array<double, 9> crdset, double dangle);

double calc_bindRadius2D(double bindRadius, Coord iFace);

void set_memProtein_sphere(Complex reactCom, Molecule& memProtein, std::vector<Molecule> moleculeList, const Membrane membraneObject);
void find_Lipid_sphere(Complex reactCom, Molecule& Lipid, std::vector<Molecule> moleculeList, const Membrane membraneObject);