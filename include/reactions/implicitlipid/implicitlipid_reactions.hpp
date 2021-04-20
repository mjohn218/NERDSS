/*! \file implicitlipid_reactions.hpp
 *
 * \brief
 *
 * ### Created on 2019-09-16 by Yiben Fu
 */
#pragma once

#include "classes/class_Membrane.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_SimulVolume.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"
#include <classes/class_copyCounters.hpp>
#include <gsl/gsl_matrix.h>

struct paramsIL {
    double R2D;
    double sigma;
    double Dtot;
    double ka;
    double kb;
    double area;
    double dt;
    int Na;
    int Nlipid;
};

/*!
 * \brief Gets the distance between one Molecule's Interface and the membrane surface, and determines if it is within Rmax, and can therefore
 * react.
 */
bool get_distance_to_surface(int pro1, int pro2, int iface1, int iface2, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    double& sep, double& R1, double Rmax, std::vector<Complex>& complexList, const ForwardRxn& currRxn,
    std::vector<Molecule>& moleculeList, const Membrane& membraneObject);

// binding probabilities for 3D and 2D
// h is the time-step dt; sigma is the bind_radius, rho is the lipid density on the surface.
double pimplicitlipid_3D(double z, paramsIL& parameters3D);
double pimplicitlipid_2D(paramsIL& parameters2D);

void determine_3D_implicitlipid_reaction_probability(int simItr, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    BiMolData& biMolData, const Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, Membrane& membraneObject, const int& relStateIndex);

void determine_2D_implicitlipid_reaction_probability(int simItr, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    std::vector<double>& ILTableIDs, BiMolData& biMolData, const Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, std::vector<double>& IL2DbindingVec, std::vector<double>& IL2DUnbindingVec, Membrane& membraneObject, const int& relStateIndex);
void check_dissociation_implicitlipid(unsigned int simItr, const Parameters& params, SimulVolume& simulVolume,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, unsigned int molItr,
    std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, const std::vector<BackRxn>& backRxns, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<CreateDestructRxn>& createDestructRxns, copyCounters& counterArrays, Membrane& membraneObject, std::vector<double>& IL2DbindingVec, std::vector<double>& IL2DUnbindingVec, std::vector<double>& ILTableIDs);

void break_interaction_implicitlipid(size_t relIface1, size_t relIface2, Molecule& reactMol1, Molecule& reactMol2,
    const BackRxn& currRxn, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, std::vector<MolTemplate>& molTemplateList);

// some functions that helps calculate the 2D binding probability.
double dissociate2D(paramsIL& parameters2D);
double integral_for_blockdistance2D(paramsIL& parameters2D);
double function2D(double u, void* parameter);
void block_distance(paramsIL& parameters2D);
double dissociate3D(double dt, double Dtot, double sigma, double ka, double kb);
