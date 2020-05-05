/*! \file bimolecular_reactions.hpp
 *
 * \brief
 *
 * ### Created on 2019-02-12 by Matthew Varga
 */
#pragma once

#include "classes/class_Rxns.hpp"
#include <classes/class_copyCounters.hpp>
#include <gsl/gsl_matrix.h>

struct BiMolData {
    int pro1Index { 0 };
    int pro2Index { 0 };
    int com1Index { 0 };
    int com2Index { 0 };
    int relIface1 { 0 };
    int relIface2 { 0 };
    int absIface1 { 0 };
    int absIface2 { 0 };
    double Dtot { 0 };
    double magMol1 { 0 };
    double magMol2 { 0 };

    BiMolData() = default;
    BiMolData(int pro1Index, int pro2Index, int com1Index, int com2Index, int relIface1, int relIface2, int absIface1,
        int absIface2, double Dtot, double magMol1, double magMol2)
        : pro1Index(pro1Index)
        , pro2Index(pro2Index)
        , com1Index(com1Index)
        , com2Index(com2Index)
        , relIface1(relIface1)
        , relIface2(relIface2)
        , absIface1(absIface1)
        , absIface2(absIface2)
        , Dtot(Dtot)
        , magMol1(magMol1)
        , magMol2(magMol2)
    {
    }
};

/*!
 * \brief Gets the distance between two Molecule's Interfaces and determines if they are within Rmax, and can therefore
 * react.
 */
bool get_distance(int pro1, int pro2, int iface1, int iface2, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    double& sep, double& R1, double Rmax, std::vector<Complex>& complexList, const ForwardRxn& currRxn,
    std::vector<Molecule>& moleculeList, Membrane& membraneObject);

double passocF(double r0, double tCurr, double Dtot, double bindRadius, double alpha, double cof);

void determine_2D_bimolecular_reaction_probability(int simItr, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    unsigned& DDTableIndex, double* tableIDs, BiMolData& biMolData, const Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, Membrane& membraneObject, std::vector<gsl_matrix*>& normMatrices,
    std::vector<gsl_matrix*>& survMatrices, std::vector<gsl_matrix*>& pirMatrices);

void determine_3D_bimolecular_reaction_probability(int simItr, int rxnIndex, int rateIndex, bool isStateChangeBackRxn,
    unsigned& DDTableIndex, double* tableIDs, BiMolData& biMolData, const Parameters& params,
    std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<BackRxn>& backRxns, Membrane& membraneObject, std::vector<gsl_matrix*>& normMatrices,
    std::vector<gsl_matrix*>& survMatrices, std::vector<gsl_matrix*>& pirMatrices);

void perform_bimolecular_state_change(int stateChangeIface, int facilitatorIface, std::array<int, 3>& rxnItr,
    Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom, Complex& facilitatorCom,
    copyCounters& counterArrays, const Parameters& params, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, Membrane& membraneObject);
void perform_bimolecular_state_change_box(int stateChangeIface, int facilitatorIface, std::array<int, 3>& rxnItr,
    Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom, Complex& facilitatorCom,
    copyCounters& counterArrays, const Parameters& params, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, Membrane& membraneObject);
void perform_bimolecular_state_change_sphere(int stateChangeIface, int facilitatorIface, std::array<int, 3>& rxnItr,
    Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom, Complex& facilitatorCom,
    copyCounters& counterArrays, const Parameters& params, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, Membrane& membraneObject);

void perform_implicitlipid_state_change(int stateChangeIface, int facilitatorIface, std::array<int, 3>& rxnItr,
    Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom, Complex& facilitatorCom,
    copyCounters& counterArrays, const Parameters& params, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, Membrane& membraneObject);
void perform_implicitlipid_state_change_box(int stateChangeIface, int facilitatorIface, std::array<int, 3>& rxnItr,
    Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom, Complex& facilitatorCom,
    copyCounters& counterArrays, const Parameters& params, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, Membrane& membraneObject);
void perform_implicitlipid_state_change_sphere(int stateChangeIface, int facilitatorIface, std::array<int, 3>& rxnItr,
    Molecule& stateChangeMol, Molecule& facilitatorMol, Complex& stateChangeCom, Complex& facilitatorCom,
    copyCounters& counterArrays, const Parameters& params, std::vector<ForwardRxn>& forwardRxns,
    std::vector<BackRxn>& backRxns, std::vector<Molecule>& moleculeList, std::vector<Complex>& complexList,
    std::vector<MolTemplate>& molTemplateList, std::map<std::string, int>& observablesList, Membrane& membraneObject);