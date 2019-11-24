/*! \file check_bimolecular_reactions_par

 * Created on 2019-03-07 by Matthew Varga
 * Purpose:
 * Notes:
 */

#pragma once

#include "reactions/bimolecular_reactions.hpp"
#include "classes/class_SimulVolume.hpp"

#include <mutex>

void search_subvolumes(unsigned pro1Index, unsigned pro1Itr, int simItr, double *tableIDs, unsigned& DDTableIndex,
                       const Parameters& params, const SimulVolume::SubVolume& subVolume, std::vector<std::mutex>&
                           molLocks,
                       std::vector<gsl_matrix *>& normMatrices, std::vector<gsl_matrix *>& survMatrices,
                       std::vector<gsl_matrix *>& pirMatrices, std::vector<Molecule>& moleculeList,
                       std::vector<Complex>& complexList, const std::vector<MolTemplate>& molTemplateList,
                       const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns);

void check_bimolecular_reactions_par(int pro1Index, int pro2Index, int simItr, double *tableIDs, unsigned& DDTableIndex,
                                     const Parameters& params, std::vector<std::mutex>& molLocks,
                                     std::vector<gsl_matrix *>& normMatrices, std::vector<gsl_matrix *>& survMatrices,
                                     std::vector<gsl_matrix *>& pirMatrices, std::vector<Molecule>& moleculeList,
                                     std::vector<Complex>& complexList, const std::vector<MolTemplate>& molTemplateList,
                                     const std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns);

void determine_2D_bimolecular_reaction_probability_par(int simItr, unsigned& DDTableIndex, double *tableIDs,
                                                       BiMolData& biMolData, std::array<int, 2> rxnItr,
                                                       const Parameters& params, std::vector<Molecule>& moleculeList,
                                                       std::vector<Complex>& complexList, std::vector<std::mutex>&
                                                           molLocks,
                                                       const std::vector<ForwardRxn>& forwardRxns,
                                                       const std::vector<BackRxn>& backRxns,
                                                       std::vector<gsl_matrix *>& normMatrices,
                                                       std::vector<gsl_matrix *>& survMatrices,
                                                       std::vector<gsl_matrix *>& pirMatrices);

void determine_3D_bimolecular_reaction_probability_par(int simItr, unsigned& DDTableIndex, double *tableIDs,
                                                       BiMolData& biMolData, std::array<int, 2> rxnItr,
                                                       const Parameters& params, std::vector<Molecule>& moleculeList,
                                                       std::vector<Complex>& complexList, std::vector<std::mutex>&
                                                           molLocks,
                                                       const std::vector<ForwardRxn>& forwardRxns,
                                                       const std::vector<BackRxn>& backRxns,
                                                       std::vector<gsl_matrix *>& normMatrices,
                                                       std::vector<gsl_matrix *>& survMatrices,
                                                       std::vector<gsl_matrix *>& pirMatrices);
