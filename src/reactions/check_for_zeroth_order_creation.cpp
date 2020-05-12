#include "math/constants.hpp"
#include "math/math_functions.hpp"
#include "math/rand_gsl.hpp"
#include "reactions/unimolecular/unimolecular_reactions.hpp"
#include "tracing.hpp"

void check_for_zeroth_order_creation(unsigned simItr, Parameters& params, std::vector<int>& emptyMolList,
    std::vector<int>& emptyComList, SimulVolume& simulVolume, const std::vector<ForwardRxn>& forwardRxns,
    const std::vector<CreateDestructRxn>& createDestructRxns, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, const std::vector<MolTemplate>& molTemplateList,
    std::map<std::string, int>& observablesList, copyCounters& counterArrays, Membrane& membraneObject)
{
    // TRACE();
    for (auto& oneRxn : createDestructRxns) {
        if (oneRxn.rxnType == ReactionType::zerothOrderCreation) {
            const MolTemplate& oneTemp = molTemplateList[oneRxn.productMolList.at(0).molTypeIndex];
            if (oneTemp.isImplicitLipid == true) {
                int indexIlState = oneRxn.productMolList.at(0).interfaceList.at(0).absIfaceIndex;
                // Calculate the probability
                // Poisson process, p(zero events) = exp(-k * V * deltaT), where k = rate (M/s), V = volume (L), deltaT=
                // timestep (s)
                // long double lambda { Constants::avogadro * oneRxn.rateList[0].rate * membraneObject.waterBox.volume
                //     * Constants::nm3ToLiters * params.timeStep * Constants::usToSeconds };
                double Volume = 0.0;
                // Change Volume depending upon box or sphere simulation
                if (membraneObject.isSphere) {
                    Volume = membraneObject.sphereVol;
                } else {
                    Volume = membraneObject.waterBox.volume;
                }
                long double lambda { Constants::avogadro * oneRxn.rateList[0].rate * Volume * Constants::nm3ToLiters * params.timeStep * Constants::usToSeconds };
                long double prob { exp(-lambda) };
                // std::cout << "create prob: " << std::setprecision(20) << prob << std::endl;
                // exit(1);
                unsigned numEvents { 0 };
                // double rNum { 1.0 * rand_gsl() };
                // double rNum2 { rNum + rand_gsl() * Constants::iRandMax }; // to get higher resolution

                while (rand_gsl() > prob) {
                    ++numEvents;
                    prob += (exp(-lambda) * pow(lambda, numEvents)) / MathFuncs::gammFactorial(numEvents);
                }

                if (numEvents > 0) {
                    std::cout << "Creating " << numEvents << " molecule(s) of type " << oneTemp.molName << " at iteration "
                              << simItr << '\n';
                }

                if (oneRxn.isObserved) {
                    auto observeItr = observablesList.find(oneRxn.observeLabel);
                    if (observeItr == observablesList.end()) {
                        std::cerr << "WARNING: Observable " << oneRxn.observeLabel << " not defined.\n";
                    } else {
                        observeItr->second += numEvents;
                    }
                }

                while (numEvents > 0) {
                    ++counterArrays.copyNumSpecies[indexIlState];
                    ++membraneObject.numberOfFreeLipidsEachState[indexIlState];
                    --numEvents;
                }
            } else {

                // Calculate the probability
                // Poisson process, p(zero events) = exp(-k * V * deltaT), where k = rate (M/s), V = volume (L), deltaT=
                // timestep (s)
                long double lambda { Constants::avogadro * oneRxn.rateList[0].rate * membraneObject.waterBox.volume
                    * Constants::nm3ToLiters * params.timeStep * Constants::usToSeconds };
                long double prob { exp(-lambda) };
                unsigned numEvents { 0 };
                // double rNum { 1.0 * rand_gsl() };
                // double rNum2 { rNum + rand_gsl() * Constants::iRandMax }; // to get higher resolution

                // std::cout << "create prob: " << std::setprecision(20) << prob << std::endl;
                // exit(1);

                while (rand_gsl() > prob) {
                    // std::cout << "create prob: " << std::setprecision(20) << prob << std::endl;
                    // std::cout << "rNum: " << std::setprecision(20) << rNum << std::endl;
                    // std::cout << "rNum2: " << std::setprecision(20) << rNum2 << std::endl;
                    ++numEvents;
                    prob += (exp(-lambda) * pow(lambda, numEvents)) / MathFuncs::gammFactorial(numEvents);
                }

                if (numEvents > 0) {
                    std::cout << "Creating " << numEvents << " molecule(s) of type " << oneTemp.molName << " at iteration "
                              << simItr << '\n';
                }

                if (oneRxn.isObserved) {
                    auto observeItr = observablesList.find(oneRxn.observeLabel);
                    if (observeItr == observablesList.end()) {
                        std::cerr << "WARNING: Observable " << oneRxn.observeLabel << " not defined.\n";
                    } else {
                        observeItr->second += numEvents;
                    }
                }

                while (numEvents > 0) {
                    // these aren't used here, just for the function
                    int newMolIndex { 0 };
                    int newComIndex { 0 };

                    create_molecule_and_complex_from_rxn(0, newMolIndex, newComIndex, false, oneTemp, params, emptyMolList,
                        emptyComList, oneRxn, simulVolume, moleculeList, complexList, molTemplateList, forwardRxns, membraneObject);

                    // update the copy number arrays
                    for (const auto& iface : moleculeList[newMolIndex].interfaceList)
                        ++counterArrays.copyNumSpecies[iface.index];

                    --numEvents;
                }
            }
        }
    }
}
