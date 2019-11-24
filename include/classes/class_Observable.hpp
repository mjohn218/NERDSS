/*! \file class_Observable.hpp
 * ### Created on 2019-02-06 by Matthew Varga
 * TODO: 02/06/2019 still being implemented
 */

/*! \defgroup SpeciesTracker
 * \brief Components of species tracking
 */

#pragma once

#include "classes/class_Rxns.hpp"

#include <string>
#include <vector>

namespace SpeciesTracker {

/*! \ingroup SpeciesTracker
 * \brief
 */
enum class ObservableType {
    none = 0, //!< only to initialize
    molecule = 1, //!< just one constituent Molecules, no bonds
    complex = 2, //!< one or two constituent Molecules, contains bonds
    reaction = 3, //!< how many times a reaction has occurred
    sum = 4, //!< perform an addition of observables. only performed once per timeWrite interval
};

/*! \ingroup SpeciesTracker
 * \brief Holds information on an Observable provided by the user in an observables block in the parameters file.
 *
 * TODO: function to check complex against Observable at points where molecules and complexes change
 */
struct Observable {
    /*!
     * \brief Holds information on a constituent of the Observable
     */
    struct Constituent {
        struct Iface {
            unsigned relIndex { 0 };
            unsigned absIndex { 0 };
            char state { '\0' };
            bool isBound { false };
            int bondIndex { -1 };

            // overloaded operators
            bool operator==(const RxnIface& rxnIface) const;

            Iface() = default;
            Iface(unsigned _relIndex, unsigned _absIndex, char _state, bool _isBound, int _bondIndex)
                : relIndex(_relIndex)
                , absIndex(_absIndex)
                , state(_state)
                , isBound(_isBound)
                , bondIndex(_bondIndex)
            {
            }
        };
        unsigned molTypeIndex {};
        std::vector<Iface> interfaceList {};
    };
    std::vector<Constituent> constituentList; //!< Molecules which compose the Observable
    unsigned currNum { 0 }; //!< current number of Observable species in the system
    std::string name; //!< name of the observable
    ObservableType observableType { ObservableType::none }; //!< type of observable
    std::array<int, 2> coupledRxn; //!< type (0 = forward, 1 = back, 2 = createdestruct) and index of reaction which
                                   //!< forms this Observable

    // functions

    // Overloaded operators
    /*!
     * \brief Checks the Molecule against the Observable
     *
     * Only used with Observables of type complex.
     */
    bool operator==(const Molecule& mol) const;

    Observable() = default;
};
}
