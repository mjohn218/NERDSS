/*! \file class_Observable.cpp
 *
 * \brief Member functions for Observable class.
 *
 * ### Created on 2019-02-06 by Matthew Varga
 * Note: Currently not fully implemented or used
 */

#include "classes/class_Observable.hpp"

using namespace SpeciesTracker;

/* OBSERVABLE::CONSTITUENT::IFACE */
bool Observable::Constituent::Iface::operator==(const RxnIface& rxnIface) const
{
    return (relIndex == rxnIface.relIfaceIndex) && (isBound == rxnIface.relIfaceIndex);
}

/* OBSERVABLE::CONSTITUENT */

/* OBSERVABLE */
bool Observable::operator==(const Molecule& mol) const
{
    if (constituentList[0].molTypeIndex != mol.molTypeIndex)
        return false;
    for (auto& obsIface : constituentList[0].interfaceList) {
        if (mol.interfaceList[obsIface.relIndex].index != obsIface.absIndex) {
            return false;
        }
    }
    return true;
}
