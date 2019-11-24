#include "reactions/association/association.hpp"

bool areSameAngle(double ang1, double ang2)
{
//    return std::abs(ang1) == std::abs(ang2);
    return std::abs(ang1 - ang2) < 1E-8;
}
