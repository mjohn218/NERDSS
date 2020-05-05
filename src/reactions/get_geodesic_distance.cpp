#include "classes/class_Coord.hpp"
#include "reactions/association/association.hpp"
#include <cmath>

double get_geodesic_distance(Coord intFace1, Coord intFace2)
{
    double r1 = intFace1.get_magnitude();
    double r2 = intFace2.get_magnitude();

    double dotProduct = intFace1.x * intFace2.x + intFace1.y * intFace2.y + intFace1.z * intFace2.z;

    double theta = std::acos(dotProduct / (r1 * r2));
    double meanR = (r1 + r2) / 2.0;

    return meanR * theta;
}