#include "reactions/association/association.hpp"

Vector create_arbitrary_vector(Vector& vec)
{
    Vector x_axis(1, 0, 0);
    Vector y_axis(0, 1, 0);
    vec.normalize();

    return (vec.dot_theta(x_axis) != 0 && (vec.dot_theta(x_axis) != M_PI)) ? Vector(vec).cross(x_axis)
                                                                           : Vector(vec).cross(y_axis);
}
