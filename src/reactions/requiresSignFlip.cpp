#include "reactions/association/association.hpp"

bool requiresSignFlip(Vector axis, Vector v1, Vector v2)
{
    Vector zAxis { 0, 0, 1 };
    Vector xAxis { 1, 0, 0};
    zAxis.magnitude = 1.0;
    xAxis.magnitude = 1.0;
    Vector u { zAxis.cross(axis) };
    u.calc_magnitude();
    double theta { zAxis.dot_theta(axis) };
    double useXAxis { false };
    if (std::abs(u.x) < 1E-8 && std::abs(u.y) < 1E-8 && std::abs(u.z) < 1E-8) {
        u  = xAxis.cross(axis);
        u.calc_magnitude();
        theta = xAxis.dot_theta(axis);
        useXAxis = true;
    }

    // rotate
    Quat rot(cos(theta / 2), sin(theta / 2) * u.x, sin(theta / 2) * u.y, sin(theta / 2) * u.z);
    rot = rot.unit();
    rot.rotate(v1);
    rot.rotate(v2);
    rot.rotate(axis);

    // if we rotated wrong way, reverse and rotate the other way
    // TODO: check this 0.01 business
    if ((zAxis.dot_theta(axis) > 0.01 && !useXAxis) || (useXAxis && xAxis.dot_theta(axis) < 0.01)) {
        rot = rot.inverse();
        rot.rotate(v1);
        rot.rotate(v2);
        rot = Quat(cos(-theta / 2), sin(-theta / 2) * u.x, sin(-theta / 2) * u.y, sin(-theta / 2) * u.z);
        rot.rotate(v1);
        rot.rotate(v2);
    }

    Vector projectedVec1;
    Vector projectedVec2;
    if (!useXAxis) {
        projectedVec1 = Vector{v1.x, v1.y, 0};
        projectedVec2 = Vector{v2.x, v2.y, 0};
        return projectedVec1.cross(projectedVec2).z > 0;
    } else {
        projectedVec1 = Vector(0, v1.y, v1.z);
        projectedVec2 = Vector(0, v2.y, v2.z);
        return projectedVec1.cross(projectedVec2).x > 0;
    }

    // if the angle sign isn't correct, return true
//    return !angleSignIsCorrect(projectedVec1, projectedVec2);
}
