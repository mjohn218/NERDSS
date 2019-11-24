#include "reactions/association/association.hpp"

void reverse_rotation(Coord& reactIface1, Molecule& reactMol1, Molecule& reactMol2, Complex& reactCom1, Complex& reactCom2, Quat rotQuatPos, Quat rotQuatNeg, std::vector<Molecule>& moleculeList)
{
    rotQuatPos = rotQuatPos.inverse();
    rotQuatNeg = rotQuatNeg.inverse();
    rotate(reactIface1, rotQuatPos, reactCom1, moleculeList);
    rotate(reactIface1, rotQuatNeg, reactCom2, moleculeList);
}
