#include "reactions/association/association.hpp"

void determine_rotation_angles(double targAngle, double currAngle, double& rotAngPos, double& rotAngNeg,
    const Complex& reactCom1, const Complex& reactCom2)
{
    double totDrz { reactCom1.Dr.z + reactCom2.Dr.z };
    double totDrx { reactCom1.Dr.x + reactCom2.Dr.x };
    double totDz { reactCom1.D.z + reactCom2.D.z };
    double tol=1e-11;
    if (totDrx < tol) {
      //no rotation, use translation in z (correct this below if Dz==0 for both molecules.
      rotAngPos = (targAngle - currAngle) * (reactCom1.D.z / totDz);
      rotAngNeg = -(targAngle - currAngle) * (reactCom2.D.z / totDz);
      
    } else {
        rotAngPos = (targAngle - currAngle) * (reactCom1.Dr.x / totDrx);
        rotAngNeg = -(targAngle - currAngle) * (reactCom2.Dr.x / totDrx);
    }
    // if the molecules are both on the membrane, use relative Dx components.
    if (reactCom1.D.z<tol && reactCom2.D.z<tol) {
      
      double totDx { reactCom1.D.x + reactCom2.D.x };
      rotAngPos = (targAngle - currAngle) * (reactCom1.D.x / totDx);
      rotAngNeg = -(targAngle - currAngle) * (reactCom2.D.x / totDx);
    } 
       
}
