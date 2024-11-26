#include "math/Faddeeva.hpp"
#include "reactions/bimolecular/bimolecular_reactions.hpp"

/**
 * @brief Association probability from Smoluchovski's reaction diffusion model
 * 
 * @param r0 original distance
 * @param tCurr current time
 * @param Dtot total diffusion constant
 * @param bindRadius sigma
 * @param alpha 
 * @param cof 
 * @return double 
 */
double passocF_1D(double r0, double tCurr, double Dtot, double bindRadius, double ka)
{
  if (Dtot==0){
    if (r0 - bindRadius > 1e-6) {
      return 0.0;
    }
    else {
      return 1.0;
    }
  }
  const double sqrtfDt{sqrt(4.0 * Dtot * tCurr)};
  const double sq_scale_t{ka * sqrt(tCurr / Dtot)};
  const double sep{(r0 - bindRadius) / sqrtfDt}; 

  const double efc1{erfc(sep)};
  const double efcx1{Faddeeva::erfcx(sep + sq_scale_t)};
  const double ep1 = exp(-sep * sep);

  return efc1 - ep1 * efcx1;
}
