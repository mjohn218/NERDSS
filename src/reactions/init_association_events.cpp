#include "boundary_conditions/reflect_functions.hpp"
#include "classes/class_Rxns.hpp"
#include "io/io.hpp"
#include "reactions/association/association.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"
#include "reactions/shared_reaction_functions.hpp"
#include "tracing.hpp"
#include <cmath>
#include <iomanip>

void init_association_events( copyCounters& counterArrays)
{
  
  /*Bundle values for n>maxSingles*/
  int spacing= 10;
  int maxSingles=10;//keep all n integers less than this value.
  int arraySize= counterArrays.eventArraySize;//final array element contains all additions above (arraySize-maxSingles)*spacing
  for(int i=0;i<arraySize;i++){
    counterArrays.events3D.push_back(0);//initialize these histograms with 0. 
    counterArrays.events3Dto2D.push_back(0);
    counterArrays.events2D.push_back(0);
  }

}
