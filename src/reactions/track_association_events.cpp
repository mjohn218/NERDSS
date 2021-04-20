#include "reactions/association/association.hpp"
#include <cmath>
#include <iomanip>

void track_association_events(Complex& reactCom1, Complex& reactCom2, bool transitionToSurface, bool isOnMembrane, copyCounters& counterArrays)
{
  int size1=reactCom1.memberList.size();
  int size2=reactCom2.memberList.size();
  int smallerSize=size1;
  bool isDimerization = false;
  int index=0;
  int arrayIndex = 0;
  int loc;
  bool exists = false;
  
  //find size of smaller complex.
  if(size2 < size1)
    smallerSize = size2;
  if(size2 == 1 && size1 ==1 ){
    //this is a dimerization event.
    isDimerization = true;
    index = 0;//assign zeroth element to 1+1
  }else
    index=smallerSize;//all indices>0 are for addition of n=1, or n>1 . 
  
  /*Bundle values for n>maxSingles*/
  int spacing= 10;
  int maxSingles=10;//keep all n integers less than this value.
  int arraySize= counterArrays.eventArraySize;//final array element contains all additions above (arraySize-maxSingles)*spacing
  if(index>maxSingles-1){
    int a=index/spacing;//automatically rounds down
    arrayIndex=maxSingles-1+a;
    if(arrayIndex>=arraySize)
      arrayIndex = arraySize-1;//max index is arraySize-1
  }else{
    arrayIndex=index;
  }
  /*Keep track of all event types, and in which dimensionality they occur*/
  if(transitionToSurface ==true){
    //3D->2D

    counterArrays.events3Dto2D[arrayIndex]++;
  }else if(isOnMembrane == true){
    //2D
    counterArrays.events2D[arrayIndex]++;

  }else{
    //3D

    counterArrays.events3D[arrayIndex]++;
  }



}
