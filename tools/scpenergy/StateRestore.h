/*  -*- c++ -*-  */
#ifndef STATERESTORE_H
#define STATERESTORE_H

#include "Real.h"
#include "Vector3DBlock.h"
#include "ScalarStructure.h"

#include <vector>
using std::vector;

namespace ProtoMol {

  class StateRestore {
  public:
    StateRestore();
    ~StateRestore();
    Vector3D myE1;
    Vector3D myE2;
    Vector3D myE3;
    Vector3D myOrigin;
    Vector3DBlock* myPositions;
    Vector3DBlock* myVelocities;
    vector<Vector3DBlock*> myForces;
    Real myTimestep;
    ScalarStructure* myEnergies;
    unsigned int mySystemsize;
    unsigned int myLevels;
    unsigned short myCRN[3];
  };
  
}




#endif
