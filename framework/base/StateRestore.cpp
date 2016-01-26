#include "StateRestore.h"
#include "Vector3DBlock.h"
#include "ScalarStructure.h"

namespace ProtoMol {

StateRestore::StateRestore() {
  myPositions = new Vector3DBlock();
  myVelocities = new Vector3DBlock();
  myEnergies = new ScalarStructure();
  mySystemsize = 0;
  myTimestep = 0.0;
  myLevels = 0;
}


StateRestore::~StateRestore() {
  if (myPositions) delete myPositions;
  if (myVelocities) delete myVelocities;
  if (myEnergies) delete myEnergies;
}

}
