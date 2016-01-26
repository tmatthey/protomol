#include "MollyForce.h"
#include "ForceGroup.h"

using std::vector;
using std::string;
using namespace ProtoMol::Report;
namespace ProtoMol {

  //_________________________________________________________________ MollyForce

  void MollyForce::addToForceGroup(ForceGroup* forceGroup){
    forceGroup->addMollyForce(this);
  }
}
