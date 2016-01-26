#include "MetaForce.h"
#include "ForceGroup.h"

using std::vector;
using std::string;
using namespace ProtoMol::Report;
namespace ProtoMol {

  //_________________________________________________________________ MetaForce

  void MetaForce::addToForceGroup(ForceGroup* forceGroup){
    forceGroup->addMetaForce(this);
  }
}
