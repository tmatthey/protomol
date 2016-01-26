#include "iSGregisterIntegratorExemplars.h"
#include "IntegratorFactory.h"
#include "iSGIntegrator.h"

#include "Vector.h"

namespace ProtoMol {

  void registerIntegratorExemplars() {
    IntegratorFactory::registerExemplar(new iSGIntegrator());
  }
  
}
