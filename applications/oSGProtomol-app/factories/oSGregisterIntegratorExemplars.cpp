#include "oSGregisterIntegratorExemplars.h"
#include "IntegratorFactory.h"
#include "muVTIntegrator.h"
#include "NfPTIntegrator.h"

#include "Vector.h"

namespace ProtoMol {

  void registerIntegratorExemplars() {
    IntegratorFactory::registerExemplar(new muVTIntegrator());
    IntegratorFactory::registerExemplar(new NfPTIntegrator());
  }

}
