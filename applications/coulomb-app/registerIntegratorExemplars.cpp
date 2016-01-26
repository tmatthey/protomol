#include "registerIntegratorExemplars.h"
#include "IntegratorFactory.h"
#include "PaulTrapIntegrator.h"
#include "Vector.h"

namespace ProtoMol {

  void registerIntegratorExemplars(){
    IntegratorFactory::registerExemplar(new PaulTrapIntegrator(),Vector<std::string>("NoseNVTLeapfrog"));
  }
}
