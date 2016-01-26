#include "BBKIntegrator.h"
#include "BSplineMOLLYIntegrator.h"
#include "CGMinimizerIntegrator.h"
#include "DLMCIntegrator.h"
#include "DMDLeapfrogIntegrator.h"
#include "DihedralHMCIntegrator.h"
#include "EquilibriumMOLLYIntegrator.h"
#include "HMCIntegrator.h"
#include "HessianInt.h"
#include "ImpulseIntegrator.h"
#include "IntegratorFactory.h"
#include "LangevinImpulseIntegrator.h"
#include "LeapfrogIntegrator.h"
#include "LeapfrogTruncatedShadow.h"
#include "NPTVerletIntegrator.h"
#include "NVTVerletIntegrator.h"
#include "NormModeInt.h"
#include "NormModeMin.h"
#include "NormModeSmplMin.h"
#include "NormModeDiag.h"
#include "NormModeVisual.h"
#include "NoseNVTLeapfrogIntegrator.h"
#include "RMTIntegrator.h"
#include "BerendsenIntegrator.h"
#include "PLeapfrogIntegrator.h"
#include "PaulTrapIntegrator.h"
#include "ShadowHMCIntegrator.h"
#include "S2HMCIntegrator.h"
#include "NormalModeLangevin.h"
#include "NormalModeMinimizer.h"
#include "NormalModeDiagonalize.h"
#include "NormalModeMori.h"
#include "NormalModeRelax.h"
#include "NormalModeBrownian.h"
#include "NumericalDifferentiation.h"

#include "registerIntegratorExemplars.h"
#include "Vector.h"

namespace ProtoMol {

  void registerIntegratorExemplars(){
    IntegratorFactory::registerExemplar(new BBKIntegrator());
    IntegratorFactory::registerExemplar(new BSplineMOLLYIntegrator(),Vector<std::string>("HBondMOLLY"));
    IntegratorFactory::registerExemplar(new CGMinimizer());
    IntegratorFactory::registerExemplar(new DLMCIntegrator());
    IntegratorFactory::registerExemplar(new DMDLeapfrogIntegrator());
    IntegratorFactory::registerExemplar(new DihedralHMCIntegrator());
    IntegratorFactory::registerExemplar(new EquilibriumMOLLYIntegrator());
    IntegratorFactory::registerExemplar(new HMCIntegrator());
    IntegratorFactory::registerExemplar(new HessianInt());
    IntegratorFactory::registerExemplar(new ImpulseIntegrator());
    IntegratorFactory::registerExemplar(new LangevinImpulseIntegrator());
    IntegratorFactory::registerExemplar(new LeapfrogIntegrator());
    IntegratorFactory::registerExemplar(new LeapfrogTruncatedShadow());
    IntegratorFactory::registerExemplar(new NPTVerletIntegrator());
    IntegratorFactory::registerExemplar(new NormModeInt());
    IntegratorFactory::registerExemplar(new NormModeMin());
    IntegratorFactory::registerExemplar(new NormModeSmplMin());
    IntegratorFactory::registerExemplar(new NormModeDiag());
    IntegratorFactory::registerExemplar(new NormModeVisual());
    IntegratorFactory::registerExemplar(new NoseNVTLeapfrogIntegrator());
    IntegratorFactory::registerExemplar(new RMTIntegrator());
    IntegratorFactory::registerExemplar(new BerendsenIntegrator());
    IntegratorFactory::registerExemplar(new PLeapfrogIntegrator());
    IntegratorFactory::registerExemplar(new PaulTrapIntegrator());
    IntegratorFactory::registerExemplar(new S2HMCIntegrator());
    IntegratorFactory::registerExemplar(new ShadowHMCIntegrator()); 
    IntegratorFactory::registerExemplar(new NormalModeLangevin());
    IntegratorFactory::registerExemplar(new NormalModeMinimizer());
    IntegratorFactory::registerExemplar(new NormalModeDiagonalize());
    IntegratorFactory::registerExemplar(new NormalModeMori());
    IntegratorFactory::registerExemplar(new NormalModeRelax());
    IntegratorFactory::registerExemplar(new NormalModeBrownian());
    IntegratorFactory::registerExemplar(new NumericalDifferentiation());

  }
}

