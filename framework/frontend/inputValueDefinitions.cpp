#include "inputValueDefinitions.h"
#include "Vector.h"
using std::string;
using std::vector;
namespace ProtoMol {

  //________________________________________________________ InputValue<>
  defineInputValue(InputOutput,"output");
  defineInputValueAndText(InputDebug,"debug","report level, suppresses all output with higher output level");
  defineInputValue(InputTemperature,"temperature");
  defineInputValue(InputSeed,"seed");
  defineInputValue(InputConfig,"config");
  defineInputValue(InputFirststep,"firststep");
  defineInputValue(InputNumsteps,"numsteps");
  defineInputValue(InputOutputfreq,"outputfreq");
  defineInputValue(InputIntegrator,"integrator");
  defineInputValueWithAliases(InputPositions,"posfile",("coords")("coordinates"));
  defineInputValue(InputVelocities,"velfile");
  defineInputValueWithAliases(InputPSF,"psffile",("structure"));
  defineInputValueWithAliases(InputPAR,"parfile",("parameters"));
  defineInputValueWithAliasesAndText(InputRemoveLinearMomentum,
				     "removeLinearMomentum",
				     ("comMotion"),
				     "removes linear momentum, where -1 for never, 0 at initialization or at STS frequency <n>");
  defineInputValueWithAliasesAndText(InputRemoveAngularMomentum,
				     "removeAngularMomentum",
				     ("angularMomentum"),
				     "removes angular momentum, where -1 for never, 0 at initialization or at STS frequency <n>");
  defineInputValue(InputAMBER,"parfile");
  defineInputValue(InputCoords,"crdfile");
  defineInputValue(InputUseBarrier,"useBarrier");
  defineInputValue(InputParallelPipe,"parallelPipe");
  defineInputValue(InputParallelMode,"parallelMode");
  defineInputValue(InputMaxPackages,"maxPackages");
  defineInputValue(InputPDBScaling,"pdbScaling");
  defineInputValue(InputBoundaryConditions,"boundaryConditions");
  defineInputValue(InputCellManager,"cellManager");
  defineInputValue(InputDihedralMultPSF,"dihedralMultPSF");
  defineInputValueAndText( InputVirialCalc,"virialCalc",
                           "Required for constant pressure simulations." );
  defineInputValueAndText( InputMolVirialCalc,"molVirialCalc",
                           "Required for constant pressure simulations." );
  defineInputValue(InputShake,"shake");
  defineInputValue(InputShakeEpsilon,"shakeEpsilon");
  defineInputValue(InputShakeMaxIter,"shakeMaxIter");
  defineInputValue(InputRattle,"rattle");
  defineInputValue(InputRattleEpsilon,"rattleEpsilon");
  defineInputValue(InputRattleMaxIter,"rattleMaxIter");
  defineInputValue(InputReducedImage,"reducedImage");  
  defineInputValueAndText(InputMinimalImage,"minimalImage","global default flag whether the coordinates should be transformed to minimal image or not");  
  defineInputValueAndText( InputShadow, "shadowEnergy",
                           "Calculate shadow Hamiltonian" );
  defineInputValueAndText( InputShadowFreq, "shadowFreq",
                           "Frequency (in # of steps) to calculate shadow" );
  defineInputValueAndText( InputShadowOrder, "shadowOrder",
                           "Order (2k) of the shadow Hamiltonian" );
  defineInputValue(InputEigenVectors,"eigfile");
  defineInputValue(InputEigenValues,"eigvalfile");
  defineInputValue(InputCheckpointFreq, "checkpointfreq");
  defineInputValue(InputCheckpointFile, "checkpointfile");
  defineInputValue(InputRestoreState, "restorestate");
  defineInputValue(InputDoSCPISM, "doscpism");
  defineInputValue(InputEigTextFile, "eigtextfile");
}
