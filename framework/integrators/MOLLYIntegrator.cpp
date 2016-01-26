#include "MOLLYIntegrator.h"
#include "ModifierAveraging.h"
#include "ModifierMollification.h"
#include "Vector3DBlock.h"
#include <algorithm>
#include "Report.h"
using namespace ProtoMol::Report;

namespace ProtoMol {


  //_________________________________________________________ MOLLYIntegrator
  MOLLYIntegrator::MOLLYIntegrator():MTSIntegrator(),mySwapPositions(NULL){}

  MOLLYIntegrator::MOLLYIntegrator (int cycles,
				    ForceGroup *overloadedForces,
				    StandardIntegrator *nextIntegrator)
    : MTSIntegrator(cycles,overloadedForces,nextIntegrator),mySwapPositions(NULL){}

  MOLLYIntegrator::~MOLLYIntegrator(){}

  void MOLLYIntegrator::initialize(GenericTopology *topo,
				   Vector3DBlock *positions,
				   Vector3DBlock *velocities,
				   ScalarStructure *energies){
    MTSIntegrator::initialize(topo,positions,velocities,energies);

  }

  void MOLLYIntegrator::addModifierBeforeInitialize(){
    Report::report <<"add ModifierAveraging"<<Report::endr;
    Report::report <<"add ModifierMollification"<<Report::endr;
    adoptPreForceModifier(new ModifierAveraging(this));
    adoptMediForceModifier(new ModifierMollification(this));
  }
  void MOLLYIntegrator::averagingPositions(){
	mySwapPositions = myPositions;
    myPositions = doAveragingPositions();

  }

  void MOLLYIntegrator::mollification(){
    myPositions = mySwapPositions;
    doMollification(myPositions);
  }
}
