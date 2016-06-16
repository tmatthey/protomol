#include "SystemCompareForce.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "GenericTopology.h"
using std::vector;
using std::string;
using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_________________________________________________________________ SystemCompareForce

	SystemCompareForce::SystemCompareForce(Force* actualForce, CompareForce* compareForce): CompareForce(actualForce, compareForce)
	{
	}

	void SystemCompareForce::evaluate(const GenericTopology* topo,
	                                  const Vector3DBlock* positions,
	                                  Vector3DBlock* forces,
	                                  ScalarStructure* energies)
	{
		preprocess(positions->size());
		(dynamic_cast<SystemForce*>(myActualForce))->evaluate(topo,
		                                                      positions,
		                                                      myForces,
		                                                      myEnergies);
		postprocess(topo, forces, energies);
	}

	void SystemCompareForce::parallelEvaluate(const GenericTopology* topo,
	                                          const Vector3DBlock* positions,
	                                          Vector3DBlock* forces,
	                                          ScalarStructure* energies)
	{
		preprocess(positions->size());
		(dynamic_cast<SystemForce*>(myActualForce))->parallelEvaluate(topo,
		                                                              positions,
		                                                              myForces,
		                                                              myEnergies);
		postprocess(topo, forces, energies);
	}
}
