#include "BerendsenIntegrator.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "ForceGroup.h"
#include "GenericTopology.h"
#include "pmconstants.h"
#include "topologyutilities.h"
#include "ModifierBerendsen.h"

using std::string;
using std::vector;

using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_________________________________________________________________ BerendsenIntegrator

	const string BerendsenIntegrator::keyword("Berendsen");

	BerendsenIntegrator::BerendsenIntegrator(): STSIntegrator(),
	                                            myTemperature(0.0),
	                                            myThermalCoupling(0.0)
	{
	}

	BerendsenIntegrator::BerendsenIntegrator(Real timestep,
	                                         Real temperature,
	                                         Real thermalCoupling,
	                                         ForceGroup* overloadedForces)

		: STSIntegrator(timestep, overloadedForces),
		  myTemperature(temperature),
		  myThermalCoupling(thermalCoupling)
	{
	}

	void BerendsenIntegrator::doBerendsen()
	{
		Real h = getTimestep() * Constant::INV_TIMEFACTOR;
		//
		Real Tsystem = 2.0 * kineticEnergy(myTopo, myVelocities) / (Constant::BOLTZMANN * myTopo->degreesOfFreedom);
		//reset energy to constant value
		Real efact;
		efact = sqrt(1.0 + (h / myThermalCoupling) * (myTemperature / Tsystem - 1.0));
		for (unsigned int k = 0; k < myVelocities->size(); k++) (*myVelocities)[k] *= efact;
	}

	void BerendsenIntegrator::initialize(GenericTopology* topo,
	                                     Vector3DBlock* positions,
	                                     Vector3DBlock* velocities,
	                                     ScalarStructure* energies)
	{
		STSIntegrator::initialize(topo, positions, velocities, energies);

		initializeForces();
	}

	void BerendsenIntegrator::addModifierAfterInitialize()
	{
		adoptPostStepModifier(new ModifierBerendsen(this));
		STSIntegrator::addModifierAfterInitialize();
	}

	void BerendsenIntegrator::getParameters(vector<Parameter>& parameters) const
	{
		STSIntegrator::getParameters(parameters);
		parameters.push_back(Parameter("temperature", Value(myTemperature, ConstraintValueType::NotNegative()), Text("preferred system temperature")));
		parameters.push_back(Parameter("thermalCoupling", Value(myThermalCoupling), Text("heat bath coupling.")));
	}

	STSIntegrator* BerendsenIntegrator::doMake(string&, const vector<Value>& values, ForceGroup* fg) const
	{
		return new BerendsenIntegrator(values[0], values[1], values[2], fg);
	}
}
