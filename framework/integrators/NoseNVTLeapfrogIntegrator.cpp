#include "NoseNVTLeapfrogIntegrator.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "ForceGroup.h"
#include "GenericTopology.h"
#include "pmconstants.h"
#include "topologyutilities.h"
#include "ModifierFriction.h"

using std::string;
using std::vector;

using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_________________________________________________________________ NoseNVTLeapfrogIntegrator

	const string NoseNVTLeapfrogIntegrator::keyword("NoseNVTLeapfrog");

	NoseNVTLeapfrogIntegrator::NoseNVTLeapfrogIntegrator(): STSIntegrator(),
	                                                        myTemperature(0.0),
	                                                        myThermalInertia(0.0),
	                                                        myBathPosition(0.0)
	{
	}

	NoseNVTLeapfrogIntegrator::NoseNVTLeapfrogIntegrator(Real timestep,
	                                                     Real temperature,
	                                                     Real thermalInertia,
	                                                     Real bathPosition,
	                                                     ForceGroup* overloadedForces)

		: STSIntegrator(timestep, overloadedForces),
		  myTemperature(temperature),
		  myThermalInertia(thermalInertia),
		  myBathPosition(bathPosition)
	{
	}

	void NoseNVTLeapfrogIntegrator::friction()
	{
		const unsigned int numberOfAtoms = myTopo->atoms.size();
		const Real h = getTimestep() * Constant::INV_TIMEFACTOR;
		Real kineticEnergy = 0.0;
		for (unsigned int i = 0; i < numberOfAtoms; ++i)
		{
			Real m = myTopo->atoms[i].scaledMass;
			kineticEnergy += ((*myVelocities)[i]).normSquared() * m;//+(*myForces)[i]*h/m
		}
		myBathPosition += (0.5 * kineticEnergy - myTargetKE) * myThermalInertia / (h * numberOfAtoms);

		for (unsigned int i = 0; i < numberOfAtoms; i++)
		{
			(*myForces)[i] -= (*myVelocities)[i] * myBathPosition * myTopo->atoms[i].scaledMass;
		}
	}

	void NoseNVTLeapfrogIntegrator::initialize(GenericTopology* topo,
	                                           Vector3DBlock* positions,
	                                           Vector3DBlock* velocities,
	                                           ScalarStructure* energies)
	{
		STSIntegrator::initialize(topo, positions, velocities, energies);

		myTargetKE = (3.0 / 2.0 * positions->size() * Constant::BOLTZMANN * myTemperature);
		mySumMass = 0.0;
		for (unsigned int i = 0; i < positions->size(); i++)
			mySumMass += myTopo->atoms[i].scaledMass;

		initializeForces();
	}

	void NoseNVTLeapfrogIntegrator::addModifierAfterInitialize()
	{
		adoptPostForceModifier(new ModifierFriction(this));
		STSIntegrator::addModifierAfterInitialize();
	}

	void NoseNVTLeapfrogIntegrator::getParameters(vector<Parameter>& parameters) const
	{
		STSIntegrator::getParameters(parameters);
		parameters.push_back(Parameter("temperature", Value(myTemperature, ConstraintValueType::NotNegative()), Text("preferred system temperature")));
		parameters.push_back(Parameter("thermal", Value(myThermalInertia), Text("heat bath coupling: 1.0 very, very strong, 0.0 none")));
		parameters.push_back(Parameter("bathPos", Value(myBathPosition), 0.0, Text("history of the difference of system and heat bath")));
	}

	STSIntegrator* NoseNVTLeapfrogIntegrator::doMake(string&, const vector<Value>& values, ForceGroup* fg) const
	{
		return new NoseNVTLeapfrogIntegrator(values[0], values[1], values[2], values[3], fg);
	}
}
