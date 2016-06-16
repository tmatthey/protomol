#include "PLeapfrogIntegrator.h"
#include "Report.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "ForceGroup.h"
#include "GenericTopology.h"
#include "topologyutilities.h"
#include "pmconstants.h"

using std::string;
using std::vector;
using namespace ProtoMol::Report;

namespace ProtoMol
{
	//__________________________________________________ PLeapfrogIntegrator

	const string PLeapfrogIntegrator::keyword("PLeapfrog");

	PLeapfrogIntegrator::PLeapfrogIntegrator() : STSIntegrator(), myTempForces(NULL)
	{
	}

	PLeapfrogIntegrator::PLeapfrogIntegrator(Real timestep,
	                                         ForceGroup* overloadedForces)
		: STSIntegrator(timestep, overloadedForces), myTempForces(new Vector3DBlock())
	{
	}


	PLeapfrogIntegrator::~PLeapfrogIntegrator()
	{
		if (myTempForces != NULL)
			delete myTempForces;
	}


	void PLeapfrogIntegrator::initialize(GenericTopology* topo,
	                                     Vector3DBlock* positions,
	                                     Vector3DBlock* velocities,
	                                     ScalarStructure* energies)
	{
		STSIntegrator::initialize(topo, positions, velocities, energies);
		initializeForces();
		myTempForces->resize(positions->size());
	}

	void PLeapfrogIntegrator::doKickdoDrift()
	{
		if (anyPreDriftOrNextModify() || anyPreStepModify() || anyPostStepModify())
		{
			doKick();
			if (anyPreStepModify() || anyPostStepModify())
			{
				preDriftOrNextModify();
				doHalfDrift();
				postStepModify();
				preStepModify();
				doHalfDrift();
				postDriftOrNextModify();
			}
			else
			{
				doDrift();
			}
		}
		else
		{
			Real h = getTimestep() * Constant::INV_TIMEFACTOR;
			const unsigned int count = myPositions->size();

			for (unsigned int i = 0; i < count; ++i)
			{
				(*myVelocities)[i] += (*myForces)[i] * h / myTopo->atoms[i].scaledMass;
				(*myPositions)[i] += (*myVelocities)[i] * h;
			}
			buildMolecularCenterOfMass(myPositions, myTopo);
			buildMolecularMomentum(myVelocities, myTopo);
			postDriftOrNextModify();
		}
	}

	void PLeapfrogIntegrator::doKickdoHalfDrift()
	{
		if (anyPreDriftOrNextModify())
		{
			doKick();
			preDriftOrNextModify();
			doHalfDrift();
		}
		else
		{
			Real h = getTimestep() * Constant::INV_TIMEFACTOR;
			const unsigned int count = myPositions->size();

			for (unsigned int i = 0; i < count; ++i)
			{
				(*myVelocities)[i] += (*myForces)[i] * h / myTopo->atoms[i].scaledMass;
				(*myPositions)[i] += (*myVelocities)[i] * 0.5 * h;
			}
			buildMolecularCenterOfMass(myPositions, myTopo);
			buildMolecularMomentum(myVelocities, myTopo);
		}
	}

	void PLeapfrogIntegrator::doHalfDrift()
	{
		Real h = getTimestep() * Constant::INV_TIMEFACTOR;
		const unsigned int count = myPositions->size();

		for (unsigned int i = 0; i < count; ++i)
		{
			(*myPositions)[i] += (*myVelocities)[i] * 0.5 * h;
		}
		buildMolecularCenterOfMass(myPositions, myTopo);
	}

	void PLeapfrogIntegrator::run(int numTimesteps)
	{
		if (numTimesteps < 1)
			return;

		preStepModify();
		doHalfDrift();
		postDriftOrNextModify();
		calculateForces();
		for (int i = 1; i < numTimesteps; ++i)
		{
			doKickdoDrift();
			calculateForces();
		}
		doKickdoHalfDrift();

		// Correction of energy ..
		Vector3DBlock* temp = myForces;
		myForces = myTempForces;
		Real t = myTopo->time;
		calculateForces();
		myTopo->time = t;
		myForces = temp;

		postStepModify();
	}

	STSIntegrator* PLeapfrogIntegrator::doMake(string&, const vector<Value>& values, ForceGroup* fg) const
	{
		return new PLeapfrogIntegrator(values[0], fg);
	}
}
