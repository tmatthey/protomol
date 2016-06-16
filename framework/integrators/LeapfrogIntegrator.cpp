#include "LeapfrogIntegrator.h"
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
	//__________________________________________________ LeapfrogIntegrator

	const string LeapfrogIntegrator::keyword("Leapfrog");

	LeapfrogIntegrator::LeapfrogIntegrator() : STSIntegrator()
	{
	}

	LeapfrogIntegrator::LeapfrogIntegrator(Real timestep,
	                                       ForceGroup* overloadedForces)
		: STSIntegrator(timestep, overloadedForces)
	{
	}


	LeapfrogIntegrator::~LeapfrogIntegrator()
	{
	}


	void LeapfrogIntegrator::initialize(GenericTopology* topo,
	                                    Vector3DBlock* positions,
	                                    Vector3DBlock* velocities,
	                                    ScalarStructure* energies)
	{
		STSIntegrator::initialize(topo, positions, velocities, energies);
		initializeForces();
	}

	void LeapfrogIntegrator::doHalfKickdoDrift()
	{
		if (anyPreDriftOrNextModify())
		{
			doHalfKick();
			doDriftOrNextIntegrator();
		}
		else
		{
			Real h = getTimestep() * Constant::INV_TIMEFACTOR;
			const unsigned int count = myPositions->size();

			//  Do a half kick on beta.
			updateBeta(0.5 * h);

			for (unsigned int i = 0; i < count; ++i)
			{
				(*myVelocities)[i] += (*myForces)[i] * h * 0.5 / myTopo->atoms[i].scaledMass;
				(*myPositions)[i] += (*myVelocities)[i] * h;
			}
			buildMolecularCenterOfMass(myPositions, myTopo);
			buildMolecularMomentum(myVelocities, myTopo);
			postDriftOrNextModify();
		}
	}

	void LeapfrogIntegrator::doKickdoDrift()
	{
		if (anyPreDriftOrNextModify() || anyPreStepModify() || anyPostStepModify())
		{
			if (anyPreStepModify() || anyPostStepModify())
			{
				doHalfKick();
				postStepModify();
				preStepModify();
				doHalfKick();
			}
			else
			{
				doKick();
			}
			doDriftOrNextIntegrator();
		}
		else
		{
			Real h = getTimestep() * Constant::INV_TIMEFACTOR;
			const unsigned int count = myPositions->size();

			updateBeta(h);

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

	void LeapfrogIntegrator::run(int numTimesteps)
	{
		if (numTimesteps < 1)
			return;
		preStepModify();
		doHalfKickdoDrift();
		calculateForces();
		for (int i = 1; i < numTimesteps; i++)
		{
			doKickdoDrift();
			calculateForces();
		}
		doHalfKick();
		postStepModify();
	}

	STSIntegrator* LeapfrogIntegrator::doMake(string&, const vector<Value>& values, ForceGroup* fg) const
	{
		return new LeapfrogIntegrator(values[0], fg);
	}


	//  --------------------------------------------------------------------  //
	//  This function is necessary to compute the shadow Hamiltonian and it   //
	//  is integrator specific.  This version is written to work with LF.     //
	//  Update beta: beta -= dt * ( q * F + 2 U )                             //
	//  --------------------------------------------------------------------  //

	void LeapfrogIntegrator::updateBeta(Real dt)
	{
		//  ----------------------------------------------------------------  //
		//  The shadow calculation is done in a postStep modifier.  If there  //
		//  aren't any, then obviously we don't need to do this calculation.  //
		//  It's possible that a different poststep modifier could make this  //
		//  execute, but no harm would be done ... only some extra cycles.    //
		//  ----------------------------------------------------------------  //

		if (! (anyPostStepModify() || top()->anyPostStepModify()))
			return;

		Real posDotF = 0.;

		for (unsigned int i = 0; i < myPositions->size(); i++)
			posDotF += (*myPositions)[i].dot((*myForces)[i]);

		myBeta -= dt * (posDotF + 2. * myPotEnergy);
	}
}
