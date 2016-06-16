#include "NormalModeLangevin.h"
#include "Report.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "ForceGroup.h"
#include "GenericTopology.h"
#include "topologyutilities.h"
#include "pmconstants.h"

#include "ModifierForceProjection.h"

using namespace std;

using namespace ProtoMol::Report;

using std::string;
using std::vector;

namespace ProtoMol
{
	//__________________________________________________ NormalModeLangevin

	const string NormalModeLangevin::keyword("NormalModeLangevin");

	NormalModeLangevin::NormalModeLangevin() : MTSIntegrator(), NormalModeUtilities()
	{
	}

	NormalModeLangevin::NormalModeLangevin(int cycles, int firstmode, int nummode, Real gamma, int seed, Real temperature, bool gencn,
	                                       ForceGroup* overloadedForces, StandardIntegrator* nextIntegrator)
		: MTSIntegrator(cycles, overloadedForces, nextIntegrator),
		  NormalModeUtilities(firstmode, nummode, gamma, seed, temperature), genCompNoise(gencn)
	{
	}


	NormalModeLangevin::~NormalModeLangevin()
	{
	}

	void NormalModeLangevin::initialize(GenericTopology* topo,
	                                    Vector3DBlock* positions,
	                                    Vector3DBlock* velocities,
	                                    ScalarStructure* energies)
	{
		MTSIntegrator::initialize(topo, positions, velocities, energies);
		//check valid eigenvectors
		//NM initialization if OK
		int nm_flags = NO_NM_FLAGS;
		if (genCompNoise) nm_flags |= GEN_COMP_NOISE;
		NormalModeUtilities::initialize((int)myPositions->size(), myTopo, myForces, nm_flags); //last int for no complimentary forces or gen noise: GEN_COMP_NOISE
		//
		//do first force calculation, and remove non sub-space part
		myEnergies->clear(); //Need this or initial error, due to inner integrator energy?
		initializeForces();
		//
		//take initial C velocites from system and remove non-subspace part
		if (*Q != NULL) subspaceVelocity(myVelocities, myVelocities);
		//
	}

	//*************************************************************************************
	//****Normal run routine***************************************************************
	//*************************************************************************************

	void NormalModeLangevin::run(int numTimesteps)
	{
		Real h = getTimestep() * Constant::INV_TIMEFACTOR;
		Real actTime;

		if (numTimesteps < 1)
			return;

		//check valid eigenvectors
		if (*Q == NULL)
			report << error << "No Eigenvectors for NormalMode integrator." << endr;
		//
		//time calculated in forces! so fix here
		actTime = myTopo->time + numTimesteps * getTimestep();
		//
		//main loop
		for (int i = 0; i < numTimesteps; i++)
		{
			//****main loop*************************************
			preStepModify();
			doHalfKick();
			//
			nmlDrift(myPositions, myVelocities, h, myTopo);
			//constraints?
			myEnergies->clear();
			//run minimizer if any remaining modes
			if (testRemainingModes()) myNextIntegrator->run(myCycleLength); //cyclelength 
			if (*Q == NULL)
			{ //rediagonalize?
				myTopo->time = actTime - (i - numTimesteps) * getTimestep();
				if (myPreviousIntegrator == NULL)
					report << error << "[NormalModeLangevin::Run] Re-diagonalization forced with NormalModeLangevin as outermost Integrator. Aborting." << endr;
				return;
			}
			//calculate sub space forces
			myEnergies->clear();
			calculateForces();
			//
			doHalfKick();
			//
			postStepModify();
		}
		//fix time
		myTopo->time = actTime;
		//
	}

	//*************************************************************************************
	//****Output int paramiters************************************************************
	//*************************************************************************************

	void NormalModeLangevin::getParameters(vector<Parameter>& parameters) const
	{
		MTSIntegrator::getParameters(parameters);
		parameters.push_back(Parameter("firstmode", Value(firstMode, ConstraintValueType::NotNegative()), 1, Text("First mode to use in set")));
		parameters.push_back(Parameter("numbermodes", Value(numMode, ConstraintValueType::NotNegative()), 1, Text("Number of modes propagated")));
		parameters.push_back(Parameter("gamma", Value(myGamma * (1000 * Constant::INV_TIMEFACTOR), ConstraintValueType::NotNegative()), 80.0, Text("Langevin Gamma")));
		parameters.push_back(Parameter("seed", Value(mySeed, ConstraintValueType::NotNegative()), 1234, Text("Langevin random seed")));
		parameters.push_back(Parameter("temperature", Value(myTemp, ConstraintValueType::NotNegative()), 300.0, Text("Langevin temperature")));
		parameters.push_back(Parameter("gencompnoise", Value(genCompNoise, ConstraintValueType::NoConstraints()), false, Text("Generate complimentary noise")));
	}

	MTSIntegrator* NormalModeLangevin::doMake(string&, const vector<Value>& values, ForceGroup* fg, StandardIntegrator* nextIntegrator) const
	{
		return new NormalModeLangevin(values[0], values[1], values[2], values[3], values[4], values[5], values[6], fg, nextIntegrator);
	}

	void NormalModeLangevin::addModifierAfterInitialize()
	{
		adoptPostForceModifier(new ModifierForceProjection(this));
		MTSIntegrator::addModifierAfterInitialize();
	}
}
