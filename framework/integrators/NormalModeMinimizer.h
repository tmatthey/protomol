/*  -*- c++ -*-  */
#ifndef NORMALMODEMINIMIZER_H
#define NORMALMODEMINIMIZER_H

#include "STSIntegrator.h"
#include "NormalModeUtilities.h"

namespace ProtoMol
{
	class ScalarStructure;
	class ForceGroup;

	//__________________________________________________ NormalModeMinimizer
	class NormalModeMinimizer : public STSIntegrator, public NormalModeUtilities
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		NormalModeMinimizer();
		NormalModeMinimizer(Real timestep, int firstmode, int nummode, Real gamma, int seed, Real temperature,
		                    Real minimlim, int rforce, bool simplemin, ForceGroup* overloadedForces);
		~NormalModeMinimizer();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class NormalModeMinimizer
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		void utilityCalculateForces();
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Makeable
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual std::string getIdNoAlias() const
		{
			return keyword;
		}

		virtual unsigned int getParameterSize() const
		{
			return 9;
		}

		virtual void getParameters(std::vector<Parameter>& parameters) const;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Integrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual void initialize(GenericTopology* topo,
		                        Vector3DBlock* positions,
		                        Vector3DBlock* velocities,
		                        ScalarStructure* energies);
		virtual void run(int numTimesteps);
	protected:
		//virtual void addModifierAfterInitialize();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class STSIntegrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		virtual STSIntegrator* doMake(std::string& errMsg, const std::vector<Value>& values, ForceGroup* fg) const;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		static const std::string keyword;
		int avItrs, itrs, avMinForceCalc, numSteps;

	private:
		int minCount, forceCalc;
		Real minLim;
		int randforce;
		NormalModeUtilities* myPreviousNormalMode;
		Real lastLambda;
		bool reDiag, simpleMin;
		Real eUFactor, randStp;
	};
}

#endif
