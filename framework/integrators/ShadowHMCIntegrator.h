#ifndef SHADOWHMCINTEGRATOR_H
#define SHADOWHMCINTEGRATOR_H

#include "MCIntegrator.h"
#include "ModifierShadow.h"
#include "Report.h"
#include "topologyutilities.h"

namespace ProtoMol
{
	class GenericTopology;
	class ScalarStructure;
	class Vector3DBlock;
	class ForceGroup;

	//  ____________________________________________________ ShadowHMCIntegrator

	class ShadowHMCIntegrator : public MCIntegrator
	{
		//  ----------------------------------------------------------------  //
		//  Constructors, destructors, assignment
		//  ----------------------------------------------------------------  //
	public:
		ShadowHMCIntegrator();
		ShadowHMCIntegrator(int cycleLen,
		                    Real initTemp,
		                    int shadowOrder,
		                    Real c,
		                    bool optimize,
		                    ForceGroup* overloadedForces,
		                    StandardIntegrator* nextIntegrator);
		virtual ~ShadowHMCIntegrator();


		//  ----------------------------------------------------------------  //
		//  New methods of class ShadowHMCIntegrator
		//  ----------------------------------------------------------------  //
	private:
		void runPreSteps(int numTimesteps);
		void calcDelGstats(Real muH, Real varH, Real& muG, Real& varG);
		void optimizeC();// Real ratio );
		void calcMGStepVar(Real& mu, Real& var, int numSteps);


		//  ----------------------------------------------------------------  //
		//  From class Makeable
		//  ----------------------------------------------------------------  //
	public:
		virtual std::string getIdNoAlias() const
		{
			return keyword;
		}

		virtual unsigned int getParameterSize() const
		{
			return 5;
		}

		virtual void getParameters(std::vector<Parameter>& parameters) const;


		//  ----------------------------------------------------------------  //
		//  From class Integrator
		//  ----------------------------------------------------------------  //
	public:
		virtual void initialize(GenericTopology* topo,
		                        Vector3DBlock* positions,
		                        Vector3DBlock* velocities,
		                        ScalarStructure* energies);

		virtual void run(int numTimesteps);


		//  ----------------------------------------------------------------  //
		//  From class MTSIntegrator
		//  ----------------------------------------------------------------  //
	private:
		virtual MTSIntegrator* doMake(std::string& errMsg,
		                              const std::vector<Value>& values, ForceGroup* fg,
		                              StandardIntegrator* nextIntegrator) const;

		//  ----------------------------------------------------------------  //
		//  From class MCIntegrator
		//  ----------------------------------------------------------------  //
	protected:
		virtual void perturbSystem();
		virtual void saveValues();
		virtual void restoreValues();

		//  ----------------------------------------------------------------  //
		//  My public data members
		//  ----------------------------------------------------------------  //
	public:
		static const std::string keyword;

		//  ----------------------------------------------------------------  //
		//  My private data members
		//  ----------------------------------------------------------------  //
	private:
		unsigned int myOrder,
		             myShadowK;

		Real myC,
		     savedPE;

		ModifierShadow* shadowMod;

		Vector3DBlock myOldForces;

		bool myOptimize;
	};
}

#endif
