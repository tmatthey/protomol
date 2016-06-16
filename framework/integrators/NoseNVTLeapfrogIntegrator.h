/*  -*- c++ -*-  */
#ifndef NOSENVTLEAPFROGINTEGRATOR_H
#define NOSENVTLEAPFROGINTEGRATOR_H


#include "STSIntegrator.h"

namespace ProtoMol
{
	//_________________________________________________________________ NoseNVTLeapfrogIntegrator

	class ScalarStructure;
	class ForceGroup;
	class ModifierFriction;

	class NoseNVTLeapfrogIntegrator : public STSIntegrator
	{
		friend class ModifierFriction;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		NoseNVTLeapfrogIntegrator();
		NoseNVTLeapfrogIntegrator(Real timestep,
		                          Real temperature,
		                          Real thermalInertia,
		                          Real bathPosition,
		                          ForceGroup* overloadedForces);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class NoseNVTLeapfrogIntegrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		void friction();


		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Makeable
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual std::string getIdNoAlias() const
		{
			return keyword;
		}

		virtual void getParameters(std::vector<Parameter>& parameters) const;

		virtual unsigned int getParameterSize() const
		{
			return 4;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Integrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual void initialize(GenericTopology* topo,
		                        Vector3DBlock* positions,
		                        Vector3DBlock* velocities,
		                        ScalarStructure* energies);

	protected:
		virtual void addModifierAfterInitialize();
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
	private:
		const Real myTemperature;
		const Real myThermalInertia;
		Real myBathPosition;
		Real myTargetKE;
		Real mySumMass;
	};

	//______________________________________________________________________ INLINES
}

#endif
