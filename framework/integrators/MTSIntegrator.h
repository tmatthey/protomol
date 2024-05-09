/*  -*- c++ -*-  */
#ifndef MTSINTEGRATOR_H
#define MTSINTEGRATOR_H

#include "StandardIntegrator.h"

namespace ProtoMol
{
	class GenericTopology;
	class ScalarStructure;
	class Vector3DBlock;
	class ForceGroup;

	//_________________________________________________________________ MTSIntegrator

	class MTSIntegrator: public StandardIntegrator
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		// totally overload the function.  Use no defaults (unless the lists passed
		//   are empty lists of forces.  Then the default forces will be assigned);
		//   however, the lists still must be passed nonetheless.
		MTSIntegrator();
		MTSIntegrator(int cycles, ForceGroup* overloadedForces, StandardIntegrator* nextIntegrator);

		virtual ~MTSIntegrator();


		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class MTSIntegrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		MTSIntegrator* make(std::string& errMsg, const std::vector<Value>& values, ForceGroup* fg, StandardIntegrator* nextIntegrator) const;
	private:
		virtual MTSIntegrator* doMake(std::string& errMsg, const std::vector<Value>& values, ForceGroup* fg, StandardIntegrator* nextIntegrator) const =0;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Makeable
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual void getParameters(std::vector<Parameter>& parameter) const;

		virtual unsigned int getParameterSize() const
		{
			return 1;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Integrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual Real getTimestep() const;
		virtual void initialize(GenericTopology* topo,
		                        Vector3DBlock* positions,
		                        Vector3DBlock* velocities,
		                        ScalarStructure* energies);

		virtual Integrator* next();
		virtual const Integrator* next() const;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class StandardIntegrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		virtual void doDriftOrNextIntegrator();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		StandardIntegrator* myNextIntegrator;
		const int myCycleLength;
	};

	//______________________________________________________________________ INLINES

	inline Real MTSIntegrator::getTimestep() const
	{
		return myCycleLength * myNextIntegrator->getTimestep();
	}

	inline Integrator* MTSIntegrator::next()
	{
		return myNextIntegrator;
	}

	inline const Integrator* MTSIntegrator::next() const
	{
		return myNextIntegrator;
	}
}
#endif
