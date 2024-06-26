/*  -*- c++ -*-  */
#ifndef NUMERICALDIFFERENTIATION_H
#define NUMERICALDIFFERENTIATION_H

#include "STSIntegrator.h"
#include "Hessian.h"

namespace ProtoMol
{
	class ScalarStructure;
	class ForceGroup;

	//__________________________________________________ NumericalDifferentiation
	class NumericalDifferentiation : public STSIntegrator
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		NumericalDifferentiation();
		NumericalDifferentiation(Real timestep, Real epsil, ForceGroup* overloadedForces);
		~NumericalDifferentiation();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class NumericalDifferentiation
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:

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
			return 2;
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
		Real epsilon;
		unsigned int _N, _3N;
		Hessian hsn;
	};
}

#endif
