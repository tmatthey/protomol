/*  -*- c++ -*-  */
#ifndef NORMMODEDIAG_H
#define NORMMODEDIAG_H

#include "MTSIntegrator.h"
#include "Vector3DBlock.h"
#include "Hessian.h"

#include "CheckpointInputStream.h"
#include "CheckpointOutputStream.h"

namespace ProtoMol
{
	class ScalarStructure;
	class ForceGroup;

	//__________________________________________________ NormModeDiag
	class NormModeDiag : public MTSIntegrator
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		NormModeDiag();
		NormModeDiag(int cycles, int avs, Real avss, int redi, int rayf, int raya, std::string ray_s, bool raysf, int mins, Real minl, ForceGroup* overloadedForces, StandardIntegrator* nextIntegrator);
		~NormModeDiag();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class NormModeDiag
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		int diagHessian(double* eigVecO, double* eigValO);
		void absSort();
		void calcRayleigh();
		int SDminimize(Real peLim, int numlp);
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
			return 10;
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
		virtual MTSIntegrator* doMake(std::string& errMsg, const std::vector<Value>& values, ForceGroup* fg, StandardIntegrator* nextIntegrator) const;
	public:
		void saveState(CheckpointOutputStream& os);
		void restoreState(CheckpointInputStream& is);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		static const std::string keyword;

	private:
		unsigned int sz;
		Vector3DBlock diagAt;
		Hessian hsn;
		double *eigVec, *eigVal;
		int* eigIndx;
		bool eigAlloc;
		int noAvStep;
		Real avStep;
		int rediagCount, nextRediag, raylFrequ, raylAverage, nextRayl, raylAvCount;
		bool raylDo;
		std::string raylFile;
		bool raylSift;
		double* mhQhQt;
		int minSteps;
		Real minLim;
		bool validMaxEigv;
	};
}

#endif
