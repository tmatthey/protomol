/*  -*- c++ -*-  */
#ifndef NORMMODEVISUAL_H
#define NORMMODEVISUAL_H

#include "STSIntegrator.h"

namespace ProtoMol
{
	class ScalarStructure;
	class ForceGroup;

	//__________________________________________________ NormModeVisual
	class NormModeVisual : public STSIntegrator
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		NormModeVisual();
		NormModeVisual(Real timestep, int fixmodes, std::string evals, int mstrt, int mend, int cyclest, Real tmpr, bool scle, std::string scndp, std::string mfil, ForceGroup* overloadedForces);
		~NormModeVisual();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class NormModeVisual
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
		virtual STSIntegrator* doMake(std::string& errMsg, const std::vector<Value>& values, ForceGroup* fg) const;
		double* vector3DBlockTOvect(Vector3DBlock* blkDat, double* vecDat);
		Vector3DBlock* vectTOvector3DBlock(double* vecDat, Vector3DBlock* blkDat);
		Vector3DBlock* subspaceProj(double* tmpC, Vector3DBlock* iPos);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		static const std::string keyword;
	private:
		Vector3DBlock* ex0;
		Vector3DBlock *eigvalp, *secondPos;
		double *tmpFX, *cPos;

		int _N, _m, _rfM, _3N, avItrs, numSteps;
		int fixMod;
		std::string evalfile;
		int modeStart, modeEnd, cycleSteps;
		Real tempT;
		bool fScale;
		std::string secondPosFile, modeFile;
		int currMode;
	};
}

#endif
