/*  -*- c++ -*-  */
#ifndef NORMALMODEBROWNIAN_H
#define NORMALMODEBROWNIAN_H

#include "STSIntegrator.h"
#include "NormalModeUtilities.h"

//####diagnostics
#include "XYZTrajectoryWriter.h"


namespace ProtoMol {

	class ScalarStructure;
	class ForceGroup;

	//__________________________________________________ NormalModeBrownian
	class NormalModeBrownian : public STSIntegrator, public NormalModeUtilities {
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		NormalModeBrownian();
		NormalModeBrownian(Real timestep, int firstmode, int nummode, Real gamma, int seed, Real temperature, 
			std::string avff, std::string inff, //####added avff, inff for diagnostics
			ForceGroup *overloadedForces);
		~NormalModeBrownian(); 

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class NormalModeBrownian
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Makeable
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual std::string getIdNoAlias() const{return keyword;}
		virtual unsigned int getParameterSize() const{return 8;}	//####6 if no diagnostics
		virtual void getParameters(std::vector<Parameter>& parameters) const;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Integrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual void initialize(GenericTopology *topo,
			Vector3DBlock   *positions,
			Vector3DBlock   *velocities, 
			ScalarStructure *energies);
		virtual void run(int numTimesteps);
	protected:
		virtual void addModifierAfterInitialize();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class STSIntegrator
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		virtual STSIntegrator* doMake(std::string& errMsg, const std::vector<Value>& values, ForceGroup* fg)const;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class NormalModeBrownian
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual void forceProjection();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		static const std::string keyword;

	private:
		Real randStp;
		int aveForceCount;
		//####diagnostics
		// mean force output:
		std::string avForceFile;
		// instantaneous force output:
		std::string inForceFile; 

		XYZTrajectoryWriter *myWriter;

		XYZTrajectoryWriter *myWriter2;


	};
}

#endif


