/*  -*- c++ -*-  */
#ifndef OUTPUTPAULTRAP_H
#define OUTPUTPAULTRAP_H

#include "OutputFile.h"
#include "Vector3DBlock.h"

namespace ProtoMol
{
	//________________________________________________________ OutputPaulTrap
	class OutputPaulTrap : public OutputFile
	{
		// Writes all output of interest for Coulomb Crystal simulations.
		//
		// Note, that you can also specify the caching, such that you only
		// write after a number of calls of run().
		//
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		OutputPaulTrap();
		OutputPaulTrap(const std::string& filename, int freq, int cacheFreq, int cacheSize,
		               Real closeTime, Real omegar, Real omegaz, const std::string& filenameLow,
		               bool doLow, bool screen);
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class OutputPaulTrap
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  From class OutputFile
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		virtual void doRunCached(int step);
		virtual void doFlushCache();
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  From class Output
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		virtual Output* doMake(std::string& errMsg, const std::vector<Value>& values) const;
		virtual void doInitialize();


		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// From class Makeable
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		virtual std::string getIdNoAlias() const
		{
			return keyword;
		}

		// Returns the identification string

		virtual unsigned int getParameterSize() const
		{
			return 10;
		}

		virtual void getParameters(std::vector<Parameter>& parameter) const;
		virtual bool adjustWithDefaultParameters(std::vector<Value>& values, const Configuration* config) const;
	private:
		void doWrite();
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		static const std::string keyword;

	private:
		Real myOmegaR;
		Real myOmegaZ;
		std::string myFilenameLow;
		bool myDoPaulLow;
		bool myScreen;
	private:
		Real myPaulOmega;
		Real myPaulQ;
		Real myPaulK;
		Real myPaulR;
		Real myPaulA;
		Real myPaulUHom;
		Real myPaulF;
		Real myPaulF2;
		Real myPaulM;
		Real myPaulLow;
		std::string myLowOut;
		Vector3DBlock myLowXYZ;
		std::string myLowComment;
		std::string myUnit;
		Real myFactor;
	};
}
#endif
