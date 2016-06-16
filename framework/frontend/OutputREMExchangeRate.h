/*  -*- c++ -*-  */
#ifndef OUTPUTREMEXCHANGERATE_H
#define OUTPUTREMEXCHANGERATE_H

#include "OutputFile.h"
#include "Parallel.h"

#include <string>
#include <vector>

namespace ProtoMol
{
	//_________________________________________________________________ OutputREMExchangeRate
	/** 
	 * Outputs the REM exchange rates at the end of the simulation. 
	 */
	class OutputREMExchangeRate : public OutputFile
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		OutputREMExchangeRate();
		OutputREMExchangeRate(const std::string& filename, int freq, int cacheFreq, int cacheSize, Real closeTime);
		~OutputREMExchangeRate();

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  From class OutputFile
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		virtual void doRunCached(int step);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  From class Output
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		virtual Output* doMake(std::string& errMesg, const std::vector<Value>& values) const;
		virtual void doInitialize();
		virtual void doFinalize(int step);

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
			return 5;
		}

		virtual bool adjustWithDefaultParameters(std::vector<Value>& values, const Configuration* config) const;
		virtual void getParameters(std::vector<Parameter>& parameter) const;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		static const std::string keyword;
	};
}

#endif
