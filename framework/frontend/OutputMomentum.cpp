#include "OutputMomentum.h"
#include "GenericTopology.h"
#include "topologyutilities.h"
#include "OutputCache.h"
#include "inputValueDefinitions.h"

#include <iomanip>

using namespace ProtoMol::Report;

using std::string;
using std::vector;
using std::setw;
using std::endl;
using std::flush;
using std::stringstream;
using std::setprecision;
using std::setiosflags;
using std::resetiosflags;
using std::ofstream;

namespace ProtoMol
{
	//________________________________________________________ Output
	const string OutputMomentum::keyword("momentumFile");

	OutputMomentum::OutputMomentum(): OutputFile()
	{
	}

	OutputMomentum::OutputMomentum(const string& filename,
	                               int freq,
	                               int cacheFreq,
	                               int cacheSize,
	                               Real closeTime): OutputFile(filename, freq, cacheFreq, cacheSize, closeTime)
	{
	}

	void OutputMomentum::doInitialize()
	{
		ofstream momentumHeaderFile(string(myFilename + ".header").c_str(), std::ios::out | std::ios::trunc);
		if (!momentumHeaderFile)
			report << error << " Can not open \'" << myFilename << ".header\' for " << getId() << "." << endr;

		momentumHeaderFile << setw(14)
			<< "Time(fs)" << " "
			<< setw(25)
			<< "10^6*LinMomentum_X" << " "
			<< setw(25)
			<< "10^6*LinMomentum_Y" << " "
			<< setw(25)
			<< "10^6*LinMomentum_Z" << " "
			<< setw(25)
			<< "10^6*AngMomentum_X" << " "
			<< setw(25)
			<< "10^6*AngMomentum_Y" << " "
			<< setw(25)
			<< "10^6*AngMomentum_Z";
		momentumHeaderFile << endl;
		momentumHeaderFile.close();
		open();
		close();
	}

	void OutputMomentum::doRunCached(int)
	{
		Real f = 1e6 * Constant::INV_TIMEFACTOR;

		myBuffer << resetiosflags(std::ios::showpoint | std::ios::fixed | std::ios::floatfield)
			<< setw(14) << setprecision(2)
			<< setiosflags(std::ios::showpoint | std::ios::fixed)
			<< myCache->time() << " "
			<< setw(25) << setprecision(5) << myCache->linearMomentum().x * f << " "
			<< setw(25) << setprecision(5) << myCache->linearMomentum().y * f << " "
			<< setw(25) << setprecision(5) << myCache->linearMomentum().z * f << " "
			<< setw(25) << setprecision(5) << myCache->angularMomentum().x * f << " "
			<< setw(25) << setprecision(5) << myCache->angularMomentum().y * f << " "
			<< setw(25) << setprecision(5) << myCache->angularMomentum().z * f;
		myBuffer << endl;
	}


	Output* OutputMomentum::doMake(string&, const vector<Value>& values) const
	{
		return (new OutputMomentum(values[0], values[1], values[2], values[3], values[4]));
	}
}
