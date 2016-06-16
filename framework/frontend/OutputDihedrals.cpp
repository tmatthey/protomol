#include "OutputDihedrals.h"
#include "GenericTopology.h"
#include "topologyutilities.h"
#include "OutputCache.h"
#include "DCDTrajectoryWriter.h"
#include "inputValueDefinitions.h"

#include <iomanip>
#include <algorithm>

using namespace ProtoMol::Report;

using std::string;
using std::vector;
using std::set;
using std::setw;
using std::endl;
using std::flush;
using std::stringstream;
using std::setprecision;
using std::setiosflags;
using std::resetiosflags;
using std::ofstream;
using std::ifstream;

#define DEBUG_OUTPUTDIHEDRALS

namespace ProtoMol
{
	//________________________________________________________ Output
	const string OutputDihedrals::keyword("dihedralsFile");

	OutputDihedrals::OutputDihedrals(): OutputFile(),
	                                    myDihedralIndex(-1),
	                                    myDihedralsSet(false),
	                                    myDCD(NULL),
	                                    myMinimalImage(false)
	{
	}

	OutputDihedrals::OutputDihedrals(const string& filename,
	                                 int freq,
	                                 int cacheFreq,
	                                 int cacheSize,
	                                 Real closeTime,
	                                 bool minimal,
	                                 int index,
	                                 bool dihset,
	                                 std::string dsetfile):
		OutputFile(filename, freq, cacheFreq, cacheSize, closeTime),
		myDihedralIndex(index),
		myDihedralsSet(dihset),
		myDihedralsSetfile(dsetfile),
		myDCD(NULL),
		myMinimalImage(minimal)
	{
		oldconformstring = "";
	}

	OutputDihedrals::~OutputDihedrals()
	{
		if (myDCD != NULL)
			delete myDCD;
		myDCD = NULL;
	}

	void OutputDihedrals::doInitialize()
	{
		if (myDihedralIndex < 0 || myDihedralIndex >= static_cast<int>(myTopology->dihedrals.size()))
		{
			report << error << "[OutputDihedrals::doInitialize] Dihedral index " << myDihedralIndex
				<< " out of range [0," << myTopology->dihedrals.size() - 1 << "]." << endr;
		}

		ofstream dihedralsHeaderFile(string(myFilename + ".header").c_str(), std::ios::out | std::ios::trunc);
		if (!dihedralsHeaderFile)
			report << error << " Can not open \'" << myFilename << ".header\' for " << getId() << "." << endr;

		dihedralsHeaderFile << setiosflags(std::ios::showpoint | std::ios::fixed)
			<< setw(14) << "Time (fs)"
			<< setw(12) << "DihIndex"
			<< setw(13) << "Val (rad)"
			<< setw(21) << "Energy (kcal/mol)"
			<< endl;

		if (myDihedralsSet)
		{
			ifstream dihedralsSetinput(string(myDihedralsSetfile).c_str(), std::ios::in);
			if (!dihedralsSetinput)
				report << error << " Could not open \'" << myDihedralsSetfile << "\' for " << getId() << "." << endr;
			int tempdihedral = 1;
			while (dihedralsSetinput >> tempdihedral)
				myDihedrals.push_back(tempdihedral);
		}
		else
		{
			myDihedrals.push_back(myDihedralIndex);
		}

		// Assigns incremental dihedral value for conformation string analysis.
		// The term maxima and minima is irrelevant in this case but remains as
		// not to conflict with the alternative dihedral well assignment below.

		myMaximas = dihInc(myDihedrals);
		myMinimas = dihInc(myDihedrals);

		// Calculate the location of energy maximas for each dihedral
		// This implementation is correct however the application is not well
		// suited for backbone phi,psi analysis as the energy function has a single well.
		// In such a case use the simple dihedral value increment assignment function above
		//
		// myMaximas = myCache->brentMaxima(myDihedrals, true);
		// myMinimas = myCache->brentMaxima(myDihedrals, false);

		dihedralsHeaderFile.close();

		//getbackbonedihedrals(); //Use this function as needed to get approximate backbone dihedral set
		//getdihedralsfromatomset(); //Use this function to get dihedrals based on an atom set

		flipcounts.resize(myMaximas.size());
		for (unsigned int i = 0; i < myMaximas.size(); ++i)
		{
			(flipcounts[i]).resize((myMaximas[i]).size());
			for (unsigned int j = 0; j < (myMaximas[i]).size(); j++)
			{
				(flipcounts[i])[j] = 0;
			}
		}
		/**
		    std::vector< std::vector< Real > >::iterator maximaset_itr = myMaximas.begin();
		    for(std::vector<Real>::iterator dihedral_itr = dihedrals.begin();
		        dihedral_itr != (dihedrals.end()); ++dihedral_itr)
		**/
		open();
		close();
	}

	// The run function outputs the dihedral values and computes the conformation strings
	// utilizing calls to Brent Maxima finding functions in the Output Cache.
	// Unique Conf Strings and their respective DCDs are stored
	void OutputDihedrals::doRunCached(int)
	{
		unsigned int old_confstrings_size = 0;

		string conformstring = "";
		//vector<int> dihStatesVec; to ID states with int instead of char
		//dihStatesVec.clear();
		vector<Real> dihedrals = myCache->dihedralPhis(myDihedrals);

		Real tempdihedral = 1.0;

		// This starts the dihedral value output file. The dihedral values output code is further
		// down the page within a loop which iterates through the dihedrals
		myBuffer << setiosflags(std::ios::showpoint | std::ios::fixed)
			<< setw(14) << setprecision(3) << myCache->time();

		//This section assigns the conformation string at each time step

		if (dihedrals.size() != myMaximas.size())
			report << error << "ERROR - SIZES DON'T MATCH" << endr;

		std::vector<std::vector<Real>>::iterator maximaset_itr = myMaximas.begin();
		int dihedralindex = 0;
		const string mycharacters("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789");

		for (std::vector<Real>::iterator dihedral_itr = dihedrals.begin();
		     dihedral_itr != (dihedrals.end()); ++dihedral_itr)
		{
			tempdihedral = *dihedral_itr;
			if (tempdihedral < 0.0)
			{
				tempdihedral = 2.0 * Constant::M_PI + tempdihedral;
			}

			std::vector<Real>::iterator maxima_itr = maximaset_itr->begin();
			int j = 0;
			while ((maxima_itr != maximaset_itr->end()) && (tempdihedral > (*maxima_itr)))
			{
				maxima_itr++;
				j++;
			}

			// if the dihedral value is greater than each of the maxima
			if (maxima_itr == maximaset_itr->end())
				j = 0;

			if (j < 36)
			{
				//dihStatesVec.pushback(j);  
				conformstring += mycharacters[j];
			}
			else
			{
				conformstring += '!';
			}

			//outputs the dihedral values if the dihedral set is small ( <=50 )
			if (dihedrals.size() <= 50)
			{
				myBuffer << setw(9) << myDihedrals[dihedralindex]
					<< setw(15) << setprecision(4) << tempdihedral
					<< setw(18) << computePhiDihedralEnergy(myTopology, myDihedrals[dihedralindex], tempdihedral);
			}

			maximaset_itr++;
			dihedralindex++;
		}

		//insert string into set to establish uniqueness
		old_confstrings_size = myConfstrings.size();
		myConfstrings.insert(conformstring);


		//if the string was unique store it's DCD
		if ((myConfstrings.size()) > old_confstrings_size)
		{
			confstrings.push_back(conformstring);
			confstringsCounter.push_back(1);
			if (myDCD == NULL)
				myDCD = new DCDTrajectoryWriter(string(myFilename + ".dcds").c_str());
			const Vector3DBlock* pos = (myMinimalImage ? myCache->minimalPositions() : myPositions);
			if (!myDCD->write(*pos))
				report << error << "Could not write " << getId() << " \'" << myDCD->getFilename() << "\'." << endr;
		}
		else
		{
			for (unsigned int i = 0; i < confstrings.size(); i++)
			{
				if (confstrings[i] == conformstring)
					confstringsCounter[i] += 1;
			}
		}

		myBuffer << endl;

		if (oldconformstring.size() > 1)
		{
			for (unsigned int i = 0; i < conformstring.size(); i++)
			{
				if (oldconformstring[i] != conformstring[i])
				{
					for (unsigned int j = 0; j < (myMaximas[i]).size(); j++)
					{
						if (conformstring[i] == mycharacters[j])
							(flipcounts[i])[j] += 1;
					}
				}
			}
		}

		oldconformstring = conformstring;
	}

	void OutputDihedrals::doFinalize(int)
	{
		const string mycharacters("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789");

		// outputs the file holding the conformation strings
		ofstream confstringsFile(string(myFilename + ".confstrings").c_str(), std::ios::out | std::ios::trunc);
		if (!confstringsFile)
			report << error << " Can not open \'" << myFilename << ".confstrings\' for " << getId() << "." << endr;

		confstringsFile << myConfstrings.size()
			<< setw(24)
			<< "unique conformation(s)"
			<< endl << endl;

		confstringsFile << "The dihedral angles for which dihedral energy maxima are found" << endl
			<< "Used to identify bounds on conformation well" << endl;
		int maxcounter = 0;
		for (std::vector<vector<Real>>::iterator myMaximaSetItr = myMaximas.begin();
		     myMaximaSetItr != myMaximas.end(); ++myMaximaSetItr)
		{
			confstringsFile << "Dihedral Index: " << myDihedrals[maxcounter] << " angles: ";
			for (std::vector<Real>::iterator myMaximasItr = (*myMaximaSetItr).begin();
			     myMaximasItr != (*myMaximaSetItr).end(); ++myMaximasItr)
			{
				confstringsFile << " " << *myMaximasItr << " ";
			}
			confstringsFile << endl << endl;
			maxcounter++;
		}

		confstringsFile << "The dihedral angles for which dihedral energy minima are found" << endl
			<< "Used to identify potential rotamer values" << endl;
		int mincounter = 0;
		for (std::vector<vector<Real>>::iterator myMinimaSetItr = myMinimas.begin();
		     myMinimaSetItr != myMinimas.end(); ++myMinimaSetItr)
		{
			confstringsFile << "Dihedral Index: " << myDihedrals[mincounter] << " angles: ";
			for (std::vector<Real>::iterator myMinimasItr = (*myMinimaSetItr).begin();
			     myMinimasItr != (*myMinimaSetItr).end(); ++myMinimasItr)
			{
				confstringsFile << " " << *myMinimasItr << " ";
			}
			confstringsFile << endl << endl;
			mincounter++;
		}

		for (unsigned int i = 0; i < confstrings.size(); i++)
			confstringsFile << "String " << confstrings[i]
				<< " occurred a total of " << confstringsCounter[i] << " times" << endl;

		confstringsFile << endl << "The conformation strings"
			<< endl;

		for (std::set<string>::iterator conf_itr = myConfstrings.begin();
		     conf_itr != (myConfstrings.end()); conf_itr++)
		{
			confstringsFile << *conf_itr << endl;
		}

		confstringsFile << endl << "The flip counts"
			<< endl;

		for (unsigned int i = 0; i < flipcounts.size(); i++)
		{
			for (unsigned int j = 0; j < flipcounts[i].size(); j++)
			{
				confstringsFile << "Dihedral " << myDihedrals[i] << " State " << mycharacters[j]
					<< " flipped a total of " << (flipcounts[i])[j] << " times" << endl;
			}
		}


		confstringsFile.close();
		close();
	}

	/** This function can be used to get an accurate set of backbone dihedrals
	 *  given an accurate set of atom indices which make up the backbone.
	 *  The function can be called during initialization by uncommenting the call and recompiling.
	 *  The name of the file with atom indices must be atomset.txt as it is currently hardcoded.
	 */
	/*void OutputDihedrals::getdihedralsfromatomset(){
  
	  std::set< int > myAtomSet; // set structure to hold the atom indices from select protein components
  
	  ifstream atomSetInput(string("atomset.txt").c_str(), std::ios::in);
	  if(!atomSetInput)
	    report << error <<" Could not open atomset.txt for "<<getId()<<"."<<endr;
	  int tempatom = 0;
	  while (atomSetInput >> tempatom)
	    myAtomSet.insert(tempatom);
	  
	  ofstream selDihedralsFile(string(myFilename + ".selDih").c_str(), std::ios::out | std::ios::trunc);
	  if(!selDihedralsFile)
	    report << error <<" Can not open \'"<<myFilename<<".selDih\' for "<<getId()<<"."<<endr;
  
	  int atom1index = 0;
	  int atom2index = 0;
	  int atom3index = 0;
	  int atom4index = 0;
  
	  int atom1type = 0;
	  int atom2type = 0;
	  int atom3type = 0;
	  int atom4type = 0;
  
	  string atom1name = "";
	  string atom2name = "";
	  string atom3name = "";
	  string atom4name = "";
  
	  std::set< int >::iterator myAtomSet_itr1 = NULL;
	  std::set< int >::iterator myAtomSet_itr2 = NULL;
	  std::set< int >::iterator myAtomSet_itr3 = NULL;
	  std::set< int >::iterator myAtomSet_itr4 = NULL;
  
	  for(unsigned int i=0; i < myTopology->dihedrals.size();i++){
	    
	    atom1index = myTopology->dihedrals[i].atom1;
	    atom2index = myTopology->dihedrals[i].atom2;
	    atom3index = myTopology->dihedrals[i].atom3;
	    atom4index = myTopology->dihedrals[i].atom4;
  
	    myAtomSet_itr1 = myAtomSet.find(atom1index);
	    myAtomSet_itr2 = myAtomSet.find(atom2index);
	    myAtomSet_itr3 = myAtomSet.find(atom3index);
	    myAtomSet_itr4 = myAtomSet.find(atom4index);
  
	    atom1type = myTopology->atoms[atom1index].type;
	    atom2type = myTopology->atoms[atom2index].type;
	    atom3type = myTopology->atoms[atom3index].type;
	    atom4type = myTopology->atoms[atom4index].type;
  
	    atom1name = myTopology->atomTypes[atom1type].name;
	    atom2name = myTopology->atomTypes[atom2type].name;
	    atom3name = myTopology->atomTypes[atom3type].name;
	    atom4name = myTopology->atomTypes[atom4type].name;
  
	    // This check makes sure that only dihedrals are selected which have atoms in atomset
	    if( (myAtomSet_itr1 != myAtomSet.end()) &&
	        (myAtomSet_itr2 != myAtomSet.end()) &&
	        (myAtomSet_itr3 != myAtomSet.end()) &&
	        (myAtomSet_itr4 != myAtomSet.end()) ){
  
	        selDihedralsFile << resetiosflags(std::ios::showpoint |  std::ios::fixed | std::ios::floatfield)
	                         << i
	                         << endl;
  
  #ifdef DEBUG_OUTPUTDIHEDRALS
	        std::cout << "dihedral index: " << setw(3) << i
	                  << "  atom1: " << setw(3) << atom1index << " " << setw(5) << atom1name
	                  << "  atom2: " << setw(3) << atom2index << " " << setw(5) << atom2name
	                  << "  atom3: " << setw(3) << atom3index << " " << setw(5) << atom3name
	                  << "  atom4: " << setw(3) << atom4index << " " << setw(5) << atom4name << endl;
  #endif
	    }
	  }
	  selDihedralsFile.close();
	}*/

	// Simple function to assign dihedral increment set instead of calculating potential Maxima 
	std::vector<std::vector<Real>> OutputDihedrals::dihInc(std::vector<int> dihedralset)
	{
		std::vector<std::vector<Real>> dihIncs;
		Real increment = 0;
		dihIncs.clear();
		dihIncs.resize(dihedralset.size());

		for (unsigned int i = 0; i < dihedralset.size(); ++i)
		{
			for (unsigned int j = 1; j <= 36; j++)
			{
				increment = 2 * Constant::M_PI * j / 36;
				// std::cout << "increment: " << increment << endl; 
				(dihIncs[i]).push_back(increment);
			}
		}
		return dihIncs;
	}

	Output* OutputDihedrals::doMake(std::string&, const std::vector<Value>& values) const
	{
		return (new OutputDihedrals(values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8]));
	}

	void OutputDihedrals::getParameters(std::vector<Parameter>& parameter) const
	{
		OutputFile::getParameters(parameter);
		parameter.push_back(Parameter(keyword + "MinimalImage", Value(myMinimalImage), Text("whether the coordinates should be transformed to minimal image or not")));
		parameter.push_back(Parameter("dihedralsIndex", Value(myDihedralIndex, ConstraintValueType::NotNegative())));
		parameter.push_back(Parameter("dihedralsSet", Value(myDihedralsSet), false));
		parameter.push_back(Parameter("dihedralsSetfile", Value(myDihedralsSetfile, ConstraintValueType::NotEmpty())));
	}


	bool OutputDihedrals::adjustWithDefaultParameters(vector<Value>& values, const Configuration* config) const
	{
		if (!checkParameterTypes(values))
			return false;
		if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
			values[1] = (*config)[InputOutputfreq::keyword];
		if (config->valid(InputMinimalImage::keyword) && !values[5].valid())
			values[5] = (*config)[InputMinimalImage::keyword];
		if (!values[6].valid())
		{
			string str;
			stringstream ss(config->get(InputIntegrator::keyword).getString());
			while (ss >> str)
			{
				int index = -1;
				if (equalNocase(str, "dihedralsIndex") && ss >> index)
				{
					if (index >= 0)
					{
						values[6] = index;
						break;
					}
				}
			}
		}

		return checkParameters(values);
	}
}
