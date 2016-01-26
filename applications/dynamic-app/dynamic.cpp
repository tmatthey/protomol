/* This application allows us to compute various dynamic information ( e. g. TCF, */
/* collective modes of vibrations, autocorrelation fn of vibrations of residues   */
/* etc. ) from MD simulation, GNM and ANM.										  */

#include<cstdlib>
#include "TrajectoryAnalyzer.h"
#include "GNM.h"
#include "ANM.h"

using namespace ProtoMol;
using namespace ProtoMol::Report;
using std::string;
using std::cout;
using std::endl;
using std::cin;

int main(int argc, char *argv[])
{

	string datatype(argv[1]);

	string analysis_type(argv[2]);


	/* 2nd parameter is the type of analysis */		

	if(datatype == "MD") {
		cout<<"MD simulation"<<endl;

		/* Input dcd filename */
		string dcdfilename(argv[3]);

		/*Input pdb filename */
		string pdbfilename(argv[4]);

		/*output filename */
		string outfilename(argv[5]);			

		TrajectoryAnalyzer mdData(dcdfilename, pdbfilename, outfilename);

		int records = atoi(argv[6]);	

		if( analysis_type == "TCF") {
			cout<<"TCF calculation"<<endl;
		
			/* N atom index */
			int natom = atoi(argv[7]);	
			int hatom = atoi(argv[8]);	
			int frequency = atoi(argv[9]);
			mdData.bond_time_correlation(natom, hatom, records, frequency);

		}
		else if (analysis_type == "CMODE") {
			cout<<"Collective mode analysis"<<endl;
			mdData.collective_mode_analysis(records);

		}
		else {
			cout<<"Fatal error:unknown analysis type"<<endl;
			return(-2);
		}
	}	
	else if ( datatype == "GNM") {
		cout<<"GNM"<<endl;

		/*Input pdb filename */
		string pdbfilename(argv[3]);

		/* outputfilename */
		string outfilename(argv[4]);

		/* cutoff */
		Real c = atof(argv[5]);

		/* gamma */
		Real gamma = atof(argv[6]);

		GNM gnm(pdbfilename, outfilename, c, gamma);
		
		gnm.connectivity_matrix();	

	}
	else if( datatype == "ANM" ) {
		cout<<"ANM"<<endl;

		/*Input pdb filename */
		string pdbfilename(argv[3]);

		/* outputfilename */
		string outfilename(argv[4]);

		/* cutoff */
		Real c = atof(argv[5]);

		/* gamma */
		Real gamma = atof(argv[6]);
		
		ANM anm(pdbfilename, outfilename, c, gamma);
		
		anm.connectivity_matrix();

	}
	else {
		cout<<"Unknown type"<<endl;
		return (-1);
	}

	/* 2nd parameter is the type of analysis */		
	
    return 0;
}
