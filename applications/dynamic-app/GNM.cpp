#include "GNM.h"


using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::fstream;


namespace ProtoMol {

	GNM::GNM(string pdbfilename, string outfilename, Real c, Real gm):pdbReader(pdbfilename), cutoff(c), gamma(gm), fname(outfilename){}

	GNM::~GNM(){}

	void GNM::connectivity_matrix()  {

		PDB pdb;

		myFile.open(fname.c_str(), std::fstream::out);	
	
		acFile.open("ac_coords.out",std::fstream::out);

		pdbReader.read(pdb);

	 	/* Now pdb.coords contain all the co-ordinates */
  		/* Need to find those which are alpha carbons */

  		unsigned int moleculeSize = pdb.atoms.size();

  		vector<int> alpha_carbons;

  		for(unsigned int i=0;i<moleculeSize;i++) {

    		if(pdb.atoms[i].elementName == "CA") {
        		std::cout<<i<<endl;
        		alpha_carbons.push_back(i);
    		}

  		} //end for..

	  unsigned int a_c_size = alpha_carbons.size();

  	  vector< vector<int> > Gamma ( a_c_size,  vector<int> (a_c_size));

	  for(unsigned int i=0;i<a_c_size;i++) {
	    int springs = 0;
    		for(unsigned int j=0;j<a_c_size;j++) {
        		if( j != i) {
        			if( (pdb.coords[alpha_carbons[i]] - pdb.coords[alpha_carbons[j]]).norm()<cutoff ) {
            		//std::cout<<"Alpha carbon "<<alpha_carbons[j]<<"connected to alpha carbon ";
            		//std::cout<<alpha_carbons[i]<<" in the Gaussian network model"<<endl;
            		Gamma[i][j] = -1;
            		springs++;
        		}
        		else {
            		Gamma[i][j] = 0;
        		}
        	}
    	}
    	Gamma[i][i] = springs;

  	}
  for(unsigned int i=0;i<a_c_size;i++) {
    for(unsigned int j=0;j<a_c_size;j++) {

        myFile << Gamma[i][j]<<" ";
    }
    myFile << endl;
  }

  for(unsigned int i=0;i<a_c_size;i++) {
  	for(unsigned int j=0;j<3;j++) {
		acFile<<pdb.coords[alpha_carbons[i]][j]<<" ";
	}
	acFile<<endl;
  }	

	myFile.close();	
	acFile.close();

	}

}
