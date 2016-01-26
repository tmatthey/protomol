#include "ANM.h"


using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::fstream;


namespace ProtoMol {

    ANM::ANM(string pdbfilename, string outfilename, Real c, Real gm):pdbReader(pdbfilename), cutoff(c), gamma(gm), fname(outfilename){}

    ANM::~ANM(){}

	void ANM::connectivity_matrix() {

  		PDB pdb;

		myFile.open(fname.c_str(), std::fstream::out);

  		pdbReader.read(pdb);

  		/* Now pdb.coords contain all the co-ordinates */
  		/* Need to find those which are alpha carbons */

  		int moleculeSize = pdb.atoms.size();

  		vector<int> alphac;

  		for(int i=0;i<moleculeSize;i++) {

    		if(pdb.atoms[i].elementName == "CA") {
        		//std::cout<<i<<endl;
        		alphac.push_back(i);
    		}

  		} //end for..
  		int a_c_size = alphac.size();

  		int a_c_size_with_anisotropy = 3*a_c_size;

  		vector< vector<Real> > H(a_c_size_with_anisotropy, vector<Real> (a_c_size_with_anisotropy));

  		for(int i=0;i<a_c_size_with_anisotropy;i++) {
    		for(int j=0;j<a_c_size_with_anisotropy;j++) {
        		H[i][j] = 0;
    		}
  		}

  		for(int i=0;i<a_c_size;i++) {
    		int i_index = i*3;
    		for(int j=0;j<a_c_size;j++) {
        		int j_index = j*3;
        		//Real ssij = (pdb.coords[alphac[i]] - pdb.coords[alphac[j]]).norm();
        		if( j != i) {
            	Vector3D vv = (pdb.coords[alphac[i]] - pdb.coords[alphac[j]]);
            	Real ssij = vv.norm();
            	if( ssij < cutoff ) {

                	Real diff_x = vv.x;
                	Real diff_y = vv.y;
                	Real diff_z = vv.z;

                	H[i_index][j_index] = (gamma*diff_x*diff_x)/(ssij*ssij);
                	H[i_index+1][j_index+1] = (gamma*diff_y*diff_y)/(ssij*ssij);
                	H[i_index+2][j_index+2] = (gamma*diff_z*diff_z)/(ssij*ssij);

                	H[i_index][j_index+1] = (-1)*((gamma*diff_x*diff_y)/(ssij*ssij));
                	H[i_index][j_index+2] = (-1)*((gamma*diff_x*diff_z)/(ssij*ssij));

                	H[i_index+1][j_index] = (-1)*((gamma*diff_y*diff_x)/(ssij*ssij));
                	H[i_index+1][j_index+2] = (-1)*((gamma*diff_y*diff_z)/(ssij*ssij));

                	H[i_index+2][j_index] = (-1)*((gamma*diff_z*diff_x)/(ssij*ssij));
                	H[i_index+2][j_index+1] = (-1)*((gamma*diff_z*diff_y)/(ssij*ssij));

            	}

        	}
        	else {
            	//i == j

        	}
    	}
   
		//Real diff_x_sum = 0, diff_y_sum = 0, diff_z_sum = 0;
    	//if i == j
    	for(int j=0;j<a_c_size;j++) {
        	Real diff_x, diff_y, diff_z;
        	Vector3D vv = (pdb.coords[alphac[i]] - pdb.coords[alphac[j]]);
        	Real ssij = vv.norm();

         	diff_x = vv.x;
         	diff_y = vv.y;
         	diff_z = vv.z;

         	//sum up the contributions for all j so that j not i

         	if(j != i) {
            	H[i_index][i_index] += (gamma*diff_x*diff_x)/(ssij*ssij);
            	H[i_index+1][i_index+1] += (gamma*diff_y*diff_y)/(ssij*ssij);
            	H[i_index+2][i_index+2] += (gamma*diff_z*diff_z)/(ssij*ssij);

            	H[i_index][i_index+1] += (gamma*diff_x*diff_y)/(ssij*ssij);
            	H[i_index][i_index+2] += (gamma*diff_x*diff_z)/(ssij*ssij);

            	H[i_index+1][i_index] += (gamma*diff_y*diff_x)/(ssij*ssij);
            	H[i_index+1][i_index+2] += (gamma*diff_y*diff_z)/(ssij*ssij);

            	H[i_index+2][i_index] += (gamma*diff_z*diff_x)/(ssij*ssij);
            	H[i_index+2][i_index+1] += (gamma*diff_z*diff_y)/(ssij*ssij);
         	}

    	}

  	}

  	for(int i=0;i<a_c_size_with_anisotropy;i++) {
    	for(int j=0;j<a_c_size_with_anisotropy;j++) {
        	myFile <<H[i][j]<<" ";

    	}
    	myFile << endl;
  	}


	myFile.close();

	}

}

