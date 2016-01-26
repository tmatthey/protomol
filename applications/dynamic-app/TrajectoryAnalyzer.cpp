#include "TrajectoryAnalyzer.h"

using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::fstream;


namespace ProtoMol {

    //Implementation of TrajectoryAnalyzer class
    TrajectoryAnalyzer::TrajectoryAnalyzer(string dcdfilenames, string
            pdbfilename, string out_filename):fname(out_filename), dcdFileSetName(dcdfilenames), dcdReader(), pdbReader(pdbfilename){}

	TrajectoryAnalyzer::~TrajectoryAnalyzer(){}

	void TrajectoryAnalyzer::getDcdFiles() {

		string s;

		dcdFileSet.open(dcdFileSetName.c_str(), std::fstream::in);

		while( dcdFileSet >> s ) {

			dcdFiles.push_back(s);
		}

		//cout<<dcdFiles.size()<<" dcd files"<<endl;
		//cout<<dcdFiles[0];	

	}

	void TrajectoryAnalyzer::copyXYZ(XYZ& to, XYZ from) {

		copy(from.coords.begin(), from.coords.end(), to.coords.begin());
		copy(from.names.begin(), from.names.end(), to.names.begin());

	}	

	void TrajectoryAnalyzer::bond_time_correlation(int natom, int hatom, int records, int fr) {
    	XYZ xyz;

    	vector<Vector3D> traj_N;
    	vector<Vector3D> traj_H;

    	int count = 0;

		getDcdFiles();
		
		unsigned int num_dcd = 0;

		while (count < records){
			if( num_dcd < dcdFiles.size() ) {
				dcdReader.open(dcdFiles[num_dcd]);
				cout<<" Opened dcd file "<<dcdFiles[num_dcd]<<endl;
				num_dcd++;
				if ( num_dcd > 1 ) dcdReader >> xyz; /* Ignore the first frame if its not 1st dcd file of the sequence */
			}
			else {
				cout<<"Fatal error: Too many records requested "<< endl;
				break;
			}
        	while ( dcdReader >> xyz ) {
            	Vector3D vv1 = xyz.coords[natom-1];
            	Vector3D vv2 = xyz.coords[hatom-1];
            	traj_N.push_back(vv1);
            	traj_H.push_back(vv2);
            	count++;
            	if(count == records) break;
        	}
			dcdReader.close();

		}

        if(traj_N.size() != traj_H.size()) {
            cout<<"Fatal error:size mismatch between N and H coordinate vectors"<<endl;
            return;
        }

        unsigned int num_of_frames = traj_N.size();

        cout<< " Total number of frames = "<<num_of_frames<<endl;

        Vector3D nh_bond_total;

        vector<Vector3D> nh_bond_over_trajectory;
       /* Calculating the N-H bond vectors from each frame and storing it in a Vector3DBlock */
        /* nh_bond_total computes the sum of the N-H bond vectors from each frame             */

        for(unsigned int i=0;i<num_of_frames;i++) {
            Vector3D Bond = traj_H[i] - traj_N[i];
            nh_bond_over_trajectory.push_back(Bond);
            nh_bond_total += Bond;
        }

        /* Calculating the average bond bector        */

        nh_bond_total.x /= num_of_frames;
        nh_bond_total.y /= num_of_frames;
        nh_bond_total.z /= num_of_frames;


        /* The difference of N-H bond vector for each frame from the average   */

        /*for(unsigned int i=0;i<num_of_frames;i++) {
            nh_bond_over_trajectory[i] -= nh_bond_total;
        }*/

        int ac_steps = (int)((num_of_frames-1)/fr);

        //cout<<"Ac steps = "<<ac_steps<<endl;

        myFile.open(fname.c_str(), std::fstream::out);

        //Real Ct;
		Real Gtao;

        for(int i=0;i<ac_steps;i++) {
            Gtao = 0;
            for(int j=0;j<static_cast<int>(num_of_frames-i);j++) {
				Vector3D vj_unit, vji_unit;
				Real vj_norm, vji_norm;

				vj_norm = nh_bond_over_trajectory[j].norm();	
				vji_norm = nh_bond_over_trajectory[j+i].norm();
				
				vj_unit = nh_bond_over_trajectory[j]/vj_norm;
				vji_unit = nh_bond_over_trajectory[j+i]/vji_norm;

				Real cosine_tao = vj_unit.dot(vji_unit)/(vj_unit.norm()*vji_unit.norm());

                /*Ct += (nh_bond_over_trajectory[j].x*nh_bond_over_trajectory[j+i].x);
                Ct += (nh_bond_over_trajectory[j].y*nh_bond_over_trajectory[j+i].y);
                Ct += (nh_bond_over_trajectory[j].z*nh_bond_over_trajectory[j+i].z);*/
				Gtao += 0.2*(0.5*(3*pow(cosine_tao,2)-1));
            }
            //Ct /= (num_of_frames-i);
			Gtao /= (num_of_frames-i);	
            //myFile<<i*1<<" "<<Ct<<endl;
            myFile<<i*1<<" "<<Gtao<<endl;

        }

        myFile.close();

	}


	void TrajectoryAnalyzer::collective_mode_analysis(int records)  {

		getDcdFiles();				
	
		/* Find out the indices of alpha carbons */
		pdbReader.read(pdb);


        vector<int> alphac;
                                                                                                                                             
        for(unsigned int i=0;i<pdb.atoms.size();i++) {
            if(pdb.atoms[i].elementName == "CA") {
                /*alpha carbon, record the index */
                alphac.push_back(i);
            }
        }
                                                                                                                                             
        cout<<"ALpha carbons "<<alphac.size()<<endl;

		int nn = alphac.size();

		myFile.open(fname.c_str(), std::fstream::out); /* opening output file */

		int count = 0;

		int num_dcd = 0;

		XYZ frame1, frame2;

		while ( count < records ) {

			if ( num_dcd >= static_cast<int>(dcdFiles.size())) {
				cout<<"Fatal error: number of records exceeded"<<endl;
				break;
			}
			else {
				dcdReader.open(dcdFiles[num_dcd]);
				num_dcd++;
				//cout<<num_dcd<<endl;
			}

			dcdReader >> frame1;
		
			if ( num_dcd == 1) count++;

			//cout<<"Read first frame"<<endl;

			while( dcdReader >> frame2) {
				//cout<<"Read next frame"<<endl;
				if(static_cast<int>(frame1.coords.size()) != static_cast<int>(frame2.coords.size())) {
					cout<<"Fatal error:corrupted frame"<<endl;
					return;
				}
				for(int i=0;i<nn;i++) {
					//cout<<"for "<<endl;
					Vector3D Z = frame2.coords[alphac[i]] - frame1.coords[alphac[i]];
					//cout<<Z.x<<" "<<Z.y<<" "<<Z.z<<endl;
					myFile<<Z.x;
					myFile<<" ";
			
					myFile<<Z.y;
					myFile<<" ";

					myFile<<Z.z;
					myFile<<" ";
				}
				myFile << endl;
				count++;
				if(count == records) break;
				copyXYZ(frame1, frame2); /* frame 2 is being copied into frame 1*/
			}

			dcdReader.close();

		}

		myFile.close();	
	}
			

} 
