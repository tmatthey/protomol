#include "Report.h"
#include "EigenvectorReader.h"
#include "typeSelection.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include "pmconstants.h"
#include "DCDTrajectoryReader.h"
#include "DCDTrajectoryWriter.h"

using namespace std;
using std::string;

using namespace ProtoMol;
using namespace ProtoMol::Report;
//_____________________________________________________________________ dcd2dcd

int main(int argc, char **argv) {
    typedef TypeSelection::Int<4>::type int32;

    ofstream myFile;
    int numSelect;

    // Parse
    if(argc < 4)
        report << quit << "usage: " << argv[0] << " <eigenvector input file> <eigenvector output file> <modes>"<< endr;

    //get number of eigs
    numSelect = atoi(argv[3]);
    if(numSelect < 1)
      report << error << "Number of modes must be greater than 0."<<endr;

    // Open
    EigenvectorInfo ei;  
    EigenvectorReader evReader;

    if(!evReader.open(argv[1]))
      report << error << "Can't open eigenvector file \'"<<argv[1]<<"\'."<<endr;
    if(!(evReader >> ei))
      report << error << "Could not parse eigenvector file \'"<<argv[1]<<endr;
    if (ei.myNumEigenvectors < 1 || ei.myNumEigenvectors < numSelect)
        report << error << "Wrong number of eigenvectors (" << ei.myNumEigenvectors << ")." << endr;

    //Process and output
    myFile.open(argv[2],ofstream::out);

	int32 vp = ei.myEigenvectorLength;
	int32 fm = numSelect;
	double ev = ei.myMaxEigenvalue;
	if (ISLITTLEENDIAN){
		swapBytes(vp);
		swapBytes(fm);
		swapBytes(ev);
	}

    myFile.write((char *)&vp, sizeof(int32));
    myFile.write((char *)&fm, sizeof(int32));
    myFile.write((char *)&ev, sizeof(double));
    int totalDoubles = numSelect * ei.myEigenvectorLength * 3;

	for(int i=0 ; i<totalDoubles; i++){
		double evec = ei.myEigenvectors[i];
		if (ISLITTLEENDIAN) swapBytes(evec);
		myFile.write((char *)&evec, sizeof(double));
	}

    //close file
    myFile.close();

    // Done
    report << hint << "Wrote eigenvector file with " << numSelect << " vector(s) of " << ei.myNumEigenvectors << ", and " << ei.myEigenvectorLength << " element(s)." << endr;  
  
    return 0;
}
