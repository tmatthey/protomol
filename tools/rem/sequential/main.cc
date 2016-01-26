#include "replica.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sys/time.h>
#include <unistd.h>  // Shell commands
#include <sys/types.h> // stat
#include <sys/stat.h> // stat
#include <unistd.h> // stat

using namespace std;

//  -----------------------------------------------------------------------  //

// #define DEBUG if(0) cout
#define DEBUG if(1) cout << "dbg: "

//  -----------------------------------------------------------------------  //

void readInput( char *, string &, string &, string &, string &, string &, int &,
        int &, vector< vector<string> > & );

vector< string > parseParameters( const string & );

//  -----------------------------------------------------------------------  //

int main(int argc, char* argv[]) {

    struct timeval tv;

    // Seed the random number generator
    gettimeofday( &tv, NULL );
    srand48( tv.tv_sec );

    // Boltmann's constant
    const double k_B = 0.001987191;

    vector< Replica > replicas;

    bool first = true;
    bool even = true;

    int numExchangeSteps;

    // If no input file name was provided, then exit
    if( argc == 1 ) {
        cout << "No REUS file specified." << endl;
        return 0;
    }

    //  Make sure we can open the given input file.

    struct stat buf;
    int result = stat(argv[1], &buf);
    if (result != 0) {
      cout << "Invalid REUS file specified.  Please check the path and try again." << endl;
      exit(-1);
    }
    
    //end replace

    //  Parse the input file
    string configFileName,
           psfFileName,
           posFileName,
           parFileName,
           rootDir;

    int numSteps;

    vector< vector< string > > parameters;

    readInput( argv[1], rootDir, configFileName, psfFileName, posFileName,
            parFileName, numSteps, numExchangeSteps, parameters );


    //  Create the output file.
    ofstream oFile;

    oFile.open( string( rootDir + "exchange.txt" ).c_str());


    vector< string > v;

    //  Read in temperatures and create replicas.
    for( unsigned int i = 0; i < parameters.size(); i++ ) {

        v.clear();

        for( unsigned int j = 1; j < parameters[i].size(); j++ )
            v.push_back( parameters[i][j] );

        Replica r( parameters[i][0], rootDir, posFileName, psfFileName,
                parFileName, configFileName );

        r.setNumSteps( numSteps );
        r.setTargetTemperature( atof( parameters[i][1].c_str() ) );
        r.setParameters( v );

        replicas.push_back( r );

    }


    //  For each number of exchanges, do:
    for (int i = 0; i < numExchangeSteps; i++) {


        // Run each replica (temperature), run the MD.
        for (unsigned int j = 0; j < replicas.size(); j++)
	  replicas[j].run(); // get rid of "first"

        first = false;


        unsigned int index = ( even ) ? 0 : 1;


        // Perform replica exchanges
        while( ( index + 1 ) < replicas.size() ) {


            // Compute delta
            int m = index,
                n = index + 1;


            // Temperature Exchange
            double betaM = 1 / ( k_B * replicas[m].getTargetTemperature() ),
                   betaN = 1 / ( k_B * replicas[n].getTargetTemperature() ),
                   UM = replicas[m].getPotentialEnergy(),
                   UN = replicas[n].getPotentialEnergy(),
                   delta = ( betaN - betaM ) * ( UM - UN );


            DEBUG << "Replica(" << m << "): targetTemp = " <<
                    replicas[m].getTargetTemperature() << endl;
            DEBUG << "Replica(" << n << "): targetTemp = " <<
                    replicas[n].getTargetTemperature() << endl;


            // Probability of exchange
            double w = ( delta < 0 ) ? 1.0 : exp( -delta ),  // XXX
                   rnd = drand48();


            oFile << "Attempting to exchange replicas " << m << " and " << n <<
                    "." << endl << "  delta = " << delta << "\tprob = " << w <<
                    endl;


            if(w > rnd) {
                // Perform the exchange
                replicas[m].exchange( replicas[n] );
                oFile << "Accept Exchange" << endl << endl;
            }
            else {
                oFile << "Deny Exchange" << endl << endl;
            }

            index += 2;

        }

        even = !(even);

    }

    return(0);

}


void readInput(char* fileName,
        string& rootDir,
        string& configFileName,
        string& psfFileName,
        string& posFileName,
        string& parFileName,
        int& numSteps,
        int& numExchSteps,
        vector< vector<string> >& parameters) {
    ifstream iFile;
    iFile.open(fileName,ios::in);

    // For file input
    char line[256];
    iFile.getline(line,256);
    rootDir = string(line);
    iFile.getline(line,256);
    configFileName = string(line);
    iFile.getline(line,256);
    psfFileName = string(line);
    iFile.getline(line,256);
    posFileName = string(line);
    iFile.getline(line,256);
    parFileName = string(line);
    iFile.getline(line,256);
    numSteps = atoi(line);
    iFile.getline(line,256);
    numExchSteps = atoi(line);

    while (! iFile.eof()) {
        iFile.getline(line,256);
        if (line[0] !=  0) {
            parameters.push_back(parseParameters(string(line)));
        }
    }
    iFile.close();
}

vector<string> parseParameters(const string& parameterString) {
    // The value to be returned
    vector<string> result;

    // For parsing
    string param;
    unsigned int start = 0,
                 end = 0;

    end = parameterString.find(',');
    while (end != string::npos) {
        param = parameterString.substr(start,end-start);
        result.push_back(param);
        start = end + 1;
        end = parameterString.find(',',start);
    }

    // Parse the final parameter
    param = parameterString.substr(start);
    result.push_back(param);

    return result;
}

