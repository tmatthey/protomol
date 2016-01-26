#include "replica.h"
#include <vector>    // Vectors
#include <string>    // Strings
#include <unistd.h>  // Shell commands unlink
#include <fstream>   // File i/o
#include <iostream>  // Standard i/o
#include <iomanip>   // i/o manipulation
#include <cmath>     // Math functions

using namespace std;


// ------------------------------------ Constructors
// Default constructor
Replica::Replica( string name, string parentDir, string pdbFile, string psfFile,
        string parFile, string configFile ) {

    myName = name;
    myParentDir = parentDir;
    myDir = parentDir + name + "/";
    myNextStep = 0;

    // Create a directory for the replica results
    string sysCommand;

    sysCommand = "/bin/rm -f -r " + myDir + "; ";
    sysCommand += "/bin/mkdir " + myDir + "; ";

    // Copy the PDB, PSF, and PAR files to the new directory
    myPDBFile = pdbFile;
    myPSFFile = psfFile;
    myPARFile = parFile;

    // Move the config file
    mySourceConfig = configFile;

    myConfig = string( "c.config" );

    sysCommand += "/bin/cp " + parentDir + pdbFile + " " + parentDir + psfFile +
            " " + parentDir + parFile + " " + parentDir + configFile + " " +
            myDir + "; ";

    // Set the path to protomol
    myProtomol = "/afs/nd.edu/user8/shampton/research/biocomp/protomol/applications/protomol-app/protomol";

    // Set the names of the output files
    myTemperatureFile = "temperature.out";
    myDihedralFile = "dihedral.out";

    sysCommand += "/bin/cp " + parentDir + "dihedral.set " + myDir;
    system(sysCommand.c_str());

    myDihedralSetFile = "dihedral.set";
    myEnergyFile = "energy.out";

    myFinalPos = "finalpos.xyz";
    myFinalVel = "finalvel.xyz";

    first_run = true;

}

// ------------------------------------ Commands
// Run the simulation
void Replica::run() {

    // Set up the configuration file
    if (first_run == true) {
        setupFirstConfig();
    }
    else {
        setupConfig();
    }

    string sysCommand;
    sysCommand = myProtomol + " " + myDir + "run.config; ";

    myNextStep = myNextStep + myNumSteps;

    if( first_run == false ) {

        // Combine the temperature output files
        sysCommand += "perl -pi -e '$_ = \"\" if ($. == 1);' " + myDir +
                myTemperatureFile + ".cur; ";

        sysCommand += "/bin/cat " + myDir + myTemperatureFile + " " + myDir +
                myTemperatureFile + ".cur > " + myDir + "cat.tmp; ";

        sysCommand += "/bin/mv " + myDir + "cat.tmp " + myDir +
                myTemperatureFile + "; ";


        // Combine the dihedral output files
        sysCommand += "perl -pi -e '$_ = \"\" if ($. == 1);' " + myDir +
                myDihedralFile + ".cur; ";

        sysCommand += "/bin/cat " + myDir + myDihedralFile + " " + myDir +
                myDihedralFile + ".cur > " + myDir + "cat.tmp; ";

        sysCommand += "/bin/mv " + myDir + "cat.tmp " + myDir + myDihedralFile + "; ";

        // Combine the energy output files
        sysCommand += "perl -pi -e '$_ = \"\" if ($. == 1);' " + myDir +
                myEnergyFile + ".cur; ";

        sysCommand += "/bin/cat " + myDir + myEnergyFile + " " + myDir +
                myEnergyFile + ".cur > " + myDir + "cat.tmp; ";

        sysCommand += "/bin/mv " + myDir + "cat.tmp " + myDir + myEnergyFile + "; ";

        //  Remove temporary files.
        sysCommand += "/bin/rm -f " + myDir + "cat.tmp " + myDir + myEnergyFile +
                ".cur* " + myDir + myDihedralFile + ".cur* " + myDir +
                myTemperatureFile + ".cur*; ";

    }
    system(sysCommand.c_str());
    first_run = false;

}

// Exchange with another replica
void Replica::exchange(Replica other) {

    // Scale the replicas to the required temperature
    scaleToTemperature(other.getTargetTemperature());
    other.scaleToTemperature(getTargetTemperature());

    string sysCommand;

    // Exchange the velocity files
    sysCommand = "/bin/mv " + myDir + myFinalVel + ".scaled " + other.myDir +
            myFinalVel + "; ";

    sysCommand += "/bin/mv " + other.myDir + myFinalVel + ".scaled " + myDir +
            myFinalVel + "; ";

    // Exchange the position files
    sysCommand += "/bin/cp " + myDir + myFinalPos + " " + myDir + "pos.tmp; ";

    sysCommand += "/bin/mv " + other.myDir + myFinalPos + " " + myDir +
            myFinalPos + "; ";

    sysCommand += "/bin/mv " + myDir + "pos.tmp " + other.myDir + myFinalPos;

    system(sysCommand.c_str());
    unlink((myDir + string("pos.tmp")).c_str());

}

// ------------------------------------ Get/Set Pairs
vector<string> Replica::getParameters() const {
    return myParameters;
}

void Replica::setParameters(const vector<string>& parameters) {

    myParameters = parameters;

    // Perform parameter replacement in the config file
    ofstream oFile;

    string fileName = myDir + "tmp.config";
    oFile.open(fileName.c_str());
    ifstream iFile;
    fileName = myDir + mySourceConfig;
    iFile.open(fileName.c_str());

    // Read the file line by line
    char line[2048];
    string sLine;
    while (! iFile.eof()) {
        iFile.getline(line,2048);
        sLine = string(line);

        int tmp1, tmp2;
        char *parSuffix;
        string par;
        for (unsigned int i = 0; i < parameters.size(); i++) {
            parSuffix = fcvt((float) i,0,&tmp1,&tmp2);
            par = "PAR" + string(parSuffix);
            strReplace(sLine,par,parameters[i]);
        }
        oFile << sLine << endl;
    }
    iFile.close();
    oFile.close();

    string sysCommand = "/bin/mv " + myDir + "tmp.config " + myDir + myConfig;
    system(sysCommand.c_str());
}

double Replica::getTargetTemperature() const {
    return myTargetTemperature;
}

void Replica::setTargetTemperature (double targetTemp) {
    myTargetTemperature = targetTemp;
}

int Replica::getNextStep() const {
    return myNextStep;
}

int Replica::getNumSteps() const {
    return myNumSteps;
}

void Replica::setNumSteps(int numSteps) {
    myNumSteps = numSteps;
}

// ------------------------------------ Accessors
// The directory in which the replica
// results are stored
string Replica::getDirectory() {
    return myDir;
}

// The value of a dihedral
double Replica::getDihedral(int dihedralIndex) const {
    // The value to be returned
    double result;

    // Use awk to extract the dihedral
    string sysCommand;
    sysCommand = "tail -n 1 " + myDir + myDihedralFile + " > " + myDir + "tail.tmp";
    system(sysCommand.c_str());

    // Read from the file created by tail
    ifstream iFile;
    string fileName = myDir + "tail.tmp";
    iFile.open(fileName.c_str());
    double time, value, energy;
    int index;
    iFile >> time;
    iFile >> index >> value >> energy;
    while ((! iFile.eof()) && (index != dihedralIndex)) {
        iFile >> index >> value >> energy;
    }
    iFile.close();

    if (index == dihedralIndex) {
        result = value;
    }
    else {
        cout << "ERROR: Requested dihedral is not included in the dihedral output file." << endl;
        result = 0.0;
    }

    unlink((myDir + string("tail.tmp")).c_str());

    return result;
}

// The value of the potential energy
double Replica::getPotentialEnergy() const {
    // The value to be returned
    double result;

    // Use awk to extract the potential
    string sysCommand;
    sysCommand = "tail -n 1 " + myDir + myEnergyFile + " > " + myDir + "tail.tmp";
    system(sysCommand.c_str());

    // Read from the file created by tail
    ifstream iFile;
    string fileName = myDir + "tail.tmp";
    iFile.open(fileName.c_str());
    double time;
    iFile >> time >> result;

    unlink((myDir + string("tail.tmp")).c_str());

    return result;
}

// The temperature
double Replica::getTemperature() {

    // The value to be returned
    double result;

    // Use awk to extract the temperature
    string sysCommand;

    sysCommand = "tail -n 1 " + myDir + myTemperatureFile + " > " + myDir +
            "tail.tmp; ";

    sysCommand += "/bin/awk '{print $2}' " + myDir + "tail.tmp > " + myDir +
            "awk.tmp";
    system(sysCommand.c_str());

    // Read from the file created by awk
    ifstream iFile;

    string fileName = myDir + "awk.tmp";
    iFile.open(fileName.c_str());

    iFile >> result;
    iFile.close();

    unlink((myDir + string("tail.tmp")).c_str());
    unlink((myDir + string("awk.tmp")).c_str());

    return result;

}

// ------------------------------------ Member Functions
// Setup the configuration file for the first run
void Replica::setupFirstConfig() {

    // Create the new configuration file
    string sysCommand;
    string fileName = myDir + "tmp.config";

    ofstream oFile;
    oFile.open(fileName.c_str());

    // Output the required lines
    // Timestep parameters
    oFile << "firstStep " << myNextStep << endl 
          << "numSteps " << myNumSteps << endl << endl

    // Initial temperature
          << "temperature " << myTargetTemperature << endl << endl

    // Temperature output
          << "tempFile " << myTemperatureFile << endl 
          << "temperatureSeparateWater true" << endl << endl

    // Dihedral output
          << "doDihedralsFile true" << endl
          << "dihedralsFile " << myDihedralFile << endl
          << "dihedralsIndex 0" << endl
          << "dihedralsSet true" << endl
          << "dihedralsSetFile " << myDihedralSetFile << endl << endl

    // Energy output
          << "allEnergiesFile " << myEnergyFile << endl << endl

    // Final position output
          << "finXYZPosFile " << myFinalPos << endl << endl

    // Final velocity output
          << "finXYZVelFile " << myFinalVel << endl << endl

    // Initial coordinates
          << "POSFile " << myPDBFile << endl << endl

    // PSF and PAR files
          << "PSFFile " << myPSFFile << endl
          << "PARFile " << myPARFile << endl << endl;

    oFile.close();

    // Create the config file
    sysCommand = "/bin/cat " + fileName + " " + myDir + myConfig + " > " + myDir
            + "run.config";
    system( sysCommand.c_str() );

    unlink(fileName.c_str());

}

// Setup the configuration file for the subsequent runs
void Replica::setupConfig() {

    // Create the new configuration file
    string sysCommand;
    string fileName = myDir + "tmp.config";

    ofstream oFile;
    oFile.open(fileName.c_str());

    // Output the required lines
    // Timestep parameters
    oFile << "firstStep " << myNextStep << endl
          << "numSteps " << myNumSteps << endl

    // Temperature output
          << "tempFile " << myTemperatureFile << ".cur" << endl
          << "temperatureSeparateWater true" << endl

    // Dihedral output
          << "doDihedralsFile true" << endl
          << "dihedralsFile " << myDihedralFile << ".cur" << endl
          << "dihedralsIndex 0" << endl
          << "dihedralsSet true" << endl
          << "dihedralsSetFile " << myDihedralSetFile << endl

    // Energy output
          << "allEnergiesFile " << myEnergyFile << ".cur" << endl

    // Final position output
          << "finXYZPosFile " << myFinalPos << endl

    // Final velocity output
          << "finXYZVelFile " << myFinalVel << endl

    // Initial coordinates
          << "POSFile " << myFinalPos << endl

    // Initial velocities
          << "VELFile " << myFinalVel << endl

    // PSF and PAR files
          << "PSFFile " << myPSFFile << endl
          << "PARFile " << myPARFile << endl;

    oFile.close();

    // Create the config file
    sysCommand = "/bin/cat " + fileName + " " + myDir + myConfig + " > " + myDir
            + "run.config";
    system(sysCommand.c_str());

    unlink(fileName.c_str());
}

// Scale the velocities to a new temperature
void Replica::scaleToTemperature(double newTemp) {
    // Determine the current temperature
    double curTemp = getTemperature();

    // The velocity scaling factor
    double factor = sqrt(newTemp / curTemp);

    // First, remove any comments from the velocity XYZ file
    string fileName = myDir + myFinalVel;
    string sysCommand;
    sysCommand = "/bin/cp " + fileName + " " + myDir + "cp.tmp; ";
    sysCommand += "sed -i -e '/\\!.*/d' " + myDir + "cp.tmp";
    system(sysCommand.c_str());

    // Read from the file
    fileName = fileName + ".scaled";
    ofstream oFile;
    oFile.open(fileName.c_str());
    fileName = myDir + "cp.tmp";
    ifstream iFile;
    iFile.open(fileName.c_str());

    int numAtoms;
    double x, y, z;
    double scaled_x, scaled_y, scaled_z;
    string groupName;
    iFile >> numAtoms;

    oFile << numAtoms << endl << "!REUS - velocities scaled from " << curTemp <<
            " to " << newTemp << "K" << endl;

    for (int i = 0; i < numAtoms; i++) {
        iFile >> groupName >> x >> y >> z;
        scaled_x = factor * x;
        scaled_y = factor * y;
        scaled_z = factor * z;
        oFile << setprecision(15) << groupName << "\t";
        oFile.width(24);
        oFile << scaled_x;
        oFile.width(24);
        oFile << scaled_y;
        oFile.width(24);
        oFile << scaled_z;
        oFile << endl;
    }

    iFile.close();
    oFile.close();

    // Remove the temporary file
    unlink((myDir + string("cp.tmp")).c_str());
}

// Perform a string replacement
bool Replica::strReplace(string &str, const string& oldStr, const string& newStr) {
    bool result = false;
    unsigned int start = 0;
    int oldLen = oldStr.length();
    int newLen = newStr.length();

    start = str.find(oldStr);
    while (start != string::npos) {
        result = true;
        str.replace(start,oldLen,newStr);
        start = str.find(oldStr,start + newLen);
    }

    return result;
}
