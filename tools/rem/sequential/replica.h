#ifndef REPLICA_H
#define REPLICA_h

#include <vector> // Vectors
#include <string> // Strings

using std::vector;
using std::string;

class Replica {
 public:
  // ------------------------------------ Constructors
  // Default constructor
  Replica(string name, string parentDir, string pdbFile, string psfFile, string parFile, string configFile);

  // ------------------------------------ Commands
  // Run the simulation
  void run();
  // Exchange with another replica
  void exchange(Replica other);
  // Scale the velocities to a new temperature
  void scaleToTemperature(double newTemp);

  // ------------------------------------ Get/Set Pairs
  vector<string> getParameters() const;
  void setParameters(const vector<string>& parameters);
  double getTargetTemperature() const;
  void setTargetTemperature(double targetTemp);
  int getNextStep() const;
  int getNumSteps() const;
  void setNumSteps(int numSteps);

  // ------------------------------------ Accessors
  // The directory in which the replica
  // results are stored
  string getDirectory();
  // The value of a dihedral
  double getDihedral(int dihedralIndex) const;
  // The value of a the potential energy
  double getPotentialEnergy() const;
  // The temperature
  double getTemperature();

 private:
  // ------------------------------------ Member Functions
  // Setup the configuration file for the first run
  void setupFirstConfig();
  // Setup the configuration file for subsequent runs
  void setupConfig();
  // Perform a string replacement
  bool strReplace(string &str, const string& oldStr, const string& newStr);

  // ------------------------------------ Member Variables
  string myName;
  string myParentDir;
  string myDir;
  double myTargetTemperature;
  vector<string> myParameters;
  string mySourceConfig;
  string myConfig;
  string myPDBFile;
  string myPSFFile;
  string myPARFile;
  int myNextStep;
  int myNumSteps;
  bool first_run;

  // File locations and names
  string myProtomol;
  string myTemperatureFile;
  string myDihedralFile;
  string myDihedralSetFile;
  string myEnergyFile;
  string myFinalPos;
  string myFinalVel;
};
#endif
