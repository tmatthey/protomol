#include "OutputState.h"
#include "GenericTopology.h"
#include "topologyutilities.h"
#include "OutputCache.h"
#include "inputValueDefinitions.h"
#include "XYZTrajectoryWriter.h"

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

namespace ProtoMol {
  //________________________________________________________ Output
  const string  OutputState::keyword("stateFile");

  OutputState::OutputState(): OutputFile(),
			      AObserved(false),
			      BObserved(false),
			      stateA(),
			      stateB(){}
  
  OutputState::OutputState(const string& filename, 
				   int freq, 
				   int cacheFreq, 
				   int cacheSize,
				   Real closeTime,
				   std::string statesInFile):
    OutputFile(filename,freq,cacheFreq,cacheSize,closeTime),
    AObserved(false),
    BObserved(false),
    myStatesFile(statesInFile),
    stateA(),
    stateB(),
    myVELA(NULL),
    myVELB(NULL),
    myPOSA(NULL),
    myPOSB(NULL){}

  OutputState::~OutputState(){
    //if(myDCD != NULL)
    //  delete myDCD;
    //myDCD = NULL;
  }

  void OutputState::doInitialize(){
    
    ofstream stateHeaderFile(string(myFilename + ".header").c_str(), std::ios::out | std::ios::trunc);
    if(!stateHeaderFile)
      report << error <<" Can not open \'"<<myFilename<<".header\' for "<<getId()<<"."<<endr;

    stateHeaderFile << setiosflags( std::ios::showpoint | std::ios::fixed ) 
                        << setw(14) << "Time (fs)"
			<< setw(12) << "DihIndex"
                        << setw(13) << "Val (rad)"
                        << endl;
            
    ifstream statesInput(string(myStatesFile).c_str(), std::ios::in);
    if(!statesInput)       
    report << error <<" Could not open \'"<<myStatesFile<<"\' for "<<getId()<<"."<<endr;

    string curLine = "";

    orderParam tempOP;

    report << debug(10) << "Before reading the StatesInput" << endr;
    while (statesInput >> curLine)
    {
      if (curLine == "stateA")
      {
        while ((statesInput >> curLine) && (curLine != "stateB"))
        {
          tempOP.dihID = atoi(curLine.c_str());
          statesInput >> tempOP.lowBound;
          statesInput >> tempOP.highBound;
          stateA.push_back(tempOP);
        }
      }
      else
      {
        report << error << "Error: stateA not defined" << endr;
      }
      if (curLine == "stateB")
      {
        while (statesInput >> curLine)
        {
          tempOP.dihID = atoi(curLine.c_str());
          statesInput >> tempOP.lowBound;
          statesInput >> tempOP.highBound;
          stateB.push_back(tempOP);
        }
      }
      else
      {
        report << error << "Error: stateB not defined" << endr;
      }
    }

    if(myVELA == NULL)
      myVELA = new XYZTrajectoryWriter((myFilename + ".velA").c_str());
    if(myVELA == NULL || !myVELA->open())
      report << error <<" Can not open \'"<<(myVELA!=NULL?myVELA->getFilename():"")
	     <<"\' for "<<getId()<<"."<<endr;

    if(myVELB == NULL)
      myVELB = new XYZTrajectoryWriter((myFilename + ".velB").c_str());
    if(myVELB == NULL || !myVELB->open())
      report << error <<" Can not open \'"<<(myVELB!=NULL?myVELB->getFilename():"")
	     <<"\' for "<<getId()<<"."<<endr;

    if(myPOSA == NULL)
      myPOSA = new XYZTrajectoryWriter((myFilename + ".posA").c_str());

    if(myPOSB == NULL)
      myPOSB = new XYZTrajectoryWriter((myFilename + ".posB").c_str());

  }

  // The run function outputs 
  void OutputState::doRunCached(int){

    vector< int > myDihedralsA;
    for(unsigned int i=0;i<stateA.size();++i)
    {
      myDihedralsA.push_back(stateA[i].dihID);
    }
    vector<Real> dihedralsA = myCache->dihedralPhis(myDihedralsA);

    vector< int > myDihedralsB;
    for(unsigned int i=0;i<stateB.size();++i)
    {
      myDihedralsB.push_back(stateB[i].dihID);
    }
    vector<Real> dihedralsB = myCache->dihedralPhis(myDihedralsB);

    Real tempdihedral = 1.0;
    
    if (!AObserved || !BObserved)
    {
      AObserved = true;
      int i = 0;
      for(std::vector<Real>::iterator dihedral_itr = dihedralsA.begin();
          dihedral_itr != (dihedralsA.end()); ++dihedral_itr)
      {
        tempdihedral = *dihedral_itr;
        if(tempdihedral < 0.0)
        {
          tempdihedral = 2.0*M_PI + tempdihedral;
        }
        if ((tempdihedral <= stateA[i].lowBound) || (tempdihedral >= stateA[i].highBound))
        {
          AObserved = false;
        }
        i = i + 1;
      }

      if (AObserved == true)
      {
        int i = 0;
        myBuffer << setiosflags( std::ios::showpoint | std::ios::fixed )
                 << setw(14) << setprecision(3) << "At time step " << myCache->time()
                 << ", State A Observed" << endl;
        for(std::vector<Real>::iterator dihedral_itr = dihedralsA.begin();
            dihedral_itr != (dihedralsA.end()); ++dihedral_itr)
        {
          tempdihedral = *dihedral_itr;
          if(tempdihedral < 0.0)
          {
            tempdihedral = 2.0*M_PI + tempdihedral;
          }
          myBuffer << " " << stateA[i].dihID << "  low bound:" << stateA[i].lowBound << " val:" << tempdihedral
                   << " high bound:" << stateA[i].highBound << endl;
          i = i+1;
        }
        myVELA->open();
        if(!myVELA->write(*myVelocities,myTopology->atoms,myTopology->atomTypes))
          report << error << "Could not write "<<getId()<<" \'"<<myVELA->getFilename()<<"\'."<<endr;
        myPOSA->open();
        if(!myPOSA->write(*myPositions,myTopology->atoms,myTopology->atomTypes))
          report << error << "Could not write "<<getId()<<" \'"<<myPOSA->getFilename()<<"\'."<<endr;
      }
    }

    if (!BObserved)
    {
      BObserved = true;
      int i = 0;
      for(std::vector<Real>::iterator dihedral_itr = dihedralsB.begin();
          dihedral_itr != (dihedralsB.end()); ++dihedral_itr)
      {
        tempdihedral = *dihedral_itr;
        if(tempdihedral < 0.0)
        {
          tempdihedral = 2.0*M_PI + tempdihedral;
        }

        if ((tempdihedral <= stateB[i].lowBound) || (tempdihedral >= stateB[i].highBound))
        {
          BObserved = false;
        }
        i = i + 1;
      }
      if (BObserved == true)
      {
        int i = 0;
        myBuffer << setiosflags( std::ios::showpoint | std::ios::fixed )
                 << setw(14) << setprecision(3) << "At time step " << myCache->time()
                 << ", State B Observed" << endl;        
        for(std::vector<Real>::iterator dihedral_itr = dihedralsB.begin();
            dihedral_itr != (dihedralsB.end()); ++dihedral_itr)
        {
          tempdihedral = *dihedral_itr;
          if(tempdihedral < 0.0)
          {
            tempdihedral = 2.0*M_PI + tempdihedral;
          }
          myBuffer << " " << stateB[i].dihID << " low bound:" << stateB[i].lowBound << " val:" << tempdihedral
                   << " high bound:" << stateB[i].highBound << endl;
          i = i + 1;
        }
        if(!myVELB->write(*myVelocities,myTopology->atoms,myTopology->atomTypes))
          report << error << "Could not write "<<getId()<<" \'"<<myVELA->getFilename()<<"\'."<<endr;
        if(!myPOSB->write(*myPositions,myTopology->atoms,myTopology->atomTypes))
          report << error << "Could not write "<<getId()<<" \'"<<myPOSB->getFilename()<<"\'."<<endr;
      }
    }
  }

  void OutputState::doFinalize(int){

    if (!AObserved)
      myBuffer << "State A was not Observed" << endl;

    if (!BObserved)
      myBuffer << "State B was not Observed" << endl;

  }
    
  Output* OutputState::doMake(std::string& , const std::vector<Value>& values) const{
    return (new OutputState(values[0],values[1],values[2],values[3],values[4],values[5]));
  }

  void OutputState::getParameters(std::vector<Parameter> &parameter) const{
    OutputFile::getParameters(parameter);
    parameter.push_back(Parameter("statesInputFile",Value(myStatesFile,ConstraintValueType::NotEmpty())));
  }


  bool OutputState::adjustWithDefaultParameters(vector<Value>& values, const Configuration* config) const{
    if(!checkParameterTypes(values))
      return false;
    if(config->valid(InputOutputfreq::keyword) && !values[1].valid())
      values[1] = (*config)[InputOutputfreq::keyword];

    return checkParameters(values);
  }

}
