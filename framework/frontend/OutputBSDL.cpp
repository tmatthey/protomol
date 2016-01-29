#include <fstream>
#include <map>
#include "OutputBSDL.h"
#include "Configuration.h"
#include "OutputCache.h"
#include "stringutilities.h"
#include "GenericTopology.h"
#include "inputValueDefinitions.h"
using namespace ProtoMol::Report;

using std::string;
using std::vector;
using std::ofstream;
using std::map;
using std::endl;
using std::stringstream;
namespace ProtoMol {
  //________________________________________________________ OutputBSDL
  const string  OutputBSDL::keyword("BSDLFile");

  OutputBSDL::OutputBSDL():Output(),myFilename(""),myMinimalImage(false),myIncludeView(true),myWater(true),myCounter(0),myScale(1.0),myMolecularDistance(false){}

  OutputBSDL::OutputBSDL(const string& filename, int freq, bool minimal, bool include, bool water):Output(freq),myFilename(filename),myMinimalImage(minimal),myIncludeView(include),myWater(water),myCounter(0),myScale(1.0),myMolecularDistance(false){}

  OutputBSDL::~OutputBSDL(){
  }

  void OutputBSDL::doInitialize(){


    // File counter
    myCounter = 0;


    // Scale to reasonable coordinates
    myScale = 1.0;
    Vector3D minb,maxb;
    if(myTopology->getVolume() < Constant::REAL_INFINITY)
      myTopology->getBoundaryConditionsBox(minb,maxb);
    else
      myTopology->getBoundingbox(*myPositions,minb,maxb);
    Real d = (maxb-minb).normSquared();
    if(d  < Constant::EPSILON){
      myScale = 1.0;      
    }
    else {
      while (d*myScale > 1e6){
	myScale /= 10.0;
      }
      while (d*myScale < 1.0){
	myScale *= 10.0;
      }
    }


    // Suppress bonds if their wrapped around 
    myMolecularDistance = myTopology->checkMoleculePairDistances(*myPositions);


    // Build a look up table to assign the colors
    myColortable.clear();
    int color = 0;
    for(unsigned long i = 0;i<myPositions->size();i++){
      if(!myWater && myTopology->molecules[myTopology->atoms[i].molecule].water)
	continue;

      string aname = myTopology->atomTypes[myTopology->atoms[i].type].symbolName;
      if(myColortable.find(aname) == myColortable.end())
	myColortable[aname] = color;
      color = (color+1)%8;
    }

    report << hint << "OutputBDSL : molecular distance : "<<(bool)myMolecularDistance<<", scaling : "<<myScale<<endr;

  }

  void OutputBSDL::doRun(int){
    const Vector3DBlock* pos = (myMinimalImage ? myCache->minimalPositions() : myPositions);

    string viewstr("");
    if(myCounter <= 0 || myIncludeView){
      viewstr = view(pos);
      if(myCounter <= 0){
	ofstream inc(string(myFilename + ".inc").c_str());
	if(!inc)
	  report << error << "Can't write include BSDL file \'"+myFilename + ".inc\'"<<endr;

	int l = static_cast<int>(pow(static_cast<Real>(myTopology->atoms.size()*1.5),1.0/3.0)+1);
	if (l < 5)
	  l = 5;

	Real u = 0.25*myScale;
	Real v = 0.1*myScale;

	inc <<"//Auto generated include BSDL file by ProtoMol.\n"
	    <<"using 3D;\n"
	    <<"pointLight (0.8, [1,1,1]) {position [10000,5000,10000];};\n"
	    <<"ambientLight (0.3, [1,1,1]);\n"
	    <<"define b whitted {diffuse [1.0, 1.0, 1.0];ambient [0.5, 0.5, 0.5];};\n"
	    <<"define d whitted {diffuse [0.7, 0.7, 0.7];ambient [0.5, 0.5, 0.5];};\n"
	    <<"define c0 whitted {diffuse [1.0, 1.0, 1.0];ambient [0.5, 0.5, 0.5];};\n"
	    <<"define c1 whitted {diffuse [1.0, 0.1, 0.1];ambient [0.5, 0.1, 0.1];};\n"
	    <<"define c2 whitted {diffuse [0.1, 1.0, 0.1];ambient [0.1, 0.5, 0.1];};\n"
	    <<"define c3 whitted {diffuse [0.1, 0.1, 1.0];ambient [0.1, 0.1, 0.5];};\n"
	    <<"define c4 whitted {diffuse [0.1, 1.0, 1.0];ambient [0.1, 0.5, 0.5];};\n"
	    <<"define c5 whitted {diffuse [1.0, 1.0, 0.1];ambient [0.5, 0.5, 0.1];};\n"
	    <<"define c6 whitted {diffuse [1.0, 0.1, 1.0];ambient [0.5, 0.1, 0.5];};\n"
	    <<"define c7 whitted {diffuse [0.1, 0.1, 0.1];ambient [0.1, 0.1, 0.1];};\n"
	    <<"#define s sphere\n"
	    <<"#define c cylinder\n"
	    <<"#define u "<<u<<"\n"
 	    <<"#define v "<<v<<"\n"	  
	    <<"grid ["<<l<<","<<l<<","<<l<<"]{\n";

	if(!myIncludeView){
	  inc << viewstr;
	  viewstr = "";
	}

      }
    } 


    char tmp[16];
    sprintf_s(tmp, "%06d", myCounter);
    ofstream out(string( myFilename +"."+tmp+".bsdl3").c_str());
    if(!out)
      report << error << "Can't write BSDL file \'"+myFilename +"."+tmp+".bsdl3\'"<<endr;
    
    out << "//Auto generated BSDL file by ProtoMol. T="<<myCache->time()<<"[fs]\n";
    out << "#include \""<<myFilename<<".inc\"\n";
    out << viewstr;
    for(unsigned long i = 0;i<myPositions->size();i++){
      if(!myWater && myTopology->molecules[myTopology->atoms[i].molecule].water)
	continue;
      Vector3D pi = (*pos)[i]*myScale;
      out <<"s(u,["<< pi.x <<","<< pi.y<<","<<pi.z<<"])"
	  <<"{c"<<myColortable[myTopology->atomTypes[myTopology->atoms[i].type].symbolName]<<";}\n";
    }
    for(unsigned long k = 0;k<myTopology->bonds.size();k++){
      unsigned int i = myTopology->bonds[k].atom1;
      unsigned int j = myTopology->bonds[k].atom2;
      if(!myWater && myTopology->molecules[myTopology->atoms[i].molecule].water)
	continue;
      Vector3D pi = (*pos)[i]*myScale;
      Vector3D pj = (*pos)[j]*myScale;
      if(myMolecularDistance || (pi-pj).norm() < myTopology->bonds[k].restLength*2.0*myScale){
	out <<"c(v,["<< pi.x <<","<< pi.y<<","<<pi.z<<"],"
	    <<"["<< pj.x <<","<< pj.y<<","<<pj.z <<"]){b;}\n";
      }

    }
    out <<"}\n";
    myCounter++;

  }

  string OutputBSDL::view(const Vector3DBlock* pos) const{
    stringstream str;
    Vector3D minb,maxb;
    if(myTopology->getVolume() < Constant::REAL_INFINITY)
      myTopology->getBoundaryConditionsBox(minb,maxb);
    else
      myTopology->getBoundingbox(*pos,minb,maxb);

    minb *= myScale;
    maxb *= myScale;

    Vector3D d(maxb-minb);
    Vector3D lookat((maxb+minb)*0.5);
    Vector3D eye(maxb+(d+Vector3D(1,1,1))*2);


    str <<"camera {\n"
	<<"  perspective {\n"
	<<"    eye ["<<eye.x<<","<<eye.y<<","<<eye.z<<"];\n"
	<<"    lookat ["<<lookat.x<<","<<lookat.y<<","<<lookat.z<<"];\n"
	<<"    up [0,0,1];\n"
	<<"    fov 40;\n"
	<<"    resolution (640, 512);\n"
	<<"    eyesep 1;\n"<<"  }\n"
	<<"  background [1,1,1];\n"
	<<"}\n"
	<<"box(["<<(lookat.x-d.x*100)<<","<<(lookat.y-d.y*100)<<","<<(minb.z - d.z/10-10)<<"],"
	<<"["<<(lookat.x+d.x*100)<<","<<(lookat.y+d.y*100)<<","<<(minb.z - d.z/10)<<"]){d;}\n";
	  

    return str.str();
  }

  Output* OutputBSDL::doMake(string&, const vector<Value>& values) const{
    return (new OutputBSDL(values[0],values[1],values[2],values[3],values[4]));
  }

  void OutputBSDL::getParameters(vector<Parameter> &parameter) const{
    parameter.push_back(Parameter(getId(),Value(myFilename,ConstraintValueType::NotEmpty())));
    parameter.push_back(Parameter(keyword+"OutputFreq",Value(myOutputFreq,ConstraintValueType::Positive())));
    parameter.push_back(Parameter(keyword+"MinimalImage",Value(myMinimalImage),true,Text("whether the coordinates should be transformed to minimal image or not")));
    parameter.push_back(Parameter("BSDLIncludeView",Value(myIncludeView),true,Text("whether the camera view should be include in each frame and adapted or not")));
    parameter.push_back(Parameter("BSDLShowWater",Value(myWater),true,Text("whether to suppress water molecules or not")));
  }

  bool OutputBSDL::adjustWithDefaultParameters(vector<Value>& values, const Configuration* config) const{
    if(!checkParameterTypes(values))
      return false;
    if(config->valid(InputOutputfreq::keyword) && !values[1].valid())
      values[1] = (*config)[InputOutputfreq::keyword];
    if(config->valid(InputMinimalImage::keyword) && !values[2].valid())
      values[2] = (*config)[InputMinimalImage::keyword];
    return checkParameters(values);
  }

}
