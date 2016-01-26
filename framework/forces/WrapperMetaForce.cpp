#include "WrapperMetaForce.h"

#include <map>
#include <set>

using std::set;
using std::map;
using std::swap;
using std::vector;
using std::string;
using namespace ProtoMol::Report;
namespace ProtoMol {

  //_________________________________________________________________ WrapperMetaForce



  WrapperMetaForce::WrapperMetaForce(string id, bool minimal, Force* force, string forceAlias):MetaForce(),keyword(id),myMinimal(minimal){
    addForce(force,forceAlias);
  }

  WrapperMetaForce::WrapperMetaForce(string id, bool minimal, Force* force1, string forceAlias1, Force* force2, string forceAlias2):MetaForce(),keyword(id),myMinimal(minimal){
    addForce(force1,forceAlias1);
    addForce(force2,forceAlias2);
  }

  WrapperMetaForce::~WrapperMetaForce(){
    for(unsigned int i=0;i<myForces.size();i++)
      delete myForces[i];
  }

  void WrapperMetaForce::addForce(Force* force, string forceAlias){
    if(myForces.empty()){
      myParameters.clear();
      myParameterSize = 0;
    }

    if(!forceAlias.empty())
      force->setAlias(forceAlias);

    std::vector<Parameter> parameters = force->getParameters();
    unsigned int count = force->getParameterSize();

    // Remove all multiple keywords
    if(myMinimal){
      set<string,ltstrNocase> unique;
      std::vector<Parameter> tmp;
      count = 0;
      for(unsigned int i=0;i<myParameters.size();i++)
	unique.insert(myParameters[i].keyword);
      for(unsigned int i=0;i<parameters.size();i++){
	if(unique.find(parameters[i].keyword) == unique.end()){
	  tmp.push_back(parameters[i]);
	  count++;
	  unique.insert(parameters[i].keyword);
	}
      }
      swap(tmp,parameters);	
    }

    // Add parameters
    myForces.push_back(force);
    myParameterSize += count;
    for(unsigned int i=0;i<parameters.size();i++)
      myParameters.push_back(parameters[i]);
    
  }

  Force* WrapperMetaForce::doMake(string& errMsg, vector<Value> val) const{
    vector<vector<Value> > values(myForces.size());

    if(myMinimal){
      map<string,Value,ltstrNocase> unique;
      for(unsigned int i=0;i<myParameters.size();i++)
	unique[myParameters[i].keyword] = val[i];
      for(unsigned int i=0;i<myForces.size();i++){
	std::vector<Parameter> p = myForces[i]->getParameters();
	for(unsigned int j=0;j<p.size();j++)
	  values[i].push_back(unique[p[j].keyword]);	
      }
    }
    else {
      unsigned int k = 0;
      for(unsigned int i=0;i<myForces.size();i++){
	for(unsigned int j=0;j<myForces[i]->getParameterSize();j++){
	  values[i].push_back(val[k++]);		
	}
      }
    }

    if(myForces.size() == 1)
      return new WrapperMetaForce(keyword,myMinimal,myForces[0]->make(errMsg,values[0]),"");
    else if(myForces.size() == 2)
      return new WrapperMetaForce(keyword,myMinimal,myForces[0]->make(errMsg,values[0]),"",myForces[1]->make(errMsg,values[1]),"");
    else
      return NULL;

  }

  void WrapperMetaForce::getDeepForces(vector<Force*>& forces) const{
    for(unsigned int i=0;i<myForces.size();i++)
       forces.push_back(myForces[i]); 
  }


}
