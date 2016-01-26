#include "Makeable.h"

#include "stringutilities.h"

using std::string;
using std::vector;

namespace ProtoMol {
  //________________________________________________________ Makeable
  
  string Makeable::getId() const{
    if(!myAlias.empty()){
      return getAlias();
    }
    else {
      return getIdNoAlias();
    }
  }

  string Makeable::getAlias() const{
    return myAlias;
  }

  string Makeable::setAlias(const string& id){
    string tmp(myAlias);
    if(equalNocase(id,getIdNoAlias())){
      myAlias = "";
    }
    else {
      myAlias = id;
    }
    return tmp;
  }

  bool Makeable::checkParameters(string& errMsg, const vector<Value>& values) const{
    vector<Parameter> tmp;
    getParameters(tmp);
    bool ok = (tmp.size() == values.size());
    if(!ok){
      errMsg += getId()+": Expected "+toString(tmp.size())+" value(s), but got "+toString(values.size())+".";
    }

    for(unsigned int i=0;i<values.size();++i){
      if(!values[i].valid()){
	if(ok)
	  errMsg += getId()+":";
	errMsg += " Parameter "+toString(i)+" \'"+tmp[i].keyword+"\' "+tmp[i].value.getDefinitionTypeString()+" undefined/missing or with non-valid value \'"+values[i].getString()+"\'.";
	ok = false;
      }
      if(!values[i].equalType(tmp[i].value)){
	if(ok)
	  errMsg += getId()+":";
	errMsg += " Expected type "+tmp[i].value.getDefinitionTypeString()+" for parameter "+toString(i)+" \'"+tmp[i].keyword+"\', but got "+values[i].getDefinitionTypeString()+".";;
	ok = false;
      }
    }
    return ok;
  }

  bool Makeable::checkParameterTypes(const vector<Value>& values) const{
    vector<Parameter> tmp;
    getParameters(tmp);
    if(tmp.size() != values.size())
      return false;

    for(unsigned int i=0;i<values.size();++i)
      if(!values[i].equalType(tmp[i].value) && values[i].defined())
	return false;

    return true;
  }

  bool Makeable::checkParameters(const vector<Value>& values) const{
    vector<Parameter> tmp;
    getParameters(tmp);
    if(tmp.size() != values.size())
      return false;

    for(unsigned int i=0;i<values.size();++i)
      if(!values[i].valid() || !values[i].equalType(tmp[i].value))
	return false;

    return true;
  }

  vector<Parameter> Makeable::getParameters() const{ 
    vector<Parameter> p;
    this->getParameters(p);
    return p;
  }

  MakeableDefinition Makeable::getDefinition() const{
    return MakeableDefinition (getId(),getParameters());
  } 



}
