#include "TopologyFactoryDetails.h"
#include "Configuration.h"
#include "HelpTextFactory.h"

using std::vector;
using std::map;
using std::set;
using std::string;
using namespace ProtoMol::Report;

namespace ProtoMol {

  //_____________________________________________________ TopologyFactoryDetails
  TopologyFactoryDetails::TopologyFactoryDetails():FactoryBase<GenericTopology>(){
    //Report::report << Report::plain <<"[TopologyFactoryDetails::TopologyFactoryDetails()]"<< Report::endr;
  }

  TopologyFactoryDetails::~TopologyFactoryDetails(){
    //Report::report << Report::plain <<"TopologyFactoryDetails::~TopologyFactoryDetails()]"<< Report::endr;
  }

  TopologyFactoryDetails::TopologyFactoryDetails(const TopologyFactoryDetails&){}

  TopologyFactoryDetails& TopologyFactoryDetails::operator=(const TopologyFactoryDetails&){return *this;}

  void TopologyFactoryDetails::registerAllExemplarsConfiguration(Configuration* config){
    TFactory::instance().doRegisterAllExemplarsConfiguration(config);
  }

  GenericTopology* TopologyFactoryDetails::make(string& errMsg, const Configuration* config){
    string id = config->get(GenericTopology::getKeyword()).getString();
    const GenericTopology* prototype = TFactory::instance().getPrototype(id);

    return TFactory::instance().doMake(errMsg,id,config->get(prototype!=NULL?prototype->getParameters():vector<Parameter>()));
  }

  GenericTopology* TopologyFactoryDetails::make(string& errMsg, const string& id, const vector<Value>& values){
    return TFactory::instance().doMake(errMsg,id,values);
  }

  GenericTopology* TopologyFactoryDetails::doMake(string& errMsg, const string& id, const vector<Value>& values) const{
    errMsg = "";
    const GenericTopology* prototype = getPrototype(id);

    if(prototype == NULL) {
      errMsg += " Could not find any match for \'"+id+"\' in "+GenericTopology::scope+"Factory.\nPossible topologies are:\n"+doPrint();
      return NULL;
    }

    // Make
    GenericTopology* newObj = prototype->make(errMsg,values);
    if(newObj == NULL)
      return NULL;

    // Adjust external alias
    newObj->setAlias(id);
    return newObj;
  }


  void TopologyFactoryDetails::doRegisterAllExemplarsConfiguration(Configuration* config) const {
    for(set<const GenericTopology*>::const_iterator i=myPointers.begin();i != myPointers.end();++i){
      vector<Parameter> parameter = (*i)->getParameters();
      for(unsigned int i=0;i<parameter.size();i++)
	config->registerKeyword(parameter[i].keyword,parameter[i].value);	
    }
    config->registerKeyword(GenericTopology::getKeyword(),Value(string(""),ConstraintValueType::NotEmpty()));
    myCache = false;
  }

  string TopologyFactoryDetails::doPrint()const{
    string res;
    
    for( map<string,const GenericTopology*,ltstrNocase>::const_iterator i=myExemplars.begin();i != myExemplars.end();++i){
      res += (i==myExemplars.begin()?"":"\n")+i->first;
      vector<Parameter> parameter(i->second->getParameters());
      for(unsigned int k=0;k<parameter.size();k++){
	if(!parameter[k].keyword.empty()){
	  res += "\n"+Constant::PRINTINDENT+Constant::PRINTINDENT+getRightFill(parameter[k].keyword,Constant::PRINTMAXWIDTH);
	}
	res +=(parameter[k].defaultValue.valid()?parameter[k].defaultValue.getDefinitionTypeString():parameter[k].value.getDefinitionTypeString());
	if(!parameter[k].text.empty())
	  res += "\t # "+parameter[k].text;
      }
    }
    res += "\nAlias:";
    for(map<string,const GenericTopology*,ltstrNocase>::const_iterator j=myAliasExemplars.begin();j != myAliasExemplars.end();++j)
      res += "\n"+j->first+" : "+j->second->getId()+" ("+j->second->getIdNoAlias()+")";
    return res;
  }

  void TopologyFactoryDetails::doRegisterHelpText() const{
    for( map<string,const GenericTopology*,ltstrNocase>::const_iterator i=myExemplars.begin();i != myExemplars.end();++i){
      HelpText helpText;
      i->second->getParameters(helpText.parameters);
      helpText.id = i->second->getIdNoAlias();
      helpText.text = i->second->getText();
      helpText.scope = i->second->getScope();
      HelpTextFactory::registerExemplar(i->second->getId(),helpText);

      HelpText alias;
      alias.text   = "alias for \'"+i->second->getId()+"\'";
      alias.scope = i->second->getScope();
      for( map<string,const GenericTopology*,ltstrNocase>::const_iterator j=myAliasExemplars.begin();j != myAliasExemplars.end();++j){
	if(j->second->getIdNoAlias() == i->second->getIdNoAlias()){
	  alias.id = j->first;
	  HelpTextFactory::registerExemplar(alias.id,alias);
	}
      }
    }

  }
}
