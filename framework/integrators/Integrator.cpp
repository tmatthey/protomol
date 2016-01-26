#include "Integrator.h"
#include "Report.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"
#include "ForceGroup.h"
#include "GenericTopology.h"
#include "topologyutilities.h"
#include "ModifierShake.h"
#include "ModifierRattle.h"
#include "ModifierShadow.h"

using namespace ProtoMol::Report;
using std::vector;
using std::string;

namespace ProtoMol {
  //_________________________________________________________________ Integrator

  //  Initialize static members.
  const string  Integrator::scope("Integrator");
  Real          Integrator::myBeta = 0.;

  Integrator::Integrator():
    Makeable(),
    myPotEnergy(0),
    myTopo(NULL),
    myPositions(NULL),
    myVelocities(NULL),
    myForces(NULL),
    myEnergies(NULL),
    myForcesToEvaluate(NULL),
    myForward(true),
    myOldForces(NULL) {}

  Integrator::Integrator(ForceGroup* forceGroup):
    Makeable(),
    myPotEnergy(0),
    myTopo(NULL),
    myPositions(NULL),
    myVelocities(NULL),
    myForces(new Vector3DBlock),
    myEnergies(NULL),
    myForcesToEvaluate(forceGroup),
    myForward(true),
    myOldForces(new Vector3DBlock) {}


  Integrator::~Integrator(){
    delete myForces;
    delete myOldForces;
    delete myForcesToEvaluate;
    for(std::set<Modifier*>::iterator i=myListModifiers.begin();i!=myListModifiers.end();++i)
      delete (*i);
  }

  void Integrator::initialize(GenericTopology *topo, 
			      Vector3DBlock   *positions, 
			      Vector3DBlock   *velocities, 
			      ScalarStructure *energies){

    //Report::report <<"[Integrator::initialize]"<<Report::endr;

    myTopo       = topo;
    myPositions  = positions;
    myVelocities = velocities;
    myEnergies   = energies;

    myForces->zero(positions->size());
    myOldForces->zero(positions->size());

    buildMolecularCenterOfMass(myPositions,myTopo);
    buildMolecularMomentum(myVelocities,myTopo);

    // Initialize only external modifiers,
    // where internal modifiers will be added
    // and initialized at appropriated time
    deleteInternalModifiers();
    initializeModifiers();
  }

  Integrator* Integrator::top(){
    Integrator* i = this;
    for(;i->previous() != NULL;i = i->previous());
    return i;
  }
  
  const Integrator* Integrator::top() const{
    const Integrator* i = this;
    for(;i->previous() != NULL;i = i->previous());
    return i;
  }
  
  Integrator* Integrator::bottom(){
    Integrator* i = this;
    for(;i->next() != NULL;i = i->next());
    return i;
  }
  
  const Integrator* Integrator::bottom() const{
    const Integrator* i = this;
    for(;i->next() != NULL;i = i->next());
    return i;
  }
  
  int Integrator::level() const{
    int n = 0;
    for(const Integrator* i=this;i->next() != NULL;i = i->next())
      n++;
    return n;    
  }

  IntegratorDefinition Integrator::getIntegratorDefinition() const{
    IntegratorDefinition tmp;
    
    // Integrator definition
    tmp.integrator.id = this->getId();
    this->getParameters(tmp.integrator.parameters);
    
    // Force definitions
    if(myForcesToEvaluate != NULL){
      myForcesToEvaluate->getDefinition(tmp.forces);
    }

    return tmp;
  }

  vector<IntegratorDefinition> Integrator::getIntegratorDefinitionAll() const{
    vector<IntegratorDefinition> res;
    for(const Integrator* i=bottom();i != NULL;i = i->previous())
      res.push_back(i->getIntegratorDefinition());
    return res;
  }

  void Integrator::uncache(){
    for(Integrator* i=top();i != NULL;i = i->next()){
      if(i->myForcesToEvaluate != NULL)
	i->myForcesToEvaluate->uncache();
      i->doUncache();
    }
  }


  void Integrator::forward(){
    for(Integrator* i=top();i != NULL;i = i->next()){
      i->myForward = true;
    }
  }
  void Integrator::backward(){
    for(Integrator* i=top();i != NULL;i = i->next()){
      i->myForward = false;
    }
  }


  Modifier* Integrator::createRattleModifier(Real eps, int maxIter){
    return (new ModifierRattle(eps, maxIter, this));
  }

  Modifier*  Integrator::createShakeModifier(Real eps, int maxIter){
    return (new ModifierShake(eps, maxIter, this));
  }

  //  ----------------------------------------------------------------------  //

  Modifier * Integrator::createShadowModifier( int order2k, int freq ) {
      return( new ModifierShadow( order2k, freq, this ) );
  }

  void Integrator::preStepModify(){
    Report::report << Report::debug(10)<<"[Integrator::preStepModify] ("<<(long)this<<") "<<(myTopo != NULL ? myTopo->time : 0.0)<<Report::endr;
    for(std::set<Modifier*>::iterator i=myPreStepModifiers.begin();i!=myPreStepModifiers.end();++i){
      (*i)->execute();
    }
  }

  void Integrator::preDriftOrNextModify(){
    Report::report << Report::debug(10)<<"[Integrator::preDriftOrNextModify] ("<<(long)this<<") "<<(myTopo != NULL ? myTopo->time : 0.0)<<Report::endr;
    for(std::set<Modifier*>::iterator i=myPreDriftOrNextModifiers.begin();i!=myPreDriftOrNextModifiers.end();++i){
      (*i)->execute();
    }
  }

  void Integrator::postDriftOrNextModify(){
    Report::report << Report::debug(10)<<"[Integrator::postDriftOrNextModify] ("<<(long)this<<") "<<(myTopo != NULL ? myTopo->time : 0.0)<<Report::endr;
    for(std::set<Modifier*>::iterator i=myPostDriftOrNextModifiers.begin();i!=myPostDriftOrNextModifiers.end();++i){
      (*i)->execute();
    }
  }

  void Integrator::preForceModify(){
    Report::report << Report::debug(10)<<"[Integrator::preForceModify] ("<<(long)this<<") "<<(myTopo != NULL ? myTopo->time : 0.0)<<Report::endr;
    for(std::set<Modifier*>::iterator i=myPreForceModifiers.begin();i!=myPreForceModifiers.end();++i){
      (*i)->execute();
    }
  }

  void Integrator::mediForceModify(){
    Report::report << Report::debug(10)<<"[Integrator::mediForceModify] ("<<(long)this<<") "<<(myTopo != NULL ? myTopo->time : 0.0)<<Report::endr;
    for(std::set<Modifier*>::iterator i=myMediForceModifiers.begin();i!=myMediForceModifiers.end();++i){
      (*i)->execute();
    }
  }

  void Integrator::postForceModify(){
    Report::report << Report::debug(10)<<"[Integrator::postForceModify] ("<<(long)this<<") "<<(myTopo != NULL ? myTopo->time : 0.0)<<Report::endr;
    for(std::set<Modifier*>::iterator i=myPostForceModifiers.begin();i!=myPostForceModifiers.end();++i){
      (*i)->execute();
    }
  }

  void Integrator::postStepModify(){
    Report::report << Report::debug(10)<<"[Integrator::postStepModify] ("<<(long)this<<") "<<(myTopo != NULL ? myTopo->time : 0.0)<<Report::endr;
    for(std::set<Modifier*>::iterator i=myPostStepModifiers.begin();i!=myPostStepModifiers.end();++i){
      (*i)->execute();
    }
  }

  //  ---------------------------------------------------------------------  //

  void Integrator::adoptPreStepModifier(Modifier* modifier){
    Report::report << Report::debug(10)<<"[Integrator::adoptPreStepModifier] "<<modifier->print()<<"("<<(long)modifier<<") "<<(myTopo != NULL ? myTopo->time : 0.0)<<Report::endr;
    if(myTopo != NULL)
      modifier->initialize(myTopo,myPositions,myVelocities,myForces,myEnergies);
    myPreStepModifiers.insert(modifier);
    addModifier(modifier);
  }

  void Integrator::adoptPreDriftOrNextModifier(Modifier* modifier){
    Report::report << Report::debug(10)<<"[Integrator::adoptPreDriftOrNextModifier] "<<modifier->print()<<"("<<(long)modifier<<") "<<(myTopo != NULL ? myTopo->time : 0.0)<<Report::endr;
    if(myTopo != NULL)
      modifier->initialize(myTopo,myPositions,myVelocities,myForces,myEnergies);
    myPreDriftOrNextModifiers.insert(modifier);
    addModifier(modifier);
  }

  void Integrator::adoptPostDriftOrNextModifier(Modifier* modifier){
    Report::report << Report::debug(10)<<"[Integrator::adoptPostDriftOrNextModifier] "<<modifier->print()<<"("<<(long)modifier<<") "<<(myTopo != NULL ? myTopo->time : 0.0)<<Report::endr;
    if(myTopo != NULL)
      modifier->initialize(myTopo,myPositions,myVelocities,myForces,myEnergies);
    myPostDriftOrNextModifiers.insert(modifier);
    addModifier(modifier);
  }

  void Integrator::adoptPreForceModifier(Modifier* modifier){
    Report::report << Report::debug(10)<<"[Integrator::adoptPreForceModifier] "<<modifier->print()<<"("<<(long)modifier<<") "<<(myTopo != NULL ? myTopo->time : 0.0)<<Report::endr;
    if(myTopo != NULL)
      modifier->initialize(myTopo,myPositions,myVelocities,myForces,myEnergies);
    myPreForceModifiers.insert(modifier);
    addModifier(modifier);
  }

  void Integrator::adoptMediForceModifier(Modifier* modifier){
    Report::report << Report::debug(10)<<"[Integrator::adoptMediForceModifier] "<<modifier->print()<<"("<<(long)modifier<<") "<<(myTopo != NULL ? myTopo->time : 0.0)<<Report::endr;
    if(myTopo != NULL)
      modifier->initialize(myTopo,myPositions,myVelocities,myForces,myEnergies);
    myMediForceModifiers.insert(modifier);
    addModifier(modifier);
  }

  void Integrator::adoptPostForceModifier(Modifier* modifier){
    Report::report << Report::debug(10)<<"[Integrator::adoptPostForceModifier] "<<modifier->print()<<"("<<(long)modifier<<") "<<(myTopo != NULL ? myTopo->time : 0.0)<<Report::endr;
    if(myTopo != NULL)
      modifier->initialize(myTopo,myPositions,myVelocities,myForces,myEnergies);
    myPostForceModifiers.insert(modifier);
    addModifier(modifier);
  }

  void Integrator::adoptPostStepModifier(Modifier* modifier){
    Report::report << Report::debug(10)<<"[Integrator::adoptPostStepModifier] "<<modifier->print()<<"("<<(long)modifier<<") "<<(myTopo != NULL ? myTopo->time : 0.0)<<Report::endr;
    if(myTopo != NULL)
      modifier->initialize(myTopo,myPositions,myVelocities,myForces,myEnergies);
    myPostStepModifiers.insert(modifier);
    addModifier(modifier);
  }

  //  ---------------------------------------------------------------------  //

  void Integrator::deleteInternalModifiers(){
    Report::report << Report::debug(10)<<"[Integrator::deleteInternalModifiers]"<<endr;
    for(std::set<Modifier*>::iterator i=myPreStepModifiers.begin();i!=myPreStepModifiers.end();++i){
      if((*i)->isInternal()){
	myPreStepModifiers.erase(i);
	deleteModifier(*i);
      }
    }
    for(std::set<Modifier*>::iterator i=myPreDriftOrNextModifiers.begin();i!=myPreDriftOrNextModifiers.end();++i){
      if((*i)->isInternal()){
	myPreDriftOrNextModifiers.erase(i);
	deleteModifier(*i);
      }
    }
    for(std::set<Modifier*>::iterator i=myPostDriftOrNextModifiers.begin();i!=myPostDriftOrNextModifiers.end();++i){
      if((*i)->isInternal()){
	myPostDriftOrNextModifiers.erase(i);
	deleteModifier(*i);
      }
    }
    for(std::set<Modifier*>::iterator i=myPreForceModifiers.begin();i!=myPreForceModifiers.end();++i){
      if((*i)->isInternal()){
	myPreForceModifiers.erase(i);
	deleteModifier(*i);
      }
    }
    for(std::set<Modifier*>::iterator i=myMediForceModifiers.begin();i!=myMediForceModifiers.end();++i){
      if((*i)->isInternal()){
	myMediForceModifiers.erase(i);
	deleteModifier(*i);
      }
    }
    for(std::set<Modifier*>::iterator i=myPostForceModifiers.begin();i!=myPostForceModifiers.end();++i){
      if((*i)->isInternal()){
	myPostForceModifiers.erase(i);
	deleteModifier(*i);
      }
    }
    for(std::set<Modifier*>::iterator i=myPostStepModifiers.begin();i!=myPostStepModifiers.end();++i){
      if((*i)->isInternal()){
	myPostStepModifiers.erase(i);
	deleteModifier(*i);
      }
    }
  }

  //  ---------------------------------------------------------------------  //

  void Integrator::deleteExternalModifiers(){
    Report::report << Report::debug(10)<<"[Integrator::deleteExternalModifiers] size="<<myListModifiers.size()<<endr;
    for(std::set<Modifier*>::iterator i=myPreStepModifiers.begin();i!=myPreStepModifiers.end();++i){
      if(!((*i)->isInternal())){
	myPreStepModifiers.erase(i);
	deleteModifier(*i);
      }
    }
    for(std::set<Modifier*>::iterator i=myPreDriftOrNextModifiers.begin();i!=myPreDriftOrNextModifiers.end();++i){
      if(!((*i)->isInternal())){
	myPreDriftOrNextModifiers.erase(i);
	deleteModifier(*i);
      }
    }
    for(std::set<Modifier*>::iterator i=myPostDriftOrNextModifiers.begin();i!=myPostDriftOrNextModifiers.end();++i){
      if(!((*i)->isInternal())){
	myPostDriftOrNextModifiers.erase(i);
	deleteModifier(*i);
      }
    }
    for(std::set<Modifier*>::iterator i=myPreForceModifiers.begin();i!=myPreForceModifiers.end();++i){
      if(!((*i)->isInternal())){
	myPreForceModifiers.erase(i);
	deleteModifier(*i);
      }
    }
    for(std::set<Modifier*>::iterator i=myMediForceModifiers.begin();i!=myMediForceModifiers.end();++i){
      if(!((*i)->isInternal())){
	myMediForceModifiers.erase(i);
	deleteModifier(*i);
      }
    }
    for(std::set<Modifier*>::iterator i=myPostForceModifiers.begin();i!=myPostForceModifiers.end();++i){
      if(!((*i)->isInternal())){
	myPostForceModifiers.erase(i);
	deleteModifier(*i);
      }
    }
    for(std::set<Modifier*>::iterator i=myPostStepModifiers.begin();i!=myPostStepModifiers.end();++i){
      if(!((*i)->isInternal())){
	myPostStepModifiers.erase(i);
	deleteModifier(*i);
      }
    }
    Report::report << Report::debug(10)<<"[Integrator::deleteExternalModifiers] end size="<<myListModifiers.size()<<endr;
  }

  //  ---------------------------------------------------------------------  //

  bool Integrator::removeModifier(const Modifier* modifier){
    Report::report << Report::debug(10)<<"[Integrator::removeModifier]"<<endr;
    bool ok = false;
    for(std::set<Modifier*>::iterator i=myPreStepModifiers.begin();i!=myPreStepModifiers.end();++i){
      if(modifier == (*i)){
	myPreStepModifiers.erase(i);
	ok = true;
      }
    }
    for(std::set<Modifier*>::iterator i=myPreDriftOrNextModifiers.begin();i!=myPreDriftOrNextModifiers.end();++i){
      if(modifier == (*i)){
	myPreDriftOrNextModifiers.erase(i);
	ok = true;
      }
    }
    for(std::set<Modifier*>::iterator i=myPostDriftOrNextModifiers.begin();i!=myPostDriftOrNextModifiers.end();++i){
      if(modifier == (*i)){
	myPostDriftOrNextModifiers.erase(i);
	ok = true;
      }
    }
    for(std::set<Modifier*>::iterator i=myPreForceModifiers.begin();i!=myPreForceModifiers.end();++i){
      if(modifier == (*i)){
	myPreForceModifiers.erase(i);
	ok = true;
      }
    }
    for(std::set<Modifier*>::iterator i=myMediForceModifiers.begin();i!=myMediForceModifiers.end();++i){
      if(modifier == (*i)){
	myMediForceModifiers.erase(i);
	ok = true;
      }
    }
    for(std::set<Modifier*>::iterator i=myPostForceModifiers.begin();i!=myPostForceModifiers.end();++i){
      if(modifier == (*i)){
	myPostForceModifiers.erase(i);
	ok = true;
      }
    }
    for(std::set<Modifier*>::iterator i=myPostStepModifiers.begin();i!=myPostStepModifiers.end();++i){
      if(modifier == (*i)){
	myPostStepModifiers.erase(i);
	ok = true;
      }
    }
    if(ok){
      std::set<Modifier*>::iterator i = myListModifiers.find(const_cast<Modifier*>(modifier));
      if(i != myListModifiers.end()){
	myListModifiers.erase(i);
      }
    }

    return ok;
  }

  //  ---------------------------------------------------------------------  //

  void Integrator::initializeModifiers(){
    Report::report << Report::debug(10)<<"[Integrator::initializeModifiers] "<<(myTopo != NULL ? myTopo->time : 0.0)<<Report::endr;
    for(std::set<Modifier*>::iterator i=myListModifiers.begin();i!=myListModifiers.end();++i){
      (*i)->initialize(myTopo,myPositions,myVelocities,myForces,myEnergies);
    }
  }

  //  ---------------------------------------------------------------------  //

  void Integrator::addModifier(Modifier* modifier){
    Report::report << Report::debug(10)<<"[Integrator::addModifier] size=";
    myListModifiers.insert(modifier);
    Report::report << myListModifiers.size()<<endr;
  }

  //  ---------------------------------------------------------------------  //

  void Integrator::deleteModifier(Modifier* modifier){
    Report::report << Report::debug(10)<<"[Integrator::deleteModifier] size="<<myListModifiers.size()<<","<<modifier->isInternal()<<endr;
    std::set<Modifier*>::iterator i = myListModifiers.find(modifier);
    if(i != myListModifiers.end()){
      Report::report << Report::debug(10) << "[Integrator::deleteModifier] delete"<<(long)(modifier)<<endr;
      delete modifier;
      myListModifiers.erase(i);
    }
    Report::report << Report::debug(10) <<"[Integrator::deleteModifier] end size="<<myListModifiers.size()<<endr;
  }
  
  
    //  --------------------------------------------------------------------  //
    //  The last modifier found with modifierName, is removed.                //
    //  --------------------------------------------------------------------  //

    bool Integrator::removeModifier( const std::string modifierName ) {

        Report::report << Report::debug(10)<<"[Integrator::removeModifier]"<<endr;

        bool found = false;
        Modifier *foundMod = NULL;

        for( std::set<Modifier*>::iterator i = myPostStepModifiers.begin();
                i != myPostStepModifiers.end();
                ++i ) {

            if( (*i)->print() == modifierName ) {
                myPostStepModifiers.erase(i);
                foundMod = *i;
                found = true;
            }

        }

        for( std::set<Modifier*>::iterator i = myPreStepModifiers.begin();
                i != myPreStepModifiers.end();
                ++i ) {

            if( (*i)->print() == modifierName ) {
                myPreStepModifiers.erase(i);
                foundMod = *i;
                found = true;
            }

        }

        for( std::set<Modifier*>::iterator i = myPreDriftOrNextModifiers.begin();
                i != myPreDriftOrNextModifiers.end();
                ++i ) {

            if( (*i)->print() == modifierName ) {
                myPreDriftOrNextModifiers.erase(i);
                foundMod = *i;
                found = true;
            }

        }

        for( std::set<Modifier*>::iterator i = myPostDriftOrNextModifiers.begin();
                i != myPostDriftOrNextModifiers.end();
                ++i ) {

            if( (*i)->print() == modifierName ) {
                myPostDriftOrNextModifiers.erase(i);
                foundMod = *i;
                found = true;
            }

        }

        for( std::set<Modifier*>::iterator i = myPreForceModifiers.begin();
                i != myPreForceModifiers.end();
                ++i ) {

            if( (*i)->print() == modifierName ) {
                myPreForceModifiers.erase(i);
                foundMod = *i;
                found = true;
            }

        }

        for( std::set<Modifier*>::iterator i = myMediForceModifiers.begin();
                i != myMediForceModifiers.end();
                ++i ) {

            if( (*i)->print() == modifierName ) {
                myMediForceModifiers.erase(i);
                foundMod = *i;
                found = true;
            }

        }

        for( std::set<Modifier*>::iterator i = myPostForceModifiers.begin();
                i != myPostForceModifiers.end();
                ++i ) {

            if( (*i)->print() == modifierName ) {
                myPostForceModifiers.erase(i);
                foundMod = *i;
                found = true;
            }

        }

        if( found ){
            std::set<Modifier*>::iterator i = myListModifiers.find(
                    const_cast<Modifier*>( foundMod ) );

            if( i != myListModifiers.end() ) {
                myListModifiers.erase(i);
            }

        }

        return found;

    }


    //  --------------------------------------------------------------------  //
    //  These methods will save/restore the forces.  This is probably only    //
    //  going to be used by *MCIntegrator methods, so it may one day be       //
    //  moved to a more appropriate position.                                 //
    //  --------------------------------------------------------------------  //

    void Integrator::saveForces() {
        (*myOldForces) = (*myForces);
    }

    void Integrator::restoreForces() {
        (*myForces) = (*myOldForces);
    }

}

