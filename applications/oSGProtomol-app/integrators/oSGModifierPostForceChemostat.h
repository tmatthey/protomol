/*  -*- c++ -*-  */
#ifndef OSGMODIFIERPOSTFORCECHEMOSTAT_H
#define OSGMODIFIERPOSTFORCECHEMOSTAT_H

#include "Modifier.h"
#include "oSGIntegrator.h"

using namespace ProtoMol::Report;

namespace ProtoMol {

  //_________________________________________________________________ ModifierPostForceChemostat
  class oSGModifierPostForceChemostat : public Modifier {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    oSGModifierPostForceChemostat(oSGIntegrator* i):Modifier(Constant::MAX_INT-400),myTheIntegrator(i){}
    oSGModifierPostForceChemostat(oSGIntegrator* i, int order):Modifier(order),myTheIntegrator(i){}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool isInternal() const {return true;}
  private:
    virtual void doExecute(){
      myTheIntegrator->PostForceChemostat();
    }
    virtual std::string doPrint()const{return std::string("PostForceChemostat");};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    oSGIntegrator* myTheIntegrator;
  };

}
#endif /* OSGMODIFIERPOSTFORCECHEMOSTAT_H */
