/*  -*- c++ -*-  */
#ifndef OSGMODIFIERPREFORCECHEMOSTAT_H
#define OSGMODIFIERPREFORCECHEMOSTAT_H

#include "Modifier.h"
#include "oSGIntegrator.h"

using namespace ProtoMol::Report;

namespace ProtoMol {

  //_________________________________________________________________ oSGModifierPreForceChemostat
  class oSGModifierPreForceChemostat : public Modifier {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    oSGModifierPreForceChemostat(oSGIntegrator* i):Modifier(Constant::MAX_INT-400),myTheIntegrator(i){}
    oSGModifierPreForceChemostat(oSGIntegrator* i, int order):Modifier(order),myTheIntegrator(i){}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool isInternal() const {return true;}
  private:
    virtual void doExecute(){
      myTheIntegrator->PreForceChemostat();
    }
    virtual std::string doPrint()const{return std::string("PreForceChemostat");}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    oSGIntegrator* myTheIntegrator;
  };

}
#endif /* OSGMODIFIERPREFORCHEMOSTAT_H */
