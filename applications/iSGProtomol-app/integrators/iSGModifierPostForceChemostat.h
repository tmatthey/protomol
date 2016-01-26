/*  -*- c++ -*-  */
#ifndef ISGMODIFIERPOSTFORCECHEMOSTAT_H
#define ISGMODIFIERPOSTFORCECHEMOSTAT_H

#include "Modifier.h"
#include "iSGIntegrator.h"

using namespace ProtoMol::Report;

namespace ProtoMol {

  //_________________________________________________________________ ModifierPostForceChemostat
  class iSGModifierPostForceChemostat : public Modifier {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    iSGModifierPostForceChemostat(iSGIntegrator* i):Modifier(Constant::MAX_INT-400),myTheIntegrator(i){}
    iSGModifierPostForceChemostat(iSGIntegrator* i, int order):Modifier(order),myTheIntegrator(i){}
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
    iSGIntegrator* myTheIntegrator;
  };

}
#endif /* ISGMODIFIERPOSTFORCECHEMOSTAT_H */
