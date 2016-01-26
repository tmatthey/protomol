/*  -*- c++ -*-  */
#ifndef MODIFIERPREFORCETHERMOSTAT_H
#define MODIFIERPREFORCETHERMOSTAT_H

#include "Modifier.h"

namespace ProtoMol {

  //_________________________________________________________________ ModifierPreForceThermostat
  template<class TIntegrator>
  class ModifierPreForceThermostat : public Modifier {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierPreForceThermostat(TIntegrator* i):Modifier(Constant::MAX_INT-400),myTheIntegrator(i){}
    ModifierPreForceThermostat(TIntegrator* i, int order):Modifier(order),myTheIntegrator(i){}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool isInternal() const {return true;}
  private:
    virtual void doExecute(){
      //Report::report <<"ModifierPreForceThermostat"<<Report::endr;
      myTheIntegrator->PreForceThermostat();
    }
    virtual std::string doPrint()const{return std::string("PreForceThermostat");};

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    TIntegrator* myTheIntegrator;
  };

}
#endif /* MODIFIERPREFORCETHERMOSTAT_H */
