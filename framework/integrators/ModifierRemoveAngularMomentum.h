/*  -*- c++ -*-  */
#ifndef MODIFIERREMOVEANGULARMOMENTUM_H
#define MODIFIERREMOVEANGULARMOMENTUM_H

#include "Modifier.h"
#include "topologyutilities.h"

namespace ProtoMol {

  //_________________________________________________________________ ModifierRemoveAngularMomentum
  class ModifierRemoveAngularMomentum : public Modifier {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierRemoveAngularMomentum(int freq):Modifier(Constant::MAX_INT-200),myStep(0),myFreq(freq){}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool isInternal() const {return false;}
  private:
    virtual void doExecute(){
      if(myFreq == 0 || 0 == (myStep = (myStep % myFreq))){
	removeAngularMomentum(myPositions,myVelocities,myTopology);
      }
      myStep++;      
    }
    virtual std::string doPrint()const{return std::string("RemoveAngularMomentum");};
  private:
    virtual void doInitialize(){myStep =0;}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    int myStep;
    int myFreq;
  };

}
#endif /* MODIFIER_H */
