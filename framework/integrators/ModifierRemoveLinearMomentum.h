/*  -*- c++ -*-  */
#ifndef MODIFIERREMOVELINEARMOMENTUM_H
#define MODIFIERREMOVELINEARMOMENTUM_H

#include "Modifier.h"
#include "topologyutilities.h"

namespace ProtoMol {

  //_________________________________________________________________ ModifierRemoveLinearMomentum
  class ModifierRemoveLinearMomentum : public Modifier {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    ModifierRemoveLinearMomentum(int freq):Modifier(Constant::MAX_INT-100),myStep(0),myFreq(freq){}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool isInternal() const {return false;}
  private:
    virtual void doExecute(){
      if(myFreq == 0 || 0 == (myStep = (myStep % myFreq))){
	removeLinearMomentum(myVelocities,myTopology);
      }
      myStep++;      
    }
    virtual std::string doPrint()const{return std::string("RemoveLinearMomentum");};

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
