/*  -*- c++ -*-  */
#ifndef iSGMODIFIERPREFORCECHEMOSTAT_H
#define iSGMODIFIERPREFORCECHEMOSTAT_H

#include "Modifier.h"
#include "iSGIntegrator.h"

using namespace ProtoMol::Report;

namespace ProtoMol {

  //_________________________________________________________________ iSGModifierPreForceChemostat
  class iSGModifierPreForceChemostat : public Modifier {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    iSGModifierPreForceChemostat(iSGIntegrator* i):Modifier(Constant::MAX_INT-400),myTheIntegrator(i){}
    iSGModifierPreForceChemostat(iSGIntegrator* i, int order):Modifier(order),myTheIntegrator(i){}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Modifier
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual bool isInternal() const {return true;}
  private:
    virtual void doExecute(){
      myTheIntegrator->PreForceChemostat();
    }
    virtual std::string doPrint()const{return std::string("PretForceChemostat");}
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    iSGIntegrator* myTheIntegrator;
  };

}
#endif /* ISGMODIFIERPREFORCHEMOSTAT_H */
