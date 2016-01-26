/* -*- c++ -*- */
#ifndef MagneticDipoleMirrorSystemForce_H
#define MagneticDipoleMirrorSystemForce_H

#include "SystemForce.h"
#include "MagneticDipoleMirrorSystemForceBase.h"

namespace ProtoMol {
  //_________________________________________________________________ MagneticDipoleMirrorSystemForce
  // This force calculates the mirror - source effect, i.e. the main-effect that positions the spheres
  // at z = 0
  template<class TBoundaryConditions>
  class MagneticDipoleMirrorSystemForce : public SystemForce, private MagneticDipoleMirrorSystemForceBase {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    MagneticDipoleMirrorSystemForce():myChi(0.0), myR(0.0), myD(0.0), myHx(0.0), myHy(0.0), myHz(0.0){}
    MagneticDipoleMirrorSystemForce(Real chi, Real r, Real Hx, Real Hy, Real Hz, Real d):
      myChi(chi), myR(r), myD(d), myHx(Hx), myHy(Hy), myHz(Hz){
      //Short cuts
      Real realChi = myChi/(1.0-2.0/3.0*myChi);
      Real volume = 4.0/3.0*Constant::M_PI*power<3>(myR);
      Real kappa = realChi/(realChi+2.0);
      Vector3D H(myHx, myHy, myHz/(1+realChi));
      Vector3D sigma = -H*volume*myChi;
      Vector3D sigma_m = sigma*kappa; sigma_m.z = -sigma_m.z;
      myF = sigma.dot(sigma_m)*(3.0);
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class MagneticDipoleMirrorSystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void evaluate(const GenericTopology* topo,
			  const Vector3DBlock* positions, 
			  Vector3DBlock* forces,
			  ScalarStructure* energies);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Force
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getKeyword() const{return keyword;}
  private:
    virtual Force* doMake(std::string&, std::vector<Value> values) const;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual void getParameters(std::vector<Parameter>& parameters) const;
    virtual unsigned int getParameterSize() const{return 6;}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Real myChi, myR, myD, myHx, myHy, myHz, myF;

  }; 

  //______________________________________________________________________ INLINES


  template<class TBoundaryConditions> 
  inline void MagneticDipoleMirrorSystemForce<TBoundaryConditions>::getParameters(std::vector<Parameter>& parameters) const {
    parameters.push_back(Parameter("-chi",Value(myChi)));
    parameters.push_back(Parameter("-r",Value(myR)));
    parameters.push_back(Parameter("-H",Value(myHx)));
    parameters.push_back(Parameter("",Value(myHy)));
    parameters.push_back(Parameter("",Value(myHz)));
    parameters.push_back(Parameter("-d",Value(myD)));

  }

  template<class TBoundaryConditions> 
  inline Force* MagneticDipoleMirrorSystemForce<TBoundaryConditions>::doMake(std::string&, std::vector<Value> values) const {
    return new MagneticDipoleMirrorSystemForce(values[0],values[1],values[2],values[3],values[4],values[5]);
  }

  template<class TBoundaryConditions> 
  inline void MagneticDipoleMirrorSystemForce<TBoundaryConditions>::evaluate(const GenericTopology* topo,
									     const Vector3DBlock* positions,
									     Vector3DBlock* forces,
									     ScalarStructure* energies){
    Real e = 0.0;
    for(unsigned int i=0;i<topo->atoms.size();i++){
      Real z_upper = myD-(*positions)[i].z*2;
      Real z_lower = myD+(*positions)[i].z*2;
      (*forces)[i].z += myF*(1.0/power<4>(z_lower)-1.0/power<4>(z_upper));
      e += myF*(1.0/power<3>(z_lower) - 1.0/power<3>(z_upper));
      // Denne energien faar du pushe tilbake ett eller annet sted hvis den skal vaere med aa bestemme totalenergien.
    }
    (*energies)[ScalarStructure::OTHER] += e;
  }
}
#endif /* MagneticDipoleMirrorSystemForce_H */
