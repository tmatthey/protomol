/* -*- c++ -*- */
#ifndef ISGCOULOMBFORCE_H
#define ISGCOULOMBFORCE_H

#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "ExclusionTable.h"
#include "Parameter.h"
#include <string>

namespace ProtoMol {

  //_________________________________________________________________ ISGCoulombForce
  class iSGCoulombForce {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // This uses the weighted charges on each atom, so the Coulomb 
    // constant here is one.
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class ISGCoulombForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void operator()(Real &energy, 
		    Real &force,
		    Real &DeltaMu,
		    Real /*distSquared*/,
		    Real rDistSquared,
		    const Vector3D &, 
		    const GenericTopology* topo, 
		    int atom1, int atom2,
		    ExclusionClass excl) const {

      // needed variables
      bool atom1_scaled, atom2_scaled;
      atom1_scaled = atom2_scaled = false;
      int myStage, OldStage, atom_stage;
      myStage = OldStage = atom_stage = 0;
      
      // Use the molecule types to determine if the atoms
      // are being transformed or not
      int M1 = topo->atoms[atom1].molecule;
      int M2 = topo->atoms[atom2].molecule;
      if (topo->molecules[M1].lambda != 0.0) {atom1_scaled = true;}
      if (topo->molecules[M2].lambda != 0.0) {atom2_scaled = true;}

      // integer used to select which interaction type to compute   
      int Choice;

      // determine the value of Choice
      // neither atom is being transformed
      if ( !(atom1_scaled) && !(atom2_scaled) ) {Choice = 0;}
      // intramolecular interaction on a transforming molecule
      else if ( (atom1_scaled) && (atom2_scaled) ) {             

        // determine the current transformation stage of the molecule
        myStage = static_cast<int>(floor(topo->molecules[M1].lambda)) + 1;

        // get the transformation stage #'s for these atoms
        int atom1_stage = topo->atoms[atom1].stageNumber;
        int atom2_stage = topo->atoms[atom2].stageNumber;

        // if the current stage is not the same as either atom1_stage or atom2_stage
        // the neither atom is currently being transformed
        if (myStage != atom1_stage && myStage != atom2_stage) Choice = 0;
        else Choice = 2;
      }
      // one of the atoms is being transformed
      else {
        // determine the current transformation stage of the molecule
        // and get the transformation stage # for this atoms
        if (atom1_scaled) {
          myStage = static_cast<int>(floor(topo->molecules[M1].lambda)) + 1;
          atom_stage = topo->atoms[atom1].stageNumber;}
        else {
          myStage = static_cast<int>(floor(topo->molecules[M2].lambda)) + 1;
          atom_stage = topo->atoms[atom2].stageNumber;}
        
        // if the current stage is not the same as atom_stage
        // then neither atom is currently being transformed
        if (myStage != atom_stage) Choice = 0;
        else Choice = 1;
      }

      // the scaled charges on each atom
      Real q1 = topo->atoms[atom1].scaledCharge;
      Real q2 = topo->atoms[atom2].scaledCharge;

      //************************************************************************
      // Compute the appropriate electrostatic interaction
      switch (Choice) {
	//----------------------------------------------------------
      case 0:
	// interaction between two non-transforming atoms
	
	//compute the interaction energy
	energy = q1 * q2 * sqrt(rDistSquared);
	
	// scale the energy and chemical potential difference if necessary
	if (topo->coulombScalingFactor != 1.0 && excl == EXCLUSION_MODIFIED)
	  energy *= topo->coulombScalingFactor;
	
	// compute the interaction force magnitude
	force = energy * rDistSquared;
	break;
	//----------------------------------------------------------
      case 1:  
	// interaction between non-transforming atom and a transforming atom
	
	// compute the interaction energy
	energy = q1 * q2 * sqrt(rDistSquared);
	
	// compute the chemical potential difference
	DeltaMu = (q1 * topo->atoms[atom2].deltaQ + q2 * topo->atoms[atom1].deltaQ)
	  * sqrt(rDistSquared);
	
	// compute the interaction force magnitude
	force = energy * rDistSquared;
	break;
	//----------------------------------------------------------
      case 2:  
	// intramolecular interaction on a transforming molecule
	
	// get the value of lambda for this molecule
	Real myLambda = topo->molecules[M1].lambda;

        // determine the current transformation stage of the molecule
        OldStage = myStage - 1;
        
	// compute the interaction energy, which is the sum of the interaction in the
	// old state scaled by (1 - lambda) plus the interaction in the new state
	// scaled by lambda 
	energy = ((myStage - myLambda) * topo->atoms[atom1].Qold * topo->atoms[atom2].Qold
	  + (myLambda - OldStage) * topo->atoms[atom1].Qnew * topo->atoms[atom2].Qnew)
	  * sqrt(rDistSquared);
	
	// compute the chemical potential difference, which is the derivative of the 
	// energy above with respect to lambda
	DeltaMu = (topo->atoms[atom1].Qnew * topo->atoms[atom2].Qnew 
          - topo->atoms[atom1].Qold * topo->atoms[atom2].Qold)
          * sqrt(rDistSquared);
	
	
	// scale the energy and chemical potential difference if necessary
	if (topo->coulombScalingFactor != 1.0 && excl == EXCLUSION_MODIFIED) {
	  energy *= topo->coulombScalingFactor;
	  DeltaMu *= topo->coulombScalingFactor;}
	
	// compute the interaction force magnitude
	force = energy * rDistSquared;
	break;
	//----------------------------------------------------------
      } // end switch structure
      
    } // end operator()

    static void accumulateEnergy(ScalarStructure* energies, Real energy, Real deltaMu) {
      (*energies)[ScalarStructure::COULOMB] += energy;
      (*energies)[ScalarStructure::COULOMB_DELTAMU] += deltaMu;
      //Report::report.precision(8);
      //Report::report << "current DeltaMu (Coulomb) = " << energies->DeltaMu << Report::endr;
      //Report::report << "current energy (Coulomb) = " << energies->coulombEnergy << Report::endr;
    }

    static Real getEnergy(const ScalarStructure* energies) {return (*energies)[ScalarStructure::COULOMB];}

    // Parsing
    static std::string getId() {return keyword;}
    static unsigned int getParameterSize() {return 0;}
    void getParameters(std::vector<Parameter>&) const{}
    static iSGCoulombForce make(std::string&, std::vector<Value>) {
      return iSGCoulombForce();
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Sub Classes
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    class C1 {
    public:
      static Real kernel(Real r)                      {return (1.0/r);}
      static Real dKernel(Real r)                     {r=1.0/r;return (-r*r);}
      static Real kernelR(Real rr)                    {return (rr);}
      static Real dKernelR(Real rr)                   {return (-rr*rr);}
      static Real smooth(Real r,Real /*c*/,Real cr)   {return (cr*(1.5-0.5*r*r*cr*cr));}
      static Real smooth0(Real /*c*/,Real cr)         {return (1.5*cr);}
      static Real dSmooth(Real r,Real /*c*/,Real cr)      {return (-r*cr*cr*cr);}
      static Real smoothKernel(Real r,Real c,Real cr) {return (r<c ?  smooth(r,c,cr) :  kernel(r));}
      static Real dSmoothKernel(Real r,Real c,Real cr){return (r<c ? dSmooth(r,c,cr) : dKernel(r));}
      static std::string getKeyword()                 {return keyword;}
      static std::string getForceKeyword()            {return iSGCoulombForce::keyword;}
    public:
      static const std::string keyword;      
    };
  public:
    class C2 {
    public:
      static Real kernel(Real r)                      {return (1.0/r);}
      static Real dKernel(Real r)                     {r=1.0/r;return (-r*r);}
      static Real kernelR(Real rr)                    {return (rr);}
      static Real dKernelR(Real rr)                   {return (-rr*rr);}
      static Real smooth(Real r,Real /*c*/,Real cr)   {r=r*r*cr*cr;return (cr*(1.875-r*(1.25-0.375*r)));}
      static Real smooth0(Real /*c*/,Real cr)         {return (1.875*cr);}
      static Real dSmooth(Real r,Real c,Real cr)      {c=r*cr*cr;return (cr*c*(1.5*r*c-2.5));}
      static Real smoothKernel(Real r,Real c,Real cr) {return (r<c ?  smooth(r,c,cr) :  kernel(r));}
      static Real dSmoothKernel(Real r,Real c,Real cr){return (r<c ? dSmooth(r,c,cr) : dKernel(r));}
      static std::string getKeyword()                 {return keyword;}
      static std::string getForceKeyword()            {return iSGCoulombForce::keyword;}
    public:
      static const std::string keyword;      
    };
  public:
    class C3 {
    public:
      static Real kernel(Real r)                      {return (1.0/r);}
      static Real dKernel(Real r)                     {r=1.0/r;return (-r*r);}
      static Real kernelR(Real rr)                    {return (rr);}
      static Real dKernelR(Real rr)                   {return (-rr*rr);}
      static Real smooth(Real r,Real /*c*/,Real cr)   {r=r*r*cr*cr;return (0.0625*cr*(35.0-r*(35.0-r*(21.0-5.0*r))));}
      static Real smooth0(Real /*c*/,Real cr)         {return (2.1875*cr);}
      static Real dSmooth(Real r,Real c,Real cr)      {c=r*r*cr*cr;return (r*cr*cr*cr*(-4.375 + c * (5.25 - 1.875 * c)));}
      static Real smoothKernel(Real r,Real c,Real cr) {return (r<c ?  smooth(r,c,cr) :  kernel(r));}
      static Real dSmoothKernel(Real r,Real c,Real cr){return (r<c ? dSmooth(r,c,cr) : dKernel(r));}
      static std::string getKeyword()                 {return keyword;}
      static std::string getForceKeyword()            {return iSGCoulombForce::keyword;}
    public:
      static const std::string keyword;      
    };
  public:
    class C4 {
    public:
      static Real kernel(Real r)                      {return (1.0/r);}
      static Real dKernel(Real r)                     {r=1.0/r;return (-r*r);}
      static Real kernelR(Real rr)                    {return (rr);}
      static Real dKernelR(Real rr)                   {return (-rr*rr);}
      static Real smooth(Real r,Real /*c*/,Real cr)   {r=r*r*cr*cr;return (0.0078125*cr*(315.0-r*(420.0-r*(378.0-r*(180.0-r*35.0)))));}
      static Real smooth0(Real /*c*/,Real cr)         {return (2.4609375*cr);}
      static Real dSmooth(Real r,Real c,Real cr)      {c=r*r*cr*cr;return(-r*cr*cr*cr*(6.5625-c*(11.8125-c*(8.4375-c*2.1875))));
      }
      static Real smoothKernel(Real r,Real c,Real cr) {return (r<c ?  smooth(r,c,cr) :  kernel(r));}
      static Real dSmoothKernel(Real r,Real c,Real cr){return (r<c ? dSmooth(r,c,cr) : dKernel(r));}
      static std::string getKeyword()                 {return keyword;}
      static std::string getForceKeyword()            {return iSGCoulombForce::keyword;}
    public:
      static const std::string keyword;      
    };
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;      
  private:

  };

  //______________________________________________________________________ INLINES
}
#endif /* ISGCOULOMBFORCE_H */
