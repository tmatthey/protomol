/*  -*- c++ -*-  */
#ifndef ISGPAR_H
#define ISGPAR_H

#include "PAR.h"

namespace ProtoMol {
  //_________________________________________________________________iSGPAR
  /**
   * Container class for Charmm19/28/XPLOR parameters
   * This struct holds all possible force constants and parameters for
   * each bond, angle, dihedral, etc...@n
   *
   * NB:
   * - angles are always kept in degrees
   * - simga's for Nonbonded are asummed to be Charmm28
   */
  class iSGPAR{
  public:

    //______________________________________________________________________Bond
    /**
     * This structure holds data for a bond, including the bond number, the two
     * atoms involved, the force constants and the distances.
     */      
    struct Bond {

      Bond(){}

      /// Constructor for iSGMD simulations that creates memory for all identities
      Bond(int a) {
	forceConstant.resize(a);
	distance.resize(a);}
      
      /// iSGbond number
      int number;
      /// atom 1 number
      std::string atom1;
      /// atom 2 number
      std::string atom2;
      /// force constants
      std::vector<Real> forceConstant;
      /// distances
      std::vector<Real> distance;
      
      friend Report::MyStreamer& operator<< (Report::MyStreamer& OS, const Bond & p);
    };

    //_____________________________________________________________________Angle
    /** This structure holds data for an angle.  There is the angle number, the three
     * atoms involved, the force constants, the actual values of the angles, and
     * the Urey-Bradley constants if they exist.  The ub_flag will be 1 if they
     * exist and 0 otherwise.
     */
    struct Angle{
      Angle(){}
    
      /// Constructor for iSGMD simulations that creates memory for all identities
      Angle(int a) {
	forceConstant.resize(a);
	angleval.resize(a);
	k_ub.resize(a);
	r_ub.resize(a);}
      
      /// angle number
      int number;
      /// atom 1 number
      std::string atom1;
      /// atom 2 number
      std::string atom2;
      /// atom 3 number
      std::string atom3;
      /// force constants
      std::vector<Real> forceConstant;
      /// angle value
      std::vector<Real> angleval;
      /**
       * Urey-Bradley flag - '1' if there are Urey-Bradley constants following
       * If '0', ignore the next two data members
       */
      bool ub_flag;
      /// Urey-Bradley force constant
      std::vector<Real> k_ub;
      /// Urey-Bradley radius
      std::vector<Real> r_ub;

      friend Report::MyStreamer& operator<< (Report::MyStreamer& OS, const Angle & p);
    };


    //__________________________________________________________________Dihedral
    /** 
     * This structure holds data for a dihedral, consisting of the number, the four
     * atoms involved, the multiplicity (default = 1), and the force constants,
     * the periodicities, and the phase shifts
     */
    struct Dihedral{

      Dihedral(){}
      
      /// Constructor for iSGMD simulations that creates memory for all identities
      Dihedral(int a) {
	multiplicity.resize(a);
	forceConstant.resize(a);
	periodicity.resize(a);
	phaseShift.resize(a);}
      
      /// dihedral number
      int number;
      /// atom 1 number
      std::string atom1;
      /// atom 2 number
      std::string atom2;
      /// atom 3 number
      std::string atom3;
      /// atom 4 number
      std::string atom4;
      /// multiplicity
      std::vector<int> multiplicity;
      /// force constant
      std::vector< std::vector<Real> > forceConstant;
      /// periodicity
      std::vector< std::vector<int> > periodicity;
      /// phase shift
      std::vector< std::vector<Real> > phaseShift;

      friend Report::MyStreamer& operator<< (Report::MyStreamer& OS, const Dihedral & p);
    };


    //__________________________________________________________________Improper
    /**
     * This structure holds data for an improper.  The data held is the same as that
     * for a dihedral - the number, four atoms involved, the force constant,
     * the periodicity, and the phase shift
     */
    struct Improper{
     
      Improper(){}
      
      /// Constructor for iSGMD simulations that creates memory for all identities
      Improper(int a){
	forceConstant.resize(a);
	periodicity.resize(a);
	phaseShift.resize(a);}
      
      /// improper number
      int number;
      /// atom 1 number
      std::string atom1;
      /// atom 2 number
      std::string atom2;
      /// atom 3 number
      std::string atom3;
      /// atom 4 number
      std::string atom4;
      /// force constant
      std::vector<Real> forceConstant;
      /// periodicity
      std::vector<int> periodicity;
      /// phase shift
      std::vector<Real> phaseShift;

      friend Report::MyStreamer& operator<< (Report::MyStreamer& OS, const Improper & p);
    };

    //_________________________________________________________________ AtomType
    /// This class contains information common to one type of atom.  
    struct AtomType {
  
      AtomType(){}
  
      /// constructor for iSGMD simulations that creates memory for all identities
      AtomType(int a, std::string b, std::string c): name(b), symbolName(c) {
	mass.resize(a);
        charge.resize(a);
	stageNumber.resize(a);
      }

      /// The name of this atom type.
      std::string name; 

      /// The masses of this atom type.
      std::vector<Real> mass;


      /// This atomtype's charge for each of its identities.
      std::vector<Real> charge;

      /// The particular transformation stage in which this atom's identity is to be transformed
      std::vector<int> stageNumber;

      /// The symbol untity name of this atom type.
      std::string symbolName;
    };
    
    //_________________________________________________________________Nonbonded
    /// This structure holds data for an atom's nonbonded parameters
    struct Nonbonded{

      Nonbonded(){}
     
      /// constructor for iSGMD simulations that creates memory for all identities 
      Nonbonded(int a) {
        polarizability.resize(a);
        epsilon.resize(a);
        sigma.resize(a);
        polarizability2.resize(a);
        epsilon14.resize(a);
        sigma14.resize(a);
	negative.resize(a);
	vdw.resize(a);
	negative2.resize(a);}
      
      /// nonbonded number
      int number;
      /// atom number
      std::string atom;
      /// polarizability or ignore (see description of negative below),  default to zero
      std::vector<Real> polarizability; 
      /// well depth or number of effective electrons (see description of negative below)
      std::vector<Real> epsilon;
      /// minimum radius divided by 2
      std::vector<Real> sigma;
      /**
       * flag for if the second term is negative - if so, second_term = epsilon or well-depth and 
       * the first term is ignored, otherwise second_term = number of effective electrons 
       * and the first term is the polarizability.@n
       *
       * default to true
       */
      std::vector<bool> negative; 
      /**
       * flag to see if there is a second set, indicating VDW parameters@n
       *
       *  default to true - likely there will be epsilon 1:4 and sigma 1:4
       */
      std::vector<bool> vdw; 
      /**
       * VDW parameter polarizability @n
       *
       * default to zero
       */
      std::vector<Real> polarizability2; 
      /// VDW parameter well depth or number of effective electrons (see above)
      std::vector<Real> epsilon14;
      /// VDW parameter minimum radius divided by 2
      std::vector<Real> sigma14;
      /**
       * flag for a negative VDW paramenter second term (see above) @n
       *
       *  default to true
       */
      std::vector<bool> negative2; 
      

      static const Real SIGMA_CHARMM19_TO_CHARMM28;
      static const Real SIGMA_CHARMM28_TO_CHARMM19;

      friend Report::MyStreamer& operator<< (Report::MyStreamer& OS, const Nonbonded & p);
    };

    //_________________________________________________________________Nbfix
    typedef PAR::Nbfix Nbfix;
    //_________________________________________________________________Hbond
    typedef PAR::Hbond Hbond;

    //_________________________________________________________________iSGPAR

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Types and Enums
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    /// Supported Charmm/XPLOR type
    enum CharmmTypeEnum {
      UNDEFINED,
      CHARMM28,
      CHARMM19
    };


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class iSGPAR
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void clear();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // iSGPAR container
    std::vector<Bond> bonds;
    std::vector<Angle> angles;
    std::vector<Dihedral> dihedrals;
    std::vector<Improper> impropers;
    std::vector<Nonbonded> nonbondeds;
    std::vector<Nbfix> nbfixs;
    std::vector<Hbond> hbonds;

  };

  //____________________________________________________________________________INLINES

}
#endif /* ISGPAR_H */
