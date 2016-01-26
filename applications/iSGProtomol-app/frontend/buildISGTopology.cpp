#include "buildISGTopology.h"

#include "ExclusionTable.h"
#include "ExclusionType.h"
#include "GenericTopology.h"
#include "iSGPAR.h"
#include "PSF.h"
#include "Report.h"
#include "pmconstants.h"
#include "mathutilities.h"
#include "stringutilities.h"
#include "topologyutilities.h"

#include <algorithm>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

using std::list;
using std::map;
using std::set;
using std::pair;
using std::string;
using std::vector;

using namespace ProtoMol::Report;

//#define DEBUG_PRINT_MOLECULETABLE

namespace ProtoMol {
 
  static void findNextNeighbor(int a, vector<int>& v, vector<PairInt>& p, vector<bool>& unused, const vector<vector<int> >& graph, set<PairInt>& pairs);
  // bool cmpSize(const vector<int>& m1, const vector<int>& m2);


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  //  buildISGTopology
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  void buildISGTopology(GenericTopology* topo,const PSF& psf, const iSGPAR& par) {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // First, generate the array of atomtypes
    // Each time a new atom comes up, we need to check if it is
    // already in the vector....
    // NOTE:  this may take a while for large systems; however, it will cut
    // down on the size of the atomTypes vector, and therefore, the amount
    // access time in the back end.
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    topo->atoms.clear();
    topo->atomTypes.clear();
    topo->bonds.clear();
    topo->angles.clear();
    topo->dihedrals.clear();
    topo->impropers.clear();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Get the atoms
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    map<string,int> atomLookUpTable;

    // loop over all atoms in the PSF object
    for (vector<PSF::Atom>::const_iterator atom = psf.atoms.begin();
	 atom != psf.atoms.end(); ++atom) {

      // store the atom_type name for this atom
      AtomType tempatomtype;
      tempatomtype.name = atom->atom_type;
      tempatomtype.symbolName = atomTypeToSymbolName(atom->atom_type);

      // Now check if this already exists (same name)    
      if (atomLookUpTable.find(tempatomtype.name) == atomLookUpTable.end()) {

	// this is a new atom type, add to the atomTypes list
	atomLookUpTable[tempatomtype.name] = topo->atomTypes.size();
	topo->atomTypes.push_back(tempatomtype);
      }

 
      Atom tempatom;
      // First, we need to find the index. (an integer corresponding
      // to the type of the atom
      tempatom.type = atomLookUpTable[tempatomtype.name];

      // Now, the mass and scaled charge.  These are straightforward.      
      tempatom.scaledCharge = (atom->charge)*Constant::SQRTCOULOMBCONSTANT;
      tempatom.scaledMass = atom->mass;

      // Now we need the size of the group for heavy atom ordering
      // We need to parse the name for any H's then any numbers following    
      // First, if the atom is an H then this is 0
      if (atom->atom_type == "H"){
	tempatom.hvyAtom = 0;
      }
      else{
	// Otherwise, we need to parse..
	// Initialize to 1
	tempatom.hvyAtom = 1;
	for (unsigned int pos = 0; pos < atom->atom_type.size(); ++pos){
	  if (atom->atom_type[pos] == 'H'){
	    string number = "";
	    while (isdigit(atom->atom_type[++pos]))   {
	      number += atom->atom_type[pos];
	    }
	    if (number == "") // never entered loop, default is 1
	      number = "1";
	    tempatom.hvyAtom += atoi(number.c_str());
	  }
	}
      }    
      // C/C++ starts at 0, where PSF/PDB at 1
      tempatom.atomNum = atom->number-1;
      // Also the molecule - using residue sequence for now
      topo->atoms.push_back(tempatom);

    } // end loop over atoms

    // calculate the # of degrees of freedom, if there are any bond constraints
    // they will be subtracted later by ModifierShake
    topo->degreesOfFreedom = 3 * topo->atoms.size() - 3;
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Get the bonds
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    // First create look-up-table (into PAR object bonds)
    map<string,vector<iSGPAR::Bond>::const_iterator> bondLookUpTable;

    // take the bonds from the PAR object and put then into the map/look-up-table
    for (vector<iSGPAR::Bond>::const_iterator bond = par.bonds.begin();
	 bond != par.bonds.end(); ++bond){
      bondLookUpTable[bond->atom1+","+bond->atom2] = bond;
      //report << (*bond)<< ", " << bond->atom1 << ", " << bond->atom2 <<endr;
    }


    // Find the parameters from PAR
    int ignoredBonds = 0;

    // loop over the bond list in the PSF object
    for (vector<PSF::Bond>::const_iterator bond = psf.bonds.begin();
	 bond != psf.bonds.end(); ++bond){
    
      // store the ID numbers of the bonded atoms
      int atom1 = bond->atom1-1;
      int atom2 = bond->atom2-1;

      // store the type names of the bonded atoms
      string bond1(topo->atomTypes[topo->atoms[atom1].type].name);
      string bond2(topo->atomTypes[topo->atoms[atom2].type].name);

      // look into the bond table and see if this bond type was stored in the PAR file
      map<string,vector<iSGPAR::Bond>::const_iterator>::const_iterator currentbond = bondLookUpTable.find(bond1+","+bond2);

      // if this bond type has not been found, try reversing the order of the atom types
      if(currentbond == bondLookUpTable.end()){
	currentbond = bondLookUpTable.find(bond2+","+bond1);}
    
      // if we still have not found this bond type in the PAR object, report an error
      if(currentbond == bondLookUpTable.end()){
	report << error << "Could not find bond \'"<<bond1<<"\'-\'"<<bond2<<"\'."<<std::endl;
	for (map<string,vector<iSGPAR::Bond>::const_iterator>::const_iterator i = bondLookUpTable.begin();
	     i != bondLookUpTable.end();
	     ++i){
	  report << plain << i->first<<std::endl;
	}            
	report << endr;
      }  // end if statement
    
      // if we have found this bond type then copy the bond parameters
      // into the topology
      Bond tempbond;

      // store the time-zero identity of this bond
      int bondIdentity = psf.atoms[atom1].identity;
      tempbond.springConstant = currentbond->second->forceConstant[bondIdentity];
      tempbond.restLength = currentbond->second->distance[bondIdentity];
      tempbond.atom1 = atom1;
      tempbond.atom2 = atom2;
      topo->bonds.push_back(tempbond);
      if(tempbond.springConstant == 0.0)
	++ignoredBonds;
    }

    if(ignoredBonds > 0)
      report << hint << "Systems contains "<<ignoredBonds<<" bonds with zero force constants."<<endr;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Get the angles
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
    // First create look-up-table (into PAR object angles)
    map<string,vector<iSGPAR::Angle>::const_iterator> angleLookUpTable;

    // take the angles from the PAR object and put them into the map
    for (vector<iSGPAR::Angle>::const_iterator angle = par.angles.begin();
	 angle != par.angles.end(); ++angle){
      angleLookUpTable[angle->atom1+","+angle->atom2+","+angle->atom3] = angle;
      //report << (*angle)<<endr;
    }
  
    // Find the parameters from PAR
    int ignoredAngles = 0;

    // loop over the angle list in the PSF object
    for (vector<PSF::Angle>::const_iterator angle = psf.angles.begin();
	 angle != psf.angles.end(); ++angle){
    
      // store the ID numbers of the atoms in this angle
      int atom1 = angle->atom1-1;
      int atom2 = angle->atom2-1;
      int atom3 = angle->atom3-1;

      // store the type names of the atoms in this angle
      string angle1(topo->atomTypes[topo->atoms[atom1].type].name);
      string angle2(topo->atomTypes[topo->atoms[atom2].type].name);
      string angle3(topo->atomTypes[topo->atoms[atom3].type].name);

      // look into the angle table and see if this angle type was stored in the PAR file
      map<string,vector<iSGPAR::Angle>::const_iterator>::const_iterator currentangle = angleLookUpTable.find(angle1+","+angle2+","+angle3);

      // if this angle type has not been found, try reversing the order of the atom types
      if(currentangle == angleLookUpTable.end())
	currentangle = angleLookUpTable.find(angle3+","+angle2+","+angle1);
    
      // if we still have not found this angle type in the PAR object, report an error
      if(currentangle == angleLookUpTable.end())
	report << error << "Could not find angle \'"<<angle1<<"\'-\'"<<angle2<<"\'-\'"<<angle3<<"\'."<<endr;

      // if we have found this angle type then copy the angle parameters
      // into the topology
      Angle tempangle;

      // store the time-zero identity of this angle
      int angleIdentity = psf.atoms[atom1].identity;
      tempangle.atom1 = atom1;
      tempangle.atom2 = atom2;	    
      tempangle.atom3 = atom3;	    
      tempangle.forceConstant = currentangle->second->forceConstant[angleIdentity];
      tempangle.restAngle = dtor(currentangle->second->angleval[angleIdentity]);

      // check to see if a Urey-Bradley term has been specified
      if (currentangle->second->ub_flag){ // do we want defaults for these	  
	tempangle.ureyBradleyConstant = currentangle->second->k_ub[angleIdentity];
	tempangle.ureyBradleyRestLength = currentangle->second->r_ub[angleIdentity];
      }
      // no Urey-Bradley term specified
      else{	  
	tempangle.ureyBradleyConstant = 0.0;
	tempangle.ureyBradleyRestLength = 0.0;
      }
      topo->angles.push_back(tempangle);
      if(tempangle.forceConstant == 0.0)
	++ignoredAngles;
    }

    if(ignoredAngles > 0)
      report << hint << "Systems contains "<<ignoredAngles<<" angles with zero force constants."<<endr;

  
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Get the dihedrals
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //                                                                   
    // One change I made was to assume that a dihedral will only appear 
    // once in the .psf file regardless of it's multiplicity.  The      
    // multiplicity should be handled in the .par file.                  
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    // First create look-up-table (into PAR object)
    map<string,vector<iSGPAR::Dihedral>::const_iterator> dihedralLookUpTable;

    // take the dihedrals from the PAR object and put them into the map
    for (vector<iSGPAR::Dihedral>::const_iterator dihedral = par.dihedrals.begin();
	 dihedral != par.dihedrals.end(); ++dihedral){
      dihedralLookUpTable[dihedral->atom1+","+dihedral->atom2+","+dihedral->atom3+","+dihedral->atom4] = dihedral;
      //report << (*dihedral)<<endr;
    }

    // Find the parameters from PAR
    // loop over the dihedral list in the PSF object
    for(vector<PSF::Dihedral>::const_iterator dihedral = psf.dihedrals.begin();
	dihedral != psf.dihedrals.end(); ++dihedral) {

      // store the ID numbers of the atoms in this dihedral
      int atom1 = dihedral->atom1 - 1;
      int atom2 = dihedral->atom2 - 1;
      int atom3 = dihedral->atom3 - 1;
      int atom4 = dihedral->atom4 - 1;

      // store the type names of the atoms in this dihedral
      string dihedral1 = topo->atomTypes[topo->atoms[atom1].type].name;
      string dihedral2 = topo->atomTypes[topo->atoms[atom2].type].name;
      string dihedral3 = topo->atomTypes[topo->atoms[atom3].type].name;
      string dihedral4 = topo->atomTypes[topo->atoms[atom4].type].name;

      // look into the dihedral table and see if this dihedral type was stored in the PAR file
      map<string,vector<iSGPAR::Dihedral>::const_iterator>::const_iterator currentdihedral = 
	dihedralLookUpTable.find(dihedral1+","+dihedral2+","+dihedral3+","+dihedral4);

      // if this dihedral type has not been found, try reversing the order of the atom types
      if(currentdihedral == dihedralLookUpTable.end())
	currentdihedral = dihedralLookUpTable.find(dihedral4+","+dihedral3+","+dihedral2+","+dihedral1);

      // Try wildcards if necessary
      if(currentdihedral == dihedralLookUpTable.end()){
	currentdihedral = dihedralLookUpTable.find("X,"+dihedral2+","+dihedral3+",X");
	if(currentdihedral == dihedralLookUpTable.end())
	  currentdihedral = dihedralLookUpTable.find("X,"+dihedral3+","+dihedral2+",X");
      }

      // if we still have not found this dihedral type in the PAR object, report an error
      if(currentdihedral == dihedralLookUpTable.end())
	report << error << "Could not find dihedral \'"<<dihedral1<<"\'-\'"<<dihedral2<<"\'-\'"<<dihedral3<<"\'-\'"<<dihedral4<<"\'."<<endr;

      // if we have found this dihedral type then copy the 
      // dihedral parameters into the topology
      Torsion torsion;

      // store the time-zero identity of this torsion
      int torsionIdentity = psf.atoms[atom1].identity;
      torsion.atom1         = atom1;
      torsion.atom2         = atom2;
      torsion.atom3         = atom3;
      torsion.atom4         = atom4;
      torsion.periodicity   = currentdihedral->second->periodicity[torsionIdentity];                
      torsion.forceConstant = currentdihedral->second->forceConstant[torsionIdentity];      
      torsion.phaseShift    = dtor(currentdihedral->second->phaseShift[torsionIdentity]);      
      torsion.multiplicity  = currentdihedral->second->multiplicity[torsionIdentity]; 
      torsion.DeltaK.resize(torsion.multiplicity);
      torsion.DeltaPhase.resize(torsion.multiplicity); 

      //report << plain
      //     << dihedral1 <<","
      //     << dihedral2 <<","
      //     << dihedral3 <<","
      //     << dihedral4 <<","
      //     << torsion.atom1 <<","
      //     << torsion.atom2 <<","
      //     << torsion.atom3 <<","
      //     << torsion.atom4 <<","
      //     << torsion.periodicity.size() <<",";
      //for (unsigned int tors=0; tors<torsion.periodicity.size(); tors++) {
      //  report << torsion.periodicity[tors] <<","
      //	 << torsion.forceConstant[tors] <<","
      //	 << torsion.phaseShift[tors] <<","; }
      //report << torsion.multiplicity << endr;
    
      if(topo->dihedrals.size() == 0 || 
	 topo->dihedrals[topo->dihedrals.size()-1].atom1 != atom1 ||
	 topo->dihedrals[topo->dihedrals.size()-1].atom2 != atom2 ||
	 topo->dihedrals[topo->dihedrals.size()-1].atom3 != atom3 ||
	 topo->dihedrals[topo->dihedrals.size()-1].atom4 != atom4){
	topo->dihedrals.push_back(torsion);
	//report << plain //<< (wildcard ? "#":"")
	//          << dihedral1 <<","
	//          << dihedral2 <<","
	//          << dihedral3 <<","
	//          << dihedral4 <<","
	//          << torsion.atom1 <<","
	//          << torsion.atom2 <<","
	//          << torsion.atom3 <<","
	//          << torsion.atom4 <<","
	//          << torsion.periodicity[0] <<","
	//          << torsion.periodicity.size() <<","
	//          << torsion.forceConstant[0] <<","
	//          << torsion.phaseShift[0] <<","
	//          << torsion.multiplicity << endr;
      }
    }


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Get the impropers
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //                                                                   
    // One change I made was to assume that a improper will only appear 
    // once in the .psf file regardless of it's multiplicity.  The      
    // multiplicity should be handled in the .par file.                  
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //         No wildcard usage is allowed for bonds and angles. For dihedrals,
    // two types are allowed; A - B - C - D (all four atoms specified) and
    // X - A - B - X (only middle two atoms specified). Double dihedral
    // specifications may be specified for the four atom type by listing a
    // given set twice. When specifying this type in the topology file, specify
    // a dihedral twice (with nothing intervening) and both forms will be used.
    //
    //         There are five choices for wildcard usage for improper dihedrals;
    // 1) A - B - C - D  (all four atoms, double specification allowed)
    // 2) A - X - X - B
    // 3) X - A - B - C
    // 4) X - A - B - X
    // 5) X - X - A - B
    // When classifying an improper dihedral, the first acceptable match (from
    // the above order) is chosen. The match may be made in either direction
    // ( A - B - C - D = D - C - B - A).
    //
    //         The periodicity value for dihedrals and improper dihedral terms
    // must be an integer. If it is positive, then a cosine functional form is used.
    // Only positive values of 1,2,3,4,5 and 6 are allowed for the vector, parallel
    // vector and cray routines. Slow and scalar routines can use any positive
    // integer and thus dihedral constrains can be of any periodicity.
    //  Reference angle 0.0 and 180.0 degree correspond to minimum in staggered 
    // and eclipsed respectively. Any reference angle is allowed. The value
    // 180 should be prefered over -180 since it is parsed faster and more
    // accuratly. When the periodicity is given as zero, for OTHER THAN THE
    // FIRST dihdral in a multiple dihedral set, then a the amplitude is a
    // constant added to the energy. This is needed to effect the
    // Ryckaert-Bellemans potential for hydrocarbons (see below). 
  

    // First create look-up-table into PAR object
    map<string,vector<iSGPAR::Improper>::const_iterator> improperLookUpTable;

    // take the impropers from the PAR object and put them into the map
    for (vector<iSGPAR::Improper>::const_iterator improper = par.impropers.begin();
	 improper != par.impropers.end(); ++improper){
      improperLookUpTable[improper->atom1+","+improper->atom2+","+improper->atom3+","+improper->atom4] = improper;
      //report << (*improper)<<endr;
    }

    // Find the parameters from PAR
    // loop over the improper list in the PSF object
    for(vector<PSF::Improper>::const_iterator improper = psf.impropers.begin();
	improper != psf.impropers.end(); ++improper) {

      // store the ID numbers of the atoms in this improper
      int atom1 = improper->atom1 - 1;
      int atom2 = improper->atom2 - 1;
      int atom3 = improper->atom3 - 1;
      int atom4 = improper->atom4 - 1;

      // store the type names of the atoms in this improper
      string improper1 = topo->atomTypes[topo->atoms[atom1].type].name;
      string improper2 = topo->atomTypes[topo->atoms[atom2].type].name;
      string improper3 = topo->atomTypes[topo->atoms[atom3].type].name;
      string improper4 = topo->atomTypes[topo->atoms[atom4].type].name;

      // look into the improper table and see if this improper type was stored in the PAR file
      map<string,vector<iSGPAR::Improper>::const_iterator>::const_iterator currentimproper = 
	improperLookUpTable.find(improper1+","+improper2+","+improper3+","+improper4);

      // if this dihedral type has not been found, try reversing the order of the atom types
      if(currentimproper == improperLookUpTable.end())
	currentimproper = improperLookUpTable.find(improper4+","+improper3+","+improper2+","+improper1);

      // Try wildcards if necessary
      // 2) A - X - X - B
      if(currentimproper == improperLookUpTable.end()){
	currentimproper = improperLookUpTable.find(improper1+",X,X,"+improper4);      
	if(currentimproper == improperLookUpTable.end())
	  currentimproper = improperLookUpTable.find(improper4+",X,X,"+improper1);      

      }
      // 3) X - A - B - C
      if(currentimproper == improperLookUpTable.end()){
	currentimproper = improperLookUpTable.find("X,"+improper2+","+improper3+","+improper4);      
	if(currentimproper == improperLookUpTable.end())
	  currentimproper = improperLookUpTable.find(improper4+","+improper3+","+improper2+",X");      
      }

      // 4) X - A - B - X
      if(currentimproper == improperLookUpTable.end()){
	currentimproper = improperLookUpTable.find("X,"+improper2+","+improper3+",X");      
	if(currentimproper == improperLookUpTable.end())
	  currentimproper = improperLookUpTable.find("X,"+improper3+","+improper2+",X");      
      }

      // 5) X - X - A - B
      if(currentimproper == improperLookUpTable.end()){
	currentimproper = improperLookUpTable.find("X,X,"+improper3+","+improper4);      
	if(currentimproper == improperLookUpTable.end())
	  currentimproper = improperLookUpTable.find(improper4+","+improper3+",X,X");      
      }

      // if we still have not found this improper type in the PAR object, report an error
      if(currentimproper == improperLookUpTable.end())
	report << error << "Could not find improper."<<endr;

      // if we have found this improper type then copy the 
      // improper parameters into the topology
      Torsion torsion;

      // store the time-zero identity of this improper
      int improperIdentity = psf.atoms[atom1].identity;
      torsion.atom1         = atom1;
      torsion.atom2         = atom2;
      torsion.atom3         = atom3;
      torsion.atom4         = atom4;
      torsion.periodicity.push_back(currentimproper->second->periodicity[improperIdentity]);                
      torsion.forceConstant.push_back(currentimproper->second->forceConstant[improperIdentity]);      
      torsion.phaseShift.push_back(dtor(currentimproper->second->phaseShift[improperIdentity]));      
      torsion.multiplicity  = 1;
      topo->impropers.push_back(torsion);
      //     report << plain<< (wildcard ? "#":"")
      //            << improper1 <<","
      //            << improper2 <<","
      //            << improper3 <<","
      //            << improper4 <<","
      //            << torsion.atom1 <<","
      //            << torsion.atom2 <<","
      //            << torsion.atom3 <<","
      //            << torsion.atom4 <<","
      //            << torsion.periodicity[0] <<","
      //            << torsion.periodicity.size() <<","
      //            << torsion.forceConstant[0] <<","
      //            << torsion.phaseShift[0] <<","
      //            << torsion.multiplicity << endr;
    }

 
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // LennardJonesParameters
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // get some array sizes
    unsigned int sizeAtomTypes  = topo->atomTypes.size();
    unsigned int sizeNonbondeds = par.nonbondeds.size();
    unsigned int sizeNbfixs     = par.nbfixs.size();
    unsigned int sizeLJbank     = par.nonbondeds[0].epsilon.size();

    // create the bank of Lennard-Jones parameter tables
    topo->isgLJParms.createTables(sizeLJbank,sizeAtomTypes); 

  
    // loop over all atom types
    for(unsigned int i=0;i<sizeAtomTypes;++i) {

      // sentinel
      int foundI = 0;
      unsigned int AtomI = sizeNonbondeds;

      // loop over all nonbonded types in the PAR object
      for(unsigned int k=0;k<sizeNonbondeds;++k) {

	// determine if this PAR nonbonded type is the same as the current atom type (i)
	int ok = equalWildcard(par.nonbondeds[k].atom,topo->atomTypes[i].name);

	// if we have found the current atom type in the PAR object the record its index (k)
	if(ok > foundI){
	  AtomI = k;
	  foundI = ok;
	}

	// if we have found the current atom type then exit the loop
	if(foundI > 1)
	  break; 
      } // end loop over k

      // report an error if we could not find atom type (i) in the PAR object
      if(foundI <= 0)
	report << error <<"Could not find matching parameter nonbonded of atom \'"<<topo->atomTypes[i].name<<"\'."<<endr;

      // loop over all other (i != j) atom types
      for(unsigned int j=i;j<sizeAtomTypes;++j) {

	// sentinel
	int foundJ = 0;
	unsigned int AtomJ = sizeNonbondeds;

	// loop over all nonbonded types in the PAR object
	for(unsigned int k=0;k<sizeNonbondeds;++k) {

	  // determine if this PAR nonbonded type is the same as the current atom type (j)
	  int ok = equalWildcard(par.nonbondeds[k].atom,topo->atomTypes[j].name);

	  // if we have found the current atom type in the PAR object the record its index (k)
	  if(ok > foundJ){
	    AtomJ = k;
	    foundJ = ok;
	  }

	  // if we have found the current atom type then exit the loop
	  if(foundJ > 1)
	    break;
	}

	// report an error if we could not find atom type (j) in the PAR object
	if(foundJ<=0)
	  report << error <<"Could not find matching parameter nonbonded of atom \'"<<topo->atomTypes[j].name<<"\'."<<endr;

	// initialize a set of LJ parameters
 	LennardJonesParameters paramsij;

	// loop over all possible identities of atom i
	for (unsigned int identity_I=0; identity_I < sizeLJbank; identity_I++) {

	  // loop over all possible identities for atom j
	  for (unsigned int identity_J = 0; identity_J < sizeLJbank; identity_J++) {

	    // calculate the (Charmm28) minimum energy distances
	    Real sigma_i = par.nonbondeds[AtomI].sigma[identity_I];
	    Real sigma_j = par.nonbondeds[AtomJ].sigma[identity_J];
	    Real sigma14_i = par.nonbondeds[AtomI].sigma14[identity_I];
	    Real sigma14_j = par.nonbondeds[AtomJ].sigma14[identity_J];
	
	    // calculate the minimum LJ pair energy
	    Real epsilon_i = par.nonbondeds[AtomI].epsilon[identity_I];
	    Real epsilon_j = par.nonbondeds[AtomJ].epsilon[identity_J];
	    Real epsilon14_i = par.nonbondeds[AtomI].epsilon14[identity_I];
	    Real epsilon14_j = par.nonbondeds[AtomJ].epsilon14[identity_J];
	
	    // calculate the ij pair minimum distances and energies
	    Real r_ij = sigma_i + sigma_j;
	    Real e_ij = sqrt(epsilon_i * epsilon_j);
	    Real r14_ij = sigma14_i + sigma14_j;
	    Real e14_ij = sqrt(epsilon14_i * epsilon14_j);
	
	    // stores these numbers into the paramsij object
	    paramsij.A = power<12>(r_ij) * e_ij;
	    paramsij.B = 2 * power<6>(r_ij) * e_ij;
	    paramsij.A14 = power<12>(r14_ij) * e14_ij;
	    paramsij.B14 = 2 * power<6>(r14_ij) * e14_ij;

	    topo->isgLJParms.set(identity_I,identity_J, i, j, paramsij);

	  } // end loop over identity_J
	} // end loop over identity_I
      } // end loop over atom types (j)
    } // end loop over atom types (i)
 

    // NbFix
    // loop over all Nbfix types
    for(unsigned int k=0;k<sizeNbfixs;++k) {

      //report << par.nbfixs[k].atom1 << ","<<par.nbfixs[k].atom2<<endr;

      int ti = 0;
      int tj = 0;
      unsigned int bi = sizeNbfixs;
      unsigned int bj = sizeNbfixs;

      for(unsigned int i=0;i<sizeAtomTypes;++i){
	int ok = equalWildcard(par.nbfixs[k].atom1,topo->atomTypes[i].name);
	if(ok > ti){
	  bi = i;
	  ti = ok;
	}
	if(ti > 1)
	  break;
      }

      if(ti <=0) // ???
	continue;

      for(unsigned int j=0;j<sizeAtomTypes;++j){
	int ok = equalWildcard(par.nbfixs[k].atom2,topo->atomTypes[j].name);
	if(ok > tj){
	  bj = j;
	  tj = ok;
	}
	if(tj > 1)
	  break;
      }
	
      if(tj<=0)
	report << error <<"Could not find matching parameter nbfix of atoms \'"<<par.nbfixs[k].atom1<<"\' - '"<<par.nbfixs[k].atom2<<"\'."<<endr;
      		
      LennardJonesParameters paramsij;

      paramsij.A = par.nbfixs[k].a;
      paramsij.B = par.nbfixs[k].b;
      paramsij.A14 = par.nbfixs[k].a14;
      paramsij.B14 = par.nbfixs[k].b14;
      topo->lennardJonesParameters.set(bi, bj, paramsij);

    } // end loop over NbFix types

    // store the molecule information
    buildMoleculeTable(topo);

    // build the nonbonded exclusion table
    buildExclusionTable(topo,topo->exclude);

    // build the list of bonds, angles, etc. belonging to each molecule
    buildMoleculeBondingLists(topo);

    // determine the time-zero identity of each molecule
    for(unsigned int i=0; i<topo->molecules.size(); i++) {
      
      // the ID# of the first atom in this molecule
      int firstAtom = topo->molecules[i][0];

      // get the identity number from the psf file
      topo->molecules[i].type = psf.atoms[firstAtom].identity;
      topo->molecules[i].newtype = psf.atoms[firstAtom].identity;
    }

  } // end build topology function


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  //  buildMoleculeTable
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  void buildMoleculeTable(GenericTopology *topo){
 
    // *** First we clear all molecules ***
    topo->molecules.clear();

    const unsigned int numAtoms = topo->atoms.size();

    // *** Collecting all possible connections, building the graph ***
    vector<vector<int> > graph(numAtoms,vector<int>());
    set<pair<int,int> > pairs;
    // *** Bonds ***
    for(unsigned int i=0;i<topo->bonds.size();++i){
      int a1 = topo->bonds[i].atom1;
      int a2 = topo->bonds[i].atom2;
      graph[a1].push_back(a2);
      graph[a2].push_back(a1);
      pairs.insert(pair<int,int>(std::min(a1,a2),std::max(a1,a2)));
    }
    unsigned int count = pairs.size();

    // *** Angles ***
    for(unsigned int i=0;i<topo->angles.size();++i){
      int a1 = topo->angles[i].atom1;
      int a2 = topo->angles[i].atom2;
      int a3 = topo->angles[i].atom3;
      graph[a1].push_back(a2);
      graph[a1].push_back(a3);
      graph[a2].push_back(a1);
      graph[a2].push_back(a3);
      graph[a3].push_back(a1);
      graph[a3].push_back(a2);
      pairs.insert(pair<int,int>(std::min(a1,a2),std::max(a1,a2)));
      pairs.insert(pair<int,int>(std::min(a3,a2),std::max(a3,a2)));
    }
    
    if(count < pairs.size())
      report << hint << "Angles added "<<pairs.size()-count<<" new bond(s)." <<endr;
    count = pairs.size();

    // *** Dihedrals ***
    for(unsigned int i=0;i<topo->dihedrals.size();++i){
      int a1 = topo->dihedrals[i].atom1;
      int a2 = topo->dihedrals[i].atom2;
      int a3 = topo->dihedrals[i].atom3;
      int a4 = topo->dihedrals[i].atom4;
      graph[a1].push_back(a2);
      graph[a1].push_back(a3);
      graph[a1].push_back(a4);
      graph[a2].push_back(a1);
      graph[a2].push_back(a3);
      graph[a2].push_back(a4);
      graph[a3].push_back(a1);
      graph[a3].push_back(a2);
      graph[a3].push_back(a4);
      graph[a4].push_back(a1);
      graph[a4].push_back(a2);
      graph[a4].push_back(a3);
      pairs.insert(pair<int,int>(std::min(a1,a2),std::max(a1,a2)));
      pairs.insert(pair<int,int>(std::min(a3,a2),std::max(a3,a2)));
      pairs.insert(pair<int,int>(std::min(a3,a4),std::max(a3,a4)));
    }
    if(count < pairs.size())
      report << hint << "Dihedrals added "<<pairs.size()-count<<" new bond(s)." <<endr;
    count = pairs.size();
  
    // *** Impropers ***
    set<pair<int,int> > pairsAddImpropers;
    // Impropers are defined over the bonds 1-2,1-3,1-4 or 4-1,4-3,4-2
    // but MTorsionSystemForce computes distances betweeen 1-2,2-3,3-4
    // we have to take care about these differences ...
    for(unsigned int i=0;i<topo->impropers.size();++i){
      int a1 = topo->impropers[i].atom1;
      int a2 = topo->impropers[i].atom2;
      int a3 = topo->impropers[i].atom3;
      int a4 = topo->impropers[i].atom4;
      graph[a1].push_back(a2);
      graph[a1].push_back(a3);
      graph[a1].push_back(a4);
      graph[a2].push_back(a1);
      graph[a2].push_back(a3);
      graph[a2].push_back(a4);
      graph[a3].push_back(a1);
      graph[a3].push_back(a2);
      graph[a3].push_back(a4);
      graph[a4].push_back(a1);
      graph[a4].push_back(a2);
      graph[a4].push_back(a3);
      pair<int,int> p0(std::min(a1,a2),std::max(a1,a2));
      pair<int,int> p1(std::min(a1,a3),std::max(a1,a3));
      pair<int,int> p2(std::min(a1,a4),std::max(a1,a4));
      pair<int,int> p3(std::min(a2,a3),std::max(a2,a3));
      pair<int,int> p4(std::min(a2,a4),std::max(a2,a4));
      pair<int,int> p5(std::min(a3,a4),std::max(a3,a4));
      int j0 = 0;
      int j1 = 0;
      int j2 = 0;
      int j3 = 0;
      int j4 = 0;
      int j5 = 0;
      if(pairs.find(p0) != pairs.end()) j0++;
      if(pairs.find(p1) != pairs.end()) j1++;
      if(pairs.find(p2) != pairs.end()) j2++;
      if(pairs.find(p3) != pairs.end()) j3++;
      if(pairs.find(p4) != pairs.end()) j4++;
      if(pairs.find(p5) != pairs.end()) j5++;
      if(j0+j1+j2+j3+j4+j5 < 3){
	pairs.insert(p0);
	pairs.insert(p1);
	pairs.insert(p2);	
      }
      pairsAddImpropers.insert(p0);
      pairsAddImpropers.insert(p3);
      pairsAddImpropers.insert(p5);
    }
    if(count < pairs.size())
      report << hint << "Impropers added "<<pairs.size()-count<<" new bond(s)." <<endr;

    // Now add the improper pairs
    for(set<pair<int,int> >::const_iterator i= pairsAddImpropers.begin();i != pairsAddImpropers.end();++i)
      pairs.insert(*i);


    count = pairs.size();
    //report << hint << count << endr;
    // To keep track which atoms already have been added 
    // to molecules.
    vector<bool> unused(numAtoms,true);


    // Recursively finding the atoms beloning to a molecule
    for(unsigned int i=0;i<numAtoms;++i){
      vector<int> v;
      vector<PairInt> p;
      findNextNeighbor(i,v,p,unused,graph,pairs);
      if(v.size()>0){
	// sort(v.begin(),v.end());
	// add this atom list to the molecules array
	Molecule mol;
	mol.atoms = v;
	mol.pairs = p;
	topo->molecules.push_back(mol);	
      }
    }

    // Uncomment to sort descending after size()
    // sort(topo->molecules.begin(),topo->molecules.end(),cmpSize);

    // Look up table for atoms
    const string h("H");
    const string o("O");
    for(unsigned int i=0;i<topo->molecules.size();++i){
      Real mass = 0.0;
      const vector<int>& mol = topo->molecules[i].atoms;
      for(unsigned int j=0;j<mol.size();++j){
	int k = mol[j];
	topo->atoms[k].molecule = i;
	mass += topo->atoms[k].scaledMass;
      }
      topo->molecules[i].mass = mass;
      topo->molecules[i].water = (mol.size() == 3 &&
				  ((topo->atomTypes[topo->atoms[mol[0]].type].symbolName == h &&
				    topo->atomTypes[topo->atoms[mol[1]].type].symbolName == h &&
				    topo->atomTypes[topo->atoms[mol[2]].type].symbolName == o) ||
				   (topo->atomTypes[topo->atoms[mol[0]].type].symbolName == h &&
				    topo->atomTypes[topo->atoms[mol[1]].type].symbolName == o &&
				    topo->atomTypes[topo->atoms[mol[2]].type].symbolName == h) ||
				   (topo->atomTypes[topo->atoms[mol[0]].type].symbolName == o &&
				    topo->atomTypes[topo->atoms[mol[1]].type].symbolName == h &&
				    topo->atomTypes[topo->atoms[mol[2]].type].symbolName == h)));

    }

#if defined(DEBUG_PRINT_MOLECULETABLE)
    report<< plain << endl << "[buildMoleculeTable]: molecule table printout:" << endl;
    for(int i=0;i<topo->molecules.size();++i){
      for(int j=0;j<topo->molecules[i].size();++j){
	report << topo->molecules[i][j]<<" ";
      }
      report << endl;
    }
    report << endr;
#endif
  }


  //____________________________________________________________findNextNeighbor
  void findNextNeighbor(int a, vector<int>& v, vector<PairInt>& p, vector<bool>& unused, 
			const vector<vector<int> >& graph, set<PairInt>& pairs){
    if(unused[a]){
      v.push_back(a);
      unused[a] = false;
      for(unsigned int i=0;i<graph[a].size();++i){
	set<PairInt>::iterator itr = pairs.find(PairInt(std::min(a,graph[a][i]),std::max(a,graph[a][i])));
	if(itr != pairs.end()){
	  p.push_back(PairInt(a,graph[a][i]));
	  pairs.erase(itr);
	}
	findNextNeighbor(graph[a][i],v,p,unused,graph,pairs);
      }
    }
  }


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  //  buildExclusionTable
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  void buildExclusionTable(GenericTopology* topo, const ExclusionType& exclusionType) {

    if(!exclusionType.valid())
      report << error <<"[buildExclusionTable()] Exclusion type not defined/valid."<<endr;

    topo->exclude = exclusionType;

    const int numBonds = topo->bonds.size();
    const int numAtoms = topo->atoms.size();
    int a1,a2,a3;

    // *** Reformatting bond list ***
    vector< set<int> > bondsParsed(numAtoms);
    for (int i = 0; i < numBonds; ++i) {
      a1=topo->bonds[i].atom1;
      a2=topo->bonds[i].atom2;
      bondsParsed[a1].insert(a2);
      bondsParsed[a2].insert(a1);
    }

    // *** Building exclusion list ***
    topo->exclusions.resize(numAtoms);
    topo->exclusions.clear();
    list<ExclusionPair> one3s;
    int exclusionsInserted = 0;
    if(exclusionType!=ExclusionType::NONE) {
      for(a2 = 0; a2 < numAtoms; ++a2) {
	for(set<int>::iterator i=bondsParsed[a2].begin();i!=bondsParsed[a2].end();++i) {
	  a1=*i;
	  if(a1<a2) {
	    topo->exclusions.add(a1,a2,EXCLUSION_FULL);
	    ++exclusionsInserted;
	  }
	  if(exclusionType==ExclusionType::ONE3||exclusionType==ExclusionType::ONE4||exclusionType==ExclusionType::ONE4MODIFIED) {
	    for(set<int>::iterator j=i; j!=bondsParsed[a2].end();++j) {
	      a3=*j;
	      topo->exclusions.add(a1,a3,EXCLUSION_FULL);
	      ++exclusionsInserted;
	      one3s.push_front(ExclusionPair(a1,a3));
	      one3s.push_front(ExclusionPair(a3,a1));
	    }
	  }
	}
      }
      if(exclusionType==ExclusionType::ONE4||exclusionType==ExclusionType::ONE4MODIFIED) {
	ExclusionClass currentType = EXCLUSION_NONE;
	if(exclusionType==ExclusionType::ONE4)
	  currentType=EXCLUSION_FULL;
	if(exclusionType==ExclusionType::ONE4MODIFIED)
	  currentType=EXCLUSION_MODIFIED;     
	int curSrc=0;
	for(list<ExclusionPair>::iterator i=one3s.begin(); i!=one3s.end(); ++i,++curSrc) {
	  a1=i->a1;
	  a2=i->a2;
	  for(set<int>::iterator j=bondsParsed[a2].lower_bound(a1);j!=bondsParsed[a2].end();++j) {
	    a3=*j;
	    if(!topo->exclusions.check(a1,a3)) {
	      topo->exclusions.add(a1,a3,currentType); 
	      ++exclusionsInserted;
	    }
	  }
	}
      }
    }
  }


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  //  buildMoleculeBondingLists
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  void buildMoleculeBondingLists(GenericTopology* topo) {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // bond lists
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all bonds in the topology
    for (unsigned int i=0; i<topo->bonds.size(); i++) {

      // get the ID#s of the two atoms in this bond
      int atom1 = topo->bonds[i].atom1;
      int atom2 = topo->bonds[i].atom2;

      // loop sentinels
      bool found1 = false;
      bool found2 = false;
      unsigned int m;

      // loop over all molecules
      for (m=0; m<topo->molecules.size(); m++) {

	// loop over all atoms in the atom list for this molecule
	for (unsigned int a=0; a<topo->molecules[m].atoms.size(); a++) {

	  if (atom1 == topo->molecules[m].atoms[a]) {found1 = true;}
	  else if (atom2 == topo->molecules[m].atoms[a]) {found2 = true;}

	  if (found1 && found2) {break;}
	} // end atom list loop

	if (found1 && found2) {break;}
      } // end loop over molecules

      // store this bond index number in the molecule's bond list
      topo->molecules[m].bondList.push_back(i);
    } // end loop over bonds


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // angle lists
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all angles in the topology
    for (unsigned int i=0; i<topo->angles.size(); i++) {

      // get the ID#s of the three atoms in this angle
      int atom1 = topo->angles[i].atom1;
      int atom2 = topo->angles[i].atom2;
      int atom3 = topo->angles[i].atom3;

      // loop sentinels
      bool found1 = false;
      bool found2 = false;
      bool found3 = false;
      unsigned int m;

      // loop over all molecules
      for (m=0; m<topo->molecules.size(); m++) {

	// loop over all atoms in the atom list for this molecule
	for (unsigned int a=0; a<topo->molecules[m].atoms.size(); a++) {

	  if (atom1 == topo->molecules[m].atoms[a]) {found1 = true;}
	  else if (atom2 == topo->molecules[m].atoms[a]) {found2 = true;}
	  else if (atom3 == topo->molecules[m].atoms[a]) {found3 = true;}

	  if (found1 && found2 && found3) {break;}
	} // end atom list loop

	if (found1 && found2 && found3) {break;}
      } // end loop over molecules

      // store this angle index number in the molecule's angle list
      topo->molecules[m].angleList.push_back(i);
    } // end loop over angles

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // dihedral lists
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all diherals in the topology
    for (unsigned int i=0; i<topo->dihedrals.size(); i++) {

      // get the ID#s of the four atoms in this dihedral
      int atom1 = topo->dihedrals[i].atom1;
      int atom2 = topo->dihedrals[i].atom2;
      int atom3 = topo->dihedrals[i].atom3;
      int atom4 = topo->dihedrals[i].atom4;

      // loop sentinels
      bool found1 = false;
      bool found2 = false;
      bool found3 = false;
      bool found4 = false;
      unsigned int m;

      // loop over all molecules
      for (m=0; m<topo->molecules.size(); m++) {

	// loop over all atoms in the atom list for this molecule
	for (unsigned int a=0; a<topo->molecules[m].atoms.size(); a++) {

	  if (atom1 == topo->molecules[m].atoms[a]) {found1 = true;}
	  else if (atom2 == topo->molecules[m].atoms[a]) {found2 = true;}
	  else if (atom3 == topo->molecules[m].atoms[a]) {found3 = true;}
	  else if (atom4 == topo->molecules[m].atoms[a]) {found4 = true;}

	  if (found1 && found2 && found3 && found4) {break;}
	} // end atom list loop

	if (found1 && found2 && found3 && found4) {break;}
      } // end loop over molecules

      // store this dihedral index number in the molecule's dihderal list
      topo->molecules[m].dihedralList.push_back(i);
    } // end loop over dihedrals


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // improper lists
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all impropers in the topology
    for (unsigned int i=0; i<topo->impropers.size(); i++) {

      // get the ID#s of the four atoms in this improper
      int atom1 = topo->impropers[i].atom1;
      int atom2 = topo->impropers[i].atom2;
      int atom3 = topo->impropers[i].atom3;
      int atom4 = topo->impropers[i].atom4;

      // loop sentinels
      bool found1 = false;
      bool found2 = false;
      bool found3 = false;
      bool found4 = false;
      unsigned int m;

      // loop over all molecules
      for (m=0; m<topo->molecules.size(); m++) {

	// loop over all atoms in the atom list for this molecule
	for (unsigned int a=0; a<topo->molecules[m].atoms.size(); a++) {

	  if (atom1 == topo->molecules[m].atoms[a]) {found1 = true;}
	  else if (atom2 == topo->molecules[m].atoms[a]) {found2 = true;}
	  else if (atom3 == topo->molecules[m].atoms[a]) {found3 = true;}
	  else if (atom4 == topo->molecules[m].atoms[a]) {found4 = true;}

	  if (found1 && found2 && found3 && found4) {break;}
	} // end atom list loop

	if (found1 && found2 && found4) {break;}
      } // end loop over molecules

      // store this bond index number in the molecule's improper list
      topo->molecules[m].improperList.push_back(i);
    } // end loop over impropers
  } // end function buildMoleculeBonding Lists
  
} // end class

