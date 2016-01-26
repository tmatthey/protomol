#include "PDB.h"

using namespace ProtoMol::Report;

namespace ProtoMol {
  //___________________________________________________________________PDB
  PDB::Atom::Atom():elementType(""),	      
		    elementNum(0),	      
		    elementName(""),	      
		    altLoc(""),	      	      
		    residueName(""),	      
		    chainID(""),	      	      
		    residueNum(0),	      
		    insertionCode(""), 
		    occupancy(0.0),	      	      
		    tempFactor(0.0),  
		    segID(""),	      	      
		    symbol(""),	      
		    charge(""),	      	      
		    hvyAtomGrpsize(0){}	      



  PDB::Atom::Atom(std::string etype,  
		  int  anum,	
		  std::string ename,  
		  std::string altloc,	
		  std::string rname,  
		  std::string chain,	
		  int  rnum,	
		  std::string insertion,
		  Real occ,	
		  Real tf,	
		  std::string segname,	
		  std::string symname,	
		  std::string c,	
		  int ha):elementType(etype),	      
			  elementNum(anum),	      
			  elementName(ename),	      
			  altLoc(altloc),	      
			  residueName(rname),	      
			  chainID(chain),	      
			  residueNum(rnum),	      
			  insertionCode(insertion),  
			  occupancy(occ),	      
			  tempFactor(tf),            
			  segID(segname),	      
			  symbol(symname),	      
			  charge(c),	      
			  hvyAtomGrpsize(ha){}	      

  PDB::Ter::Ter():elementType(""),	      
		  elementNum(0),	      
		  residueName(""),	      
		  chainID(""),	      	      
		  residueNum(0),	      
		  insertionCode(""){}	      



  PDB::Ter::Ter(std::string etype,  
		int  anum,	
		std::string rname,  
		std::string chain,	
		int  rnum,	
		std::string insertion):elementType(etype),	      
				       elementNum(anum),	      
				       residueName(rname),	      
				       chainID(chain),	      
				       residueNum(rnum),	      
				       insertionCode(insertion){}	      

  void PDB::clear(){
    coords.clear();
    atoms.clear();
    ters.clear();
  }

  // Atom
  //
  // COLUMNS   DATA TYPE     FIELD       DEFINITION
  // ------------------------------------------------------------------------
  //  1 -  6   Record name   "ATOM  "
  //  7 - 11   Integer       serial      Atom serial number.
  // 13 - 16   Atom          name        Atom name.
  // 17 - 17   Character     altLoc      Alternate location indicator.
  // 18 - 20   Residue name  resName     Residue name.
  // 22 - 22   Character     chainID     Chain identifier.
  // 23 - 26   Integer       resSeq      Residue sequence number.
  // 27 - 27   AChar         iCode       Code for insertion of residues.
  // 31 - 38   Real(8.3)     x           Orthogonal coordinates for X in
  //                                     Angstroms.
  // 39 - 46   Real(8.3)     y           Orthogonal coordinates for Y in
  //                                     Angstroms.
  // 47 - 54   Real(8.3)     z           Orthogonal coordinates for Z in
  //                                     Angstroms.
  // 55 - 60   Real(6.2)     occupancy   Occupancy.
  // 61 - 66   Real(6.2)     tempFactor  Temperature factor.
  // 73 - 76   LString(4)    segID       Segment identifier, left-justified.
  // 77 - 78   LString(2)    element     Element symbol, right-justified.
  // 79 - 80   LString(2)    charge      Charge on the atom.
  // NOTE: The PDB says the length of the residue name is only 3 characters
  //  whereas XPLOR allows 4 character names.  We choose 4 for compatability
  //  with both systems (since we never change the length, we you give us is
  //  what we use)


  MyStreamer& operator<< (MyStreamer& OS, const PDB::Atom & p) {
    OS <<p.elementType<<","<<p.elementNum<<","<<p.elementName<<","<<p.altLoc<<","<<p.residueName<<","<<p.chainID<<","<<p.residueNum<<","<<p.insertionCode<<","<<p.occupancy<<","<<p.tempFactor<<","<<p.segID<<","<<p.symbol<<","<<p.charge;
    return OS;
  }

}
