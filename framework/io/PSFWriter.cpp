#include "PSFWriter.h"

#include "Report.h"
#include "stringutilities.h"
#include "systemutilities.h"

#include <iomanip>

using std::string;
using std::endl;
using std::setprecision;
using std::setw;
using std::left;
using std::right;
using std::setiosflags;
using namespace ProtoMol::Report;

namespace ProtoMol {

  //_________________________________________________________________PSFWriter

  PSFWriter::PSFWriter():Writer(){}

  PSFWriter::PSFWriter(const std::string& filename):Writer(filename){}

  PSFWriter::~PSFWriter(){}

  bool PSFWriter::write(const PSF& psf){
    if(!open())
      return false;

    // Create some comments at the beginning
    myFile << "PSF" << endl
	   << endl
	   << "       5 !NTITLE" << endl
	   << " REMARKS FILENAME=" << myFilename <<" by "<<getUserName()<< endl
	   << " REMARKS ProtoMol (built on "<< __DATE__ << " at " << __TIME__<< ")" << endl
	   << " REMARKS This .psf file was created by PSFWriter" << endl
	   << " REMARKS It was not manually assembled" << endl
	   << " REMARKS "<<myComment<< endl
	   << endl;

    unsigned int count = psf.atoms.size();
    myFile << setw(8) << right << count << " !NATOM" << endl;
    for(unsigned int i=0;i<count;++i){
      myFile << setw(8) << right
	     << i+1
	     << " " << setw(4) << left
	     << (psf.atoms[i].seg_id+"    ").substr(0,4)
	     << " " << setw(4)
	     << psf.atoms[i].residue_sequence
	     << " " << setw(4)
	     << psf.atoms[i].residue_name
	     << " " << setw(4)
	     << psf.atoms[i].atom_name
	     << " " << setw(5)
	     << psf.atoms[i].atom_type
	     << " " << setw(15) << setprecision(7) 
	     << psf.atoms[i].charge
	     << " " << setw(7) << setprecision(6) 
	     << psf.atoms[i].mass
	     << " " << setw(11) << right
	     << psf.atoms[i].identity
	     << endl;
    }
    count = psf.bonds.size();
    myFile << endl << setw(8) << right << count << " !NBOND: bonds" ;
    for(unsigned int i=0;i<count;++i){
      myFile << ((i % 4 == 0)?"\n":"") << setw(8) 
	     << psf.bonds[i].atom1
	     << setw(8)
	     << psf.bonds[i].atom2;
    }
    count = psf.angles.size();
    myFile << endl << endl << setw(8) << right << count << " !NTHETA: angles" ;
    for(unsigned int i=0;i<count;++i){
      myFile << ((i % 3 == 0)?"\n":"") << setw(8) 
	     << psf.angles[i].atom1
	     << setw(8)
	     << psf.angles[i].atom2
	     << setw(8)
	     << psf.angles[i].atom3;
    }
 
    count = psf.dihedrals.size();
    myFile << endl << endl << setw(8) << right << count << " !NPHI: dihedrals" ;
    for(unsigned int i=0;i<count;++i){
      myFile << ((i % 2 == 0)?"\n":"") << setw(8) 
	     << psf.dihedrals[i].atom1
	     << setw(8)
	     << psf.dihedrals[i].atom2
	     << setw(8)
	     << psf.dihedrals[i].atom3
	     << setw(8)
	     << psf.dihedrals[i].atom4;
    }
    count = psf.impropers.size();
    myFile << endl << endl << setw(8) << right << count << " !NIMPHI: impropers" ;
    for(unsigned int i=0;i<count;++i){
      myFile << ((i % 2 == 0)?"\n":"") << setw(8) 
	     << psf.impropers[i].atom1
	     << setw(8)
	     << psf.impropers[i].atom2
	     << setw(8)
	     << psf.impropers[i].atom3
	     << setw(8)
	     << psf.impropers[i].atom4;
    }
    count = psf.donors.size();
    myFile << endl << endl << setw(8) << right << count << " !NDON: donors" ;
    for(unsigned int i=0;i<count;++i){
      myFile << ((i % 4 == 0)?"\n":"") << setw(8) 
	     << psf.donors[i].atom1
	     << setw(8)
	     << psf.donors[i].atom2;
    }
    count = psf.acceptors.size();
    myFile << endl << endl << setw(8) << right << count << " !NACC: acceptors" ;
    for(unsigned int i=0;i<count;++i){
      myFile << ((i % 4 == 0)?"\n":"") << setw(8) 
	     << psf.acceptors[i].atom1
	     << setw(8)
	     << psf.acceptors[i].atom2;
    }
    myFile << endl << endl << setw(8) << right 
	   << 0  // 05/22/01 - not sure what the first number means here
	   << " !NNB" ;
    count = psf.nonbondeds.size();
    for(unsigned int i=0;i<count;++i){
      myFile << ((i % 8 == 0)?"\n":"") << setw(8) 
	     << 0; // psf.nonbondeds[i].atom1
    }
    count = psf.ngrp.size();
    myFile << endl << endl << setw(8) << right << count << setw(8) << 0 << " !NGRP";
    for(unsigned int i=0;i<count;++i){
      myFile << ((i % 3 == 0)?"\n":"") << setw(8) 
	     << psf.ngrp[i].atom1
	     << setw(8)
	     << psf.ngrp[i].atom2
	     << setw(8)
	     << psf.ngrp[i].atom3;
    }
    myFile << endl;
    close();
    return !myFile.fail();
  }

  PSFWriter& operator<<(PSFWriter& psfWriter, const PSF& psf){
    psfWriter.write(psf);
    return psfWriter;
  }
}
