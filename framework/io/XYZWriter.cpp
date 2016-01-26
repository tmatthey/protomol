#include "XYZWriter.h"

#include <iomanip>

#include "Report.h"
#include "stringutilities.h"
#include "systemutilities.h"

using std::string;
using std::endl;
using std::setprecision;
using namespace ProtoMol::Report;

namespace ProtoMol {

  //_________________________________________________________________XYZWriter

  XYZWriter::XYZWriter():Writer(){}

  XYZWriter::XYZWriter(const std::string& filename):Writer(filename){}

  bool XYZWriter::write(const XYZ& xyz){
    return write(xyz.coords,xyz.names);
  }


  bool XYZWriter::write(const Vector3DBlock& coords, 
			const std::vector<Atom>& atoms, 
			const std::vector<AtomType>& atomTypes){
    const unsigned int count = atoms.size();
    std::vector<std::string> names(count);
    for(unsigned int i=0;i<count;++i)
      names[i] = atomTypes[atoms[i].type].name;
    return write(coords,names);
  }

  bool XYZWriter::write(const Vector3DBlock& coords, const std::vector<std::string>& names){
    if (!open())
      return false;

    const unsigned int count = coords.size();
    if(names.size() != count)
      report << error << "[XYZWriter::write]"
	     << " Coorindate and atom name size are not equal."<< endr;
    
    // First, write the number of atoms
    myFile << count << endl;

    // Comment
    myFile << "!ProtoMol (built on "<< __DATE__ << " at " << __TIME__<< ") generated this XYZ file by "<<getUserName()<<". " <<myComment<< endl;

    // Write atoms
    myFile << setprecision(15); // This should be some FLT_DIG or DBL_DIG ...
    for(unsigned int i=0;i<count;++i){
      myFile << names[i] << "\t";
      myFile.width(24);
      myFile << coords[i].x;
      myFile.width(24);
      myFile << coords[i].y;
      myFile.width(24);
      myFile << coords[i].z;
      myFile << endl;
    }
    close();
    return !myFile.fail();
  }

  XYZWriter& operator<<(XYZWriter& xyzWriter, const XYZ& xyz){
    xyzWriter.write(xyz.coords,xyz.names);
    return xyzWriter;
  }
}
