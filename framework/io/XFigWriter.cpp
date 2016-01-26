#include "XFigWriter.h"

#include <iomanip>

#include "Report.h"
#include "pmconstants.h"
#include "stringutilities.h"
#include "systemutilities.h"

using std::string;
using std::endl;
using std::setprecision;
using namespace ProtoMol::Report;

namespace ProtoMol {

  static const Real defaultScale = 1000.0;
  static const int xFigXMax = 12500;
  static const int xFigYMax = 9700;
  //static const int xFigXYMin = (xFigXMax < xFigYMax ? xFigXMax : xFigYMax);
  static const int xFigColorsN = 32;
  static const int xFigColors[xFigColorsN] = {0,20,13,10,31,1,2,3,4,5,6,7,8,9,11,12,14,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30};
  static const int xFigBWN = 4;
  static const std::string xFigBW1[xFigBWN] ={"1 3 0 1 7 0",
					      "1 3 0 2 0 0",
					      "1 3 1 2 0 0",
					      "1 3 2 2 0 0"};
  static const std::string xFigBW2[xFigBWN] ={"0 20 0.000 1 0.0000",
					      "0 -1 0.000 1 0.0000",
					      "0 -1 4.500 1 0.0000",
					      "0 -1 4.500 1 0.0000"};

  //_________________________________________________________________XFigWriter

  XFigWriter::XFigWriter():Writer(),
			   myRadius(1.0),
			   myFromZ(-Constant::MAXREAL),
			   myToZ(Constant::MAXREAL),
			   myColor(true),
			   myAxes(true),
			   myLegend(true),
			   myAutofit(true),
			   myMat(Matrix3by3().getIdentity()){}

  XFigWriter::XFigWriter(const std::string& filename):Writer(filename),
						      myRadius(1.0),
						      myFromZ(-Constant::MAXREAL),
						      myToZ(Constant::MAXREAL),
						      myColor(true),
						      myAxes(true),
						      myLegend(true),
						      myAutofit(true),
						      myMat(Matrix3by3().getIdentity()){}
  
  bool XFigWriter::write(const XYZ& xyz){
    return write(xyz.coords,xyz.names);
  }


  bool XFigWriter::write(const Vector3DBlock& coords, const std::vector<std::string>& names){
    if (!open())
      return false;

    const unsigned int count = coords.size();
    if(names.size() != count)
      report << error << "[XFigWriter::write]"
	     << " Coorindate and atom name size are not equal."<< endr;
    
    XYZ xyz;
    Real minx =  Constant::MAXREAL;
    Real maxx = -Constant::MAXREAL;
    Real miny =  Constant::MAXREAL;
    Real maxy = -Constant::MAXREAL;
    Real minz =  Constant::MAXREAL;
    Real maxz = -Constant::MAXREAL;
    std::map<std::string,int> atomTypesMap;

    for(unsigned int i=0;i<coords.size();i++){
      Vector3D v(myMat*coords[i]);
      v.y = -v.y;
      if( v.z >= myFromZ && v.z <= myToZ){
	minx = std::min(minx,v.x);
	maxx = std::max(maxx,v.x);
	miny = std::min(miny,v.y);
	maxy = std::max(maxy,v.y);
	minz = std::min(minz,v.z);
	maxz = std::max(maxz,v.z);
	xyz.coords.push_back(v);
	xyz.names.push_back(names[i]);
	if(atomTypesMap.find(names[i]) == atomTypesMap.end()){
	  int size = atomTypesMap.size();
	  atomTypesMap[names[i]] = size;
	}
      }
    }
    Real dx = std::max(maxx-minx,Constant::EPSILON);
    Real dy = std::max(maxy-miny,Constant::EPSILON);
    Real dz = std::max(maxz-minz,Constant::EPSILON);

    Real scale = defaultScale;
    int r = static_cast<int>(myRadius*scale/10);
    if(myAutofit){
      scale = std::min(xFigXMax/dx,xFigYMax/dy);
      r = std::max(static_cast<int>(myRadius*scale/10),25);
    }
    myFile << "#FIG 3.2\nLandscape\nCenter\nMetric\nLetter  \n100.00\nSingle\n-2\n";
    myFile << "#ProtoMol (built on "<< __DATE__ << " at " << __TIME__<< ") generated this XFig file by "<<getUserName()<<".\n";    
    myFile << "# "<<myComment << "\n";
    myFile << "1200 2\n";

    for(unsigned int i=0;i<xyz.coords.size();i++){
      int x = static_cast<int>((xyz.coords[i].x-minx)*scale);
      int y = static_cast<int>((xyz.coords[i].y-miny)*scale);
      int z = static_cast<int>(1+(xyz.coords[i].z-minz)*998/dz);
      int col = (atomTypesMap[xyz.names[i]]%(myColor?xFigColorsN:xFigBWN));
      if(myColor){
	myFile <<"1 3 0 1 "<<xFigColors[col]<<" "<<xFigColors[col]<<" "<<z<<" 0 41 0.000 1 0.0000 "
	       << x <<" "
	       << y <<" "
	       << r<<" "
	       << r <<" "
	       << x <<" "
	       << y <<" "
	       << x+r <<" "
	       << y << std::endl;
      }
      else {
	myFile <<xFigBW1[col]<<" "<<z<<" "<<xFigBW2[col]<<" "
	       << x <<" "
	       << y <<" "
	       << r <<" "
	       << r <<" "
	       << x <<" "
	       << y <<" "
	       << x+r <<" "
	       << y << std::endl;
      }
     
    }
    if(myAxes){
      int x0 = static_cast<int>(-minx*scale);
      int y0 = static_cast<int>(-miny*scale);
      if(y0 >= 0 && y0<= xFigYMax){
	myFile << "2 1 0 1 0 7 50 0 -1 0.000 0 0 -1 1 0 2\n        0 0 1.00 60.00 120.00\n         0 "<<y0<<" "<<xFigXMax<<" "<<y0<<"\n";
      }
      if(x0 >= 0 && x0<= xFigXMax){
	myFile << "2 1 0 1 0 7 50 0 -1 0.000 0 0 -1 1 0 2\n        0 0 1.00 60.00 120.00\n         "<<x0<<" "<<xFigYMax<<" "<<x0<<" 0\n";
      }

    }


    close();
    return !myFile.fail();
  }

  Real XFigWriter::setRadius(Real radius){
    Real tmp = myRadius;
    myRadius = radius;
    return tmp;
  }

  Real XFigWriter::setFromZ(Real fromZ){
    Real tmp = myFromZ;
    myFromZ = fromZ;
    return tmp;
  }

  Real XFigWriter::setToZ(Real toZ){
    Real tmp = myToZ;
    myToZ = toZ;
    return tmp;
  }

  bool XFigWriter::setColor(bool color){
    bool tmp = myColor;
    myColor = color;
    return tmp;
  }

  bool XFigWriter::setAxes(bool axes){
    bool tmp = myAxes;
    myAxes = axes;
    return tmp;
  }

  bool XFigWriter::setLegend(bool legend){
    bool tmp = myLegend;
    myLegend = legend;
    return tmp;
  }

  bool XFigWriter::setAutofit(bool autofit){
    bool tmp = myAutofit;
    myAutofit = autofit;
    return tmp;
  }

  Matrix3by3 XFigWriter::setTransformation(const Matrix3by3& mat){
    Matrix3by3 tmp = myMat;
    myMat = mat;
    return tmp;
  }


  XFigWriter& operator<<(XFigWriter& xyzWriter, const XYZ& xyz){
    xyzWriter.write(xyz.coords,xyz.names);
    return xyzWriter;
  }

}
