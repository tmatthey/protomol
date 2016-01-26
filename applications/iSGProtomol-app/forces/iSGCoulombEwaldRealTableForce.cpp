#include "iSGCoulombEwaldRealTableForce.h"
using std::string;
using namespace ProtoMol::Report;
namespace ProtoMol {
  //_________________________________________________________________ iSGCoulombEwaldRealTableForce
  const string iSGCoulombEwaldRealTableForce::keyword("iSGCoulombEwaldRealTable");

  //_________________________________________________________________ empty constructor
  iSGCoulombEwaldRealTableForce::iSGCoulombEwaldRealTableForce():myTable(NULL),myERFTable(NULL),myAlpha(-1.0),myRc(-1.0),myN(0){}

  //_________________________________________________________________ constructor
  iSGCoulombEwaldRealTableForce::iSGCoulombEwaldRealTableForce(Real a, Real rc, unsigned int n):myTable(new Real[2*n]), myERFTable(new Real[n]),
											  myAlpha(a),
											  myRc(rc),
											  myAlphaSquared(a*a),
											  my2AlphaPI(2.0*a/sqrt(M_PI)),
											  myFac(((Real)n-2.0)/sqrt(rc)),
											  myN(n){
    // precompute values of the complementary error function for various r values
    for(unsigned int i=0;i<myN;++i){
      Real r = power<2>((i+(i==0?0.01:0.0))/myFac);
      
      // complementary error function      
      Real b = erfc(myAlpha*r)/r;
      myTable[i*2  ] = b;
      myTable[i*2+1] = (b+my2AlphaPI*exp(-myAlphaSquared*r*r))/(r*r); //Coulomb force contribution
      
      // error function, needed for the intramolecular correction terms
      Real c = erf(myAlpha*r)/r;
      myERFTable[i] = c;
    }
    report <<hint<<"Using iSGEwald real approximation with "<<n<<" points, range ["<<power<2>(1.0/myFac)<<","<<power<2>((Real)n/myFac)<<"[, table[0,1,"<<myN-1<<"] = ("<<myTable[0]<<","<<myTable[1]<<"),("<<myTable[2]<<","<<myTable[3]<<"),("<<myTable[(int)myN*2-2]<<","<<myTable[(int)myN*2-1]<<")."<<endr;
  }

  //_________________________________________________________________ destructor
  iSGCoulombEwaldRealTableForce::~iSGCoulombEwaldRealTableForce(){
    if(myTable != NULL)
      delete [] myTable;
    if(myERFTable != NULL)
      delete [] myERFTable;
    myTable = NULL;
    myERFTable = NULL;
  }

  //_________________________________________________________________ copy constructor
  iSGCoulombEwaldRealTableForce::iSGCoulombEwaldRealTableForce(iSGCoulombEwaldRealTableForce const& other):myTable(new Real[2*other.myN]), myERFTable(new Real[other.myN]),
												  myAlpha(other.myAlpha),
												  myRc(other.myRc),
												  myAlphaSquared(other.myAlpha*other.myAlpha),
												  my2AlphaPI(2.0*other.myAlpha/sqrt(M_PI)),
												  myFac(((Real)other.myN-2.0)/sqrt(other.myRc)),
												  myN(other.myN){
    for(unsigned int i=0;i<myN*2;++i) myTable[i] = other.myTable[i];
    for(unsigned int i=0;i<myN;++i) myERFTable[i] = other.myERFTable[i];
  }

  //_________________________________________________________________ copy operator
  iSGCoulombEwaldRealTableForce& iSGCoulombEwaldRealTableForce::operator=(iSGCoulombEwaldRealTableForce const& other){
    if(&other != this){
      if(myTable != NULL)
	delete [] myTable;
      if(myERFTable != NULL) 
        delete [] myERFTable;
      myAlpha        = other.myAlpha;
      myRc           = other.myRc;
      myAlphaSquared = other.myAlphaSquared;
      my2AlphaPI     = other.my2AlphaPI;
      myFac          = other.myFac;
      myN            = other.myN;
      myTable        = new Real[2*myN];
      myERFTable     = new Real[myN];
      for(unsigned int i=0;i<myN*2;++i) myTable[i] = other.myTable[i];
      for(unsigned int i=0;i<myN;++i) myERFTable[i] = other.myERFTable[i];
 
    }
    return *this;
  }

  //_________________________________________________________________ getParameters
  void iSGCoulombEwaldRealTableForce::getParameters(std::vector<Parameter>& parameters) const{
    parameters.push_back(Parameter("-alpha",Value(myAlpha,ConstraintValueType::Positive()),Text("splitting")));
    parameters.push_back(Parameter("-tableCutoff",Value(myRc,ConstraintValueType::Positive()),Text("cutoff for table look up")));
    parameters.push_back(Parameter("-tableSize",Value(myN,ConstraintValueType::Positive()),static_cast<unsigned int>(DEFAULT_TABLE_SIZE),Text("size of table look up")));
  }

  //_________________________________________________________________ make
  iSGCoulombEwaldRealTableForce iSGCoulombEwaldRealTableForce::make(std::string& , const std::vector<Value>& values) {
    return iSGCoulombEwaldRealTableForce(values[0],values[1],values[2]);
  }

}
