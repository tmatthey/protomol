#include "iSGLennardJonesTableForce.h"
using std::string;
using namespace ProtoMol::Report;
namespace ProtoMol {
  //_________________________________________________________________ iSGLennardJonesTableForce

  const string iSGLennardJonesTableForce::keyword("iSGLennardJonesTable");

  //_________________________________________________________________ empty constructor
  iSGLennardJonesTableForce::iSGLennardJonesTableForce():myTable(NULL),myRc(-1.0),myN(0),alpha(0.0) {}

  //_________________________________________________________________ constructor
  iSGLennardJonesTableForce::iSGLennardJonesTableForce(Real rc, unsigned int n, Real myAlpha):myTable(new Real[4*n]),
											  myRc(rc),
											  myFac(((Real)n-2.0)/rc),
											  myN(n),
                                                                                          alpha(myAlpha) {
    for(unsigned int i=0;i<myN;++i){
      Real r = (i+(i==0?0.01:0.0))/myFac;      

      Real r2   = 1.0/(r*r);
      Real r6   = r2*r2*r2;
      Real r12  = r6*r6;
      //Real r6B  = r6;
      //Real r12A = r12;
      //energy = r12A - r6B;
      //force  = 12.0*r12A*r2 - 6.0*r6B*r2;

      myTable[i*4  ] = r12;
      myTable[i*4+1] = -r6;
      myTable[i*4+2] = 12.0*r12*r2;
      myTable[i*4+3] = -6.0*r6*r2;
    }
    report <<hint<<"Using iSGLennardJones approximation with "<<n<<" points, range ["<<(1.0/myFac)<<","<<((Real)n/myFac)<<"[, table[0,1,"<<myN-1<<"] = ("<<myTable[0]<<","<<myTable[1]<<"),("<<myTable[2]<<","<<myTable[3]<<"),("<<myTable[(int)myN*4-2]<<","<<myTable[(int)myN*4-1]<<")."<<endr;
  }

  //_________________________________________________________________ destructor
  iSGLennardJonesTableForce::~iSGLennardJonesTableForce(){
    if(myTable != NULL)
      delete [] myTable;
    myTable = NULL;
  }

  //_________________________________________________________________ copy constructor
  iSGLennardJonesTableForce::iSGLennardJonesTableForce(iSGLennardJonesTableForce const& other):myTable(new Real[4*other.myN]),
												  myRc(other.myRc),
												  myFac(((Real)other.myN-2.0)/other.myRc),
												  myN(other.myN),
                                                                                                  alpha(other.alpha) {
    for(unsigned int i=0;i<myN*4;++i)
      myTable[i] = other.myTable[i];    
  }

  //_________________________________________________________________ copy operator
  iSGLennardJonesTableForce& iSGLennardJonesTableForce::operator=(iSGLennardJonesTableForce const& other){
    if(&other != this){
      if(myTable != NULL)
	delete [] myTable;
      myRc           = other.myRc;
      myFac          = other.myFac;
      myN            = other.myN;
      alpha          = other.alpha;
      myTable        = new Real[4*myN];
      for(unsigned int i=0;i<myN*4;++i)
	myTable[i] = other.myTable[i];    
    }
    return *this;
  }

  //_________________________________________________________________ get parameters())
  void iSGLennardJonesTableForce::getParameters(std::vector<Parameter>& parameters) const{
    parameters.push_back(Parameter("-tableCutoff",Value(myRc,ConstraintValueType::Positive()),Text("cutoff for table look up")));
    parameters.push_back(Parameter("-tableSize",Value(myN,ConstraintValueType::Positive()),static_cast<unsigned int>(DEFAULT_TABLE_SIZE),Text("size of table look up")));
    parameters.push_back(Parameter("-alpha", Value(alpha)));
  }

  //_________________________________________________________________ make())
  iSGLennardJonesTableForce iSGLennardJonesTableForce::make(std::string& , const std::vector<Value>& values) {
    return iSGLennardJonesTableForce(values[0],values[1],values[2]);
  }

}
