#include "Lagrange.h"

#include "mathutilities.h"
#include "Report.h"
using namespace ProtoMol::Report;
using std::string;

//#define DEBUG_LAGRANGE

namespace ProtoMol {
  //_________________________________________________________________ Lagrange

  const string Lagrange::keyword("Lagrange");

  Lagrange::Lagrange():myInterOrder(0),
		       theta(NULL),
		       dTheta(NULL){
  }

  Lagrange::Lagrange(unsigned int order):myInterOrder(order),
					 theta(new Real[order]),
					 dTheta(new Real[order]){
  }

  Lagrange::Lagrange(unsigned int order, Real w):myInterOrder(order),
						 theta(new Real[order]),
						 dTheta(new Real[order]){
    set(w);
  }

  Lagrange::~Lagrange(){
    if(theta != NULL){ 
      delete [] theta;
      delete [] dTheta;
    }
  }

  Lagrange::Lagrange(const Lagrange& Lagrange){
    myInterOrder = Lagrange.myInterOrder;
    theta  = new Real[myInterOrder];
    dTheta = new Real[myInterOrder];
    for(unsigned int k=0;k<myInterOrder;k++){
      theta[k] = Lagrange.theta[k];
      dTheta[k] = Lagrange.dTheta[k];      
    }
  }

  void Lagrange::setOrder(unsigned int order){
    if(order == myInterOrder && theta != NULL) 
      return;
    delete [] theta;
    delete [] dTheta;
    theta = new Real[order];
    dTheta = new Real[order];
    myInterOrder = order;
  }

  void Lagrange::set(Real w){
    switch(myInterOrder){
    case 0:
      report << error << "[Lagrange::set] Interpolation order is zero!"<<endr;
      break;
    case 1:
      theta[0]  = 1;
      dTheta[0] = 0;
      break;
    case 2:
      theta[1]  = w;
      theta[0]  = 1.0 - w;
      dTheta[1] = 1;
      dTheta[0] = -1;
      break;
    default:
      report << error << "[Lagrange::set] Order "<<myInterOrder<<" not supported."<<endr;
      break;
    }
#if defined(DEBUG_LAGRANGE)
    report <<plain << "Lagrange: order="<<myInterOrder<<", w="<<w<<endr;
    report << plain << "theta: ";
    for(unsigned int i=0;i<myInterOrder;i++)
      report << theta[i]<<" ";
    report << endr;
    report << plain << "dTheta: ";
    for(unsigned int i=0;i<myInterOrder;i++)
      report << dTheta[i]<<" ";
    report << endr;
#endif
  }

}  
