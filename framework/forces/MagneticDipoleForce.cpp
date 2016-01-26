#include "MagneticDipoleForce.h"
#include "mathutilities.h"
using std::string;
using std::vector;

namespace ProtoMol {
  //_________________________________________________________________ MagneticDipoleForce

  const string MagneticDipoleForce::keyword("MagneticDipole");
  MagneticDipoleForce::MagneticDipoleForce():myChi(0.0),myR(0.0),myOmega(0.0),myPhi(0.0),myHx(0.0),myHy(0.0),myHz(0.0),myD(0.0){}

  MagneticDipoleForce::MagneticDipoleForce(Real chi, Real radius, Real omega, Real phi, Real Hx,Real Hy, Real Hz, Real D):
    myChi(chi),myR(radius),myOmega(omega),myPhi(phi),myHx(Hx),myHy(Hy),myHz(Hz),myD(D){
    // Shortcuts
    volum     = 4.0/3.0*Constant::M_PI*power<3>(myR);
    expfactor = -1000.0*myChi*myChi*myR*myHx*myHx;
    realChi   = myChi/(1-2.0/3.0*myChi);
    kappa     = realChi/(realChi+2.0);
  }

  void MagneticDipoleForce::operator()(Real &energy, 
				       Real &force, 
				       Real /*distSquared*/,
				       Real rDistSquared,
				       Vector3D &diff, 
				       const GenericTopology* topo, 
				       int /*atom1*/, int /*atom2*/,
				       ExclusionClass) const {

    //defining variables used to calculate the force
    Real t = topo->time;
    Real inv_r = sqrt(rDistSquared);
    Real r = 1.0/inv_r;
    Real inv_r2 = rDistSquared;
    Real inv_r3 = inv_r2*inv_r; 
    Real inv_r5 = inv_r3*inv_r2; 
    Real inv_r7 = inv_r5*inv_r2;
    Vector3D H(myHx*cos(t*myOmega+myPhi),myHy*sin(t*myOmega+myPhi),myHz/(1.0+realChi));
    Vector3D sigma(-H*volum*myChi);
    Real sigmaDotSigma = sigma.dot(sigma);
    Real sigmaDotDiff = sigma.dot(diff);
    //Dipole-dipole effects:
    Vector3D F(diff*sigmaDotSigma*3.0*inv_r5
	       - diff*power<2>(sigmaDotDiff)*15.0*inv_r7
	       + sigma*6.0*sigmaDotDiff*inv_r5);
    energy = sigmaDotSigma*inv_r3-3*sigmaDotDiff*sigmaDotDiff*inv_r5;
      
    //Mirror-effect 1st approx (myD=0 turns mirroreffect off)
    if (myD != 0.0){
      // creating mirror-dipoles
      Vector3D sigma_m(sigma.x,sigma.y,-sigma.z); 
      sigma_m *= kappa;
      Vector3D diff_m_lower(diff.x,diff.y, myD);
      Vector3D diff_m_upper(diff.x,diff.y,-myD);

      // mirror-self-coupling
      // Not possible without z-coordinates
      // Done in MagneticDipoleMirrorSystemForce*.[Ch]
	  
      // mirror-dipole-coupling (assumes that the dipoles are all sentered at z=0)
      Real inv_r_m2 = 1.0/diff_m_upper.normSquared();
      Real inv_r_m3 = inv_r_m2*sqrt(inv_r_m2);
      Real inv_r_m5 = inv_r_m2*inv_r_m3;
      Real inv_r_m7 = inv_r_m5*inv_r_m2;
	  
      // now the dipole - mirror-force is added to F
      Real sigma_mDotsigma = sigma_m.dot(sigma);
      Real sigma_mDotdiff_m = sigma_m.dot(diff_m_upper);
      Real sigmaDotdiff_m = sigma.dot(diff_m_upper);
      F += diff_m_upper*(sigma_mDotsigma*3.0*inv_r_m5
			 - sigma_mDotdiff_m*sigmaDotdiff_m*15.0*inv_r_m7)
	+ sigma_m*sigmaDotdiff_m*6.0*inv_r_m5;
      energy += sigma_mDotsigma*inv_r_m3 - 3.0*sigma_mDotdiff_m*sigmaDotdiff_m*inv_r_m5;

      sigma_mDotdiff_m = sigma_m.dot(diff_m_lower);
      sigmaDotdiff_m = sigma.dot(diff_m_lower);
      F += diff_m_lower*(sigma_mDotsigma*3.0*inv_r_m5
			 - sigma_mDotdiff_m*sigmaDotdiff_m*15.0*inv_r_m7)
	+ sigma_m*sigmaDotdiff_m*6.0*inv_r_m5;
      energy += sigma_mDotsigma*inv_r_m3 - 3.0*sigma_mDotdiff_m*sigmaDotdiff_m*inv_r_m5;

    }
    //Peak-potetial at the edge of the spheres:
    Real f = exp(expfactor*(r-2*myR))*inv_r;
    F += diff*f;

    force = F.norm();
    diff = F/force;
  }


  void MagneticDipoleForce::getParameters(vector<Parameter>& parameters) const{
    parameters.push_back(Parameter("-chi",Value(myChi)));
    parameters.push_back(Parameter("-r",Value(myR)));
    parameters.push_back(Parameter("-omega",Value(myOmega)));
    parameters.push_back(Parameter("-phi",Value(myPhi)));
    parameters.push_back(Parameter("-H",Value(myHx)));
    parameters.push_back(Parameter("",Value(myHy)));
    parameters.push_back(Parameter("",Value(myHz)));
    parameters.push_back(Parameter("-d",Value(myD)));
  }

  MagneticDipoleForce MagneticDipoleForce::make(string& , const vector<Value>& values) {
    Real chi;    
    Real radius; 
    Real omega; 
    Real phi;    
    Real Hx;    
    Real Hy;    
    Real Hz;    
    Real D;     
    values[0].get(chi);    
    values[1].get(radius);
    values[2].get(omega);
    values[3].get(phi);
    values[4].get(Hx);
    values[5].get(Hy);
    values[6].get(Hz);
    values[7].get(D);
    return MagneticDipoleForce(chi,radius,omega,phi,Hx,Hy,Hz,D);
  }

}
