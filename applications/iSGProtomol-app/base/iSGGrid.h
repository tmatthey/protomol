/*  -*- c++ -*-  */
#ifndef ISGGRID_H
#define ISGGRID_H

#include "Array.h"
#include "Vector3D.h"
#include "mathutilities.h"
#include "FFTComplex.h"
#include "Report.h"
#include "ScalarStructure.h"

namespace ProtoMol {
  //_________________________________________________________________ iSGGrid
  // A simple Grid class using T as interpolation scheme, assuming periodic 
  // boundary conditions.  This modification is for iSGMD simulations -- TIM


  template<class TInterpolation>
  class iSGGrid {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Typedef & const
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    struct Int3D {int x; int y; int z;};
    struct Interpolation3D {
      TInterpolation x;
      TInterpolation y;
      TInterpolation z;
      Interpolation3D(){};
      Interpolation3D(unsigned int order):x(order),y(order),z(order){}
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    iSGGrid();
    ~iSGGrid();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Grid
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    // interpolate the charge and deltaQ onto separate FFT grids
    void anterpolateCharge(Real q, Real Delta_q, const Vector3D& coord, unsigned int index);

    // forward and backward FFTs
    void fftBack(){
      myFFT.backward();
      myDMU_FFT.backward();
    }
    void fftForward(){myFFT.forward();}

    // calculation of energy and chemical potential difference
    void scalarSum(ScalarStructure* energies, Real& energy, Real& deltaMu);
    void scalarSum(ScalarStructure* energies, Real& energy, Real& deltaMu, unsigned int block, unsigned int nBlocks);

    // calculation of the force and pressure
    void interpolateForce(Real q, unsigned int index, Vector3D& force);

    // function that initializes the grids and precomputes certain terms
    void initialize(Real width, Real length, Real height, Real alpha,
		    unsigned int nx, unsigned int ny, unsigned int nz, 
		    unsigned int interOrder,
		    unsigned int atomCount);
    // The order of arguments seems not nice, but g++ 2.95.2 on SUN 5.7 does
    // segfault otherwise ... do not ask why!

    // function that clears the charge distribution grids
    void clear();

    // for parallel operations
    void getQ(Real*& begin, Real*& end) {begin=&(myQ.begin()->re);end=&(myQ.end()->re);}
    void print();
    
  private:
    // computes the array B(m1,m2,m3)
    void dftmod(unsigned int order, unsigned int n, Real* interpolation, Real* interpolationMod);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Array<zomplex,3> myQ;           // the atomic charge distribution on a grid
    Array<zomplex,3> myQTmp;        // the atomic charge distribution on a grid
    Array<zomplex,3> myDMU_Q;       // the atomic charge difference (Qnew - Qold) distribution on a grid
    Array<zomplex,3> myDMU_QTmp;    // the atomic charge difference (Qnew - Qold) distribution on a grid

    unsigned int myNX;              // the number of grid points in the x-direction
    unsigned int myNY;              // the number of grid points in the y-direction
    unsigned int myNZ;              // the number of grid points in the z-direction
    
    int myNXOffset;
    int myNYOffset;
    int myNZOffset;
    
    Real myWidth;                   // cell basis vector
    Real myLength;                  // cell basis vector 
    Real myHeight;                  // cell basis vector
    Real myWidthr;                  // reciprocal cell basis vector
    Real myLengthr;                 // reciprocal cell basis vector
    Real myHeightr;                 // reciprocal cell basis vector

    Real myV;                       // cell volume
    Real myHX;                      // grid point spacing in x-direction
    Real myHY;                      // grid point spacing in y-direction
    Real myHZ;                      // grid point spacing in z-direction
    Real myHXr;                     // reciprocal grid point spacing in x-direction
    Real myHYr;                     // reciprocal grid point spacing in y-direction
    Real myHZr;                     // reciprocal grid point spacing in z-direction

    Real myAlpha;
    std::vector<Int3D> myScaledParticleIntPositions;
    std::vector<Interpolation3D> myInterpolations;

    Real* myInerpolationModX;       // bx(mx)
    Real* myInerpolationModY;       // by(my)
    Real* myInerpolationModZ;       // bz(mz)

    Real* myExpX;                   // Cx(mx)
    Real* myExpY;                   // Cy(my)
    Real* myExpZ;                   // Cz(mz)

    unsigned int myInterOrder;
    FFTComplex myFFT;               // object used to perform the FFT operations on the charge (Q) grid
    FFTComplex myDMU_FFT;           // object used to perform the FFT operations on the deltaQ (DMU_Q) grid
    unsigned int myAtomCount;       // # of atoms in the system
    Real myFac;                     // (pi^2 / alpha^2)
  };
  
  //______________________________________________________________________ INLINES

  //______________________________________________________________________ anterpolateCharge
  template<class TInterpolation>
  inline void iSGGrid<TInterpolation>::anterpolateCharge(Real q, Real Delta_q, const Vector3D& coord, unsigned int index)  {

    // scale the coordinates
    Real x = coord.x*myHXr;
    Real y = coord.y*myHYr;
    Real z = coord.z*myHZr;

    // compute the scaled fractional coordinates
    while(x < 0.0) x += myNX;
    while(x >= myNX) x -= myNX;
    while(y < 0.0) y += myNY;
    while(y >= myNY) y -= myNY;
    while(z < 0.0) z += myNZ;
    while(z >= myNZ) z -= myNZ;

    // get the integer part of the scaled fractional coordinates
    int intX = (int)x;
    int intY = (int)y;
    int intZ = (int)z;
    int i0 = intX+myNXOffset;
    int j0 = intY+myNYOffset;
    int k0 = intZ+myNZOffset;
    myScaledParticleIntPositions[index].x = i0;
    myScaledParticleIntPositions[index].y = j0;
    myScaledParticleIntPositions[index].z = k0;

    // decimal (or non-integer) part of the scaled fractional coordinates
    // this is the quantity (u - [u]) of eq. (3.2) of Essmann et. al.,
    // "A Smooth Particle Mesh Ewald Method", J. Chem. Phys. 103(19), 1995, p.8577-8593.
    myInterpolations[index].x.set(x-intX);
    myInterpolations[index].y.set(y-intY);
    myInterpolations[index].z.set(z-intZ);

    // theta[1] = (u - [u])
    // theta[0] = 1 - (u - [u])
    Real*  thetaX = myInterpolations[index].x.theta;
    Real*  thetaY = myInterpolations[index].y.theta;
    Real*  thetaZ = myInterpolations[index].z.theta;

    // loop over the interpolation order in x
    for(unsigned int i=0;i<myInterOrder;i++){

      // add the charge and deltaQ to the spline
      Real a = q*thetaX[i];
      Real b = Delta_q*thetaX[i];
      
      int i1 = (i+i0) % myNX;

      // create references to the Q and DeltaQ grids
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
      Array<zomplex,3>::RefArray<2> rQX = myQ[i1];
      Array<zomplex,3>::RefArray<2> rDMU_QX = myDMU_Q[i1];
#else
      RefArray<zomplex,2> rQX = myQ[i1];
      RefArray<zomplex,2> rDMU_QX = myDMU_Q[i1];
#endif
      
      // loop over the interpolation order in y
      for(unsigned int j=0;j<myInterOrder;j++){
	
	// add the charge and deltaQ to the spline
	Real ab = a*thetaY[j];
	Real bb = b*thetaY[j];
	
	int j1 = (j+j0) % myNY;
	
	// create references to the Q and DeltaQ grids
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	Array<zomplex,3>::RefArray<1> rQXY = rQX[j1];
	Array<zomplex,3>::RefArray<1> rDMU_QXY = rDMU_QX[j1];
#else
	RefArray<zomplex,1> rQXY = rQX[j1];
	RefArray<zomplex,1> rDMU_QXY = rDMU_QX[j1];
#endif
	
	// loop over the interpolation order in z
	for(unsigned int k=0;k<myInterOrder;k++){
	  int k1 = (k+k0) % myNZ;
	  
	  // add the charge and deltaQ to the spline
	  rQXY[k1].re += ab*thetaZ[k];
	  rDMU_QXY[k1].re += bb*thetaZ[k];

	}  // end loop over k
      }  // end loop over j
    }  // end loop over i
  }  // end function


  //______________________________________________________________________ interpolateForce
  template<class TInterpolation>
  inline void iSGGrid<TInterpolation>::interpolateForce(Real q, unsigned int index, Vector3D& force){
    int i0 = myScaledParticleIntPositions[index].x;
    int j0 = myScaledParticleIntPositions[index].y;
    int k0 = myScaledParticleIntPositions[index].z;
    Real fx = 0.0;
    Real fy = 0.0;
    Real fz = 0.0;
    Real*  thetaX = myInterpolations[index].x.theta;
    Real*  thetaY = myInterpolations[index].y.theta;
    Real*  thetaZ = myInterpolations[index].z.theta;
    Real*  dThetaX = myInterpolations[index].x.dTheta;
    Real*  dThetaY = myInterpolations[index].y.dTheta;
    Real*  dThetaZ = myInterpolations[index].z.dTheta;
    for(unsigned int i=0;i<myInterOrder;i++){
      int i1 = (i+i0)%myNX; 
      for(unsigned int j=0;j<myInterOrder;j++){
	int j1 = (j+j0)%myNY;
	Real xij = dThetaX[i]*thetaY[j];
	Real yij = thetaX[i]*dThetaY[j];
	Real zij = thetaX[i]*thetaY[j];
	for(unsigned int k=0;k<myInterOrder;k++){
	  int k1 = (k+k0)%myNZ;
	  Real term = myQ[i1][j1][k1].re;
	  fx -= term*xij*thetaZ[k];
	  fy -= term*yij*thetaZ[k];
	  fz -= term*zij*dThetaZ[k];
	}
      }
    }
    force += Vector3D(fx*myHXr*q,fy*myHYr*q,fz*myHZr*q);
  }


  //______________________________________________________________________ empty constructor
  template<class TInterpolation>
  iSGGrid<TInterpolation>::iSGGrid():
    myNX(0), myNY(0), myNZ(0),
    myNXOffset(0), myNYOffset(0), myNZOffset(0),
    myWidth(0.0), myLength(0.0), myHeight(0.0),
    myWidthr(0.0), myLengthr(0.0), myHeightr(0.0),
    myV(0.0),
    myHX(0.0), myHY(0.0), myHZ(0.0),
    myHXr(0.0), myHYr(0.0), myHZr(0.0),
    myAlpha(0.0),
    myInerpolationModX(NULL), myInerpolationModY(NULL), myInerpolationModZ(NULL),
    myExpX(NULL), myExpY(NULL), myExpZ(NULL),
    myInterOrder(0),
    myAtomCount(0),
    myFac(0.0){
  }
  

  //______________________________________________________________________ destructor
  template<class TInterpolation>
  iSGGrid<TInterpolation>::~iSGGrid(){
    if(myInerpolationModX != NULL) delete [] myInerpolationModX;
    if(myInerpolationModY != NULL) delete [] myInerpolationModY;
    if(myInerpolationModZ != NULL) delete [] myInerpolationModZ;
    if(myExpX != NULL) delete [] myExpX;
    if(myExpY != NULL) delete [] myExpY;
    if(myExpZ != NULL) delete [] myExpZ;
  }


  //______________________________________________________________________ scalarSum function
  template<class TInterpolation>
  void iSGGrid<TInterpolation>::scalarSum(ScalarStructure* energies, Real& energy, Real& deltaMu){
    
    // initialize the energy, deltaMu, and pressure
    energy = deltaMu = 0.0;
    Real virialxx = 0.0;
    Real virialxy = 0.0;
    Real virialxz = 0.0;
    Real virialyy = 0.0;
    Real virialyz = 0.0;
    Real virialzz = 0.0;

    // find out if we need to calculate the virial
    bool doVirial = energies->virial();
    bool doMolVirial = energies->molecularVirial();

    // 1 / (pi * Volume)
    Real piVr = 1.0/(M_PI*myV);
			
    int count = 0;

    // loop over all grid points in the x-direction
    for (unsigned int i = 0; i < myNX; i++){
      
      // determine the grid x-index
      int i0 = i <= myNX/2 ? i : i-myNX;
      
      // compute the x-component of the m vector
      Real mi = i0*myWidthr;
      
      // compute 1/(pi*V) * exp(-pi^2 * mx^2 / alpha^2)
      Real ex = myExpX[i]*piVr;
      
      // create references to the inverse FFT'd charge and deltaQ grids
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
      Array<zomplex,3>::RefArray<2> rQX = myQ[i];
      Array<zomplex,3>::RefArray<2> rDMU_QX = myDMU_Q[i];
#else
      RefArray<zomplex,2> rQX = myQ[i];
      RefArray<zomplex,2> rDMU_QX = myDMU_Q[i];
#endif
      
      // loop over all grid points in the y-direction
      for (unsigned int j = 0 ; j < myNY; j++){
	
	// determine the grid y-index
	int j0 = j <= myNY/2 ? j : j-myNY;
	
	// compute |bx(mx)|^2 * |by(my)|^2
	Real interpolationModXY = myInerpolationModX[i]*myInerpolationModY[j];

	// compute the y-component of the m vector
	Real mj = j0*myLengthr;
	
	// compute mx^2 + my^2	
	Real mij = mi*mi + mj*mj;
	
	// compute 1/(pi*V) * exp[-pi^2 * (mx^2 + my^2) / alpha^2]
	Real exy = ex*myExpY[j];
				
	// create references to the inverse FFT'd charge and deltaQ grids
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	Array<zomplex,3>::RefArray<1> rQXY = rQX[j];
	Array<zomplex,3>::RefArray<1> rDMU_QXY = rDMU_QX[j];
#else
	RefArray<zomplex,1> rQXY = rQX[j];
	RefArray<zomplex,1> rDMU_QXY = rDMU_QX[j];
#endif
	
	// loop over all grid points in the z-direction
	for (unsigned int k = (i!=0 || j!=0 ? 0:1); k <= myNZ/2; k++){
	  count++;
	  
	  // determine the grid z-index
	  int k0 = k <= myNZ/2 ? k : k-myNZ;
					
	  // compute the z-component of the m vector
	  Real mk = k0*myHeightr;
	  
	  // compute (mx^2 + my^2 + mz^2) = m^2
	  Real mHatSquared = mij+mk*mk;
	  
	  // compute the inverse fourier transform of theta
	  // this is equal to the product of arrays B and C, B*C(mx,my,mz) of
	  // Essmann et. al., "A Smooth Particle Mesh Ewald Method",
	  // J. Chem. Phys. 103(19), 1995, p.8577-8593. {see eq. (4.7)}
	  Real theta = interpolationModXY*myInerpolationModZ[k]*exy*myExpZ[k]/mHatSquared;
	  
	  // Energy and chemical potential difference
	  // q = power<2>(myQ[i][j][k].re)+power<2>(myQ[i][j][k].im);
	  Real q = power<2>(rQXY[k].re)+power<2>(rQXY[k].im);
	  Real q_Dq = (rQXY[k].re * rDMU_QXY[k].re) + (rQXY[k].im * rDMU_QXY[k].im);
	  
	  // the energy is equal to B*C * F(Q)[mx,my,mz] * F(Q)[-mx,-my-mz]
	  // which is q * theta.  See eq. (4.7) of Essmann et. al.
	  // the chemical potential difference is equal to
	  // B*C * {F(DeltaQ)[mx,my,mz] * F(Q)[-mx,-my-mz] + F(Q)[mx,my,mz] * F(DeltaQ)[-mx,-my-mz]}
	  Real e = q*theta;
	  Real dmu = 2.0*q_Dq*theta;
	  
	  // this is part of the Kroenecker delta term of the virial,
	  // see eq. (2.7) of Essmann et al.
	  Real v = 2.0*(1.0/mHatSquared + myFac);
	  
	  // Symmetric 
	  if(k > 0 && ((k != myNZ/2) || (myNZ & 1))){
	    e *= 2.0;
	    dmu *= 2.0;
	    
	    // overwrite the charge grid Q with the convolution (theta * Q)
	    zomplex& w = myQ[(myNX-i)%myNX][(myNY-j)%myNY][(myNZ-k)%myNZ];
	    w.re *= theta;	
	    w.im *= theta;    	     
	  }
	  
	  // sum over all the grid points to get the total energy and
	  // chemical potential difference
	  energy += e;
	  deltaMu += dmu;
	  
	  // Virial
	  if(doVirial){
	    virialxx += e *(1.0 - v * mi*mi);
	    virialxy -= e *(v * mi*mj);
	    virialxz -= e *(v * mi*mk);
	    virialyy += e *(1.0 - v * mj*mj);
	    virialyz -= e *(v * mj*mk);
	    virialzz += e *(1.0 - v * mk*mk);
	  }
	  
	  // convolute theta and Q
	  rQXY[k].re *= theta;
	  rQXY[k].im *= theta;
	  //myQ[i][j][k].re *= theta;
	  //myQ[i][j][k].im *= theta;
	}  // end loop over z grid points
      }  // end loop over y grid points
    }  // end loop over x grid points
  
    // Just clear (0,0,0) since we did this not the nested loop.
    myQ[0][0][0].re = 0.0;
    myQ[0][0][0].im = 0.0;

    // add to the atomic and molecular virials
    if(doVirial){
      (*energies)[ScalarStructure::VIRIALXX] += virialxx*0.5;
      (*energies)[ScalarStructure::VIRIALXY] += virialxy*0.5;
      (*energies)[ScalarStructure::VIRIALXZ] += virialxz*0.5;
      (*energies)[ScalarStructure::VIRIALYX] += virialxy*0.5;
      (*energies)[ScalarStructure::VIRIALYY] += virialyy*0.5;
      (*energies)[ScalarStructure::VIRIALYZ] += virialyz*0.5;
      (*energies)[ScalarStructure::VIRIALZX] += virialxz*0.5;
      (*energies)[ScalarStructure::VIRIALZY] += virialyz*0.5;
      (*energies)[ScalarStructure::VIRIALZZ] += virialzz*0.5;
    }
    if(doMolVirial) {
      (*energies)[ScalarStructure::MOLVIRIALXX] += virialxx*0.5;
      (*energies)[ScalarStructure::MOLVIRIALXY] += virialxy*0.5;
      (*energies)[ScalarStructure::MOLVIRIALXZ] += virialxz*0.5;
      (*energies)[ScalarStructure::MOLVIRIALYX] += virialxy*0.5;
      (*energies)[ScalarStructure::MOLVIRIALYY] += virialyy*0.5;
      (*energies)[ScalarStructure::MOLVIRIALYZ] += virialyz*0.5;
      (*energies)[ScalarStructure::MOLVIRIALZX] += virialxz*0.5;
      (*energies)[ScalarStructure::MOLVIRIALZY] += virialyz*0.5;
      (*energies)[ScalarStructure::MOLVIRIALZZ] += virialzz*0.5;
    }

    // finally, multiply in the factor of 1/2 to the
    // energy and chemical potential
    energy *= 0.5;
    deltaMu *= 0.5;
    
  }


  //______________________________________________________________________ scalarSum function for parallel simulations
  template<class TInterpolation>
  void iSGGrid<TInterpolation>::scalarSum(ScalarStructure* energies, Real& energy, Real& deltaMu, unsigned int block, unsigned int nBlocks){

    // initialize the energy, deltaMu, and pressure
    energy = deltaMu = 0.0;
    Real virialxx = 0.0;
    Real virialxy = 0.0;
    Real virialxz = 0.0;
    Real virialyy = 0.0;
    Real virialyz = 0.0;
    Real virialzz = 0.0;

    // find out if we need to calculate the virial
    bool doVirial = energies->virial();
    bool doMolVirial = energies->molecularVirial();

    // 1 / (pi*Volume)
    Real piVr = 1.0/(M_PI*myV);

    // allocate memory for the charge distribution arrays
    myQTmp.resize(ArraySizes(myNX)(myNY)(myNZ));
    myDMU_QTmp.resize(ArraySizes(myNX)(myNY)(myNZ));


    int m = myQ.size();
    zomplex *q = myQ.begin();
    zomplex *t = myQTmp.begin();
    for(int i=0;i<m;i++){
      t[i].re = q[i].re;
      t[i].im = q[i].im;      
      q[i].re = 0.0;
      q[i].im = 0.0;
    }

    int m_DMU = myDMU_Q.size();
    zomplex *q_DMU = myDMU_Q.begin();
    zomplex *t_DMU = myDMU_QTmp.begin();
    for(int i=0;i<m_DMU;i++){
      t_DMU[i].re = q_DMU[i].re;
      t_DMU[i].im = q_DMU[i].im;      
      q_DMU[i].re = 0.0;
      q_DMU[i].im = 0.0;
    }

    // # of grid points
    int nx = myNX;
    int ny = myNY;
    int nz = myNZ/2+1;

    int nyz = ny*nz;
    int n  = nx*nyz;
    int sn = (n*block)/nBlocks + (block==0?1:0); // Add 1 to skip i,j,k == 0
    int en = (n*(block+1))/nBlocks - 1;

    int count = 0;
    int size = en-sn+1;
    if(size == 0) return;

    int k = sn % nz;
    int j = (sn / nz) % ny;
    int i = (sn / nyz);
    int ez = (en % nz)+1;
    int ey = ((en / nz) % ny)+1;
    int ex = (en / nyz)+1;
    if(j < ey-1)
      ez = nz;
    if(i < ex-1){
      ey = ny;
      ez = nz;
    }


    // loop over all grid points in the x-direction 
    for (; i < ex; i++,j=0){

      // determine the grid x-index
      int i0 = i <= static_cast<int>(myNX/2) ? i : i-myNX;

      // compute the x-component of the m vector
      Real mi = i0*myWidthr;

      // compute 1/(pi*V) * exp(-pi^2 * mx^2 / alpha^2)
      Real ex = myExpX[i]*piVr;

      // create references to the inverse FFT'd charge and deltaQ grids
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
      Array<zomplex,3>::RefArray<2> rQX  = myQ[i];
      Array<zomplex,3>::RefArray<2> rTQX = myQTmp[i];
      //Array<zomplex,3>::RefArray<2> rDMU_QX  = myDMU_Q[i];
      Array<zomplex,3>::RefArray<2> rDMU_TQX = myDMU_QTmp[i];
#else
      RefArray<zomplex,2> rQX  = myQ[i];
      RefArray<zomplex,2> rTQX = myQTmp[i];
      //RefArray<zomplex,2> rDMU_QX  = myDMU_Q[i];
      RefArray<zomplex,2> rDMU_TQX = myDMU_QTmp[i];
#endif

      // loop over all grid points in the y-direction 
      for (; j < ey; j++,k=0){

	// determine the grid y-index
	int j0 = j <= static_cast<int>(myNY/2) ? j : j-myNY;

	// compute |bx(mx)|^2 * |by(my)|^2
	Real interpolationModXY = myInerpolationModX[i]*myInerpolationModY[j];

	// compute the y-component of the m vector
	Real mj = j0*myLengthr;

	// compute mx^2 + my^2
	Real mij = mi*mi + mj*mj;

	// compute 1/(pi*V) * exp[-pi^2 * (mx^2 + my^2) / alpha^2]
	Real exy = ex*myExpY[j];	

	// create references to the inverse FFT'd charge and deltaQ grids
#ifdef NO_PARTIAL_TEMPLATE_SPECIALIZATION
	Array<zomplex,3>::RefArray<1> rQXY  = rQX[j];
	Array<zomplex,3>::RefArray<1> rTQXY = rTQX[j];
	//Array<zomplex,3>::RefArray<1> rDMU_QXY  = rDMU_QX[j];
	Array<zomplex,3>::RefArray<1> rDMU_TQXY = rDMU_TQX[j];
#else
	RefArray<zomplex,1> rQXY  = rQX[j];
	RefArray<zomplex,1> rTQXY = rTQX[j];
	//RefArray<zomplex,1> rDMU_QXY  = rDMU_QX[j];
	RefArray<zomplex,1> rDMU_TQXY = rDMU_TQX[j];
#endif
	
	// loop over all grid points in the z-direction 
	for (; k < ez; k++){

	  // determine the grid z-index
	  int k0 = k <= static_cast<int>(myNZ/2) ? k : k-myNZ;

	  // compute the z-component of the m vector
	  Real mk = k0*myHeightr;

	  // compute (mx^2 + my^2 + mz^2) = m^2
	  Real mHatSquared = mij+mk*mk;

	  // compute the inverse fourier transform of theta
	  // this is equal to the product of arrays B and C, B*C(mx,my,mz) of
	  // Essmann et. al., "A Smooth Particle Mesh Ewald Method",
	  // J. Chem. Phys. 103(19), 1995, p.8577-8593. {see eq. (4.7)}
	  Real theta = interpolationModXY*myInerpolationModZ[k]*exy*myExpZ[k]/mHatSquared;
	  
	  // Energy and chemical potential difference
	  //Real q = power<2>(myQ[i][j][k].re)+power<2>(myQ[i][j][k].im);
	  Real q = power<2>(rTQXY[k].re)+power<2>(rTQXY[k].im);
	  Real q_Dq = (rDMU_TQXY[k].re * rTQXY[k].re) + (rDMU_TQXY[k].im * rTQXY[k].im);
	  
	  // the energy is equal to B*C * F(Q)[mx,my,mz] * F(Q)[-mx,-my-mz]
	  // which is q * theta.  See eq. (4.7) of Essmann et. al.
	  // the chemical potential difference is equal to
	  // B*C * {F(DeltaQ)[mx,my,mz] * F(Q)[-mx,-my-mz] + F(Q)[mx,my,mz] * F(DeltaQ)[-mx,-my-mz]}
	  Real e = q*theta;
	  Real dmu = 2.0*q_Dq*theta;

	  // this is part of the Kroenecker delta term of the virial
	  // see eq. (2.7) of Essmann et. al.
	  Real v = 2.0*(1.0/mHatSquared + myFac);
	  
	  // Symmetric 
	  if(k > 0 && ((k != static_cast<int>(myNZ/2)) || (myNZ & 1))){
	    e *= 2.0;
	    dmu *= 2.0;
	    myQ[(myNX-i)%myNX][(myNY-j)%myNY][(myNZ-k)%myNZ].re = myQTmp[(myNX-i)%myNX][(myNY-j)%myNY][(myNZ-k)%myNZ].re * theta;
	    myQ[(myNX-i)%myNX][(myNY-j)%myNY][(myNZ-k)%myNZ].im = myQTmp[(myNX-i)%myNX][(myNY-j)%myNY][(myNZ-k)%myNZ].im * theta;         
	  }

	  // sum over all grid points to get the total energy and chemical potential difference
	  energy += e;
	  deltaMu += dmu;

	  // Virial
	  if(doVirial){
	    virialxx += e *(1.0 - v * mi*mi);
	    virialxy -= e *(v * mi*mj);
	    virialxz -= e *(v * mi*mk);
	    virialyy += e *(1.0 - v * mj*mj);
	    virialyz -= e *(v * mj*mk);
	    virialzz += e *(1.0 - v * mk*mk);
	  }
	  
	  // convolute theta and Q
	  rQXY[k].re = rTQXY[k].re * theta;
	  rQXY[k].im = rTQXY[k].im * theta;
	  //myQ[i][j][k].re *= theta;
	  //myQ[i][j][k].im *= theta;
	  
	  count++;
	  if(count >= size){
	    i = myNX;
	    j = myNY;		
	    k = myNZ;	
	  }  // end loop over z grid points
	}  // end loop over y grid points
      }  // end loop over x grid points
    }
    
    // Just clear (0,0,0) since we did this not the nested loop.
    myQ[0][0][0].re = 0.0;
    myQ[0][0][0].im = 0.0;

    // add to the atomic and molecular virials
    if(doVirial){
      (*energies)[ScalarStructure::VIRIALXX] += virialxx*0.5;
      (*energies)[ScalarStructure::VIRIALXY] += virialxy*0.5;
      (*energies)[ScalarStructure::VIRIALXZ] += virialxz*0.5;
      (*energies)[ScalarStructure::VIRIALYX] += virialxy*0.5;
      (*energies)[ScalarStructure::VIRIALYY] += virialyy*0.5;
      (*energies)[ScalarStructure::VIRIALYZ] += virialyz*0.5;
      (*energies)[ScalarStructure::VIRIALZX] += virialxz*0.5;
      (*energies)[ScalarStructure::VIRIALZY] += virialyz*0.5;
      (*energies)[ScalarStructure::VIRIALZZ] += virialzz*0.5;
    }
    if(doMolVirial) {
      (*energies)[ScalarStructure::MOLVIRIALXX] += virialxx*0.5;
      (*energies)[ScalarStructure::MOLVIRIALXY] += virialxy*0.5;
      (*energies)[ScalarStructure::MOLVIRIALXZ] += virialxz*0.5;
      (*energies)[ScalarStructure::MOLVIRIALYX] += virialxy*0.5;
      (*energies)[ScalarStructure::MOLVIRIALYY] += virialyy*0.5;
      (*energies)[ScalarStructure::MOLVIRIALYZ] += virialyz*0.5;
      (*energies)[ScalarStructure::MOLVIRIALZX] += virialxz*0.5;
      (*energies)[ScalarStructure::MOLVIRIALZY] += virialyz*0.5;
      (*energies)[ScalarStructure::MOLVIRIALZZ] += virialzz*0.5;
    }

    // finally, multiply in the factor of 1/2 to the energy and chemical potential difference
    energy *= 0.5;
    deltaMu *= 0.5;

  }


  //______________________________________________________________________ initialize function  
  template<class TInterpolation>
  void iSGGrid<TInterpolation>::initialize(Real width, Real length, Real height, Real alpha,
					   unsigned int nx, unsigned int ny, unsigned int nz, 
					   unsigned int interOrder, unsigned int atomCount){

    // allocate memory for the charge distribution arrays
    if(!myQ.resize(ArraySizes(nx)(ny)(nz)))
      report << error <<"[Grid<>::initialize] Could not allocate memory for Q["
	     <<nx<<"]["<<ny<<"]["<<nz<<"]."<<endr;
    if(!myDMU_Q.resize(ArraySizes(nx)(ny)(nz)))
      report << error <<"[Grid<>::initialize] Could not allocate memory for deltaQ["
	     <<nx<<"]["<<ny<<"]["<<nz<<"]."<<endr;
    
    // # of grid points in each direction   
    myNX = nx; 
    myNY = ny;
    myNZ = nz;
    
    myNXOffset = -((int)interOrder-1)/2 + (int)nx;
    myNYOffset = -((int)interOrder-1)/2 + (int)ny;
    myNZOffset = -((int)interOrder-1)/2 + (int)nz;
    
    // the current cell basis vectors and
    // reciprocal basis vectors
    myWidth = width; 
    myLength = length; 
    myHeight = height;
    myWidthr = 1.0/width; 
    myLengthr = 1.0/length; 
    myHeightr = 1.0/height;
    
    // volume
    myV = width*length*height;
    
    // spacing and reciprocal spacing between grid points in each direction
    myHX = width/nx;
    myHY = length/ny;
    myHZ = height/nz;
    myHXr = nx/width;
    myHYr = ny/length;
    myHZr = nz/height;
    
    myAlpha = alpha;
    myInterOrder = interOrder;
    myAtomCount = atomCount;
    myFac =  M_PI*M_PI/(myAlpha*myAlpha);

    // Precompute exp(-pi^2m^2/alpha^2), this is the
    // array C(mx,my,mz) in Essmann, et. al.
    if(myExpX != NULL) 
      delete [] myExpX;
    myExpX = new Real[nx];
    if(myExpY != NULL)
      delete [] myExpY;
    myExpY = new Real[ny];
    if(myExpZ != NULL) 
      delete [] myExpZ;
    myExpZ = new Real[nz];

    for (unsigned int i = 0; i < myNX; i++){
      int i0 = i <= myNX/2 ? i : i-myNX;
      myExpX[i] = exp(-myFac*power<2>(i0*myWidthr));
    }
    for (unsigned int j = 0 ; j < myNY; j++){
      int j0 = j <= myNY/2 ? j : j-myNY;
      myExpY[j] = exp(-myFac*power<2>(j0*myLengthr));
    }

    //for (unsigned int k = 0; k < myNZ; k++){
    for (unsigned int k = 0; k<= myNZ/2 ; k++){
      int k0 = k <= myNZ/2 ? k : k-myNZ;
      myExpZ[k] = exp(-myFac*power<2>(k0*myHeightr));
    }


    // Precompute the mod TInterpolation, this is the
    // array B(mx,my,mz) n Essmann, et. al.
    if(myInerpolationModX != NULL) 
      delete [] myInerpolationModX;
    myInerpolationModX = new Real[nx];
    if(myInerpolationModY != NULL)
      delete [] myInerpolationModY;
    myInerpolationModY = new Real[ny];
    if(myInerpolationModZ != NULL) 
      delete [] myInerpolationModZ;
    myInerpolationModZ = new Real[nz];

    TInterpolation interpolation = TInterpolation(myInterOrder,0.0);
    dftmod(myInterOrder,nx,interpolation.theta,myInerpolationModX);
    dftmod(myInterOrder,ny,interpolation.theta,myInerpolationModY);
    dftmod(myInterOrder,nz,interpolation.theta,myInerpolationModZ);
    //for(unsigned int i=0;i<nx;i++)
    //  report << plain <<"myInerpolationModX["<<i<<"]:"<<1.0/myInerpolationModX[i]<<endr;
    //for(unsigned int i=0;i<ny;i++)
    //  report << plain <<"myInerpolationModY["<<i<<"]:"<<1.0/myInerpolationModY[i]<<endr;
    //for(unsigned int i=0;i<nz;i++)
    //  report << plain <<"myInerpolationModZ["<<i<<"]:"<<1.0/myInerpolationModZ[i]<<endr;



    // Resize the vector data members
    myScaledParticleIntPositions.resize(myAtomCount);
    myInterpolations.resize(myAtomCount,Interpolation3D(myInterOrder));

    // initialize the FFT routines for the charge and deltaQ grids
    myFFT.initialize(myNX,myNY,myNZ,&myQ[0][0][0]);
    myDMU_FFT.initialize(myNX,myNY,myNZ,&myDMU_Q[0][0][0]);
    
  }


  //______________________________________________________________________ dftmod function  
  template<class TInterpolation>
  void iSGGrid<TInterpolation>::dftmod(unsigned int order, unsigned int n, Real* interpolation, Real* interpolationMod){

    // compute the denominator of the bi(m)^2 coefficients,
    // see eq. (4.4) of Essmann, et. al.
    for(unsigned int i=0;i<n;i++){
      Real sumCos = 0.0;
      Real sumSin = 0.0;
      for(unsigned int j=0;j<order;j++){
	Real x = M_PI*2.0*i*j/(Real)n;
	sumCos += interpolation[j]*cos(x);
	sumSin += interpolation[j]*sin(x);
      }
      interpolationMod[i] = 1.0/(sumCos*sumCos + sumSin*sumSin);
    }
  }


  //______________________________________________________________________ clear function
  template<class TInterpolation>
  void iSGGrid<TInterpolation>::clear(){

    // clear off the charge grid Q
    int n = myQ.size();
    zomplex *q = myQ.begin();
    for(int i=0;i<n;i++){
      q[i].re = 0.0;
      q[i].im = 0.0;
    }

    // clear off the deltaQ grid DMU_Q
    int m = myDMU_Q.size();
    zomplex *q_Dq = myDMU_Q.begin();
    for(int i=0;i<m;i++){
      q_Dq[i].re = 0.0;
      q_Dq[i].im = 0.0;
    }
  }


  //______________________________________________________________________ print function
  template<class TInterpolation>
  void iSGGrid<TInterpolation>::print(){
    Real q = 0.0;
    report << plain;
    for (unsigned int i = 0; i < myNX; i++){
      for (unsigned int j = 0 ; j < myNY; j++){
	report << "Q["<<i<<"]["<<j<<"][0-"<<myNZ-1<<"] : ";
	for (unsigned int k = 0; k < myNZ; k++){
	  report << "("<<myQ[i][j][k].re <<","<<myQ[i][j][k].im<<")";
	  q += myQ[i][j][k].re;
	}
	report <<std::endl;
      }
    }
    report <<"Sum Q :"<<q<<endr;
  }
}
#endif
