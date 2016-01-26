#include "BSplineMOLLYIntegrator.h"
#include "ReducedHessAngle.h"
#include "reducedHessBond.h"
#include "ForceGroup.h"
#include "SystemForce.h"
#include "mathutilities.h"
#include "Report.h"
#include "pmconstants.h"
#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "Vector3DBlock.h"

using namespace ProtoMol::Report;
using std::vector;
using std::string;
using std::map;

namespace ProtoMol {

  //_____________________________________________________BSplineMOLLYIntegrator
  const string BSplineMOLLYIntegrator::keyword("BSplineMOLLY");
  
  BSplineMOLLYIntegrator::BSplineMOLLYIntegrator() : MOLLYIntegrator(),
						     myCached(false),
						     myMOLLYForcesBonded(NULL),
						     myMOLLYForcesHBonded(NULL),
						     myHBondForces(NULL),
						     myMOLLYPositions(NULL),
						     myMOLLYVelocities(NULL),
						     myMOLLYForces(NULL),
						     myMOLLYEnergies(NULL),
						     myB(NULL),
						     myAngleFilter(NULL),  
						     myPxAngle(NULL),
						     myXxAngle(NULL),
						     myBxAngle(NULL),
						     myTypeOfBSpline(),
						     myMOLLYStepsize(0.0),
						     myMOLLYStepsizeP(0.0),
						     myNumIter(0),
						     myBond(false),
						     myAngle(false),
						     myCoulomb(false),
						     myLennardJones(false){}
  
  
  BSplineMOLLYIntegrator::BSplineMOLLYIntegrator(int cycles,
						 const BSplineType& typeOfBSpline,
						 Real mollyStepsize,
						 ForceGroup *overloadedForces,
						 StandardIntegrator *nextIntegrator)
    : MOLLYIntegrator(cycles,overloadedForces,nextIntegrator),
      myCached(false),
      myMOLLYForcesBonded(new ForceGroup),
      myMOLLYForcesHBonded(new ForceGroup),
      myHBondForces(new ForceGroup),
      myMOLLYPositions(new Vector3DBlock),
      myMOLLYVelocities(new Vector3DBlock),
      myMOLLYForces(new Vector3DBlock),
      myMOLLYEnergies(new ScalarStructure),
      myB(new Vector3DBlock),
      myAngleFilter(new vector< ReducedHessAngle >()),  
      myPxAngle(new vector< ReducedHessAngle >()),
      myXxAngle(new vector< ReducedHessAngle >()),
      myBxAngle(new vector< ReducedHessAngle >()),
      myTypeOfBSpline(typeOfBSpline),
      myMOLLYStepsize(mollyStepsize/Constant::TIMEFACTOR),
      myMOLLYStepsizeP(mollyStepsize),
      myNumIter(0),
      myBond(false),
      myAngle(false),
      myCoulomb(false),
      myLennardJones(false){

    // Retrieve all meta forces and assign them to the foruce groups
    vector<Force*> mollyForces = overloadedForces->getDeepMetaForces();
    for(unsigned int i=0;i<mollyForces.size();i++){
      if(equalNocase(mollyForces[i]->getId(),"MollyBond")){
	myMOLLYForcesBonded->addForce(Makeable::copy(mollyForces[i]));
	myBond = true;
      }
      else if(equalNocase(mollyForces[i]->getId(),"MollyAngle")){
	myMOLLYForcesBonded->addForce(Makeable::copy(mollyForces[i]));
	myAngle = true;
      }
      else if(equalNocase(mollyForces[i]->getId(),"MollyCoulomb")){
	myHBondForces->addForce(Makeable::copy(mollyForces[i]));
	myCoulomb = true;
      }
      else if(equalNocase(mollyForces[i]->getId(),"HBondCoulomb")){
	myMOLLYForcesHBonded->addForce(Makeable::copy(mollyForces[i]));
	myCoulomb = true;
      }
      else if(equalNocase(mollyForces[i]->getId(),"MollyLennardJones")){
	myHBondForces->addForce(Makeable::copy(mollyForces[i]));
	myLennardJones = true;
      }
      else if(equalNocase(mollyForces[i]->getId(),"HBondLennardJones")){
	myMOLLYForcesHBonded->addForce(Makeable::copy(mollyForces[i]));
	myLennardJones = true;
      }
      else if(equalNocase(mollyForces[i]->getId(),"MollyLennardJonesCoulomb")){
	myHBondForces->addForce(Makeable::copy(mollyForces[i]));
	myLennardJones = true;
	myCoulomb = true;
      }
      else if(equalNocase(mollyForces[i]->getId(),"HBondLennardJonesCoulomb")){
	myMOLLYForcesHBonded->addForce(Makeable::copy(mollyForces[i]));
	myLennardJones = true;
	myCoulomb = true;
      }
      else {
	report << recoverable << getId() <<" could not recognize molly force \'"
	       <<mollyForces[i]->getId()<<"\'." <<endr;
      }
    }
    setMyNumIterMOLLYStepsize();
    report << hint << "iteration MOLLY step size="<<myMOLLYStepsize*Constant::TIMEFACTOR<<endr;
  }


  BSplineMOLLYIntegrator::~BSplineMOLLYIntegrator() {
    if(myMOLLYForcesBonded != NULL)     delete myMOLLYForcesBonded;
    if(myMOLLYForcesHBonded != NULL)    delete myMOLLYForcesHBonded; 
    if(myHBondForces != NULL)           delete myHBondForces;
    if(myMOLLYPositions != NULL)        delete myMOLLYPositions;
    if(myMOLLYVelocities != NULL)       delete myMOLLYVelocities;
    if(myMOLLYForces != NULL)           delete myMOLLYForces;
    if(myMOLLYEnergies != NULL)         delete myMOLLYEnergies;
    if(myB != NULL)                     delete myB;
    if(myAngleFilter != NULL)           delete myAngleFilter;  
    if(myPxAngle != NULL)               delete myPxAngle;  
    if(myXxAngle != NULL)               delete myXxAngle;  
    if(myBxAngle != NULL)               delete myBxAngle;  
  }

  void BSplineMOLLYIntegrator::initialize(GenericTopology *topo,
					  Vector3DBlock *positions,
					  Vector3DBlock *velocities,
					  ScalarStructure *energies){
    
    MOLLYIntegrator::initialize(topo,positions,velocities,energies);

    if(!myCached)
      init();
    
    myMOLLYPositions->resize(myPositions->size());
    myMOLLYPositions->zero();
    myMOLLYVelocities->resize(myPositions->size());
    myMOLLYVelocities->zero();
    myMOLLYForces->resize(myPositions->size());
    myMOLLYForces->zero();
    myB->resize(myPositions->size());
    myB->zero();
    myAngleFilter->resize(myAngleIndexes.size());
    myPxAngle->resize(myAngleIndexes.size());
    myXxAngle->resize(myAngleIndexes.size());
    myBxAngle->resize(myAngleIndexes.size());
    
    initializeForces();
  }



  Vector3DBlock* BSplineMOLLYIntegrator::doAveragingPositions() {
    if(!myCached)
      init();
    myMOLLYVelocities->zero();// init myMOLLYVelocities
    myMOLLYPositions->intoAssign(*myPositions);// init myMOLLYPositions
    
    calcHessiansBondsAnglesHBonds();
    updateB_Bx_Px_for1stHalfKick();
    updateXx();
    
    calcMOLLYForcesHalfKickOneDrift();
    for(unsigned int i = 1; i < myNumIter; i++){
      calcHessiansBondsAnglesHBonds();
      updateB_Bx_Px_for1Kick();
      updateXx();
      calcMOLLYForcesOneKickOneDrift();
    }
    updateB_Bx();
    // last half step, no need to update Px.
    
    myMOLLYPositions->intoWeighted(1./myNumIter, *myB);
    
    for (unsigned int i = 0; i < myAngleIndexes.size(); i++)
      (*myAngleFilter)[i] = (*myBxAngle)[i] / myNumIter;
    
    return myMOLLYPositions;
  }


  // Mollify the slow force. Note that the anglefilter is assembled from bond 
  // Hessian, angle Hessian and hbond Hessians. 
  void BSplineMOLLYIntegrator::doMollification(Vector3DBlock*){
    if(!myCached)
      init();

    //  (*myAngleFilter) needes to be transposed.
    for (unsigned int i = 0; i < myAngleIndexes.size(); i++){
      int a1 = myTopo->angles[myAngleIndexes[i].angle].atom1;
      int a2 = myTopo->angles[myAngleIndexes[i].angle].atom2;
      int a3 = myTopo->angles[myAngleIndexes[i].angle].atom3;
      
      //    Matrix3by3 I(1,0,0, 0,1,0, 0,0,1);
      ReducedHessAngle tm = (*myAngleFilter)[i].transposed();
      
      Vector3D ftemp1=tm(0,0)*(*myForces)[a1] + tm(0,1)*(*myForces)[a2] + tm(0,2)*(*myForces)[a3];
      Vector3D ftemp2=tm(1,0)*(*myForces)[a1] + tm(1,1)*(*myForces)[a2] + tm(1,2)*(*myForces)[a3];
      Vector3D ftemp3=tm(2,0)*(*myForces)[a1] + tm(2,1)*(*myForces)[a2] + tm(2,2)*(*myForces)[a3];
      
      (*myForces)[a1] = ftemp1; 
      (*myForces)[a2] = ftemp2; 
      (*myForces)[a3] = ftemp3;
    } 
  }
  
  
  void BSplineMOLLYIntegrator::calcMOLLYForcesOneKickOneDrift(){
    calculateMOLLYForcesBonded();
    for (unsigned int i=0;i<myMOLLYVelocities->size();i++){
      Real mass = myTopo->atoms[i].scaledMass;
      (*myMOLLYVelocities)[i] += (*myMOLLYForces)[i]*(myMOLLYStepsize/mass); 
      (*myMOLLYPositions)[i] += (*myMOLLYVelocities)[i]*myMOLLYStepsize;
    }
    if(myMOLLYForcesHBonded->anyForces()){
      calculateMOLLYForcesHBonded();
      for (unsigned int i=0;i<myMOLLYVelocities->size();i++){
	Real mass = myTopo->atoms[i].scaledMass;
	if(mass < 2){// so that the heavy atoms do not move. Hydrogen has mass 1.0008
	  (*myMOLLYVelocities)[i] += (*myMOLLYForces)[i]*(myMOLLYStepsize/mass); 
	  (*myMOLLYPositions)[i] += (*myMOLLYVelocities)[i]*myMOLLYStepsize;
	}
      }
    }
  }
  

  void BSplineMOLLYIntegrator::calcMOLLYForcesHalfKickOneDrift(){
    calculateMOLLYForcesBonded();
    for (unsigned int i=0;i<myMOLLYVelocities->size();i++){
      Real mass = myTopo->atoms[i].scaledMass;
      (*myMOLLYVelocities)[i] += (*myMOLLYForces)[i]*(myMOLLYStepsize/mass/2); 
      (*myMOLLYPositions)[i] += (*myMOLLYVelocities)[i]*myMOLLYStepsize;
    }
    if(myMOLLYForcesHBonded->anyForces()){
      calculateMOLLYForcesHBonded();
      for (unsigned int i=0;i<myMOLLYVelocities->size();i++){
	Real mass = myTopo->atoms[i].scaledMass;
	if(mass < 2){// so that the heavy atoms do not move. Hydrogen has mass 1.0008
	  (*myMOLLYVelocities)[i] += (*myMOLLYForces)[i]*(myMOLLYStepsize/mass/2); 
	  (*myMOLLYPositions)[i] += (*myMOLLYVelocities)[i]*myMOLLYStepsize;
	}
      }
    }
  }
  
  
  // Energies have to be cleared for STS integrator
  // Potential energy is U^{fastest}(X1), forces are gradient 
  // of potential, thus they are set to zero here also
  void BSplineMOLLYIntegrator::calculateMOLLYForcesBonded() {
    myMOLLYEnergies->clear();
    myMOLLYForces->zero();
    myMOLLYForcesBonded->evaluateSystemForces(myTopo,
					      myMOLLYPositions,
					      myMOLLYForces,
					      myMOLLYEnergies);
  }

  void BSplineMOLLYIntegrator::calculateMOLLYForcesHBonded() {
    myMOLLYEnergies->clear();
    myMOLLYForces->zero();
    myMOLLYForcesHBonded->evaluateSystemForces(myTopo,
					       myMOLLYPositions,
					       myMOLLYForces,
					       myMOLLYEnergies);
  }
  

  // calculate the Hessian matrices for bond and/or angle. After this function call,
  // we have the ReducedHessBond for every bond pair, saved in *myBondFilter and 
  // all the ReducedHessAngle for every angle triplet, in *myAngleFilter.
  void BSplineMOLLYIntegrator::calcHessiansBondsAnglesHBonds() {
    for (unsigned int i = 0; i < myAngleIndexes.size(); i++){
      int a1 = myTopo->angles[myAngleIndexes[i].angle].atom1;
      int a2 = myTopo->angles[myAngleIndexes[i].angle].atom2;
      int a3 = myTopo->angles[myAngleIndexes[i].angle].atom3; 
      if(myAngle){
	Real theta0  = myTopo->angles[myAngleIndexes[i].angle].restAngle;
	Real k_t     = myTopo->angles[myAngleIndexes[i].angle].forceConstant;
	(*myAngleFilter)[i].evaluate((*myMOLLYPositions)[a1],(*myMOLLYPositions)[a2],
				     (*myMOLLYPositions)[a3],k_t,theta0);
	// ReducedHessAngle for atoms a1, a2 and a3
      }else
	(*myAngleFilter)[i].clear();
      
      if(myBond){
	Real r_0  = myTopo->bonds[myAngleIndexes[i].bond1].restLength;
	Real k    = myTopo->bonds[myAngleIndexes[i].bond1].springConstant;
	Matrix3by3 bondHess12 = reducedHessBond((*myMOLLYPositions)[a1],(*myMOLLYPositions)[a2],k,r_0);
	
	r_0  = myTopo->bonds[myAngleIndexes[i].bond2].restLength;
	k    = myTopo->bonds[myAngleIndexes[i].bond2].springConstant;
	Matrix3by3 bondHess23 = reducedHessBond((*myMOLLYPositions)[a2],(*myMOLLYPositions)[a3],k,r_0);
	
	(*myAngleFilter)[i].accumulateTo(0,0,bondHess12);
	(*myAngleFilter)[i].accumulateTo(1,1,bondHess12);
	(*myAngleFilter)[i].accumulateNegTo(1,0,bondHess12);
	(*myAngleFilter)[i].accumulateNegTo(0,1,bondHess12);
	
	(*myAngleFilter)[i].accumulateTo(1,1,bondHess23);
	(*myAngleFilter)[i].accumulateTo(2,2,bondHess23);
	(*myAngleFilter)[i].accumulateNegTo(1,2,bondHess23);
	(*myAngleFilter)[i].accumulateNegTo(2,1,bondHess23);
      }
    }
    
    // if only bonds and angles are included, then it is a B-spline MOLLY integrator
    myHBondForces->evaluateMollyForces(myTopo,myMOLLYPositions,myAngleFilter);
    //here the anglefilter is updated to incorporate the hbond hessian.
  }
  
  
  // X: = x, Xx = I
  // P: = 0, Px = 0
  // B: = 0, Bx = 0
  // Initiate these matrices before integrate equations of motion and 
  // auxiliary equations for mollification matrix, Ax(x).
  // for the first half kick, the update is easily done:
  void BSplineMOLLYIntegrator::updateB_Bx_Px_for1stHalfKick(){

    myB->intoWeighted(0.5,*myMOLLYPositions);
    // update B for the first half kick
    
    for (unsigned int i = 0; i < myAngleIndexes.size(); i++){
      (*myPxAngle)[i] = (*myAngleFilter)[i]*(-0.5);
      (*myXxAngle)[i].identity();
      (*myBxAngle)[i] = (*myXxAngle)[i]*0.5;
    }// init the matrices to identity or zero.
  }
  
  // note that the Px for angle should be negatively accumulated because
  // Fx = -H
  void BSplineMOLLYIntegrator::updateB_Bx_Px_for1Kick() {
    myB->intoAdd(*myMOLLYPositions);
    // update B, getting rid of the dt factor (later we would not bother with dividing
    // B by Dt, instead, we divide B with myNumIter.
    
    for (unsigned int i = 0; i < myAngleIndexes.size(); i++){
      (*myPxAngle)[i]-=(*myAngleFilter)[i]*(*myXxAngle)[i];
      (*myBxAngle)[i]+=(*myXxAngle)[i];
    }// update Bx and Px
  }
  
  
  void BSplineMOLLYIntegrator::updateB_Bx() {
    myB->intoWeightedAdd(0.5,*myMOLLYPositions);
    // update B for the last half step
    
    for (unsigned int i = 0; i < myAngleIndexes.size(); i++){
      (*myBxAngle)[i]+=(*myXxAngle)[i]*0.5;
    }// update Bx for the last half step
  }
  
  // update Xx
  void BSplineMOLLYIntegrator::updateXx() {
    for (unsigned int i = 0; i < myAngleIndexes.size(); i++){
      int a1 = myTopo->angles[myAngleIndexes[i].angle].atom1;
      int a2 = myTopo->angles[myAngleIndexes[i].angle].atom2;
      int a3 = myTopo->angles[myAngleIndexes[i].angle].atom3;
      
      Real M1 = myTopo->atoms[a1].scaledMass;
      Real M2 = myTopo->atoms[a2].scaledMass;
      Real M3 = myTopo->atoms[a3].scaledMass;
      // we know for water M3 = M1;
      
      Real t0 = (myMOLLYStepsizeSquare/M1);
      Real t1 = (myMOLLYStepsizeSquare/M2);
      Real t2 = (myMOLLYStepsizeSquare/M3);
      
      for(int j=0;j<3;j++){
	(*myXxAngle)[i].accumulateTo(0,j,(*myPxAngle)[i](0,j) * t0);
	(*myXxAngle)[i].accumulateTo(1,j,(*myPxAngle)[i](1,j) * t1);
	(*myXxAngle)[i].accumulateTo(2,j,(*myPxAngle)[i](2,j) * t2);
      }
    }
  }
  




  //myNumIter is constant throughout the simulation. This function is called
  //only once in the constructor. For BSplineMOLLY, only SHORT and LONG averaging
  //are supported. SHORT is not as stable as LONG. LONG is prefered for 
  //long-time-step simulations.
  void BSplineMOLLYIntegrator::setMyNumIterMOLLYStepsize(){

    Real thisLevelStepsize = getTimestep()/Constant::TIMEFACTOR;
    Real temp = thisLevelStepsize/myMOLLYStepsize-0.01;
    myNumIter = static_cast<unsigned int>(ceil(temp));
    myMOLLYStepsize = thisLevelStepsize/myNumIter; 
    
    myMOLLYStepsizeSquare = myMOLLYStepsize * myMOLLYStepsize;
    // to be used in the Xx updates. In Px updates, the dt factor is removed
    // for efficiency, (later compensated in the Xx updates). They are equivalent.
    // In terms of dynamics: displacement = velocity * t, .or.
    //                       displacement = 0.5 * accelerating rate * t * t.
    
    if(myTypeOfBSpline == BSplineType::LONG){
      return;
    }else if(myTypeOfBSpline == BSplineType::SHORT){
      if(myNumIter%2 == 1){
	myMOLLYStepsize =  thisLevelStepsize/(myNumIter+1); 
	myNumIter = (myNumIter+1)/2;
      }else
	myNumIter =  myNumIter/2;
    }
    
    // the above change is necessary, e.g., if myThisLevelStepsize == 5, and STS 
    // level step size equals 1, then myNumIter = ceil(5/2) = 3, and 
    // myMOLLYStepsize = 5.0/2/3 = 5.0/6. This adjustment makes the SHORT method
    // always consistent (integral of the weight function from -inf to inf is 1).
    else{
      myTypeOfBSpline = BSplineType::LONG;
      report << recoverable << "[BSplineMOLLYIntegrator::setMyNumIterMOLLYStepsize] "
	     << "BSplineMOLLY only supports long and short as averaging B-splines."   << endr 
	     << "Integration proceeds using long as averaging B-spline." << endr;
    }
  }
  
  void BSplineMOLLYIntegrator::getParameters(vector<Parameter>& parameters) const {
    MTSIntegrator::getParameters(parameters);
    parameters.push_back(Parameter("BSplineType",Value((string)myTypeOfBSpline,ConstraintValueType::NotEmpty()),Text(std::string("interpolation scheme (")+BSplineType::getPossibleValues()+std::string(")"))));
    parameters.push_back(Parameter("mollyStepsize",Value(myMOLLYStepsizeP,ConstraintValueType::Positive())));
  }
 
  MTSIntegrator* BSplineMOLLYIntegrator::doMake(std::string& errMsg, const std::vector<Value>& values, ForceGroup* fg, StandardIntegrator *nextIntegrator)const{
    BSplineType bsplineType(values[1].getString());
    if(!bsplineType.valid()){
      errMsg += " BSplineType \'"+values[1].getString()+"\' not valid.";
      return NULL;
    }
    return (new BSplineMOLLYIntegrator(values[0],bsplineType,values[2],fg,nextIntegrator));
  }


  void BSplineMOLLYIntegrator::init(){
    map<BondIndex,int> bonds;
    for (unsigned int i = 0; i < myTopo->bonds.size(); i++){
      BondIndex b(myTopo->bonds[i].atom1,myTopo->bonds[i].atom2);
      if(bonds.find(b) == bonds.end()){
	bonds[b] = i;
      }
      else {
	int j = bonds[b];
	report <<recoverable<<"[BSplineMOLLYIntegrator::initialize] bond["<<i+1<<"]("<<b.b1+1<<","<<b.b2+1<<") already found at bond["<<j+1<<"]!"<<endr;
      }
    }

    myAngleIndexes.clear();
    for (unsigned int i = 0; i < myTopo->angles.size(); i++){
      BondIndex a1(myTopo->angles[i].atom1,myTopo->angles[i].atom2);
      BondIndex a2(myTopo->angles[i].atom2,myTopo->angles[i].atom3);
      if(bonds.find(a1) != bonds.end() && bonds.find(a2) != bonds.end()){
	myAngleIndexes.push_back(AngleIndex(i,bonds[a1],bonds[a2]));
      }
      else {
	report <<recoverable<<"[BSplineMOLLYIntegrator::initialize] angle["<<i+1<<"]("
	       <<myTopo->angles[i].atom1+1<<","<<myTopo->angles[i].atom2+1<<","<<myTopo->angles[i].atom3+1<<" has not two bonds"<<endr;
      }
    }
    myCached = true;
  }

  void BSplineMOLLYIntegrator::doUncache(){
    myCached = false;
  }
 
}
