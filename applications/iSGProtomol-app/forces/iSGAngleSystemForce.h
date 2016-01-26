/* -*- c++ -*- */
#ifndef ISGANGLESYSTEMFORCE_H
#define ISGANGLESYSTEMFORCE_H

#include "SystemForce.h"
#include "iSGAngleSystemForceBase.h"
#include "ScalarStructure.h"
#include "Parallel.h"

namespace ProtoMol {
  //_________________________________________________________________ iSGAngleSystemForce
  template<class TBoundaryConditions>
  class iSGAngleSystemForce : public SystemForce, private iSGAngleSystemForceBase {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class iSGAngleSystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void calcAngle(const TBoundaryConditions &boundary,
		   const Angle& currentAngle,
		   const Vector3DBlock* positions, 
		   Vector3DBlock* forces,
		   ScalarStructure* energies);

    Real calcAngleEnergy(const TBoundaryConditions &boundary,
			 const Angle& currentAngle,
			 const Vector3DBlock* positions);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class SystemForce
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void evaluate(const GenericTopology* topo,
			  const Vector3DBlock* positions, 
			  Vector3DBlock* forces,
 			  ScalarStructure* energies);

    virtual void parallelEvaluate(const GenericTopology* topo,
				  const Vector3DBlock* positions, 
				  Vector3DBlock* forces,
				  ScalarStructure* energies);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Force
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getKeyword() const{return keyword;}
    virtual unsigned int numberOfBlocks(const GenericTopology* topo,
					const Vector3DBlock* pos);
  private:
    virtual Force* doMake(std::string&, std::vector<Value>) const {
      return (new iSGAngleSystemForce());
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const{return keyword;}
    virtual unsigned int getParameterSize() const{return 0;}
    virtual void getParameters(std::vector<Parameter>&) const {}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:

  }; 

  //______________________________________________________________________ INLINES

  template<class TBoundaryConditions> 
  inline void iSGAngleSystemForce<TBoundaryConditions>::evaluate(const GenericTopology* topo,
                                                                 const Vector3DBlock* positions,
							         Vector3DBlock* forces,
							         ScalarStructure* energies) {
  
    const TBoundaryConditions &boundary = 
      (dynamic_cast<const SemiGenericTopology<TBoundaryConditions>& >(*topo)).boundaryConditions;
    for (unsigned int i = 0; i < topo->angles.size(); i++) 
      calcAngle(boundary, topo->angles[i], positions, forces, energies);
  }

  template<class TBoundaryConditions> 
  inline void iSGAngleSystemForce<TBoundaryConditions>::calcAngle(const TBoundaryConditions &boundary,
		                                                  const Angle& currentAngle,
							          const Vector3DBlock* positions,
							          Vector3DBlock* forces,
							          ScalarStructure* energies){
  
    int a1 = currentAngle.atom1;
    int a2 = currentAngle.atom2;
    int a3 = currentAngle.atom3;

    Real restAngle     = currentAngle.restAngle;
    Real forceConstant = currentAngle.forceConstant;
    Real ureyBradleyConstant = currentAngle.ureyBradleyConstant;
    Real ureyBradleyRestLength = currentAngle.ureyBradleyRestLength;

    // Variables for computation of the chemical potential difference
    Real DeltaK        = currentAngle.DeltaK;
    Real DeltaTheta0   = currentAngle.DeltaTheta0;
    Real Delta_ubK     = currentAngle.Delta_ubK;
    Real Delta_ubR0    = currentAngle.Delta_ubR0;
  
    Vector3D atom1((*positions)[a1]);
    Vector3D atom2((*positions)[a2]);
    Vector3D atom3((*positions)[a3]);
  
    //Vector3D r12 = atom1 - atom2;                     // Vector from atom 1 to atom 2.
    Vector3D r12(boundary.minimalDifference(atom2,atom1));
    //Vector3D r32 = atom3 - atom2;                     // Vector from atom 3 to atom 2.
    Vector3D r32(boundary.minimalDifference(atom2,atom3));
    //Vector3D r13 = atom1 - atom3;                     // Vector from atom 1 to atom 3.
    Vector3D r13(boundary.minimalDifference(atom3,atom1));
    Real d12 = r12.norm();                        // Distance between atom 1 and 2.
    Real d32 = r32.norm();                        // Distance between atom 3 and 2.
    Real d13 = r13.norm();                        // Distance between atom 1 and 3.
  
  
    // Calculate theta. 
    Real theta = atan2((r12.cross(r32)).norm(),r12.dot(r32));
    Real sinTheta = sin(theta);
    Real cosTheta = cos(theta);
  
    // Calculate dpot/dtheta
    Real dpotdtheta = 2.0 * forceConstant * (theta - restAngle);  
 
    // Calculate dr/dx, dr/dy, dr/dz.
    Vector3D dr12(r12/d12);
    Vector3D dr32(r32/d32);
    Vector3D dr13(r13/d13);
    // Calulate dtheta/dx, dtheta/dy, dtheta/dz.
    Vector3D dtheta1((dr12 * cosTheta - dr32)/(sinTheta*d12));     // atom1
    Vector3D dtheta3((dr32 * cosTheta - dr12)/(sinTheta*d32));     // atom3
  
    // Calculate Urey Bradley force.
    Vector3D ureyBradleyforce1(dr13 * (2.0 * ureyBradleyConstant * (d13 - ureyBradleyRestLength)));
    Vector3D ureyBradleyforce3(-ureyBradleyforce1);
  
    // Calculate force on atom1 due to atom 2 and 3.
    Vector3D force1(-dtheta1 * dpotdtheta - ureyBradleyforce1);
  
    // Calculate force on atom3 due to atom 1 and 2.
    Vector3D force3(-dtheta3 * dpotdtheta - ureyBradleyforce3);
  
    // Calculate force on atom2 due to atom 1 and 3.
    Vector3D force2(-force1 - force3);
     
    // Add to the total force.
    (*forces)[a1] += force1;
    (*forces)[a2] += force2;
    (*forces)[a3] += force3;
  
    // Calculate Energy. 
    Real eHarmonic = forceConstant * (theta - restAngle)*(theta - restAngle); 
    Real eUreyBradley = ureyBradleyConstant * (d13 - ureyBradleyRestLength)*(d13 - ureyBradleyRestLength);

    // Calculate Chemical potential difference
    if (DeltaK != 0.0 && DeltaTheta0 != 0.0) {
      Real DeltaMuHarmonic = eHarmonic * DeltaK / forceConstant - dpotdtheta * DeltaTheta0;
      Real DeltaMuUreyBradley = 0.0;
      if (ureyBradleyConstant != 0) DeltaMuUreyBradley = eUreyBradley * Delta_ubK / ureyBradleyConstant
        - 2.0 * ureyBradleyConstant * (d13 - ureyBradleyRestLength) * Delta_ubR0;

      // add the chemical potential difference from this angle to the total
      (*energies)[ScalarStructure::ANGLE_DELTAMU] += DeltaMuHarmonic + DeltaMuUreyBradley;
    }
       
    // Add Energy 
    (*energies)[ScalarStructure::ANGLE] += eHarmonic + eUreyBradley;
    
//#define ISG_DEBUG_ANGLE
#ifdef ISG_DEBUG_ANGLE
    Report::report.precision(8);
    Report::report << "Angle atoms: " << a1 << ", " << a2 << ", " << a3 << Report::endr;
    Report::report << "k = " << forceConstant <<  ", theta0 = " << restAngle << Report::endr;
    Report::report << "DeltaK = " << DeltaK <<  ", DeltaTheta0 = " << DeltaTheta0 << Report::endr;
    Report::report << "theta = " << theta << Report::endr;
    Report::report << "current DeltaMu (angle) = " << (*energies)[ScalarStructure::ANGLE_DELTAMU] << Report::endr;
    Report::report << "current energy (angle) = " << (*energies)[ScalarStructure::ANGLE] << Report::endr;
#endif

    // Add virial
    if(energies->virial()) {
      Real xy = force1.x * r12.y + force3.x * r32.y;
      Real xz = force1.x * r12.z + force3.x * r32.z;
      Real yz = force1.y * r12.z + force3.y * r32.z;
      (*energies)[ScalarStructure::VIRIALXX] += force1.x * r12.x + force3.x * r32.x;
      (*energies)[ScalarStructure::VIRIALXY] += xy;
      (*energies)[ScalarStructure::VIRIALXZ] += xz;
      (*energies)[ScalarStructure::VIRIALYX] += xy;
      (*energies)[ScalarStructure::VIRIALYY] += force1.y * r12.y + force3.y * r32.y;
      (*energies)[ScalarStructure::VIRIALYZ] += yz;
      (*energies)[ScalarStructure::VIRIALZX] += xz;
      (*energies)[ScalarStructure::VIRIALZY] += yz;
      (*energies)[ScalarStructure::VIRIALZZ] += force1.z * r12.z + force3.z * r32.z;
    }    
  }

  template<class TBoundaryConditions> 
  inline Real iSGAngleSystemForce<TBoundaryConditions>::calcAngleEnergy(const TBoundaryConditions &boundary,
								        const Angle& currentAngle,
								        const Vector3DBlock* positions)
  {
    int a1, a2, a3;
    Real restAngle, forceConstant, ureyBradleyConstant, ureyBradleyRestLength;
    Vector3D atom1, atom2, atom3, r12, r32, r13;
    Real d13, theta, eHarmonic, eUreyBradley;

    a1 = currentAngle.atom1;
    a2 = currentAngle.atom2;
    a3 = currentAngle.atom3;
    restAngle     = currentAngle.restAngle;
    forceConstant = currentAngle.forceConstant;
    ureyBradleyConstant = currentAngle.ureyBradleyConstant;
    ureyBradleyRestLength = currentAngle.ureyBradleyRestLength;

    atom1 = (*positions)[a1];
    atom2 = (*positions)[a2];
    atom3 = (*positions)[a3];

    //r12 = atom1 - atom2;                     // Vector from atom 1 to atom 2.
    r12 = boundary.minimalDifference(atom2,atom1);   
    //r32 = atom3 - atom2;                     // Vector from atom 3 to atom 2.
    r32 = boundary.minimalDifference(atom2,atom3);    
    //r13 = atom1 - atom3;                     // Vector from atom 1 to atom 3.
    r13 = boundary.minimalDifference(atom3,atom1);    
    d13 = r13.norm();                        // Distance between atom 1 and 3.

    // Calculate theta.
    theta = atan2((r12.cross(r32)).norm(),r12.dot(r32));

    // Calculate Energy. 
    eHarmonic = forceConstant * (theta - restAngle)*(theta - restAngle); 
    eUreyBradley = ureyBradleyConstant * (d13 - ureyBradleyRestLength)
      *(d13 - ureyBradleyRestLength);

    return (eHarmonic + eUreyBradley);
  }
  
  template<class TBoundaryConditions> 
  inline void iSGAngleSystemForce<TBoundaryConditions>::parallelEvaluate(const GenericTopology* topo,
								         const Vector3DBlock* positions, 
								         Vector3DBlock* forces,
								         ScalarStructure* energies)
  {
    const TBoundaryConditions &boundary = 
      (dynamic_cast<const SemiGenericTopology<TBoundaryConditions>& >(*topo)).boundaryConditions;

    unsigned int n = topo->angles.size();
    unsigned int count = numberOfBlocks(topo,positions);
  
    for(unsigned int i = 0;i<count;i++){
      if(Parallel::next()){
	int to = (n*(i+1))/count;
	if(to > static_cast<int>(n))
	  to = n;
	int from = (n*i)/count;
	for (int j = from; j < to; j++)
	  calcAngle(boundary, topo->angles[j], positions, forces, energies);
      }
    }
  }

  template<class TBoundaryConditions> 
  inline unsigned int iSGAngleSystemForce<TBoundaryConditions>::numberOfBlocks(const GenericTopology* topo,
									       const Vector3DBlock*){    
    return std::min(Parallel::getAvailableNum(),static_cast<int>(topo->angles.size()));
  }

}
#endif /* ISGANGLESYSTEMFORCE_H */
