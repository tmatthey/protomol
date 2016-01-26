#include "CheckpointInputStream.h"

void CheckpointInputStream::readArchive(StateRestore& rst)
{
    (*this) >> rst.mySystemsize;
    // Integrator levels
    (*this) >> rst.myLevels;
    // Then the timestep
    (*this) >> rst.myTimestep;
    // Random number generator state
    (*this) >> rst.myCRN[0];
    (*this) >> rst.myCRN[1];
    (*this) >> rst.myCRN[2];
    // Simulation box (All are 0 for VBC)
    (*this) >> rst.myE1.x;
    (*this) >> rst.myE1.y;
    (*this) >> rst.myE1.z;
    (*this) >> rst.myE2.x;
    (*this) >> rst.myE2.y;
    (*this) >> rst.myE2.z;
    (*this) >> rst.myE3.x;
    (*this) >> rst.myE3.y;
    (*this) >> rst.myE3.z;
    (*this) >> rst.myOrigin.x;
    (*this) >> rst.myOrigin.y;
    (*this) >> rst.myOrigin.z;
//     // Atom information
     rst.myPositions->resize(rst.mySystemsize);
     rst.myVelocities->resize(rst.mySystemsize);
     rst.myForces.resize(rst.myLevels+1);
     for (unsigned int i = 0; i <= rst.myLevels; i++) {
       rst.myForces[i] = new Vector3DBlock();
       rst.myForces[i]->resize(rst.mySystemsize);
     }
     for (unsigned int i = 0; i < rst.mySystemsize; i++) {
       // Position
       (*this) >> (*(rst.myPositions))[i].x;
       (*this) >> (*(rst.myPositions))[i].y;
       (*this) >> (*(rst.myPositions))[i].z;
       // Velocity
       (*this) >> (*(rst.myVelocities))[i].x;
       (*this) >> (*(rst.myVelocities))[i].y;
       (*this) >> (*(rst.myVelocities))[i].z;
       // Force
       for (unsigned int j = 0; j <= rst.myLevels; j++) {
	 (*this) >> (*((rst.myForces)[j]))[i].x;
	 (*this) >> (*((rst.myForces)[j]))[i].y;
	 (*this) >> (*((rst.myForces)[j]))[i].z;
       }
     }
//     // Energies
     unsigned int first = (rst.myEnergies)->FIRST;
     unsigned int llast = (rst.myEnergies)->LAST;
     for (unsigned int i = first; i < llast; i++) {
       Real r;
       (*this) >> r;
       (*(rst.myEnergies))[(ScalarStructure::Index)i] = r;
     }
}

