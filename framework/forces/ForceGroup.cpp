#include "ForceGroup.h"

#include "SystemForce.h"
#include "Parallel.h"
#include "Vector3DBlock.h"
#include "GenericTopology.h"
#include "ScalarStructure.h"
#include "ExtendedForce.h"
#include "MollyForce.h"
#include "MetaForce.h"
#include "TimerStatistic.h"

//#define DEBUG_OUTSTANDING_MSG

#ifdef DEBUG_OUTSTANDING_MSG
#include <mpi.h>
#endif
using namespace ProtoMol::Report;
using std::string;
using std::vector;
using std::list;

namespace ProtoMol
{
	//_________________________________________________________________ ForceGroup
	ForceGroup::ForceGroup()
	{
	}


	ForceGroup::~ForceGroup()
	{
		for (list<SystemForce*>::iterator currentForce = mySystemForcesList.begin(); currentForce != mySystemForcesList.end(); ++currentForce)
			delete (*currentForce);
		for (list<ExtendedForce*>::iterator currentForce = myExtendedForcesList.begin(); currentForce != myExtendedForcesList.end(); ++currentForce)
			delete (*currentForce);
		for (list<MollyForce*>::iterator currentForce = myMollyForcesList.begin(); currentForce != myMollyForcesList.end(); ++currentForce)
			delete (*currentForce);
		for (list<MetaForce*>::iterator currentForce = myMetaForcesList.begin(); currentForce != myMetaForcesList.end(); ++currentForce)
			delete (*currentForce);
	}

	void ForceGroup::evaluateSystemForces(GenericTopology* topo,
	                                      const Vector3DBlock* positions,
	                                      Vector3DBlock* forces,
	                                      ScalarStructure* energies) const
	{
		if (mySystemForcesList.empty())
			return;
		topo->uncacheCellList();
		Parallel::distribute(energies, forces);

		if (Parallel::isDynamic())
		{
			// Collecting the number of blocks of each force.
			vector<int> blocks;
			for (list<SystemForce*>::const_iterator currentForce = mySystemForcesList.begin();
			     currentForce != mySystemForcesList.end();
			     ++currentForce)
			{
				blocks.push_back((*currentForce)->numberOfBlocks(topo, positions));
			}
			Parallel::resetNext(blocks);
		}

		if (Parallel::iAmSlave())
		{
			Parallel::resetNext();
			TimerStatistic::timer[TimerStatistic::FORCES].start();
			for (list<SystemForce*>::const_iterator currentForce = mySystemForcesList.begin();
			     currentForce != mySystemForcesList.end();
			     ++currentForce)
			{
				if (Parallel::isParallel())
				{
					(*currentForce)->parallelEvaluate(topo, positions, forces, energies);
				}
				else
				{
					(*currentForce)->evaluate(topo, positions, forces, energies);
				}
			}
			TimerStatistic::timer[TimerStatistic::FORCES].stop();
		}
#ifdef DEBUG_OUTSTANDING_MSG
    report << allnodes <<plain <<"Node "<<Parallel::getId()<<" done."<<endr;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Status status;
    int test =0;
    MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD, &test, &status);
    if(test!= 0)
      report << plain <<allnodes<<"Node "<<Parallel::getId()<<" outstanding msg from "<<status.MPI_SOURCE<<endr;
#endif    

		Parallel::reduce(energies, forces);
	}

	// Evaluate all system forces in this group.

	void ForceGroup::evaluateExtendedForces(GenericTopology* topo,
	                                        const Vector3DBlock* positions,
	                                        const Vector3DBlock* velocities,
	                                        Vector3DBlock* forces,
	                                        ScalarStructure* energies) const
	{
		if (myExtendedForcesList.empty())
			return;

		topo->uncacheCellList();

		Parallel::distribute(energies, forces);

		if (Parallel::isDynamic())
		{
			// Collecting the number of blocks of each force.
			vector<int> blocks;
			for (list<ExtendedForce*>::const_iterator currentForce = myExtendedForcesList.begin();
			     currentForce != myExtendedForcesList.end();
			     ++currentForce)
			{
				blocks.push_back((*currentForce)->numberOfBlocks(topo, positions));
			}
			Parallel::resetNext(blocks);
		}

		if (Parallel::iAmSlave())
		{
			Parallel::resetNext();
			TimerStatistic::timer[TimerStatistic::FORCES].start();
			for (list<ExtendedForce*>::const_iterator currentForce = myExtendedForcesList.begin();
			     currentForce != myExtendedForcesList.end();
			     ++currentForce)
			{
				if (Parallel::isParallel())
				{
					(*currentForce)->parallelEvaluate(topo, positions, velocities, forces, energies);
				}
				else
				{
					(*currentForce)->evaluate(topo, positions, velocities, forces, energies);
				}
			}
			TimerStatistic::timer[TimerStatistic::FORCES].stop();
		}

#ifdef DEBUG_OUTSTANDING_MSG
    report << allnodes <<plain <<"Node "<<Parallel::getId()<<" done."<<endr;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Status status;
    int test =0;
    MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD, &test, &status);
    if(test!= 0)
      report << plain <<allnodes<<"Node "<<Parallel::getId()<<" outstanding msg from "<<status.MPI_SOURCE<<endr;
#endif    

		Parallel::reduce(energies, forces);
	}

	void ForceGroup::evaluateMollyForces(GenericTopology* topo,
	                                     const Vector3DBlock* positions,
	                                     vector<ReducedHessAngle>* angleFilter) const
	{
		if (myMollyForcesList.empty())
			return;

		topo->uncacheCellList();

		TimerStatistic::timer[TimerStatistic::FORCES].start();

		for (list<MollyForce*>::const_iterator currentForce = myMollyForcesList.begin();
		     currentForce != myMollyForcesList.end();
		     ++currentForce)
			(*currentForce)->evaluate(topo, positions, angleFilter);

		TimerStatistic::timer[TimerStatistic::FORCES].stop();
	}

	void ForceGroup::addSystemForce(SystemForce* force)
	{
		if (force != NULL)
			mySystemForcesList.push_back(force);
	}


	void ForceGroup::addExtendedForce(ExtendedForce* force)
	{
		if (force != NULL)
			myExtendedForcesList.push_back(force);
	}


	void ForceGroup::addMollyForce(MollyForce* force)
	{
		if (force != NULL)
			myMollyForcesList.push_back(force);
	}


	void ForceGroup::addMetaForce(MetaForce* force)
	{
		if (force != NULL)
			myMetaForcesList.push_back(force);
	}


	void ForceGroup::addForce(Force* force)
	{
		force->addToForceGroup(this);
	}


	void ForceGroup::getDefinition(vector<MakeableDefinition>& forces) const
	{
		for (list<SystemForce*>::const_iterator currentForce = mySystemForcesList.begin();
		     currentForce != mySystemForcesList.end();
		     ++currentForce)
			forces.push_back((*currentForce)->getDefinition());

		for (list<ExtendedForce*>::const_iterator currentForce = myExtendedForcesList.begin();
		     currentForce != myExtendedForcesList.end();
		     ++currentForce)
			forces.push_back((*currentForce)->getDefinition());

		for (list<MollyForce*>::const_iterator currentForce = myMollyForcesList.begin();
		     currentForce != myMollyForcesList.end();
		     ++currentForce)
			forces.push_back((*currentForce)->getDefinition());

		for (list<MetaForce*>::const_iterator currentForce = myMetaForcesList.begin();
		     currentForce != myMetaForcesList.end();
		     ++currentForce)
			forces.push_back((*currentForce)->getDefinition());
	}

	void ForceGroup::uncache()
	{
		for (list<SystemForce*>::iterator currentForce = mySystemForcesList.begin();
		     currentForce != mySystemForcesList.end();
		     ++currentForce)
			(*currentForce)->uncache();

		for (list<ExtendedForce*>::iterator currentForce = myExtendedForcesList.begin();
		     currentForce != myExtendedForcesList.end();
		     ++currentForce)
			(*currentForce)->uncache();

		for (list<ExtendedForce*>::iterator currentForce = myExtendedForcesList.begin();
		     currentForce != myExtendedForcesList.end();
		     ++currentForce)
			(*currentForce)->uncache();

		for (list<MetaForce*>::iterator currentForce = myMetaForcesList.begin();
		     currentForce != myMetaForcesList.end();
		     ++currentForce)
			(*currentForce)->uncache();
	}


	vector<Force*> ForceGroup::getForces() const
	{
		vector<Force*> res;
		for (list<SystemForce*>::const_iterator currentForce = mySystemForcesList.begin();
		     currentForce != mySystemForcesList.end();
		     ++currentForce)
			res.push_back(*currentForce);

		for (list<ExtendedForce*>::const_iterator currentForce = myExtendedForcesList.begin();
		     currentForce != myExtendedForcesList.end();
		     ++currentForce)
			res.push_back(*currentForce);

		for (list<ExtendedForce*>::const_iterator currentForce = myExtendedForcesList.begin();
		     currentForce != myExtendedForcesList.end();
		     ++currentForce)
			res.push_back(*currentForce);

		for (list<MetaForce*>::const_iterator currentForce = myMetaForcesList.begin();
		     currentForce != myMetaForcesList.end();
		     ++currentForce)
			res.push_back(*currentForce);
		return res;
	}

	vector<Force*> ForceGroup::getDeepMetaForces() const
	{
		vector<Force*> res;
		for (list<MetaForce*>::const_iterator currentForce = myMetaForcesList.begin();
		     currentForce != myMetaForcesList.end();
		     ++currentForce)
			(*currentForce)->getDeepForces(res);
		return res;
	}
}
