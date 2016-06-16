/*  -*- c++ -*-  */
#ifndef INPUTPOSVEL_H
#define INPUTPOSVEL_H

#include "InputPosVelType.h"
#include "XYZ.h"
#include "PDB.h"

namespace ProtoMol
{
	class Configuration;

	//________________________________________________________ InputPosVel
	class InputPosVel
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors (both default here), assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		InputPosVel();
		explicit InputPosVel(const std::string& filename);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class InputPosVel
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		void setFilename(const std::string& filename);

		bool open();
		bool open(const std::string& filename);
		bool tryFormat(InputPosVelType::Enum type);

		operator void*() const;
		bool operator!() const;
		// enable expression testing

		InputPosVelType getType() const;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Friends
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		friend InputPosVel& operator>>(InputPosVel& posReader, PDB& pdb);
		friend InputPosVel& operator>>(InputPosVel& posReader, XYZ& xyz);
		friend InputPosVel& operator>>(InputPosVel& posReader, Vector3DBlock& coords);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		std::string myFilename;
		bool myOk;
		InputPosVelType myType;
	};
}
#endif
