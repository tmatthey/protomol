/* -*- c++ -*- */
#ifndef INPUTPOSVELTYPE_H
#define INPUTPOSVELTYPE_H

#include "AbstractEnumType.h"

namespace ProtoMol
{
	//_____________________________________________________ InputPosVelEnum

	class InputPosVelEnum
	{
	public:
		virtual ~InputPosVelEnum()
		{
		}

		enum Enum
		{
			FIRST = 0, // Used internally only
			UNDEFINED = 0, // InputPosVel returned when no string matches
			PDB,
			XYZ,
			XYZBIN,
			LAST // Used internally only
		};

	protected:
		static const std::string str[];
	};

	//_____________________________________________________ InputPosVelType
	typedef AbstractEnumType<InputPosVelEnum> InputPosVelType;
}
#endif /* INPUTPOSVELTYPE_H */
