/* -*- c++ -*- */
#ifndef CONSTRAINTVALUETYPE_H
#define CONSTRAINTVALUETYPE_H

#include "AbstractEnumType.h"

namespace ProtoMol
{
	//_____________________________________________________ ConstraintValueEnum

	/// Map of Value constraints
	class ConstraintValueEnum
	{
	public:
		virtual ~ConstraintValueEnum()
		{
		}

		enum Enum
		{
			FIRST = 0, // Used internally only
			UNDEFINED = 0, // ConstraintValue returned when no string matches
			NOCONSTRAINTS,
			EMPTY,
			NOTEMPTY,
			ZERO,
			NOTZERO,
			POSITIVE,
			NEGATIVE,
			NOTPOSITIVE,
			NOTNEGATIVE,
			LAST // Used internally only
		};

	protected:
		static const std::string str[];

	public:
		template <Enum e>
		struct Enum2Type
		{
			operator Enum() const
			{
				return e;
			}

			//operator std::string() const;
			enum
			{
				value = e
			};
		};

	public:
		// Define types from the enum values
		typedef Enum2Type<NOCONSTRAINTS> NoConstraints;
		typedef Enum2Type<EMPTY> Empty;
		typedef Enum2Type<NOTEMPTY> NotEmpty;
		typedef Enum2Type<ZERO> Zero;
		typedef Enum2Type<NOTZERO> NotZero;
		typedef Enum2Type<POSITIVE> Positive;
		typedef Enum2Type<NEGATIVE> Negative;
		typedef Enum2Type<NOTPOSITIVE> NotPositive;
		typedef Enum2Type<NOTNEGATIVE> NotNegative;
	};

	//_____________________________________________________ ConstraintValueType
	typedef AbstractEnumType<ConstraintValueEnum> ConstraintValueType;

	//template<ConstraintValueEnum::Enum e>
	//inline ConstraintValueEnum::Enum2Type<e>::operator std::string() const{
	//  return ConstraintValueType::getString(e);
	//}
}
#endif /* CONSTRAINTVALUETYPE_H */
