/*  -*- c++ -*-  */
#ifndef UNIVERSALSWITCHINGFUNCTIONBASE_H
#define UNIVERSALSWITCHINGFUNCTIONBASE_H

#include <string>

namespace ProtoMol
{
	//_________________________________________________________________ UniversalSwitchingFunctionBase

	class UniversalSwitchingFunctionBase
	{
	public:
		virtual ~UniversalSwitchingFunctionBase()
		{
		}

		static const std::string keyword;
	};
}
#endif /* COMPLEMENTSWITCHINGFUNCTIONBASE_H */
