/*  -*- c++ -*-  */
#ifndef DIHEDRALSYSTEMFORCEBASE_H
#define DIHEDRALSYSTEMFORCEBASE_H

#include<string>

namespace ProtoMol
{
	//_________________________________________________________________ DihedralSystemForceBase

	class DihedralSystemForceBase
	{
	public:
		virtual ~DihedralSystemForceBase()
		{
		}

		static const std::string keyword;
	};
}
#endif /* DIHEDRALSYSTEMFORCEBASE_H */
