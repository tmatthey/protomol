/*  -*- c++ -*-  */
#ifndef TYPESELECTION_H
#define TYPESELECTION_H

#include "Real.h"

namespace ProtoMol
{
	namespace Private
	{
		template <bool cmp, class A, class B>
		struct SelectTypeHelper
		{
			typedef B type;
		};

		template <class A, class B>
		struct SelectTypeHelper<true, A, B>
		{
			typedef A type;
		};

		template <unsigned int size, class A, class B>
		struct SelectType
		{
			typedef typename SelectTypeHelper<sizeof(A) == size, A, B>::type type;
		};

		template <bool cmp, class A>
		struct SelectTypeCheckHelper
		{
		};

		template <class A>
		struct SelectTypeCheckHelper<true, A>
		{
			typedef A type;
		};

		template <class A, unsigned int size>
		struct SelectTypeCheck
		{
			typedef typename SelectTypeCheckHelper<sizeof(A) == size, A>::type type;
		};
	}

	/**
	 * Enables to select the right type of an int or float 
	 * with a given sizeof, bails out if there is not adequate
	 * type. @n
	 * 
	 * Usage:@n
	 * typedef TypeSelection::Int<4>::type int32;@n
	 * typedef TypeSelection::Float<4>::type float32;@n 
	 */
	namespace TypeSelection
	{
		/**
		 * Select the right type among short, int, long or ong long according the given sizeof.
		 */
		template <unsigned int size>
		struct Int
		{
			typedef typename Private::SelectTypeCheck<typename Private::SelectType<size,
			                                                                       typename Private::SelectType<size, int, short>::type,
			                                                                       typename Private::SelectType<size, long, long long>::type>::type, size>::type type;
		};

		/**
		 * Select the right type among float or double according the given sizeof.
		 */
		template <unsigned int size>
		struct Float
		{
			typedef typename Private::SelectTypeCheck<typename Private::SelectType<size,
			                                                                       typename Private::SelectType<size, float, double>::type,
			                                                                       Real>::type, size>::type type;
		};
	}
}
#endif /* TYPESELECTION_H */ 
