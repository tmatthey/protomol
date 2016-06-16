/* -*- c++ -*- */
#ifndef FACTORY_H
#define FACTORY_H

#include "stringutilities.h"
//#include "systemutilities.h"
#include "Report.h"

#include <vector>
#include <map>
#include <set>

namespace ProtoMol
{
	//_____________________________________________________ FactoryBase
	/**
	   Base class of all factories templated with the family type.
	   Container to keep pointers for each prototype exemplar and their aliases, where
	   the real prototype is keep in a separate set. 
	*/
	template <typename Type>
	class FactoryBase
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		FactoryBase(): myCache(false)
		{
			//Report::report << Report::plain <<"[FactoryBase<"<<Type::scope<<">::FactoryBase()]"<< Report::endr;
		}

		virtual ~FactoryBase()
		{
			//Report::report << Report::plain <<"[FactoryBase<"<<Type::scope<<">::~FactoryBase()]"<< Report::endr;
			clear();
		}

	private:
		FactoryBase(const FactoryBase&)
		{
		}

		FactoryBase& operator=(const FactoryBase&)
		{
			return *this;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class FactoryBase
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		const Type* getPrototype(const std::string& id) const
		{
			const Type* prototype = NULL;

			if (myExemplars.find(id) != myExemplars.end())
			{
				prototype = myExemplars.find(id)->second;
			}
			else if (myAliasExemplars.find(id) != myAliasExemplars.end())
			{
				prototype = myAliasExemplars.find(id)->second;
			}

			return prototype;
		}

		void clear()
		{
			for (typename std::set<const Type*>::iterator i = myPointers.begin(); i != myPointers.end(); ++i)
				delete (*i);
			myExemplars.clear();
			//realclear(myExemplars);
			myAliasExemplars.clear();
			//realclear(myAliasExemplars);
			myPointers.clear();
			//realclear(myPointers);
			myCache = false;
		}

	protected:
		virtual std::string doPrint() const =0; ///< Hook method called from static method print 
		virtual void doRegisterHelpText() const =0; ///< Hook method called from static method registerHelpText 

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// private data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		std::map<std::string, const Type*, ltstrNocase> myExemplars;
		std::map<std::string, const Type*, ltstrNocase> myAliasExemplars;
		std::set<const Type*> myPointers;
		mutable bool myCache;
	};


	//_____________________________________________________ FactoryTraits<>
	template <class T>
	class FactoryTraits;

	//_____________________________________________________ Factory
	/** 
	    Concrete factory implementing the singleton pattern and 
	    inheriting the implementation details by traits. The concept of three
	    levels enables to combine the singleton pattern and implementation details, 
	    where the base class FactoryBase provides the essential methods and
	    data members for the implementation details. 
	*/

	template <typename Type>
	class Factory : public FactoryTraits<Type>::Details
	{
	public:
		typedef typename std::set<const Type*>::const_iterator const_iterator;
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	protected:
		Factory(): FactoryTraits<Type>::Details()
		{
			//Report::report << Report::plain <<"[Factory<"<<Type::scope<<">::Factory()]"<< Report::endr;
		}

		virtual ~Factory()
		{
			//Report::report << Report::plain <<"[Factory<"<<Type::scope<<">::~Factory()]"<< Report::endr;
		}

	private:
		Factory(const Factory&)
		{
		}

		Factory& operator=(const Factory&)
		{
			return *this;
		}

	private:
		/// Call by atexit() to clean up.
		static void kill()
		{
			Factory* p = obj;
			obj = NULL;
			p->~Factory();
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class Factory
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:

		static void registerExemplar(const Type* exemplar, std::string id = "")
		{
			if (id.empty())
				id = exemplar->getId();
			Factory::instance().doRegisterExemplar(id, exemplar);
		}

		static void registerExemplar(const Type* exemplar, const std::vector<std::string>& aliases, std::string id = "")
		{
			if (id.empty())
				id = exemplar->getId();
			Factory::instance().doRegisterExemplar(id, exemplar, aliases);
		}

		static bool unregisterExemplar(const std::string& id)
		{
			return Factory::instance().doUnregisterExemplar(id);
		}

		static void unregisterAllExemplars()
		{
			Factory::instance().clear();
		}

		static std::string print()
		{
			return Factory::instance().doPrint();
		}

		static void registerHelpText()
		{
			Factory::instance().doRegisterHelpText();
		}

		static bool empty()
		{
			return Factory::instance().myPointers.empty();
		}

		static Factory& instance()
		{
			//static Factory obj;
			//return obj;
			// We have to do it ourself ... M$ problem ...
			if (obj == NULL)
			{
				obj = new Factory();
				std::atexit(kill);
			}
			return *obj;
		}

		static const_iterator begin()
		{
			return Factory::instance().myPointers.begin();
		}

		static const_iterator end()
		{
			return Factory::instance().myPointers.end();
		}


		static const Type* find(const std::string& id)
		{
			return Factory::instance().getPrototype(id);
		}

	private:
		void doRegisterExemplar(const std::string& id, const Type* exemplar)
		{
			if (this->myExemplars.find(id) != this->myExemplars.end())
			{
				Report::report << Report::hint << "Prototype \'" << id << "\' already registered in " << Type::scope << "Factory" << Type::scope << ", overwriting." << Report::endr;
			}
			if (exemplar->getParameterSize() != exemplar->getParameters().size())
			{
				Report::report << Report::error << Type::scope << "Factory" << Type::scope << " prototype \'" << id << "\'  has different parameter size definitions " << exemplar->getParameterSize() << " != " << exemplar->getParameters().size() << ", fix it." << Report::endr;
			}
			this->myExemplars[id] = exemplar;
			this->myPointers.insert(exemplar);
			this->myCache = false;
		}

		void doRegisterExemplar(const std::string& id, const Type* exemplar, const std::vector<std::string>& aliases)
		{
			doRegisterExemplar(id, exemplar);
			for (unsigned int i = 0; i < aliases.size(); i++)
			{
				this->myAliasExemplars[aliases[i]] = exemplar;
			}
			this->myCache = false;
		}

		bool doUnregisterExemplar(const std::string& id)
		{
			// Get object pointer
			const Type* p = this->getPrototype(id);
			if (p == NULL)
				return false;

			// Remove pointers
			for (typename std::map<std::string, const Type*, ltstrNocase>::iterator i = this->myExemplars.begin(); i != this->myExemplars.end();)
			{
				if (i->second == p)
					this->myExemplars.erase(i++);
				else
					++i;
			}
			for (typename std::map<std::string, const Type*, ltstrNocase>::iterator i = this->myAliasExemplars.begin(); i != this->myAliasExemplars.end();)
			{
				if (i->second == p)
					this->myAliasExemplars.erase(i++);
				else
					++i;
			}
			this->myPointers.erase(p);

			// ... and delete the object
			delete p;

			this->myCache = false;
			return true;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// private data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		static Factory* obj;
	};

	template <typename Type>
	Factory<Type>* Factory<Type>::obj = NULL;
}
#endif /* FACTORY_H */
