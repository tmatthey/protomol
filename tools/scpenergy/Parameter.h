/*  -*- c++ -*-  */
#ifndef PARAMETER_H
#define PARAMETER_H

#include "Value.h"
#include "simpleTypes.h"

namespace ProtoMol {
  //________________________________________________________ Parameter
  struct Parameter {
    /**
     * Container struct for parameters providing wide range
     * of constructors.
     */

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Parameter();

    Parameter(const std::string& k, const Value& val);
    Parameter(const std::string& k, const Value& val, const Value& def);
    template<typename T>
    Parameter(const std::string& k, const Value& val, T def):keyword(k),value(val),defaultValue(val) {defaultValue.set(def);}
    Parameter(const char* k, const Value& val);
    Parameter(const char* k, const Value& val, const Value& def);
    template<typename T>
    Parameter(const char* k, const Value& val, T def):keyword(std::string(k)),value(val),defaultValue(val) {defaultValue.set(def);}

    Parameter(const std::string& k, const Value& val, const Text& t);
    Parameter(const std::string& k, const Value& val, const Value& def, const Text& t);
    template<typename T>
    Parameter(const std::string& k, const Value& val, T def, const Text& t):keyword(k),value(val),defaultValue(val),text(t.text) {defaultValue.set(def);}
    Parameter(const char* k, const Value& val, const Text& t);
    Parameter(const char* k, const Value& val, const Value& def, const Text& t);
    template<typename T>
    Parameter(const char* k, const Value& val, T def, const Text& t):keyword(std::string(k)),value(val),defaultValue(val),text(t.text) {defaultValue.set(def);}
//      template<typename T>
//      Parameter(const std::string& k, T val):keyword(k),value(val),defaultValue(val,Value::undefined) {}
//      template<typename T>
//      Parameter(const std::string& k, T val, T def):keyword(k),value(val),defaultValue(def) {}
//      template<typename T, typename C>
//      Parameter(const std::string& k, T val, const C& con):keyword(k),value(val,con),defaultValue(val,con,Value::undefined) {}
//      template<typename T, typename C>
//      Parameter(const std::string& k, T val, T def, const C& con):keyword(k),value(val,con),defaultValue(def,con) {}
  
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// the keyword of the parameter
    std::string keyword;
    /// the value of the parameter
    Value       value;
    /// optional default value of the parameter
    Value       defaultValue;
    /// optional help text of the parameter
    std::string text;
  };

}
#endif /* PARAMETER_H */
