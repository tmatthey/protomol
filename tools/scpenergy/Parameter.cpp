#include "Parameter.h"

using std::vector;
using std::string;

namespace ProtoMol {
  //________________________________________________________ Parameter
  Parameter::Parameter():keyword(""),value(Value()),defaultValue(Value()) {}
  Parameter::Parameter(const string& k, const Value& val):keyword(k),value(val),defaultValue(val,Value::undefined) {}
  Parameter::Parameter(const string& k, const Value& val, const Value& def):keyword(k),value(val),defaultValue(def) {}
  Parameter::Parameter(const char* k, const Value& val):keyword(string(k)),value(val),defaultValue(val,Value::undefined) {}
  Parameter::Parameter(const char* k, const Value& val, const Value& def):keyword(string(k)),value(val),defaultValue(def) {}

  Parameter::Parameter(const string& k, const Value& val, const Text& t):keyword(k),value(val),defaultValue(val,Value::undefined),text(t.text) {}
  Parameter::Parameter(const string& k, const Value& val, const Value& def, const Text& t):keyword(k),value(val),defaultValue(def),text(t.text) {}
  Parameter::Parameter(const char* k, const Value& val, const Text& t):keyword(string(k)),value(val),defaultValue(val,Value::undefined),text(t.text) {}
  Parameter::Parameter(const char* k, const Value& val, const Value& def, const Text& t):keyword(string(k)),value(val),defaultValue(def),text(t.text) {}  
}
