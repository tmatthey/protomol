#include "systemutilities.h"

#include "stringutilities.h"

#ifdef WIN32

// Define the missing symbols from <unistd.h> for M$ ....
#include <direct.h>
#define CHDIR _chdir
#define PATHSEP '\\'
#define PATHSEPSTR "\\"
#define access _access
#include <fcntl.h>
#include <io.h>
#define F_OK 0
#define W_OK 2
#define R_OK 4

#else

#include <unistd.h>
#define CHDIR chdir
#define PATHSEP '/'
#define PATHSEPSTR "/"
#include <pwd.h>
#endif

#include <sys/stat.h>

using std::string;

namespace ProtoMol {

  //_____________________________________________________________________ changeDirectory
  bool changeDirectory(const string& fileName){
    char *confFile = (char*)fileName.c_str();
    char *currentdir = confFile;
    char *tmp = NULL;

#ifdef WIN32
    // Replace all '/' by '\'
    for(tmp=confFile;*tmp;++tmp){
      if(*tmp == '/')
	*tmp = '\\';
    }
#endif

    for(tmp=confFile;*tmp;++tmp); // find final null
    for(;tmp != confFile && *tmp != PATHSEP; --tmp); // find last '/'
    if (tmp != confFile ){
      *tmp = 0; 
      confFile = tmp + 1;
      if (CHDIR(currentdir))
	return false;
    }
    else if (*tmp == PATHSEP) // config file in / is odd, but it might happen
      if (CHDIR(PATHSEPSTR)){
	return false;
      }

    return true;
  }

  //_____________________________________________________________________ isAccessible
  bool isAccessible(const string& fileName){
    return (::access(fileName.c_str(), F_OK) == 0);
  }


  //_____________________________________________________________________ protomolAbort

  static void (*myAbortFunction)() = NULL;

  void protomolAbort(){
    if(myAbortFunction != NULL){
      (*myAbortFunction)();    
    }
    exit(EXIT_FAILURE);
  }

  //_____________________________________________________________________ setProtomolAbort
  void setProtomolAbort(void (*abortFunction)()){
    myAbortFunction = abortFunction;
  }

  //_____________________________________________________________________ protomolExit

  static void (*myExitFunction)() = NULL;

  void protomolExit(){
    if(myExitFunction != NULL){
      (*myExitFunction)();    
    }
    exit(EXIT_SUCCESS);
  }

  //_____________________________________________________________________ setProtomolExit
  void setProtomolExit(void (*exitFunction)()){
    myExitFunction = exitFunction;
  }


  //_____________________________________________________________________ protomolStartSerial
  static void (*myStartSerial)(bool) = NULL;

  void protomolStartSerial(bool exludeMaster){
    if(myStartSerial != NULL){
      (*myStartSerial)(exludeMaster);    
    }
  }

  //_____________________________________________________________________ setProtomolExit
  void setProtomolStartSerial(void (*startSerialFunction)(bool)){
    myStartSerial = startSerialFunction;
  }

  //_____________________________________________________________________ protomolEndSerial
  static void (*myEndSerial)(bool) = NULL;

  void protomolEndSerial(bool exludeMaster){
    if(myEndSerial != NULL){
      (*myEndSerial)(exludeMaster);    
    }
  }

  //_____________________________________________________________________ setProtomolExit
  void setProtomolEndSerial(void (*endSerialFunction)(bool)){
    myEndSerial = endSerialFunction;
  }

  //_____________________________________________________________________ ISLITTLEENDIAN
  struct Endian {
    // Helper class to make sure that we get endianess correct ... M$
    static bool isLittleEndian(){
      unsigned int tmp = 1;
      return (0 != *(reinterpret_cast<const char*>(&tmp)));
    }
  };
  const bool ISLITTLEENDIAN = Endian::isLittleEndian();

  //_____________________________________________________________________ getUserName
  string getUserName(){
#ifdef WIN32
    return "Win32";
#else
    if (getpwuid(getuid()) != NULL)
      return string(getpwuid(getuid())->pw_name);
    else 
      return toString(getuid());
#endif
  }
}
