#include "TRANSReader.h"

#include "stringutilities.h"
#include "Report.h"
#include <algorithm>

using std::string;
using std::vector;
using std::endl;
using std::find;
using std::stringstream;
using namespace ProtoMol::Report;

namespace ProtoMol {

  //_________________________________________________________________TRANSReader

  TRANSReader::TRANSReader():
    Reader(),
    myTRANS(NULL){}

  TRANSReader::TRANSReader(const std::string& filename):
    Reader(filename),
    myTRANS(NULL){}

  TRANSReader::TRANSReader(const char* filename):
    Reader(string(filename)),
    myTRANS(NULL){}

  TRANSReader::~TRANSReader(){
    if(myTRANS != NULL)
      delete myTRANS;
  }


  TRANS* TRANSReader::orphanTRANS(){
    TRANS* tmp = myTRANS;
    myTRANS = NULL;
    return tmp;
  }

  bool TRANSReader::tryFormat(){
    if (!open())
      return false;
    while(!myFile.eof()){
      string str;

      // read in the first line
      str = getline();
      stringstream ss(str);
      str.resize(find(str.begin(),str.end(),'!') - str.begin());
      str.resize(find(str.begin(),str.end(),'*') - str.begin());
      str.resize(find(str.begin(),str.end(),'{') - str.begin());
      str.resize(find(str.begin(),str.end(),'}') - str.begin());
      string i;
      ss >> i;
      if(i.size() > 0  && !equalStartNocase("END"  ,i)){
	close();
	return true;
      }
    }    
    myFile.setstate(std::ios::failbit);
    close();
    return false;
  }


  bool TRANSReader::read() {
    if(myTRANS == NULL)
      myTRANS = new TRANS();
    return read(*myTRANS);
  }

  bool TRANSReader::read(TRANS& trans) {
    // exit if the format is incorrect
    if(!tryFormat())
      return false;
    // exit if the file could not be opened
    if(!open())
      return false;

    // initialize the TRANS data object
    trans.clear();

    vector<vector<string> > input; // Stripped TRANS file input
    vector<string> signatures;     // 'd' for Real or int, 'w' for a word
    vector<int> lines;             // line number in TRANS file
    int count = 0;                 // TRANS file line counter
    
    int comment = 0;

    // make sure we are not at the end of the file
    while(!myFile.eof()){
      
      // create a data string
      vector<string> data;
      string numbers = "";

      // read in the next line
      string str(getline());
      ++count;

      // Remove {} comments
      if(find(str.begin(),str.end(),'}') == str.end() && find(str.begin(),str.end(),'{') == str.end()){
	if(comment > 0)
	  continue;
      }
      else {
	string tmp = "";
	for(unsigned int i=0;i<str.size();++i){
	  if(str[i] == '{')
	    ++comment;
	  else if(str[i] == '}'){
	    --comment;          
	    tmp += " ";
	    // Force to close comments if } found at the end of a line
	    if(comment > 0){
	      stringstream ss(str.substr(i,str.size()-i));
	      string rest;
	      ss >> rest;
	      if(rest == "}" && ss.eof()){
		comment = 0;
		report << warning << "The comments at line "<<count<<" in TRANS file are not properly closed."<<endr;
	      }
	    }
	    if(comment < 0){
	      report << warning << "The comments at line "<<count<<" in TRANS file has more closing \'}\'."<<endr;
	      comment = 0;
	    }
	  }
	  else if(comment == 0)
	    tmp += str[i];
	}
	str = tmp;
      }

      // Remove ! comments
      str.resize(find(str.begin(),str.end(),'!') - str.begin());

      stringstream ss(str);
      while(!ss.eof()){
	string i;
	ss >> i;
	// Skip * comments
	if(i[0] == '*')
	  break;
	if(i.size() > 0){        
	  data.push_back(i);
	  numbers += (isReal(i)?"d":"w");
	}
      }

      // Remove REMARK and SET
      if(data.size()>0 &&(equalStartNocase("REMARK",data[0]) || equalStartNocase("SET",data[0]))){
	data.clear();
	numbers = "";
      }

      // Store line
      if(data.size()>0){
	input.push_back(data);
	signatures.push_back(numbers);
	lines.push_back(count);

      }
    }


    // Parse
    TRANSRecordTypeEnum type = UNDEFINED;
    sized = got_alphaLJ = false;
    numTypes = 0;
    vector<string>::const_iterator sig  = signatures.begin();
    vector<int>::const_iterator    line = lines.begin();

    for(vector<vector<string> >::const_iterator i=input.begin();i != input.end();++i,++sig,++line){
      vector<string>::const_iterator j = i->begin();
      string s(*sig);

      if(equalStartNocase("iden",(*j))){
	type = IDENTITIES;
	continue;
      }
      else if(equalStartNocase("ideal",(*j))){
        type = IDEAL_GAS_DELTAMU;
        continue;                                                                     
      }
      else if(equalStartNocase("type",(*j))){
	type = ATOMTYPE;		 	    		     	  	     	      
	continue;		 	    		     	  	     	      
      }
      else if(equalStartNocase("stag",(*j))){
        type = STAGES;
        continue;
      }
      else if(equalStartNocase("charg",(*j))){
        type = ATOMCHARGE;
        continue;
      } 		 	    			  	        
      else if(equalStartNocase("end" ,(*j))){ 
	if(type == IDENTITIES || type == ATOMTYPE || type == ATOMCHARGE || type == IDEAL_GAS_DELTAMU) continue;
	report << warning <<"Unknown parameter \'"<<(*j)<<"\' in TRANS file at line "<<(*line)<<"."<<endr;
	continue;
      }
      
      
      // The definition as one string
      string definition;
      for(vector<string>::const_iterator l=j;l != i->end();++l){
	definition += (definition.size() > 0 ? " ":"") + (*l);
      }    

      // Add the type to the corresponding container ...
      switch(type){

      case IDENTITIES: {
	if(s == "d") 
	  numIdentities = toInt(j[0]);
	else
	  report << warning <<"Unknown # of identities definition \'"<<definition<<"\' ("<<s<<") in TRANS file at line "<<(*line)<<"."<<endr;
	break;
      }

      case STAGES: {
        if(s == "d") {
          numStages = toInt(j[0]);
          trans.NumberOfStages = numStages;}
        else
          report << warning <<"Unknown # of transformation stages definition \'"<<definition<<"\' ("<<s<<") in TRANS file at line "<<(*line)<<"."<<endr;
        break;
      }

      case IDEAL_GAS_DELTAMU: {
        // resize the chemical potential difference matrix if necessary
        if (!sized) {
          trans.DeltaMuIG.resize(ArraySizes(numIdentities)(numIdentities)(numStages));
          // zero out all the elements
          for (int o=0; o<numIdentities; o++) {
            for (int n=0; n<numIdentities; n++) {
              for (int s=0; s<numStages; s++) {
                trans.DeltaMuIG[o][n][s] = 0.0;
              }
            }
          }
          sized = true;
        }
        if(s == "dddd") { // stage(Uint)...oldType(Uint)...newType(Uint)...deltaMu(Real)

          // store the ideal gas state chemical potential difference for this stage
          trans.DeltaMuIG[toInt(j[1]) - 1][toInt(j[2]) - 1][toInt(j[0]) - 1] = toReal(j[3]);
        } // end if (s == "dddd") statement

        else
          report << warning <<"Unknown ideal gas deltaMu definition \'"<<definition<<"\' ("<<s<<") in TRANS file at line "<<(*line)<<"."<<endr;
        break;
      }
            
      case ATOMTYPE: {
	if(s == "wdddd") { // atomtype(string)...mass(Real)...charge(Real)...identity(int)...stage(Uint)
        
	  // loop over all previously stored atomtype names and see if this is a new atomtype
	  bool NewType = true;
	  unsigned int Index = myTypes.size();
	  for (unsigned int T=0; T<myTypes.size(); T++) {
	    if (myTypes[T] == j[0]) {
	      NewType = false;
	      Index = T;}
	  }

	  // if this is a new atomtype then store the atomtype name and number
	  if (NewType) {
	    myTypes.push_back(j[0]);
	    trans.atomTypes.push_back(TRANS::AtomType(numIdentities, numStages, myTypes[Index]));
	  }  // end if (NewType) statement

	  // now store the mass and stage # for this atomtype's particular identity
	  trans.atomTypes[Index].mass[toInt(j[3]) - 1] = toReal(j[1]);
	  trans.atomTypes[Index].charge[toInt(j[3]) - 1] = toReal(j[2]);
          trans.atomTypes[Index].stage[toInt(j[3]) - 1] = toInt(j[4]);
	} // end if (s == "wdddd") statement

	else
	  report << warning <<"Unknown mass & charge definition \'"<<definition<<"\' ("<<s<<") in TRANS file at line "<<(*line)<<"."<<endr;
	break;
      }

      case ATOMCHARGE: {

        if(s == "dwdddd") { // stage(int)...atomtype(string)...oldidentity(int)...
                            // newidentity(int)...oldcharge(Real)...newcharge(Real)
        
          // loop over all previously stored atomtype names and find the right one
          bool NewType = true;
          unsigned int Index = myTypes.size();
          for (unsigned int T=0; T<myTypes.size(); T++) {
            if (myTypes[T] == j[1]) {
              NewType = false;
              // we have found the right atomtype, set Index = T
              Index = T;}
          }

          // if this is a new atomtype then store the atomtype name and number
          if (NewType) {
            myTypes.push_back(j[1]);
            trans.atomTypes.push_back(TRANS::AtomType(numIdentities, numStages, myTypes[Index]));
          }  // end if (NewType) statement

          // now store this atomtype's charges for this particular stage
          int o = toInt(j[2]) - 1;
          int n = toInt(j[3]) - 1;
          int s = toInt(j[0]) - 1;
          trans.atomTypes[Index].old_charge[o][n][s] = toReal(j[4]);
          trans.atomTypes[Index].new_charge[o][n][s] = toReal(j[5]);
        } // end if (s == "dwdddd") statement

        else if(s == "dwddddd") { // stage(int)...atomtype(string)...oldidentity(int)...
                             // newidentity(int)...oldcharge(Real)...newcharge(Real)...alphaLJ(Real)

          // the user specified their own alphaLJs for this atomtype
          if (!got_alphaLJ) got_alphaLJ = true;
          
          // loop over all previously stored atomtype names and find the right one
          bool NewType = true;
          unsigned int Index = myTypes.size();
          for (unsigned int T=0; T<myTypes.size(); T++) {
            if (myTypes[T] == j[1]) {
              NewType = false;
              // we have found the right atomtype, set Index = T
              Index = T;}
          }

          // if this is a new atomtype then store the atomtype name and number
          if (NewType) {
            myTypes.push_back(j[1]);
            trans.atomTypes.push_back(TRANS::AtomType(numIdentities, numStages, myTypes[Index]));
          }  // end if (NewType) statement

          // now store this atomtype's charges for this particular stage and its alphaLJ for this
          // particular type of transformation
          int o = toInt(j[2]) - 1;
          int n = toInt(j[3]) - 1;
          int s = toInt(j[0]) - 1;
          trans.atomTypes[Index].old_charge[o][n][s] = toReal(j[4]);
          trans.atomTypes[Index].new_charge[o][n][s] = toReal(j[5]);
          trans.atomTypes[Index].alphaLJ[o][n]       = toReal(j[6]);
        } // end if (s == "dwddddd") statement
        
        else
          report << warning <<"Unknown set of charges definition \'"<<definition<<"\' ("<<s<<") in TRANS file at line "<<(*line)<<"."<<endr;
        break;
      }
      
      case UNDEFINED:
      default: {
	report << warning <<"Unknown definition \'"<<definition<<"\' ("<<s<<") in TRANS file at line "<<(*line)<<"."<<endr;
	break;
      }
      }    
    }
    close();
    if (!got_alphaLJ) report << warning << "Using alphaLJ = 0.5 for all atomtypes" << endr;
    return !myFile.fail();
  }


  TRANSReader& operator>>(TRANSReader& transReader, TRANS& trans){
    transReader.read(trans);
    return transReader;
  }

}
