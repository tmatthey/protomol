#include "STAGEReader.h"

#include "InputPosVel.h"
#include "PSFReader.h"
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

  //_________________________________________________________________STAGEReader

  STAGEReader::STAGEReader():
    Reader(),
    mySTAGE(NULL){}

  STAGEReader::STAGEReader(const std::string& filename):
    Reader(filename),
    mySTAGE(NULL){}

  STAGEReader::STAGEReader(const char* filename):
    Reader(string(filename)),
    mySTAGE(NULL){}

  STAGEReader::~STAGEReader(){
    if(mySTAGE != NULL)
      delete mySTAGE;
  }


  STAGE* STAGEReader::orphanSTAGE(){
    STAGE* tmp = mySTAGE;
    mySTAGE = NULL;
    return tmp;
  }

  bool STAGEReader::tryFormat(){
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


  bool STAGEReader::read() {
    if(mySTAGE == NULL)
      mySTAGE = new STAGE();
    return read(*mySTAGE);
  }

  bool STAGEReader::read(STAGE& stage) {
    // exit if the format is incorrect
    if(!tryFormat())
      return false;
    // exit if the file could not be opened
    if(!open())
      return false;

    // initialize the STAGE data object
    stage.clear();

    vector<vector<string> > input; // Stripped STAGE file input
    vector<string> signatures;     // 'd' for Real or int, 'w' for a word
    vector<int> lines;             // line number in STAGE file
    int count = 0;                 // STAGE file line counter

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
		report << warning << "The comments at line "<<count<<" in STAGE file are not properly closed."<<endr;
	      }
	    }
	    if(comment < 0){
	      report << warning << "The comments at line "<<count<<" in STAGE file has more closing \'}\'."<<endr;
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
    STAGERecordTypeEnum type = UNDEFINED;
    got_alphaLJ = false;
    numOsmotics = 0;
    vector<string>::const_iterator sig  = signatures.begin();
    vector<int>::const_iterator    line = lines.begin();

    for(vector<vector<string> >::const_iterator i=input.begin();i != input.end();++i,++sig,++line){
      vector<string>::const_iterator j = i->begin();
      string s(*sig);

      if(equalStartNocase("osmot",(*j))){
	type = OSMOTICS;
	continue;
      }
      else if(equalStartNocase("type",(*j))){
	type = ATOMTYPE;
	continue;
      }
      else if(equalStartNocase("charg",(*j))){
        type = ATOMCHARGE;
        continue;
      }
      else if(equalStartNocase("end" ,(*j))){
	if(type == ATOMTYPE || type == ATOMCHARGE || type == OSMOTICS) continue;
	report << warning <<"Unknown parameter \'"<<(*j)<<"\' in STAGE file at line "<<(*line)<<"."<<endr;
	continue;
      }


      // The definition as one string
      string definition;
      for(vector<string>::const_iterator l=j;l != i->end();++l){
	definition += (definition.size() > 0 ? " ":"") + (*l);
      }

      // Add the type to the corresponding container ...
      switch(type){

      case OSMOTICS: {
	if(s == "wdddww") { // molecule name(string)...target fugacity(Real)...# of stages(Uint)...component ID#(Uint)
                            //...PSFname(string)...PDBname(string)
          numOsmotics++;
          Real tempMu = toReal(j[1]);
          int tempStgs = toInt(j[2]);
          if (tempMu < 0) report << error << "Negative fugacity specified for " << j[0] << " ("
                                 << tempMu << ") in STAGE file at line " << (*line) << "." << endr;
          if (tempStgs < 1) report << error << "Invalid # of transformation stages specified for " << j[0]
                                   << "( " << tempStgs << ") in STAGE file at line " << (*line) << "." << endr;
          stage.components.push_back(STAGE::osmoticType(j[0], tempMu, toUInt(j[2])));
          stage.components[numOsmotics-1].myMolecule.type = toUInt(j[3]) - 1;
          stage.components[numOsmotics-1].myMolecule.newtype = toUInt(j[3]) - 1;
          //report << plain << "Molec. name: " << stage.components[numOsmotics-1].myMolecule.name << endr;
          PSFnames.push_back(j[4]);
          Coordnames.push_back(j[5]);
        }
	else
	  report << warning <<"Unknown osmotic molecule definition \'"<<definition<<"\' ("<<s<<") in STAGE file at line "<<(*line)<<"."<<endr;
	break;
      }

      case ATOMTYPE: {
	if(s == "wwddd") { // moleculetype(string)...atomtype(string)...mass(Real)...charge(Real)...stage(Uint)

          // identify the osmotic molecule type that this atomtype belongs to
          bool FoundMolecule = false;
          unsigned int moleculeIndex = 0;
          for (unsigned int M=0; M<stage.components.size(); M++) {
            if (j[0] == stage.components[M].myMolecule.name) {
              FoundMolecule = true;
              moleculeIndex = M;}
          }

          if (!FoundMolecule) report << error << "Unknown osmotic molecule name \'" << j[0]
                                     << "\' in STAGE file at line " << (*line) << "." << endr;

	  // loop over all previously stored atomtype names and see if this is a new atomtype
	  bool NewType = true;
	  for (unsigned int T=0; T<stage.components[moleculeIndex].atomTypes.size(); T++) {
	    if (stage.components[moleculeIndex].atomTypes[T].type_name == j[1]) {
	      NewType = false;}
	  }

	  // report an error if this is a duplicate atom type
	  if (!NewType)
            report << error << "Duplicate atomtype found in STAGE file at line" << (*line) << "." << endr;

          // now store the mass, charge, and stage # for this atomtype
          stage.components[moleculeIndex].atomTypes.push_back(STAGE::AtomType(j[1], toReal(j[2]),
                                                             toReal(j[3]), toUInt(j[4]),
                                                             stage.components[moleculeIndex].NumberOfStages));
          stage.components[moleculeIndex].myMolecule.mass += toReal(j[2]);
          stage.components[moleculeIndex].myMolecule.atoms.push_back(0);

        } // end if (s == "wdddd") statement

        else if(s == "wwdddd") { // moleculetype(string)...atomtype(string)...mass(Real)...charge(Real)
                                 //...stage(int)...alphaLJ(Real)

          // the user specified their own alphaLJs for this atomtype
          if (!got_alphaLJ) got_alphaLJ = true;

          if (toReal(j[5]) < 0) report << error << "Negative alphaLJ specified for " << j[1] << " ("
                                       << toReal(j[5]) << ") in STAGE file at line " << (*line) << "." << endr;

          // identify the osmotic molecule type that this atomtype belongs to
          bool FoundMolecule = false;
          unsigned int moleculeIndex = 0;
          for (unsigned int M=0; M<stage.components.size(); M++) {
            if (j[0] == stage.components[M].myMolecule.name) {
              FoundMolecule = true;
              moleculeIndex = M;}
          }

          if (!FoundMolecule) report << error << "Unknown osmotic molecule name \'" << j[0]
                                     << "\' in STAGE file at line " << (*line) << "." << endr;

	  // loop over all previously stored atomtype names and see if this is a new atomtype
	  bool NewType = true;
	  for (unsigned int T=0; T<stage.components[moleculeIndex].atomTypes.size(); T++) {
	    if (stage.components[moleculeIndex].atomTypes[T].type_name == j[1]) {
	      NewType = false;}
	  }

          // report and error if this is a duplicate atom type
	  if (!NewType)
           report << error << "Duplicate atomtype found in STAGE file at line" << (*line) << "." << endr;

	  // now store the mass, charge, and stage # for this atomtype
	  stage.components[moleculeIndex].atomTypes.push_back(STAGE::AtomType(j[1], toReal(j[2]),
                                                              toReal(j[3]), toUInt(j[4]), toReal(j[5]),
                                                              stage.components[moleculeIndex].NumberOfStages));
          stage.components[moleculeIndex].myMolecule.mass += toReal(j[2]);
          stage.components[moleculeIndex].myMolecule.atoms.push_back(0);

	} // end if (s == "wwdddd") statement

	else
	  report << warning <<"Unknown mass/charge definition \'"<<definition<<"\' ("<<s<<") in STAGE file at line "<<(*line)<<"."<<endr;
	break;
      }

      case ATOMCHARGE: {

        if(s == "dwwdd") { // stage(int)...moleculename(string)...atomtype(string)...oldcharge(Real)...newcharge(Real)

          // identify the osmotic molecule type that this atomtype belongs to
          bool FoundMolecule = false;
          unsigned int moleculeIndex = 0;
          for (unsigned int M=0; M<stage.components.size(); M++) {
            if (j[1] == stage.components[M].myMolecule.name) {
              FoundMolecule = true;
              moleculeIndex = M;}
          }

          if (!FoundMolecule) report << error << "Unknown osmotic molecule name \'" << j[1]
                                     << "\' in STAGE file at line " << (*line) << "." << endr;

          // loop over the atomtype names and find the right one
	  bool NewType = true;
          unsigned int Index = stage.components[moleculeIndex].atomTypes.size();
	  for (unsigned int T=0; T<stage.components[moleculeIndex].atomTypes.size(); T++) {
	    if (stage.components[moleculeIndex].atomTypes[T].type_name == j[2]) {
	      NewType = false;
              // we have found the right atomtype, set Index = T
              Index = T;}
	  } // end loop over T

          if (NewType) report << error << "Atom type found in CHARGES section in STAGE file at line " << (*line)
                              << " that was not specified in the TYPES section (" << j[2] << "). " << endr;

          // now store this atomtype's charges for this particular stage
          unsigned int s = toUInt(j[0]) - 1;
          if (s < 0) report << error << "Invalid transformation stage # specified for " << j[2]
                            << "( " << s << ") in STAGE file at line " << (*line) << "." << endr;
          stage.components[moleculeIndex].atomTypes[Index].delete_charge[s] = toReal(j[3]);
          stage.components[moleculeIndex].atomTypes[Index].insert_charge[s] = toReal(j[4]);
        } // end if (s == "dwddd") statement

        else
          report << warning <<"Unknown set of charges definition \'"<<definition<<"\' ("<<s<<") in STAGE file at line "<<(*line)<<"."<<endr;
        break;
      }

      case UNDEFINED:
      default: {
	report << warning <<"Unknown definition \'"<<definition<<"\' ("<<s<<") in STAGE file at line "<<(*line)<<"."<<endr;
	break;
      }
      }
    }
    close();
    if (!got_alphaLJ) report << warning << "Using alphaLJ = 0.5 for all atomtypes" << endr;

    // next, read in the structure and coordinates for each osmotic component
    for (unsigned int i=0; i<numOsmotics; i++) {
      // start with the coordinates
      InputPosVel coordsreader;

      // open the positions file
      if(!coordsreader.open(Coordnames[i]))
        report << error << "Can't open position file \'" << Coordnames[i] << "\'." << endr;

      // read in the positions file information
      if(coordsreader.tryFormat(InputPosVelType::PDB)){
        PDB pdb;
        if(!(coordsreader >> pdb))
          report << error << "Could not parse PDB position file \'" << Coordnames[i] << "\'." << endr;
        // store the xyz coordinates
        swap(stage.components[i].myCoordinates,pdb.coords);
      }
      else if(!(coordsreader >> stage.components[i].myCoordinates))
        report << error << "Could not parse position file \'" << Coordnames[i]
               << "\'. Supported formats are : "<< InputPosVelType::getPossibleValues(", ") << "." << endr;
      report << plain << "Using "<< coordsreader.getType() << " posfile \'" << Coordnames[i]
             << "\' (" << stage.components[i].myCoordinates.size() << ")." << endr;

      // create the PSF reader object
      PSFReader psfReader;
      if(!psfReader.open(PSFnames[i]))
        report << error << "Can't open PSF file \'"<< PSFnames[i] << "\'." << endr;

      // read in the PSF file, which will contain the current mass and charge of each atom
      if(!(psfReader >> stage.components[i].myStructure))
        report << error << "Could not parse PSF file \'" << PSFnames[i] << "\'." << endr;
    } // end loop over osmotic components
          //report << plain << "# of osmotics: " << stage.components.size() << endr;

    return !myFile.fail();
  }


  STAGEReader& operator>>(STAGEReader& stageReader, STAGE& stage){
    stageReader.read(stage);
    return stageReader;
  }

}
