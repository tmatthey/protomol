#include "iSGPARReader.h"

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

  //_________________________________________________________________PARReader
  iSGPARReader::iSGPARReader(iSGPAR::CharmmTypeEnum charmmType):
    Reader(),
    myPAR(NULL),myCharmmType(charmmType),myCharmmTypeDetected(iSGPAR::UNDEFINED){}

  iSGPARReader::iSGPARReader(const std::string& filename, iSGPAR::CharmmTypeEnum charmmType):
    Reader(filename),
    myPAR(NULL),myCharmmType(charmmType),myCharmmTypeDetected(iSGPAR::UNDEFINED){}

  iSGPARReader::iSGPARReader(const char* filename, iSGPAR::CharmmTypeEnum charmmType):
    Reader(string(filename)),
    myPAR(NULL),myCharmmType(charmmType),myCharmmTypeDetected(iSGPAR::UNDEFINED){}

  iSGPARReader::~iSGPARReader(){
    if(myPAR != NULL)
      delete myPAR;
  }

  bool iSGPARReader::open(const std::string& filename, iSGPAR::CharmmTypeEnum charmmType){
    myCharmmType = charmmType;
    myCharmmTypeDetected = iSGPAR::UNDEFINED;
    return open(filename);
  }

  bool iSGPARReader::open(iSGPAR::CharmmTypeEnum charmmType){
    myCharmmType = charmmType;
    myCharmmTypeDetected = iSGPAR::UNDEFINED;
    return open();
  }

  iSGPAR* iSGPARReader::orphanPAR(){
    iSGPAR* tmp = myPAR;
    myPAR = NULL;
    return tmp;
  }

  bool iSGPARReader::tryFormat(){
    if (!open())
      return false;
    while(!myFile.eof()){
      string str;
      str = getline();
      stringstream ss(str);
      str.resize(find(str.begin(),str.end(),'!') - str.begin());
      str.resize(find(str.begin(),str.end(),'*') - str.begin());
      str.resize(find(str.begin(),str.end(),'{') - str.begin());
      str.resize(find(str.begin(),str.end(),'}') - str.begin());
      string i;
      ss >> i;
      if(!i.empty() && (isKeywordCharmm19(i) || isKeywordCharmm28(i)) && !equalStartNocase("END"  ,i)){
	close();
	return true;
      }
    }    
    myFile.setstate(std::ios::failbit);
    close();
    return false;
  }


  bool iSGPARReader::read() {
    if(myPAR == NULL)
      myPAR = new iSGPAR();
    return read(*myPAR);
  }

  bool iSGPARReader::read(iSGPAR& par) {
    if(!tryFormat())
      return false;
    if(!open())
      return false;
    par.clear();

    vector<vector<string> > input; // Stripped PAR file input 
    vector<string> signatures;     // 'd' for Real or int, 'w' for a word
    vector<int> lines;             // line number in PAR file
    int charmm19 = 0;
    int charmm28 = 0;
    int count = 0;                 // PAR file line counter
    
    int comment = 0;
    while(!myFile.eof()){
      
      vector<string> data;
      string numbers = "";
      string str(removeBeginEndBlanks(getline()));
      ++count;

      // Remove {} comments
      if(!str.empty() && str[0] == '*')
	continue;

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
	    // Froce to close comments if } found at the end of a line
	    if(comment > 0){
	      stringstream ss(str.substr(i,str.size()-i));
	      string rest;
	      ss >> rest;
	      if(rest == "}" && ss.eof()){
		comment = 0;
		report << warning << "The comments at line "<<count<<" in PAR file are not properly closed."<<endr;
	      }
	    }
	    if(comment < 0){
	      report << warning << "The comments at line "<<count<<" in PAR file has more closing \'}\'."<<endr;
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
	if(!i.empty()){        
	  data.push_back(i);
	  numbers += (isReal(i)?"d":"w");
	}
      }

      // Remove REMARK and SET
      if(!data.empty() &&(equalStartNocase("REMARK",data[0]) || equalStartNocase("SET",data[0]))){
	data.clear();
	numbers = "";
      }

      // Store line
      if(!data.empty()){
	input.push_back(data);
	signatures.push_back(numbers);
	lines.push_back(count);

	// Statistics for format selection
	if(numbers.size() == 1 && isKeywordCharmm28(data[0])){
	  ++charmm28;
	}
	else if(numbers.size() > 1 && isKeywordCharmm19(data[0]) || find(numbers.begin(),numbers.end(),'w') == numbers.end()){
	  ++charmm19;
	}
	else if(numbers.size() > 1  && !isKeywordCharmm19(data[0]) && !isKeywordCharmm28(data[0])){
	  ++charmm28;
	}
      }
    }

    // Selection of format
    myCharmmTypeDetected = (charmm19 > charmm28?iSGPAR::CHARMM19:iSGPAR::CHARMM28);
    bool isCharmm19 = (myCharmmType == iSGPAR::UNDEFINED ?charmm19 > charmm28:(myCharmmType == iSGPAR::CHARMM19));

    // Parse
    PARRecordTypeEnum type = UNDEFINED;
    vector<string>::const_iterator sig  = signatures.begin();
    vector<int>::const_iterator    line = lines.begin();

    for(vector<vector<string> >::const_iterator i=input.begin();i != input.end();++i,++sig,++line){
      vector<string>::const_iterator j = i->begin();
      string s(*sig);
      if(isCharmm19){ // Charmm19
	if(equalStartNocase("AEXP" ,(*j)) ||
	   equalStartNocase("REXP" ,(*j))||
	   equalStartNocase("HAEX" ,(*j))||
	   equalStartNocase("AAEX" ,(*j))||
	   equalStartNocase("NBOND",(*j))||
	   equalStartNocase("CUTNB",(*j))||
	   equalStartNocase("END"  ,(*j))||
	   equalStartNocase("CTONN",(*j))||
	   equalStartNocase("EPS"  ,(*j))||
	   equalStartNocase("VSWI" ,(*j))||
	   equalStartNocase("NBXM" ,(*j))||
	   equalStartNocase("INHI" ,(*j)))
	  continue;
	if(equalStartNocase("bond",(*j))){
	  type = BOND;
	}
	else if(equalStartNocase("angl",(*j))){
	  type = ANGLE;
	}
	else if(equalStartNocase("dihe",(*j))){
	  type = DIHEDRAL;
	}
	else if(equalStartNocase("impr",(*j))){
	  type = IMPROPER;
	}
	else if(equalStartNocase("nonb",(*j))){
	  type = NONBONDED;
	}
	else if(equalStartNocase("nbfi",(*j))){
	  type = NBFIX;
	}
	else if(equalStartNocase("hbon",(*j))){
	  type = HBOND;
	}
	else {
	  report << warning <<"Unknown parameter \'"<<(*j)<<"\' in X-Plor parameter file at line "<<(*line)<<"."<<endr;
	  continue;
	}
	// Increment array pointer and remove first signature since the first entry is a Charmm19 keyword
	++j;
	s = s.substr(1,s.size()-1);      
      }
      else { // Charmm28
	if(equalStartNocase("bond",(*j))){
	  type = BOND;
	  continue;
	}
	else if(equalStartNocase("angl",(*j)) || equalStartNocase("thet",(*j))){
	  type = ANGLE;		 	    		     	  	     	      
	  continue;		 	    		     	  	     	      
	}			 	    			  	        
	else if(equalStartNocase("dihe",(*j)) || equalStartNocase("phi" ,(*j))){
	  type = DIHEDRAL;	 	    	     		  	     		      
	  continue;		 	    		     	  	     	      
	}			 	    			  	        
	else if(equalStartNocase("impr",(*j)) || equalStartNocase("imph",(*j))){
	  type = IMPROPER;	 	    	     		  	     		      
	  continue;		 	    		     	  	     	      
	}			 	    			  	        
	else if(equalStartNocase("nonb",(*j)) || equalStartNocase("nbon",(*j))){
	  type = NONBONDED;	 	    	     
	  continue;		 	    		     
	}			 	    			     
	else if(equalStartNocase("nbfi",(*j))){
	  type = NBFIX;		 	    		     
	  continue;		 	    		     
	}			 	    			     
	else if(equalStartNocase("hbon",(*j))){
	  type = HBOND;
	  continue;
	}    
	else if(equalStartNocase("end" ,(*j)) || equalStartNocase("nbxm",(*j)) ||
		equalStartNocase("grou",(*j)) || equalStartNocase("cdie",(*j)) ||
		equalStartNocase("shif",(*j)) || equalStartNocase("vgro",(*j)) ||
		equalStartNocase("vdis",(*j)) || equalStartNocase("vswi",(*j)) ||
		equalStartNocase("cutn",(*j)) || equalStartNocase("ctof",(*j)) ||
		equalStartNocase("cton",(*j)) || equalStartNocase("eps" ,(*j)) || 
		equalStartNocase("e14f",(*j)) || equalStartNocase("wmin",(*j)) ||
		equalStartNocase("aexp",(*j)) || equalStartNocase("rexp",(*j)) ||
		equalStartNocase("haex",(*j)) || equalStartNocase("aaex",(*j)) ||
		equalStartNocase("noac",(*j)) || equalStartNocase("hbno",(*j)) ||
		equalStartNocase("cuth",(*j)) || equalStartNocase("ctof",(*j)) ||
		equalStartNocase("cton",(*j)) || equalStartNocase("cuth",(*j)) ||
		equalStartNocase("ctof",(*j)) || equalStartNocase("cton",(*j))){ 
	  if(type == NONBONDED || type == NBFIX || type == HBOND)
	    continue;
	  report << warning <<"Unknown parameter \'"<<(*j)<<"\' in Charmm28 parameter file at line "<<(*line)<<"."<<endr;
	  continue;
	}
      }
      
      // The definition as one string
      string definition;
      for(vector<string>::const_iterator l=j;l != i->end();++l){
	definition += (definition.empty() ? "":" ") + (*l);
      }    
      // Add the type to the corresponding container ...
      switch(type){

      case BOND: {
	// if (s == "wwdd") // this is not a parameter for an iSG atom type, so ignore it       
	if (s == "wwddd") {
	  string atom1(j[0]);
	  string atom2(j[1]);
	  int l = par.bonds.size()-1;
	  if ( (l < 0 ||
		((par.bonds[l].atom1 != atom1 || par.bonds[l].atom2 != atom2) && 
		 (par.bonds[l].atom1 != atom2 || par.bonds[l].atom2 != atom1)))) {

	    // this is a bond type that we have not seen yet
	    par.bonds.push_back(iSGPAR::Bond(numComp));
	    l = par.bonds.size() - 1;
	    par.bonds[l].number = par.bonds.size();
	    par.bonds[l].atom1 = atom1;
	    par.bonds[l].atom2 = atom2;
	  }  // end if NOT new bond type statement

	  par.bonds[l].forceConstant[toInt(j[4])] = toReal(j[2]);
	  par.bonds[l].distance[toInt(j[4])] = toReal(j[3]);

	}  // end if (s == "wwddd") statement
	else if (s != "wwdd" ) {
	  report << warning <<"Unknown bond definition \'"<<definition<<"\' ("<<s<<") in PAR file at line "<<(*line)<<"."<<endr;
	}
	break; 
      }  // end case BOND statement

      case ANGLE: {
	if ( (s != "wwwdd") && (s != "wwwdddd") && (s != "wwwddwdd") &&
	     (s != "wwwddd") && (s != "wwwddddd") && (s != "wwwddwddd") )
	  report << warning <<"Unknown angle definition \'"<<definition<<"\' ("<<s<<") in PAR file at line "<<(*line)<<"."<<endr;
	else if ( (s == "wwwddd") || (s == "wwwddddd") || (s == "wwwddwddd") ) {
	  string atom1(j[0]);
	  string atom2(j[1]);
	  string atom3(j[2]);
	  int l = par.angles.size() - 1;

	  if ( (l < 0 ||
		((par.angles[l].atom1 != atom1 || par.angles[l].atom2 != atom2 || par.angles[l].atom3 != atom3) && 
		 (par.angles[l].atom1 != atom3 || par.angles[l].atom2 != atom2 || par.angles[l].atom3 != atom1)))) {

	    // this is an angle type that we have not seen yet
	    par.angles.push_back(iSGPAR::Angle(numComp));
	    l = par.angles.size() - 1;
	    par.angles[l].number = par.angles.size();
	    par.angles[l].atom1 = atom1;
	    par.angles[l].atom2 = atom2;
	    par.angles[l].atom3 = atom3;
	  }  // end if NOT new angle type statement

	  if (s == "wwwddd") {
	    par.angles[l].ub_flag = false;
	    par.angles[l].forceConstant[toInt(j[5])] = toReal(j[3]);
	    par.angles[l].angleval[toInt(j[5])] = toReal(j[4]);
	    par.angles[l].k_ub[toInt(j[5])] = 0.0;
	    par.angles[l].r_ub[toInt(j[5])] = 0.0;}
	  else if(s == "wwwddddd") {
	    par.angles[l].ub_flag = true;
	    par.angles[l].forceConstant[toInt(j[5])] = toReal(j[3]);
	    par.angles[l].angleval[toInt(j[5])] = toReal(j[4]);
	    par.angles[l].k_ub[toInt(j[5])] = toReal(j[5]);
	    par.angles[l].r_ub[toInt(j[5])] = toReal(j[6]);}
	  else if(s == "wwwddwdd") {
	    par.angles[l].ub_flag = true;
	    par.angles[l].forceConstant[toInt(j[5])] = toReal(j[3]);
	    par.angles[l].angleval[toInt(j[5])] = toReal(j[4]);
	    par.angles[l].k_ub[toInt(j[5])] = toReal(j[6]);
	    par.angles[l].r_ub[toInt(j[5])] = toReal(j[7]);
	  }  // end inner else-if statements
	}  // end outer else-if statement
	break;  
      }  // end case ANGLE statement

      case DIHEDRAL: {
	if ( (s != "wwwwddd") && (s != "wwwwwdddd") &&
	     (s != "wwwwdddd") && (s != "wwwwwddddd") )
	  report << warning <<"Unknown dihedral definition \'"<<definition<<"\' ("<<s<<") in PAR file at line "<<(*line)<<"."<<endr;
	else if(s == "wwwwdddd"){
	  // Charmm19 with multiplicity == 1
	  // Charmm28 with multiplicity >= 1
	  string atom1(j[0]);
	  string atom2(j[1]);
	  string atom3(j[2]);
	  string atom4(j[3]);
	  int l = par.dihedrals.size()-1;
	  int Identity = toInt(j[7]);
	  if (l >= 0 && !isCharmm19 &&
	      ((par.dihedrals[l].atom1 == atom1 && par.dihedrals[l].atom2 == atom2 &&
		par.dihedrals[l].atom3 == atom3 && par.dihedrals[l].atom4 == atom4) || 
	       (par.dihedrals[l].atom1 == atom4 && par.dihedrals[l].atom2 == atom3 &&
		par.dihedrals[l].atom3 == atom2 && par.dihedrals[l].atom4 == atom1))) {

	    // this is a dihedral type that we have already seen
	    // Charmm28 and previous is the same, increase multiplicity
	    par.dihedrals[l].multiplicity[Identity]++;
	    par.dihedrals[l].forceConstant[Identity].push_back(toReal(j[4]));
	    par.dihedrals[l].periodicity[Identity].push_back(toInt(j[5]));                                                        
	    par.dihedrals[l].phaseShift[Identity].push_back(toReal(j[6]));}
	  else {
	    // this is a dihedral type that we have not seen yet
	    par.dihedrals.push_back(iSGPAR::Dihedral(numComp));
	    l = par.dihedrals.size() - 1;
	    par.dihedrals[l].number = par.dihedrals.size();
	    par.dihedrals[l].atom1 = atom1;
	    par.dihedrals[l].atom2 = atom2;
	    par.dihedrals[l].atom3 = atom3;
	    par.dihedrals[l].atom4 = atom4;
	    par.dihedrals[l].multiplicity[Identity] = 1;
	    par.dihedrals[l].forceConstant[Identity].push_back(toReal(j[4]));
	    par.dihedrals[l].periodicity[Identity].push_back(toInt(j[5]));                                                        
	    par.dihedrals[l].phaseShift[Identity].push_back(toReal(j[6]));     
	  } // end inner if-else statement
	}
	else if (s == "wwwwwddddd") {

	  /*
	  // Charmm19 with multiplicity > 1        
	  // Add the dihedral ...
	  iSGPAR::Dihedral dihedral(par.dihedrals.size()+1,j[0],j[1],j[2],j[3],toReal(j[6]),toInt(j[7]),toReal(j[8]));
	  par.dihedrals.push_back(dihedral);
	  int multiplicity = toInt(j[5]);
	  // ... and the multiplicity
	  if(multiplicity > 1){
	  for(int m = 1;m<multiplicity;++m){
	  // Next line
	  ++i;
	  ++sig;
	  ++line;
	  if(i == input.end()){
	  report << warning <<"Uncomplete dihedral definition in PAR file at line "<<(*line)<<"."<<endr;
	  break;
	  }
	  if(!equalStart("ddd",(*sig))){
	  report << warning <<"Wrong multiple dihedral definition in PAR file at line "<<(*line)<<"."<<endr;
	  --i;
	  --sig;
	  --line;
	  break;
	  }
	  j = i->begin();            
	  par.dihedrals[par.dihedrals.size()-1].multiplicity++;
	  par.dihedrals[par.dihedrals.size()-1].forceConstant.push_back(toReal(j[0]));
	  par.dihedrals[par.dihedrals.size()-1].periodicity.push_back(toInt(j[1]));
	  par.dihedrals[par.dihedrals.size()-1].phaseShift.push_back(toReal(j[2]));
	  }
	  }
	  */
	  report << warning << "Charm19 dihedrals not currently accepted.  Ignoring the dihedral definition \'"
		 << definition <<"\' ("<<s<<")in PAR file at line " <<(*line)<<"."<<endr;
	}    
	break;
      }  // end case DIHEDRAL statement

      case IMPROPER: {
	if (s == "wwwwdddd") {
	  string atom1(j[0]);
	  string atom2(j[1]);
	  string atom3(j[2]);
	  string atom4(j[3]);
	  int l = par.impropers.size() - 1;
	  if ( (l < 0 ||
		((par.impropers[l].atom1 != atom1 && par.impropers[l].atom2 != atom2 && 
		  par.impropers[l].atom3 != atom3 && par.impropers[l].atom4 != atom4) &&
		 (par.impropers[l].atom1 != atom4 && par.impropers[l].atom2 != atom3 &&
		  par.impropers[l].atom3 != atom2 && par.impropers[l].atom4 != atom1)))) {

	    // this is an improper type that we have not seen yet
	    par.impropers.push_back(iSGPAR::Improper(numComp));
	    l = par.impropers.size() - 1;
	    par.impropers[l].number = par.impropers.size();
	    par.impropers[l].atom1 = atom1;
	    par.impropers[l].atom2 = atom2;
	    par.impropers[l].atom3 = atom3;
	    par.impropers[l].atom4 = atom4;
	  }  // end if NOT new bond type statement

	  par.impropers[l].forceConstant[toInt(j[7])] = toReal(j[4]);
	  par.impropers[l].periodicity[toInt(j[7])] = toInt(j[5]);
	  par.impropers[l].phaseShift[toInt(j[7])] = toReal(j[6]);
	  
	}  // end if (s = "wwwwdddd") statement     
	else if (s != "wwwwddd")
	  report << warning <<"Unknown improper definition \'"<<definition<<"\' ("<<s<<") in PAR file at line "<<(*line)<<"."<<endr;
	break;
      }  // end case IMPROPER statement

      case NONBONDED: {     
 	if((!isCharmm19 && s != "wdddddd" && s != "wddd" && s != "wdddd" && s != "wddddddd") || 
	   (isCharmm19 && s != "wdddd" && s != "wddddd")){
	  report << warning <<"Unknown nonbonded definition \'"<<definition<<"\' ("<<s<<") in PAR file at line "<<(*line)<<"."<<endr;
	  break;
	}
	else if (s == "wdddd" || s == "wddddd" || s == "wddddddd") {
	  string atom1(j[0]);
	  int l = par.nonbondeds.size() - 1;
	  if ( l < 0 || par.nonbondeds[l].atom != atom1) {
	    // this is an atom type that we have not seen yet
	    par.nonbondeds.push_back(iSGPAR::Nonbonded(numComp));
	    l = par.nonbondeds.size() - 1;
	    par.nonbondeds[l].number   = par.nonbondeds.size();
	    par.nonbondeds[l].atom     = j[0]; 
	  }

	  int Identity;
	  if (s == "wdddd") {Identity = toInt(j[4]);}
	  else {Identity = (isCharmm19 ? toInt(j[5]) : toInt(j[7]));}

	  par.nonbondeds[l].polarizability[Identity]  = (isCharmm19?0.0:toReal(j[1]));
	  par.nonbondeds[l].epsilon[Identity]         = toReal(j[1+(isCharmm19?0:1)]);
	  par.nonbondeds[l].sigma[Identity]           = toReal(j[2+(isCharmm19?0:1)])*(isCharmm19?iSGPAR::Nonbonded::SIGMA_CHARMM19_TO_CHARMM28:1.0);

	  par.nonbondeds[l].polarizability2[Identity] = (s.size()>=8?toReal(j[4]):0.0);
	  par.nonbondeds[l].epsilon14[Identity]       = (s.size()>=6?toReal(j[3+(isCharmm19?0:2)]):0.0);
	  par.nonbondeds[l].sigma14[Identity]         = (s.size()>=6?toReal(j[4+(isCharmm19?0:2)]):0.0)*(isCharmm19?iSGPAR::Nonbonded::SIGMA_CHARMM19_TO_CHARMM28:1.0);

	  par.nonbondeds[l].negative[Identity]        = (par.nonbondeds[l].epsilon[Identity] < 0.0);
	  par.nonbondeds[l].vdw[Identity]             = (s.size()>=5);
	  par.nonbondeds[l].negative2[Identity]       = (par.nonbondeds[l].epsilon14[Identity] < 0.0);
	}  // end else-if statement
	break;
      }

      case NBFIX: {
	if(s == "wwdddd")
	  par.nbfixs.push_back(iSGPAR::Nbfix(par.nbfixs.size()+1,j[0],j[1],toReal(j[2]),toReal(j[3]),toReal(j[4]),toReal(j[5])));
	else if(s == "wwdd")
	  par.nbfixs.push_back(iSGPAR::Nbfix(par.nbfixs.size()+1,j[0],j[1],toReal(j[2]),toReal(j[3]),toReal(j[2]),toReal(j[3])));
	else
	  report << warning <<"Unknown nbfix definition \'"<<definition<<"\' ("<<s<<") in PAR file at line "<<(*line)<<"."<<endr;

	break;
      }

      case HBOND: {
	if(s == "wwdd")
	  par.hbonds.push_back(iSGPAR::Hbond(par.hbonds.size()+1,j[0],j[1],toReal(j[2]),toReal(j[3])));
	else
	  report << warning <<"Unknown hbond definition \'"<<definition<<"\' ("<<s<<") in PAR file at line "<<(*line)<<"."<<endr;
	break;
      }

      case UNDEFINED:
      default: {
	report << warning <<"Unknown definition \'"<<definition<<"\' ("<<s<<") in PAR file at line "<<(*line)<<"."<<endr;
	break;
      }
      }    
    }
    close();
    return !myFile.fail();
  }

  //_________________________________________________________________isKeywordCharmm19
  bool iSGPARReader::isKeywordCharmm19(const string& word){
    return (equalStartNocase("bond",word) || equalStartNocase("angl",word) ||
	    equalStartNocase("dihe",word) || equalStartNocase("impr",word) ||
	    equalStartNocase("nonb",word) || equalStartNocase("nbfi",word) ||
	    equalStartNocase("hbon",word) );
  }

  //_________________________________________________________________isKeywordCharmm28
  bool iSGPARReader::isKeywordCharmm28(const string& word) {
    return (equalStartNocase("BOND"     ,word) || equalStartNocase("ANGL" ,word) ||
	    equalStartNocase("THET"     ,word) || equalStartNocase("DIHE" ,word) ||
	    equalStartNocase("PHI"      ,word) || equalStartNocase("IMPH" ,word) ||
	    equalStartNocase("IMPROPER" ,word) || equalStartNocase("NBOND",word) ||
	    equalStartNocase("NONBONDED",word) || equalStartNocase("NBFIX",word) ||
	    equalStartNocase("HBOND"    ,word) || equalStartNocase("END"  ,word));
  }

  iSGPARReader& operator>>(iSGPARReader& parReader, iSGPAR& par){
    parReader.read(par);
    return parReader;
  }

}
