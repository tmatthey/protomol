#include "XSCReader.h"

#include "Report.h"
#include "stringutilities.h"

using std::string;
using std::vector;
using std::stringstream;
using namespace ProtoMol::Report;

namespace ProtoMol {

  //_________________________________________________________________XSCReader

  XSCReader::XSCReader():Reader(),myXSCinfo(NULL){}

  XSCReader::XSCReader(const std::string& filename):Reader(filename),myXSCinfo(NULL){}

  XSCReader::~XSCReader(){
    if(myXSCinfo != NULL)
      delete myXSCinfo;
  }

  bool XSCReader::tryFormat() {
    if (!open())
      return false;
    return (myFile.good() );
  }

  bool XSCReader::read() {
    if(myXSCinfo == NULL)
      myXSCinfo = new XSC();
    return read(*myXSCinfo);
  }

  bool XSCReader::read(XSC& xsc) {

    if(!tryFormat())
      return false;
    if (!open())
      return false;

    // Remove comments and reformat
    stringstream all;
    while(!myFile.eof() && !myFile.fail()){

        //  Read a line of data and remove surrounding white space.  Ignore line
        //  if it begins with a '#'.
        string line(getline());      
        stringstream ss(string(line.begin(),std::find(line.begin(),line.end(),'#')));
        string str;
        //  Read line from stream and print a space if the string is not the empty
        //  line.
        while(ss >> str){
            all << ( all.str().empty() ? "" : " " ) << str;
        }
    }

    string str;
    XSCRecordTypeEnum type = UNDEFINED;

    //  Read configuration file and store any keyword pairs.

    while( all >> str ) {

        if( equalNocase("integrator", str) ) {

            all >> str;

            if(equalNocase("isgverlet", str)){
                type = ISG;
            }
            else if(equalNocase("nptleapfrog", str) || equalNocase("nptverlet", str)){
                type = NPT;
            }
            else if(equalNocase("nvtleapfrog", str)){
                type = NVTLEAPFROG;
            }

            //report << "integrator == " << str << endr;

        }

        else if( equalNocase( "molecule", str ) ) {

            all >> str;
            
            if( !toUInt( str, xsc.myMolecule ) )
                report << warning << "XSCReader: Problem with myMolecule." <<
                        endr;

            //report << "molecule == " << xsc.myMolecule << endr;

        }

        else if( equalNocase( "oldType", str ) ) {

            all >> str;
            
            if( !toUInt( str, xsc.old_type ) )
                report << warning << "XSCReader: Problem with oldType." <<
                        endr;

            //report << "oldType == " << xsc.old_type << endr;

        }

        else if( equalNocase( "newType", str ) ) {

            all >> str;
            
            if( !toUInt( str, xsc.new_type ) )
                report << warning << "XSCReader: Problem with newType." <<
                        endr;

            //report << "newType == " << xsc.new_type << endr;

        }

        else if( equalNocase( "Lambda", str ) ) {

            all >> str;
            
            if( !toReal( str, xsc.Lambda ) )
                report << warning << "XSCReader: Problem with Lambda." <<
                        endr;

            //report << "Lambda == " << xsc.Lambda << endr;

        }

        else if( equalNocase( "LambdaVel", str ) ) {

            all >> str;
            
            if( !toReal( str, xsc.Lambda_vel ) )
                report << warning << "XSCReader: Problem with LambdaVel." <<
                        endr;

            //report << "LambdaVel == " << xsc.Lambda_vel << endr;

        }

        else if( equalNocase( "Eta", str ) ) {

            all >> str;
            
            if( !toReal( str, xsc.Eta ) )
                report << warning << "XSCReader: Problem with Eta." <<
                        endr;

            //report << "Eta == " << xsc.Eta << endr;

        }

        else if( equalNocase( "EtaVel", str ) ) {

            all >> str;
            
            if( !toReal( str, xsc.Eta_vel ) )
                report << warning << "XSCReader: Problem with EtaVel." <<
                        endr;

            //report << "EtaVel == " << xsc.Eta_vel << endr;

        }

        else if( equalNocase( "EtaVol", str ) ) {

            all >> str;
            
            if( !toReal( str, xsc.EtaVol ) )
                report << warning << "XSCReader: Problem with EtaVol." <<
                        endr;

            //report << "EtaVol == " << xsc.EtaVol << endr;

        }

        else if( equalNocase( "EtaVolVel", str ) ) {

            all >> str;
            
            if( !toReal( str, xsc.EtaVol_vel ) )
                report << warning << "XSCReader: Problem with EtaVolVel." <<
                        endr;

            //report << "EtaVolVel == " << xsc.EtaVol_vel << endr;

        }

        else if( equalNocase( "EtaLambda", str ) ) {

            all >> str;
            
            if( !toReal( str, xsc.EtaLambda ) )
                report << warning << "XSCReader: Problem with EtaLambda." <<
                        endr;

            //report << "EtaLambda == " << xsc.EtaLambda << endr;

        }

        else if( equalNocase( "EtaLambdaVel", str ) ) {

            all >> str;
            
            if( !toReal( str, xsc.EtaLambda_vel ) )
                report << warning << "XSCReader: Problem with EtaLambdaVel." <<
                        endr;

            //report << "EtaLambdaVel == " << xsc.EtaLambda_vel << endr;

        }
        
        else if( equalNocase( "Volume", str ) ) {

            all >> str;
            
            if( !toReal( str, xsc.Vol ) )
                report << warning << "XSCReader: Problem with Volume." <<
                        endr;

            //report << "Volume == " << xsc.Vol << endr;

        }

        else if( equalNocase( "EpsilonVel", str ) ) {

            all >> str;
            
            if( !toReal( str, xsc.Epsilon_vel ) )
                report << warning << "XSCReader: Problem with EpsilonVel." <<
                        endr;

            //report << "EpsilonVel == " << xsc.Epsilon_vel << endr;

        }

        else {

            report << warning << "XSCReader - unknown keyword: " << str << endr;

        }

    }

    // Read in the appropriate XSC inputs for the particular integrator type...
    switch(type) {
        //-------------------------------------------------------------------
        case ISG: 
        case NPT: 
        case NVTLEAPFROG: 
            break;
        default:
            report << error << "XSCReader - unknown integrator." << endr;

            //*******************************************************************
    } // end switch structure

    return !myFile.fail();
  }

  XSC* XSCReader::orphanXSC(){
    XSC* tmp = myXSCinfo;
    myXSCinfo = NULL;
    return tmp;
  }

  XSCReader& operator>>(XSCReader& xscReader, XSC& xsc){
    xscReader.read(xsc);
    return xscReader;
  }

}
