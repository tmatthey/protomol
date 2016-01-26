#include "buildOSGTopology.h"

#include "ExclusionTable.h"
#include "ExclusionType.h"
#include "GenericTopology.h"
#include "LennardJonesParameterTable.h"
#include "PAR.h"
#include "PSF.h"
#include "Report.h"
#include "pmconstants.h"
#include "mathutilities.h"
#include "stringutilities.h"
#include "topologyutilities.h"
//#include "Molecule.h"
//#include "Vector3D.h"

#include <algorithm>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

using std::list;
using std::map;
using std::set;
using std::pair;
using std::string;
using std::vector;
using std::swap;

using namespace ProtoMol::Report;

//#define DEBUG_PRINT_MOLECULETABLE

namespace ProtoMol {

    static void findNextNeighbor(int a, vector<int>& v, vector<PairInt>& p, vector<char>& unused, const vector<vector<int> >& graph, set<PairInt>& pairs);
    // bool cmpSize(const vector<int>& m1, const vector<int>& m2);


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //  buildTopology
    //
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    void buildOSGTopology(GenericTopology* topo,const PSF& psf, const PAR& par, bool dihedralMultPSF){
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // First, generate the array of atomtypes
        // Each time a new atom comes up, we need to check if it is
        // already in the vector....
        // NOTE:  this may take a while for large systems; however, it will cut
        // down on the size of the atomTypes vector, and therefore, the amount
        // access time in the back end.
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        topo->atoms.clear();
        topo->atomTypes.clear();
        topo->bonds.clear();
        topo->angles.clear();
        topo->dihedrals.clear();
        topo->impropers.clear();

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Get the atoms
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        map<string,int> atomLookUpTable;

        // loop over all atoms in the PSF object
        for (vector<PSF::Atom>::const_iterator atom = psf.atoms.begin();
                atom != psf.atoms.end();
                ++atom) {
            // Two data members for AtomType, name and mass
            AtomType tempatomtype;
            tempatomtype.name = atom->atom_type;
            tempatomtype.mass = atom->mass;
            // tempatomtype.symbolName = charmmDefaults->getNonbondedData(atom->atom_type).mySymbolName;
            tempatomtype.symbolName = atomTypeToSymbolName(atom->atom_type);
            tempatomtype.charge = atom->charge;
            // Now check if this already exists (same name)
            if (atomLookUpTable.find(tempatomtype.name) == atomLookUpTable.end()) {
                atomLookUpTable[tempatomtype.name] = topo->atomTypes.size();
                topo->atomTypes.push_back(tempatomtype);
            }

            Atom tempatom;
            // First, we need to find the index. (an integer corresponding
            // to the type of the atom
            tempatom.type = atomLookUpTable[tempatomtype.name];
            tempatom.name = atom->atom_name;
            // Now, the scaled charge.  This is straightforward.
            tempatom.scaledCharge = (atom->charge)*Constant::SQRTCOULOMBCONSTANT;
            tempatom.scaledMass = atom->mass;

            // Now we need the size of the group for heavy atom ordering
            // We need to parse the name for any H's then any numbers following
            // First, if the atom is an H then this is 0
            if (atom->atom_type == "H"){
                tempatom.hvyAtom = 0;
            }
            else{
                // Otherwise, we need to parse..
                // Initialize to 1
                tempatom.hvyAtom = 1;
                for (unsigned int pos = 0; pos < atom->atom_type.size(); ++pos){
                    if (atom->atom_type[pos] == 'H'){
                        string number = "";
                        while (isdigit(atom->atom_type[++pos]))   {
                            number += atom->atom_type[pos];
                        }
                        if (number == "") // never entered loop, default is 1
                            number = "1";
                        tempatom.hvyAtom += atoi(number.c_str());
                    }
                }
            }
            // C/C++ starts at 0, where PSF/PDB at 1
            tempatom.atomNum = atom->number-1;
            // Also the molecule - using residue sequence for now
            topo->atoms.push_back(tempatom);


        }

        // calculate the # of degrees of freedom, if there are any bond constraints
        // they will be subtracted later by ModifierShake
        topo->degreesOfFreedom = 3 * topo->atoms.size() - 3;

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Get the bonds
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // First create look-up-table
        map<string,vector<PAR::Bond>::const_iterator> bondLookUpTable;
        for (vector<PAR::Bond>::const_iterator bond = par.bonds.begin();
                bond != par.bonds.end();
                ++bond){
            bondLookUpTable[bond->atom1+","+bond->atom2] = bond;
            //report << (*bond)<< ", " << bond->atom1 << ", " << bond->atom2 <<endr;
        }


        // Find the parameters from PAR
        int ignoredBonds = 0;
        for (vector<PSF::Bond>::const_iterator bond = psf.bonds.begin();
                bond != psf.bonds.end();
                ++bond){

            // store the ID numbers of the bonded atoms
            int atom1 = bond->atom1-1;
            int atom2 = bond->atom2-1;

            // store the type names of the bonded atoms
            string bond1(topo->atomTypes[topo->atoms[atom1].type].name);
            string bond2(topo->atomTypes[topo->atoms[atom2].type].name);

            map<string,vector<PAR::Bond>::const_iterator>::const_iterator currentbond = bondLookUpTable.find(bond1+","+bond2);
            if(currentbond == bondLookUpTable.end()){
                currentbond = bondLookUpTable.find(bond2+","+bond1);
            }

            // if we still have not found this bond type in the PAR object, report an error
            if(currentbond == bondLookUpTable.end()){
                report << error << "Could not find bond \'"<<bond1<<"\'-\'"<<bond2<<"\' ("<<bond->atom1 <<","<< bond->atom2<<")"<<std::endl;
                for (map<string,vector<PAR::Bond>::const_iterator>::const_iterator i = bondLookUpTable.begin();
                        i != bondLookUpTable.end();
                        i++){
                    report << plain << i->first<<std::endl;
                }
                report << endr;
            }

            // if we have found this bond type then copy the bond parameters
            // into the topology
            Bond tempbond;
            tempbond.springConstant = currentbond->second->forceConstant;
            tempbond.restLength = currentbond->second->distance;
            tempbond.atom1 = atom1;
            tempbond.atom2 = atom2;
            topo->bonds.push_back(tempbond);

            // populate the vector of bonds maintained at each atom
            //report << std::endl <<"Size of Bonds Vector "<<topo->bonds.size() <<endr;
            topo->atoms[atom1].mybonds.push_back((topo->bonds.size())-1);
            topo->atoms[atom2].mybonds.push_back((topo->bonds.size())-1);
            // output to screen for testing purposes
            //for (int j = 0; j < topo->atoms[atom1].mybonds.size(); j++){
            //report <<"Atom " << atom1 << " bond index = "<< topo->atoms[atom1].mybonds[j] <<endr;
            //}
            //for (int k = 0; k < topo->atoms[atom2].mybonds.size(); k++){
            //report <<"Atom " << atom2 << " bond index = "<< topo->atoms[atom2].mybonds[k] <<endr;
            //}

            if(tempbond.springConstant == 0.0)
                ignoredBonds++;
        }

        if(ignoredBonds > 0)
            report << hint << "Systems contains "<<ignoredBonds<<" bonds with zero force constants."<<endr;

        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Get the angles
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // First create look-up-table
        map<string,vector<PAR::Angle>::const_iterator> angleLookUpTable;
        for (vector<PAR::Angle>::const_iterator angle = par.angles.begin();
                angle != par.angles.end();
                ++angle){
            angleLookUpTable[angle->atom1+","+angle->atom2+","+angle->atom3] = angle;
            //report << (*angle)<<endr;
        }

        // Find the parameters from PAR
        int ignoredAngles = 0;

        // loop over the angle list in the PSF object
        for (vector<PSF::Angle>::const_iterator angle = psf.angles.begin();
                angle != psf.angles.end(); ++angle){

            // store the ID numbers of the atoms in this angle
            int atom1 = angle->atom1-1;
            int atom2 = angle->atom2-1;
            int atom3 = angle->atom3-1;

            // store the type names of the atoms in this angle
            string angle1(topo->atomTypes[topo->atoms[atom1].type].name);
            string angle2(topo->atomTypes[topo->atoms[atom2].type].name);
            string angle3(topo->atomTypes[topo->atoms[atom3].type].name);

            map<string,vector<PAR::Angle>::const_iterator>::const_iterator currentangle = angleLookUpTable.find(angle1+","+angle2+","+angle3);
            if(currentangle == angleLookUpTable.end())
                currentangle = angleLookUpTable.find(angle3+","+angle2+","+angle1);

            // if we still have not found this angle type in the PAR object, report an error
            if(currentangle == angleLookUpTable.end())
                report << error << "Could not find angle \'"<<angle1<<"\'-\'"<<angle2<<"\'-\'"<<angle3<<"\'."<<endr;

            // if we have found this angle type then copy the angle parameters
            // into the topology
            Angle tempangle;
            tempangle.atom1 = atom1;
            tempangle.atom2 = atom2;
            tempangle.atom3 = atom3;
            tempangle.forceConstant = currentangle->second->forceConstant;
            tempangle.restAngle = dtor(currentangle->second->angleval);
            if (currentangle->second->ub_flag){ // do we want defaults for these
                tempangle.ureyBradleyConstant = currentangle->second->k_ub;
                tempangle.ureyBradleyRestLength = currentangle->second->r_ub;
            }
            // no Urey-Bradley term specified
            else{
                tempangle.ureyBradleyConstant = 0.0;
                tempangle.ureyBradleyRestLength = 0.0;
            }
            topo->angles.push_back(tempangle);
            if(tempangle.forceConstant == 0.0)
                ignoredAngles++;
        }

        if(ignoredAngles > 0)
            report << hint << "Systems contains "<<ignoredAngles<<" angles with zero force constants."<<endr;


        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Get the dihedrals
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        // One change I made was to assume that a dihedral will only appear
        // once in the .psf file regardless of it's multiplicity.  The
        // multiplicity should be handled in the .par file.
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


        // First create look-up-table
        map<string,vector<PAR::Dihedral>::const_iterator> dihedralLookUpTable;
        for (vector<PAR::Dihedral>::const_iterator dihedral = par.dihedrals.begin();
                dihedral != par.dihedrals.end();
                ++dihedral){
            dihedralLookUpTable[dihedral->atom1+","+dihedral->atom2+","+dihedral->atom3+","+dihedral->atom4] = dihedral;
            //report << (*dihedral)<<endr;
        }

        // Find the parameters from PAR
        // loop over the dihedral list in the PSF object
        for(vector<PSF::Dihedral>::const_iterator dihedral = psf.dihedrals.begin();
                dihedral != psf.dihedrals.end(); ++dihedral) {

            // store the ID numbers of the atoms in this dihedral
            int atom1 = dihedral->atom1 - 1;
            int atom2 = dihedral->atom2 - 1;
            int atom3 = dihedral->atom3 - 1;
            int atom4 = dihedral->atom4 - 1;

            // store the type names of the atoms in this dihedral
            string dihedral1 = topo->atomTypes[topo->atoms[atom1].type].name;
            string dihedral2 = topo->atomTypes[topo->atoms[atom2].type].name;
            string dihedral3 = topo->atomTypes[topo->atoms[atom3].type].name;
            string dihedral4 = topo->atomTypes[topo->atoms[atom4].type].name;

            map<string,vector<PAR::Dihedral>::const_iterator>::const_iterator currentdihedral =
                    dihedralLookUpTable.find(dihedral1+","+dihedral2+","+dihedral3+","+dihedral4);

            // if this dihedral type has not been found, try reversing the order of the atom types
            if(currentdihedral == dihedralLookUpTable.end())
                currentdihedral = dihedralLookUpTable.find(dihedral4+","+dihedral3+","+dihedral2+","+dihedral1);

            // Try wildcards if necessary
            if(currentdihedral == dihedralLookUpTable.end()){
                currentdihedral = dihedralLookUpTable.find("X,"+dihedral2+","+dihedral3+",X");
                if(currentdihedral == dihedralLookUpTable.end())
                    currentdihedral = dihedralLookUpTable.find("X,"+dihedral3+","+dihedral2+",X");
            }

            // if we still have not found this dihedral type in the PAR object, report an error
            if(currentdihedral == dihedralLookUpTable.end())
                report << error << "Could not find dihedral \'"<<dihedral1<<"\'-\'"<<dihedral2<<"\'-\'"<<dihedral3<<"\'-\'"<<dihedral4<<"\'."<<endr;

            // if we have found this dihedral type then copy the
            // dihedral parameters into the topology
            Torsion torsion;
            torsion.atom1         = atom1;
            torsion.atom2         = atom2;
            torsion.atom3         = atom3;
            torsion.atom4         = atom4;

            torsion.periodicity   = currentdihedral->second->periodicity;
            torsion.forceConstant = currentdihedral->second->forceConstant;
            torsion.phaseShift    = dtor(currentdihedral->second->phaseShift);
            torsion.multiplicity  = currentdihedral->second->multiplicity;
            if(topo->dihedrals.empty() ||
                    topo->dihedrals[topo->dihedrals.size()-1].atom1 != atom1 ||
                    topo->dihedrals[topo->dihedrals.size()-1].atom2 != atom2 ||
                    topo->dihedrals[topo->dihedrals.size()-1].atom3 != atom3 ||
                    topo->dihedrals[topo->dihedrals.size()-1].atom4 != atom4){
                if(dihedralMultPSF){
                    torsion.periodicity.resize(1);
                    torsion.forceConstant.resize(1);
                    torsion.phaseShift.resize(1);
                    torsion.multiplicity = 1;
                }
                topo->dihedrals.push_back(torsion);
            }
            else {
                if(dihedralMultPSF){
                    Torsion& tmp = topo->dihedrals[topo->dihedrals.size()-1];
                    if(tmp.multiplicity > torsion.multiplicity )
                        report << error<< "PSF multiplicity definition of dihedral ("
                                << dihedral1 <<","
                                << dihedral2 <<","
                                << dihedral3 <<","
                                << dihedral4 <<") exceeded PAR definition.";
                    tmp.periodicity.push_back(torsion.periodicity[tmp.multiplicity]);
                    tmp.forceConstant.push_back(torsion.forceConstant[tmp.multiplicity]);
                    tmp.phaseShift.push_back(torsion.phaseShift[tmp.multiplicity]);
                    tmp.multiplicity++;
                }
                else {
                    report << error << "Unexpected PSF multiplicity definition of dihedral ("
                            << dihedral1 <<","
                            << dihedral2 <<","
                            << dihedral3 <<","
                            << dihedral4 <<") occurred, use dihedral multiplicity PSF.";
                }
            }
        }


        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Get the impropers
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        // One change I made was to assume that a improper will only appear
        // once in the .psf file regardless of it's multiplicity.  The
        // multiplicity should be handled in the .par file.
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        //         No wildcard usage is allowed for bonds and angles. For dihedrals,
        // two types are allowed; A - B - C - D (all four atoms specified) and
        // X - A - B - X (only middle two atoms specified). Double dihedral
        // specifications may be specified for the four atom type by listing a
        // given set twice. When specifying this type in the topology file, specify
        // a dihedral twice (with nothing intervening) and both forms will be used.
        //
        //         There are five choices for wildcard usage for improper dihedrals;
        // 1) A - B - C - D  (all four atoms, double specification allowed)
        // 2) A - X - X - B
        // 3) X - A - B - C
        // 4) X - A - B - X
        // 5) X - X - A - B
        // When classifying an improper dihedral, the first acceptable match (from
        // the above order) is chosen. The match may be made in either direction
        // ( A - B - C - D = D - C - B - A).
        //
        //         The periodicity value for dihedrals and improper dihedral terms
        // must be an integer. If it is positive, then a cosine functional form is used.
        // Only positive values of 1,2,3,4,5 and 6 are allowed for the vector, parallel
        // vector and cray routines. Slow and scalar routines can use any positive
        // integer and thus dihedral constrains can be of any periodicity.
        //  Reference angle 0.0 and 180.0 degree correspond to minimum in staggered
        // and eclipsed respectively. Any reference angle is allowed. The value
        // 180 should be prefered over -180 since it is parsed faster and more
        // accuratly. When the periodicity is given as zero, for OTHER THAN THE
        // FIRST dihdral in a multiple dihedral set, then a the amplitude is a
        // constant added to the energy. This is needed to effect the
        // Ryckaert-Bellemans potential for hydrocarbons (see below).




        // First create look-up-table
        map<string,vector<PAR::Improper>::const_iterator> improperLookUpTable;
        for (vector<PAR::Improper>::const_iterator improper = par.impropers.begin();
                improper != par.impropers.end();
                improper++ ){
            improperLookUpTable[improper->atom1+","+improper->atom2+","+improper->atom3+","+improper->atom4] = improper;
            //report << (*improper)<<endr;
        }

        // Find the parameters from PAR
        // loop over the improper list in the PSF object
        for(vector<PSF::Improper>::const_iterator improper = psf.impropers.begin();
                improper != psf.impropers.end(); improper++ ) {

            // store the ID numbers of the atoms in this improper
            int atom1 = improper->atom1 - 1;
            int atom2 = improper->atom2 - 1;
            int atom3 = improper->atom3 - 1;
            int atom4 = improper->atom4 - 1;

            // store the type names of the atoms in this improper
            string improper1 = topo->atomTypes[topo->atoms[atom1].type].name;
            string improper2 = topo->atomTypes[topo->atoms[atom2].type].name;
            string improper3 = topo->atomTypes[topo->atoms[atom3].type].name;
            string improper4 = topo->atomTypes[topo->atoms[atom4].type].name;


            map<string,vector<PAR::Improper>::const_iterator>::const_iterator currentimproper =
                    improperLookUpTable.find(improper1+","+improper2+","+improper3+","+improper4);
            if(currentimproper == improperLookUpTable.end())
                currentimproper = improperLookUpTable.find(improper4+","+improper3+","+improper2+","+improper1);

            // Try wildcards if necessary
            // 2) A - X - X - B
            if(currentimproper == improperLookUpTable.end()){
                currentimproper = improperLookUpTable.find(improper1+",X,X,"+improper4);
                if(currentimproper == improperLookUpTable.end())
                    currentimproper = improperLookUpTable.find(improper4+",X,X,"+improper1);

            }
            // 3) X - A - B - C
            if(currentimproper == improperLookUpTable.end()){
                currentimproper = improperLookUpTable.find("X,"+improper2+","+improper3+","+improper4);
                if(currentimproper == improperLookUpTable.end())
                    currentimproper = improperLookUpTable.find(improper4+","+improper3+","+improper2+",X");
            }

            // 4) X - A - B - X
            if(currentimproper == improperLookUpTable.end()){
                currentimproper = improperLookUpTable.find("X,"+improper2+","+improper3+",X");
                if(currentimproper == improperLookUpTable.end())
                    currentimproper = improperLookUpTable.find("X,"+improper3+","+improper2+",X");
            }

            // 5) X - X - A - B
            if(currentimproper == improperLookUpTable.end()){
                currentimproper = improperLookUpTable.find("X,X,"+improper3+","+improper4);
                if(currentimproper == improperLookUpTable.end())
                    currentimproper = improperLookUpTable.find(improper4+","+improper3+",X,X");
            }

            // if we still have not found this improper type in the PAR object, report an error
            if(currentimproper == improperLookUpTable.end())
                report << error << "Could not find improper."<<endr;

            // if we have found this improper type then copy the
            // improper parameters into the topology
            Torsion torsion;
            torsion.atom1         = atom1;
            torsion.atom2         = atom2;
            torsion.atom3         = atom3;
            torsion.atom4         = atom4;
            torsion.periodicity.push_back(currentimproper->second->periodicity);
            torsion.forceConstant.push_back(currentimproper->second->forceConstant);
            torsion.phaseShift.push_back(dtor(currentimproper->second->phaseShift));
            torsion.multiplicity  = 1;
            topo->impropers.push_back(torsion);
            //     report << plain<< (wildcard ? "#":"")
            //            << improper1 <<","
            //            << improper2 <<","
            //            << improper3 <<","
            //            << improper4 <<","
            //            << torsion.atom1 <<","
            //            << torsion.atom2 <<","
            //            << torsion.atom3 <<","
            //            << torsion.atom4 <<","
            //            << torsion.periodicity[0] <<","
            //            << torsion.periodicity.size() <<","
            //            << torsion.forceConstant[0] <<","
            //            << torsion.phaseShift[0] <<","
            //            << torsion.multiplicity << endr;
        }

        // store the molecule information
        buildMoleculeTable(topo);
        buildExclusionTable(topo,topo->exclude);  /// WHY DID I COMMENT THIS OUT BEFORE???

	// build the list of bonds, angles, etc. belonging to each molecule
    	buildMoleculeBondingLists(topo);

        // determine the time-zero identity of each molecule
        for(unsigned int i=0; i<topo->molecules.size(); i++) {

          // the ID# of the first atom in this molecule
          int firstAtom = topo->molecules[i][0];

          // get the identity number from the psf file
          topo->molecules[i].type = psf.atoms[firstAtom].identity;
          topo->molecules[i].newtype = psf.atoms[firstAtom].identity;
          topo->molecules[i].name = psf.atoms[firstAtom].seg_id;
        }

    }


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //  buildMoleculeTable
    //
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    void buildMoleculeTable(GenericTopology *topo){

        // *** First we clear all molecules ***
        topo->molecules.clear();

        const unsigned int numAtoms = topo->atoms.size();

        // *** Collecting all possible connections, building the graph ***
        vector<vector<int> > graph(numAtoms,vector<int>());
        set<pair<int,int> > pairs;
        // *** Bonds ***
        for(unsigned int i=0;i<topo->bonds.size();i++){
            int a1 = topo->bonds[i].atom1;
            int a2 = topo->bonds[i].atom2;
            graph[a1].push_back(a2);
            graph[a2].push_back(a1);
            pairs.insert(pair<int,int>(std::min(a1,a2),std::max(a1,a2)));
        }
        unsigned int count = pairs.size();

        // *** Angles ***
        for(unsigned int i=0;i<topo->angles.size();i++){
            int a1 = topo->angles[i].atom1;
            int a2 = topo->angles[i].atom2;
            int a3 = topo->angles[i].atom3;
            graph[a1].push_back(a2);
            graph[a1].push_back(a3);
            graph[a2].push_back(a1);
            graph[a2].push_back(a3);
            graph[a3].push_back(a1);
            graph[a3].push_back(a2);
            pairs.insert(pair<int,int>(std::min(a1,a2),std::max(a1,a2)));
            pairs.insert(pair<int,int>(std::min(a3,a2),std::max(a3,a2)));
        }

        if(count < pairs.size())
            report << hint << "Angles added "<<pairs.size()-count<<" new bond(s)." <<endr;
        count = pairs.size();

        // *** Dihedrals ***
        for(unsigned int i=0;i<topo->dihedrals.size();i++){
            int a1 = topo->dihedrals[i].atom1;
            int a2 = topo->dihedrals[i].atom2;
            int a3 = topo->dihedrals[i].atom3;
            int a4 = topo->dihedrals[i].atom4;
            graph[a1].push_back(a2);
            graph[a1].push_back(a3);
            graph[a1].push_back(a4);
            graph[a2].push_back(a1);
            graph[a2].push_back(a3);
            graph[a2].push_back(a4);
            graph[a3].push_back(a1);
            graph[a3].push_back(a2);
            graph[a3].push_back(a4);
            graph[a4].push_back(a1);
            graph[a4].push_back(a2);
            graph[a4].push_back(a3);
            pairs.insert(pair<int,int>(std::min(a1,a2),std::max(a1,a2)));
            pairs.insert(pair<int,int>(std::min(a3,a2),std::max(a3,a2)));
            pairs.insert(pair<int,int>(std::min(a3,a4),std::max(a3,a4)));
        }
        if(count < pairs.size())
            report << hint << "Dihedrals added "<<pairs.size()-count<<" new bond(s)." <<endr;
        count = pairs.size();

        // *** Impropers ***
        set<pair<int,int> > pairsAddImpropers;
        // Impropers are defined over the bonds 1-2,1-3,1-4 or 4-1,4-3,4-2
        // but MTorsionSystemForce computes distances betweeen 1-2,2-3,3-4
        // we have to take care about these differences ...
        for(unsigned int i=0;i<topo->impropers.size();i++){
            int a1 = topo->impropers[i].atom1;
            int a2 = topo->impropers[i].atom2;
            int a3 = topo->impropers[i].atom3;
            int a4 = topo->impropers[i].atom4;
            graph[a1].push_back(a2);
            graph[a1].push_back(a3);
            graph[a1].push_back(a4);
            graph[a2].push_back(a1);
            graph[a2].push_back(a3);
            graph[a2].push_back(a4);
            graph[a3].push_back(a1);
            graph[a3].push_back(a2);
            graph[a3].push_back(a4);
            graph[a4].push_back(a1);
            graph[a4].push_back(a2);
            graph[a4].push_back(a3);
            pair<int,int> p0(std::min(a1,a2),std::max(a1,a2));
            pair<int,int> p1(std::min(a1,a3),std::max(a1,a3));
            pair<int,int> p2(std::min(a1,a4),std::max(a1,a4));
            pair<int,int> p3(std::min(a2,a3),std::max(a2,a3));
            pair<int,int> p4(std::min(a2,a4),std::max(a2,a4));
            pair<int,int> p5(std::min(a3,a4),std::max(a3,a4));
            int j0 = 0;
            int j1 = 0;
            int j2 = 0;
            int j3 = 0;
            int j4 = 0;
            int j5 = 0;
            if(pairs.find(p0) != pairs.end()) j0++;
            if(pairs.find(p1) != pairs.end()) j1++;
            if(pairs.find(p2) != pairs.end()) j2++;
            if(pairs.find(p3) != pairs.end()) j3++;
            if(pairs.find(p4) != pairs.end()) j4++;
            if(pairs.find(p5) != pairs.end()) j5++;
            if(j0+j1+j2+j3+j4+j5 < 3){
                pairs.insert(p0);
                pairs.insert(p1);
                pairs.insert(p2);
            }
            pairsAddImpropers.insert(p0);
            pairsAddImpropers.insert(p3);
            pairsAddImpropers.insert(p5);
        }
        if(count < pairs.size())
            report << hint << "Impropers added "<<pairs.size()-count<<" new bond(s)." <<endr;

        // Now add the improper pairs
        for(set<pair<int,int> >::const_iterator i= pairsAddImpropers.begin();i != pairsAddImpropers.end();i++)
            pairs.insert(*i);


        count = pairs.size();
        //report << hint << count << endr;
        // To keep track which atoms already have been added
        // to molecules.
        vector<char> unused(numAtoms,1);


        // Recursively finding the atoms beloning to a molecule
        for(unsigned int i=0;i<numAtoms;i++){
            vector<int> v;
            vector<PairInt> p;
            findNextNeighbor(i,v,p,unused,graph,pairs);
            if(!v.empty()){
                std::sort(v.begin(),v.end());
                // add this atom list to the molecules array
                Molecule mol;
                mol.atoms = v;
                for(unsigned int j=0;j<p.size();j++)
                    if(p[j].first > p[j].second)
                        swap(p[j].first,p[j].second);
                std::sort(p.begin(),p.end());
                mol.pairs = p;
                topo->molecules.push_back(mol);
            }
        }

        // Uncomment to sort descending after size()
        // std::sort(topo->molecules.begin(),topo->molecules.end(),cmpSize);

        // Look up table for atoms
        const string h("H");
        const string o("O");
        for(unsigned int i=0;i<topo->molecules.size();i++){
            Real mass = 0.0;
            const vector<int>& mol = topo->molecules[i].atoms;
            for(unsigned int j=0;j<mol.size();j++){
                int k = mol[j];
                topo->atoms[k].molecule = i;
                mass += topo->atoms[k].scaledMass;
            }
            topo->molecules[i].mass = mass;
            topo->molecules[i].water = (mol.size() == 3 &&
                    ((topo->atomTypes[topo->atoms[mol[0]].type].symbolName == h &&
                      topo->atomTypes[topo->atoms[mol[1]].type].symbolName == h &&
                      topo->atomTypes[topo->atoms[mol[2]].type].symbolName == o) ||
                     (topo->atomTypes[topo->atoms[mol[0]].type].symbolName == h &&
                      topo->atomTypes[topo->atoms[mol[1]].type].symbolName == o &&
                      topo->atomTypes[topo->atoms[mol[2]].type].symbolName == h) ||
                     (topo->atomTypes[topo->atoms[mol[0]].type].symbolName == o &&
                      topo->atomTypes[topo->atoms[mol[1]].type].symbolName == h &&
                      topo->atomTypes[topo->atoms[mol[2]].type].symbolName == h)));

        }

#if defined(DEBUG_PRINT_MOLECULETABLE)
        report<< plain << endl << "[buildMoleculeTable]: molecule table printout:" << endl;
        for(int i=0;i<topo->molecules.size();i++){
            for(int j=0;j<topo->molecules[i].size();j++){
                report << topo->molecules[i][j]<<" ";
            }
            report << endl;
        }
        report << endr;
#endif
    }


    //____________________________________________________________findNextNeighbor
    void findNextNeighbor(int a, vector<int>& v, vector<PairInt>& p, vector<char>& unused,
            const vector<vector<int> >& graph, set<PairInt>& pairs){
        if(unused[a]>0){
            v.push_back(a);
            unused[a] = 0;
            for(unsigned int i=0;i<graph[a].size();i++){
                set<PairInt>::iterator itr = pairs.find(PairInt(std::min(a,graph[a][i]),std::max(a,graph[a][i])));
                if(itr != pairs.end()){
                    p.push_back(PairInt(a,graph[a][i]));
                    pairs.erase(itr);
                }
                findNextNeighbor(graph[a][i],v,p,unused,graph,pairs);
            }
        }
    }

    //_____________________________________________________________________ cmpSize
    // bool cmpSize(const vector<int>& m1, const vector<int>& m2){
    //   return (m1.size() > m2.size());
    // }



    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //  buildExclusionTable
    //
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    void buildExclusionTable(GenericTopology* topo, const ExclusionType& exclusionType) {


        if(!exclusionType.valid())
            report << error <<"[buildExclusionTable()] Exclusion type not defined/valid."<<endr;


        topo->exclude = exclusionType;

        //  Resize array.
	topo->exclusions.resize(topo->atoms.size());

        //  If exclusionType is equal to NONE, return.
        if( exclusionType == ExclusionType::NONE )
            return;


        const int numBonds     = topo->bonds.size(),
                  numAngles    = topo->angles.size(),
                  numDihedrals = topo->dihedrals.size();



        //  Add excluded bonds.
        for( int i = 0; i < numBonds; i++ )
            topo->exclusions.add( topo->bonds[i].atom1, topo->bonds[i].atom2,
                    EXCLUSION_FULL );

        if( exclusionType != ExclusionType::ONE2 ) {

            //  Add excluded angles.
            for( int i = 0; i < numAngles; i++ )
                topo->exclusions.add( topo->angles[i].atom1,
                        topo->angles[i].atom3, EXCLUSION_FULL );

            if( exclusionType != ExclusionType::ONE3 ) {

                //  Add excluded dihedrals.
                for( int i = 0; i < numDihedrals; i++ ) {

                    if( exclusionType == ExclusionType::ONE4 )
                        topo->exclusions.add( topo->dihedrals[i].atom1,
                                topo->dihedrals[i].atom4, EXCLUSION_FULL );
                    else
                        topo->exclusions.add( topo->dihedrals[i].atom1,
                                topo->dihedrals[i].atom4, EXCLUSION_MODIFIED );

                }

            }


        }

        topo->exclusions.optimize();

    }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  //  buildMoleculeBondingLists
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  void buildMoleculeBondingLists(GenericTopology* topo) {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // bond lists
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all bonds in the topology
    for (unsigned int i=0; i<topo->bonds.size(); i++) {

      // get the ID#s of the two atoms in this bond
      int atom1 = topo->bonds[i].atom1;
      int atom2 = topo->bonds[i].atom2;

      // loop sentinels
      bool found1 = false;
      bool found2 = false;
      unsigned int m;

      // loop over all molecules
      for (m=0; m<topo->molecules.size(); m++) {

	// loop over all atoms in the atom list for this molecule
	for (unsigned int a=0; a<topo->molecules[m].atoms.size(); a++) {

	  if (atom1 == topo->molecules[m].atoms[a]) {found1 = true;}
	  else if (atom2 == topo->molecules[m].atoms[a]) {found2 = true;}

	  if (found1 && found2) {break;}
	} // end atom list loop

	if (found1 && found2) {break;}
      } // end loop over molecules

      // store this bond index number in the molecule's bond list
      topo->molecules[m].bondList.push_back(i);
    } // end loop over bonds


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // angle lists
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all angles in the topology
    for (unsigned int i=0; i<topo->angles.size(); i++) {

      // get the ID#s of the three atoms in this angle
      int atom1 = topo->angles[i].atom1;
      int atom2 = topo->angles[i].atom2;
      int atom3 = topo->angles[i].atom3;

      // loop sentinels
      bool found1 = false;
      bool found2 = false;
      bool found3 = false;
      unsigned int m;

      // loop over all molecules
      for (m=0; m<topo->molecules.size(); m++) {

	// loop over all atoms in the atom list for this molecule
	for (unsigned int a=0; a<topo->molecules[m].atoms.size(); a++) {

	  if (atom1 == topo->molecules[m].atoms[a]) {found1 = true;}
	  else if (atom2 == topo->molecules[m].atoms[a]) {found2 = true;}
	  else if (atom3 == topo->molecules[m].atoms[a]) {found3 = true;}

	  if (found1 && found2 && found3) {break;}
	} // end atom list loop

	if (found1 && found2 && found3) {break;}
      } // end loop over molecules

      // store this angle index number in the molecule's angle list
      topo->molecules[m].angleList.push_back(i);
    } // end loop over angles

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // dihedral lists
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all diherals in the topology
    for (unsigned int i=0; i<topo->dihedrals.size(); i++) {

      // get the ID#s of the four atoms in this dihedral
      int atom1 = topo->dihedrals[i].atom1;
      int atom2 = topo->dihedrals[i].atom2;
      int atom3 = topo->dihedrals[i].atom3;
      int atom4 = topo->dihedrals[i].atom4;

      // loop sentinels
      bool found1 = false;
      bool found2 = false;
      bool found3 = false;
      bool found4 = false;
      unsigned int m;

      // loop over all molecules
      for (m=0; m<topo->molecules.size(); m++) {

	// loop over all atoms in the atom list for this molecule
	for (unsigned int a=0; a<topo->molecules[m].atoms.size(); a++) {

	  if (atom1 == topo->molecules[m].atoms[a]) {found1 = true;}
	  else if (atom2 == topo->molecules[m].atoms[a]) {found2 = true;}
	  else if (atom3 == topo->molecules[m].atoms[a]) {found3 = true;}
	  else if (atom4 == topo->molecules[m].atoms[a]) {found4 = true;}

	  if (found1 && found2 && found3 && found4) {break;}
	} // end atom list loop

	if (found1 && found2 && found3 && found4) {break;}
      } // end loop over molecules

      // store this dihedral index number in the molecule's dihderal list
      topo->molecules[m].dihedralList.push_back(i);
    } // end loop over dihedrals


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // improper lists
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // loop over all impropers in the topology
    for (unsigned int i=0; i<topo->impropers.size(); i++) {

      // get the ID#s of the four atoms in this improper
      int atom1 = topo->impropers[i].atom1;
      int atom2 = topo->impropers[i].atom2;
      int atom3 = topo->impropers[i].atom3;
      int atom4 = topo->impropers[i].atom4;

      // loop sentinels
      bool found1 = false;
      bool found2 = false;
      bool found3 = false;
      bool found4 = false;
      unsigned int m;

      // loop over all molecules
      for (m=0; m<topo->molecules.size(); m++) {

	// loop over all atoms in the atom list for this molecule
	for (unsigned int a=0; a<topo->molecules[m].atoms.size(); a++) {

	  if (atom1 == topo->molecules[m].atoms[a]) {found1 = true;}
	  else if (atom2 == topo->molecules[m].atoms[a]) {found2 = true;}
	  else if (atom3 == topo->molecules[m].atoms[a]) {found3 = true;}
	  else if (atom4 == topo->molecules[m].atoms[a]) {found4 = true;}

	  if (found1 && found2 && found3 && found4) {break;}
	} // end atom list loop

	if (found1 && found2 && found4) {break;}
      } // end loop over molecules

      // store this bond index number in the molecule's improper list
      topo->molecules[m].improperList.push_back(i);
    } // end loop over impropers
  } // end function buildMoleculeBonding Lists

}

