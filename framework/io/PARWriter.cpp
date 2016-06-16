#include "PARWriter.h"

#include <iomanip>

#include "Report.h"
#include "stringutilities.h"
#include "systemutilities.h"

using std::string;
using std::endl;
using std::setprecision;
using std::setw;
using std::left;
using std::right;
using std::setiosflags;

using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_________________________________________________________________PARWriter
	const PAR::CharmmTypeEnum PARWriter::defaultCharmmType = PAR::CHARMM28;

	PARWriter::PARWriter(PAR::CharmmTypeEnum charmmType):
		Writer(), myCharmmType(charmmType)
	{
	}

	PARWriter::PARWriter(const string& filename, PAR::CharmmTypeEnum charmmType):
		Writer(filename), myCharmmType(charmmType)
	{
	}

	PARWriter::PARWriter(const char* filename, PAR::CharmmTypeEnum charmmType):
		Writer(string(filename)), myCharmmType(charmmType)
	{
	}

	bool PARWriter::open(const std::string& filename, PAR::CharmmTypeEnum charmmType)
	{
		myCharmmType = charmmType;
		return open(filename);
	}

	bool PARWriter::open(PAR::CharmmTypeEnum charmmType)
	{
		myCharmmType = charmmType;
		return open();
	}

	bool PARWriter::open(const char* filename, PAR::CharmmTypeEnum charmmType)
	{
		return open(std::string(filename), charmmType);
	}

	void PARWriter::setCharmmType(PAR::CharmmTypeEnum charmmType)
	{
		myCharmmType = charmmType;
	}

	PAR::CharmmTypeEnum PARWriter::getCharmmType() const
	{
		return myCharmmType;
	}


	bool PARWriter::write(const PAR& par)
	{
		if (!open())
			return false;

		if (myCharmmType == PAR::UNDEFINED)
		{
			report << hint << "[PARWriter::write] Charmm type undefined, using default Charmm type " << ((defaultCharmmType == PAR::CHARMM28) ? "28" : "19") << "." << endr;
			myCharmmType = defaultCharmmType;
		}

		// Comment
		myFile << "!ProtoMol (built on " << __DATE__ << " at " << __TIME__ << ") generated this PAR file by " << getUserName() << "." << endl;
		if (!myComment.empty())
			myFile << "!" << myComment << endl;
		myFile << endl;

		if (myCharmmType == PAR::CHARMM28)
		{
			// Charmm28
			myFile << setiosflags(std::ios::showpoint | std::ios::fixed);
			unsigned int count = par.bonds.size();
			if (count > 0)
			{
				myFile << "BONDS" << endl;
				for (unsigned int i = 0; i < count; ++i)
				{
					myFile << setw(5) << left
						<< par.bonds[i].atom1
						<< " " << setw(5)
						<< par.bonds[i].atom2
						<< " " << setw(10) << setprecision(4) << right
						<< par.bonds[i].forceConstant
						<< " " << setw(10) << setprecision(4)
						<< par.bonds[i].distance
						<< endl;
				}
				myFile << endl;
			}
			count = par.angles.size();
			if (count > 0)
			{
				myFile << "ANGLES" << endl;
				for (unsigned int i = 0; i < count; ++i)
				{
					myFile << setw(5) << left
						<< par.angles[i].atom1
						<< " " << setw(5)
						<< par.angles[i].atom2
						<< " " << setw(5)
						<< par.angles[i].atom3
						<< " " << setw(10) << setprecision(4) << right
						<< par.angles[i].forceConstant
						<< " " << setw(10) << setprecision(4)
						<< par.angles[i].angleval;
					if (par.angles[i].ub_flag)
					{
						myFile << " " << setw(10) << setprecision(4)
							<< par.angles[i].k_ub
							<< " " << setw(10) << setprecision(6)
							<< par.angles[i].r_ub;
					}
					myFile << endl;
				}
				myFile << endl;
			}
			count = par.dihedrals.size();
			if (count > 0)
			{
				myFile << "DIHEDRALS" << endl;
				for (unsigned int i = 0; i < count; ++i)
				{
					for (unsigned int j = 0; j < par.dihedrals[i].forceConstant.size(); ++j)
					{
						myFile << setw(5) << left
							<< par.dihedrals[i].atom1
							<< " " << setw(5)
							<< par.dihedrals[i].atom2
							<< " " << setw(5)
							<< par.dihedrals[i].atom3
							<< " " << setw(5)
							<< par.dihedrals[i].atom4
							<< " " << setw(10) << setprecision(4) << right
							<< par.dihedrals[i].forceConstant[j]
							<< " " << setw(3) << setprecision(0)
							<< par.dihedrals[i].periodicity[j]
							<< " " << setw(10) << setprecision(4)
							<< par.dihedrals[i].phaseShift[j]
							<< endl;
					}
				}
				myFile << endl;
			}
			count = par.impropers.size();
			if (count > 0)
			{
				myFile << "IMPROPER" << endl;
				for (unsigned int i = 0; i < count; ++i)
				{
					myFile << setw(5) << left
						<< par.impropers[i].atom1
						<< " " << setw(5)
						<< par.impropers[i].atom2
						<< " " << setw(5)
						<< par.impropers[i].atom3
						<< " " << setw(5)
						<< par.impropers[i].atom4
						<< " " << setw(10) << setprecision(4) << right
						<< par.impropers[i].forceConstant
						<< " " << setw(3) << setprecision(0)
						<< par.impropers[i].periodicity
						<< " " << setw(10) << setprecision(4)
						<< par.impropers[i].phaseShift
						<< endl;
				}
				myFile << endl;
			}
			count = par.nonbondeds.size();
			if (count > 0)
			{
				myFile << "NONBONDED" << endl;
				for (unsigned int i = 0; i < count; ++i)
				{
					myFile << setw(5) << left
						<< par.nonbondeds[i].atom
						<< " " << setw(5) << setprecision(6) << right
						<< par.nonbondeds[i].polarizability
						<< " " << setw(10) << setprecision(6)
						<< par.nonbondeds[i].epsilon
						<< " " << setw(10) << setprecision(6)
						<< par.nonbondeds[i].sigma;
					if (par.nonbondeds[i].vdw)
					{
						myFile << " " << setw(10) << setprecision(6)
							<< par.nonbondeds[i].polarizability2
							<< " " << setw(10) << setprecision(6)
							<< par.nonbondeds[i].epsilon14
							<< " " << setw(10) << setprecision(6)
							<< par.nonbondeds[i].sigma14;
					}
					myFile << endl;
				}
				myFile << endl;
			}
			count = par.nbfixs.size();
			if (count > 0)
			{
				myFile << "NBFIX" << endl;
				for (unsigned int i = 0; i < count; ++i)
				{
					myFile << setw(5) << left
						<< par.nbfixs[i].atom1
						<< " " << setw(5)
						<< par.nbfixs[i].atom2
						<< " " << setw(10) << setprecision(4) << right
						<< par.nbfixs[i].a
						<< " " << setw(10) << setprecision(4)
						<< par.nbfixs[i].b;
					if (par.nbfixs[i].a != par.nbfixs[i].a14 || par.nbfixs[i].b != par.nbfixs[i].b14)
					{
						myFile << " " << setw(10) << setprecision(4)
							<< par.nbfixs[i].a14
							<< " " << setw(10) << setprecision(4)
							<< par.nbfixs[i].b14;
					}
					myFile << endl;
				}
				myFile << endl;
			}
			count = par.hbonds.size();
			if (count > 0)
			{
				myFile << "HBOND" << endl;
				for (unsigned int i = 0; i < count; ++i)
				{
					myFile << setw(5) << left
						<< par.hbonds[i].atom1
						<< " " << setw(5)
						<< par.hbonds[i].atom2
						<< " " << setw(10) << setprecision(4) << right
						<< par.hbonds[i].emin
						<< " " << setw(10) << setprecision(4)
						<< par.hbonds[i].rmin
						<< endl;
				}
				myFile << endl;
			}
		}
		else
		{
			// Charmm19
			myFile << setiosflags(std::ios::showpoint | std::ios::fixed);
			unsigned int count = par.bonds.size();
			if (count > 0)
			{
				for (unsigned int i = 0; i < count; ++i)
				{
					myFile << "BOND " << setw(5) << left
						<< par.bonds[i].atom1
						<< " " << setw(5)
						<< par.bonds[i].atom2
						<< " " << setw(10) << setprecision(4) << right
						<< par.bonds[i].forceConstant
						<< " " << setw(10) << setprecision(4)
						<< par.bonds[i].distance
						<< endl;
				}
				myFile << endl;
			}
			count = par.angles.size();
			if (count > 0)
			{
				for (unsigned int i = 0; i < count; ++i)
				{
					myFile << "ANGLE " << setw(5) << left
						<< par.angles[i].atom1
						<< " " << setw(5)
						<< par.angles[i].atom2
						<< " " << setw(5)
						<< par.angles[i].atom3
						<< " " << setw(10) << setprecision(4) << right
						<< par.angles[i].forceConstant
						<< " " << setw(10) << setprecision(4)
						<< par.angles[i].angleval;
					if (par.angles[i].ub_flag)
					{
						myFile << " UB " << setw(10) << setprecision(4)
							<< par.angles[i].k_ub
							<< " " << setw(10) << setprecision(6)
							<< par.angles[i].r_ub;
					}
					myFile << endl;
				}
				myFile << endl;
			}
			count = par.dihedrals.size();
			if (count > 0)
			{
				for (unsigned int i = 0; i < count; ++i)
				{
					myFile << "DIHEDRAL " << setw(5) << left
						<< par.dihedrals[i].atom1
						<< " " << setw(5)
						<< par.dihedrals[i].atom2
						<< " " << setw(5)
						<< par.dihedrals[i].atom3
						<< " " << setw(5)
						<< par.dihedrals[i].atom4;
					if (par.dihedrals[i].forceConstant.size() > 1)
					{
						myFile << " MULTIPLE= " << setw(2) << setprecision(0) << right
							<< par.dihedrals[i].forceConstant.size();
					}
					myFile << " " << setw(10) << setprecision(4) << right
						<< par.dihedrals[i].forceConstant[0]
						<< " " << setw(3) << setprecision(0)
						<< par.dihedrals[i].periodicity[0]
						<< " " << setw(10) << setprecision(4)
						<< par.dihedrals[i].phaseShift[0]
						<< endl;
					for (unsigned int j = 1; j < par.dihedrals[i].forceConstant.size(); ++j)
					{
						myFile << " " << setw(55) << setprecision(4) << right
							<< par.dihedrals[i].forceConstant[j]
							<< " " << setw(3) << setprecision(0)
							<< par.dihedrals[i].periodicity[j]
							<< " " << setw(10) << setprecision(4)
							<< par.dihedrals[i].phaseShift[j]
							<< endl;
					}
				}
				myFile << endl;
			}
			count = par.impropers.size();
			if (count > 0)
			{
				for (unsigned int i = 0; i < count; ++i)
				{
					myFile << "IMPROPER " << setw(5) << left
						<< par.impropers[i].atom1
						<< " " << setw(5)
						<< par.impropers[i].atom2
						<< " " << setw(5)
						<< par.impropers[i].atom3
						<< " " << setw(5)
						<< par.impropers[i].atom4
						<< " " << setw(10) << setprecision(4) << right
						<< par.impropers[i].forceConstant
						<< " " << setw(3) << setprecision(0)
						<< par.impropers[i].periodicity
						<< " " << setw(10) << setprecision(4)
						<< par.impropers[i].phaseShift
						<< endl;
				}
				myFile << endl;
			}
			count = par.nonbondeds.size();
			if (count > 0)
			{
				for (unsigned int i = 0; i < count; ++i)
				{
					myFile << "NONBONDED " << setw(5) << left
						<< par.nonbondeds[i].atom
						<< " " << setw(5) << setprecision(6) << right
						<< par.nonbondeds[i].epsilon
						<< " " << setw(10) << setprecision(6)
						<< par.nonbondeds[i].sigma * PAR::Nonbonded::SIGMA_CHARMM28_TO_CHARMM19
						<< " " << setw(10) << setprecision(6)
						<< (par.nonbondeds[i].vdw ? par.nonbondeds[i].epsilon14 : par.nonbondeds[i].epsilon)
						<< " " << setw(10) << setprecision(6)
						<< (par.nonbondeds[i].vdw ? par.nonbondeds[i].sigma14 : par.nonbondeds[i].sigma) * PAR::Nonbonded::SIGMA_CHARMM28_TO_CHARMM19
						<< endl;
				}
				myFile << endl;
			}
			count = par.nbfixs.size();
			if (count > 0)
			{
				for (unsigned int i = 0; i < count; ++i)
				{
					myFile << "NBFIX " << setw(5) << left
						<< par.nbfixs[i].atom1
						<< " " << setw(5)
						<< par.nbfixs[i].atom2
						<< " " << setw(10) << setprecision(4) << right
						<< par.nbfixs[i].a
						<< " " << setw(10) << setprecision(4)
						<< par.nbfixs[i].b
						<< " " << setw(10) << setprecision(4)
						<< par.nbfixs[i].a14
						<< " " << setw(10) << setprecision(4)
						<< par.nbfixs[i].b14
						<< endl;
				}
				myFile << endl;
			}
			count = par.hbonds.size();
			if (count > 0)
			{
				myFile << "HBOND" << endl;
				for (unsigned int i = 0; i < count; ++i)
				{
					myFile << setw(5) << left
						<< par.hbonds[i].atom1
						<< " " << setw(5)
						<< par.hbonds[i].atom2
						<< " " << setw(10) << setprecision(4) << right
						<< par.hbonds[i].emin
						<< " " << setw(10) << setprecision(4)
						<< par.hbonds[i].rmin
						<< endl;
				}
				myFile << endl;
			}
		}

		close();
		return !myFile.fail();
	}

	PARWriter& operator<<(PARWriter& parWriter, const PAR& par)
	{
		parWriter.write(par);
		return parWriter;
	}
}
