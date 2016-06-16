#include "PSFReader.h"

#include "stringutilities.h"
#include "Report.h"
#include "mathutilities.h"

using std::string;
using std::vector;
using std::endl;
using std::find;
using std::stringstream;
using namespace ProtoMol::Report;

//#define DEBUG_PSF

namespace ProtoMol
{
	//_________________________________________________________________PSFReader

	PSFReader::PSFReader():
		Reader(), myPSF(NULL)
	{
	}

	PSFReader::PSFReader(const std::string& filename):
		Reader(filename), myPSF(NULL)
	{
	}


	PSFReader::~PSFReader()
	{
		if (myPSF != NULL)
			delete myPSF;
	}

	bool PSFReader::tryFormat()
	{
		if (!open())
			return false;
		string psfHead;
		myFile >> psfHead;
		return (myFile.good() && equalNocase("PSF", psfHead));
	}

	bool PSFReader::read()
	{
		if (myPSF == NULL)
			myPSF = new PSF();
		return read(*myPSF);
	}

	bool PSFReader::read(PSF& psf)
	{
		if (!tryFormat())
			return false;
		if (!open())
			return false;
		psf.clear();

		// Find header
		string psfHead;
		myFile >> psfHead;
		if (!equalNocase("PSF", psfHead))
		{
			myFile.setstate(std::ios::failbit);
			close();
			return false;
		}

		while (!myFile.eof())
		{
			string line(removeBeginEndBlanks(getline()));

			// Exit if nothing more to read
			if (line.empty() && myFile.eof())
			{
				close();
				return true;
			}

			// Find '!', otherwise continue with next line
			if (std::find(line.begin(), line.end(), '!') == line.end())
				continue;

			// Parse header
			stringstream ss(line);
			vector<string> header;
			int index = -1;
			string str;
			while (ss >> str)
			{
				if (str[0] == '!' && index < 0)
					index = header.size();
				header.push_back(str);
			}
			if (myFile.fail() || index < 1 || index > 2 || header.size() < 2)
			{
				report << recoverable
					<< " PSF file \'" << myFilename
					<< "\' has corrupt record header." << endl
					<< "\'" << line << "\'"
					<< endr;
				myFile.setstate(std::ios::failbit);
				return false;
			}
			string keyword = header[index];
			int numrecords = toInt(header[0]);

			// Branch
			if (index == 1)
			{
				// Normal case: int,keyword

				if (!isInt(header[0]))
				{
					report << recoverable
						<< " PSF file \'" << myFilename
						<< "\' has corrupt record header, number of records should be an integer." << endl
						<< "\'" << line << "\'"
						<< endr;
					myFile.setstate(std::ios::failbit);
					close();
					return false;
				}

				if (equalStartNocase("!NTITLE", keyword))
				{ // title
					myComment = "";
					for (int counter = 0; counter < numrecords; ++counter)
					{
						line = removeBeginEndBlanks(getline()); // lines of title + one empty line
						if (!line.empty())
							myComment += (myComment.empty() ? "" : "\n") + line;
					}
					continue;
				}
				else if (equalStartNocase("!NATOM", keyword))
				{
					for (int counter = 1; counter <= numrecords; ++counter)
					{
						PSF::Atom temp_atom;
						myFile >> temp_atom.number // read atom number
							>> temp_atom.seg_id // read segment identifier
							>> str // read in residue sequence
							>> temp_atom.residue_name // read in residue name
							>> temp_atom.atom_name // read in atom name
							>> temp_atom.atom_type // read in name of the second atom
							>> temp_atom.charge // read in charge
							>> temp_atom.mass // read in mass
							>> temp_atom.identity; // atom's chemical identity (used only by iSGProtomol)
						temp_atom.residue_sequence = toInt(str);
						if (!isInt(str))
						{
							report << recoverable
								<< "[PSF::read] Expecting a number for residue sequence. "
								<< "I do not know what to do with \'" << str << "\'." << endr;
						}
						psf.atoms.push_back(temp_atom);
						if (myFile.fail() || (myFile.eof() && counter < numrecords))
						{
							myFile.setstate(std::ios::failbit);
							close();
							return false;
						}
					}
					continue;
				}
				else if (equalStartNocase("!NBOND", keyword))
				{
					for (int counter = 1; counter <= numrecords; ++counter)
					{
						PSF::Bond temp_bond;
						temp_bond.number = counter;
						myFile >> temp_bond.atom1;
						myFile >> temp_bond.atom2;
						psf.bonds.push_back(temp_bond);
						if (myFile.fail() || (myFile.eof() && counter < numrecords))
						{
							myFile.setstate(std::ios::failbit);
							close();
							return false;
						}
					}
					continue;
				}
				else if (equalStartNocase("!NTHETA", keyword))
				{
					for (int counter = 1; counter <= numrecords; ++counter)
					{
						PSF::Angle temp_angle;
						temp_angle.number = counter;
						myFile >> temp_angle.atom1;
						myFile >> temp_angle.atom2;
						myFile >> temp_angle.atom3;
						psf.angles.push_back(temp_angle);
						if (myFile.fail() || (myFile.eof() && counter < numrecords))
						{
							myFile.setstate(std::ios::failbit);
							close();
							return false;
						}
					}
					continue;
				}
				else if (equalStartNocase("!NPHI", keyword))
				{
					for (int counter = 1; counter <= numrecords; ++counter)
					{
						PSF::Dihedral temp_dihedral;
						temp_dihedral.number = counter;
						myFile >> temp_dihedral.atom1;
						myFile >> temp_dihedral.atom2;
						myFile >> temp_dihedral.atom3;
						myFile >> temp_dihedral.atom4;
						psf.dihedrals.push_back(temp_dihedral);
						if (myFile.fail() || (myFile.eof() && counter < numrecords))
						{
							myFile.setstate(std::ios::failbit);
							close();
							return false;
						}
					}
					continue;
				}
				else if (equalStartNocase("!NIMPHI", keyword))
				{
					for (int counter = 1; counter <= numrecords; ++counter)
					{
						PSF::Improper temp_improper;
						temp_improper.number = counter;
						myFile >> temp_improper.atom1;
						myFile >> temp_improper.atom2;
						myFile >> temp_improper.atom3;
						myFile >> temp_improper.atom4;
						psf.impropers.push_back(temp_improper);
						if (myFile.fail() || (myFile.eof() && counter < numrecords))
						{
							myFile.setstate(std::ios::failbit);
							close();
							return false;
						}
					}
					continue;
				}
				else if (equalStartNocase("!NDON", keyword))
				{
					for (int counter = 1; counter <= numrecords; ++counter)
					{
						PSF::Donor temp_donor;
						temp_donor.number = counter;
						myFile >> temp_donor.atom1;
						myFile >> temp_donor.atom2;
						psf.donors.push_back(temp_donor);
						if (myFile.fail() || (myFile.eof() && counter < numrecords))
						{
							myFile.setstate(std::ios::failbit);
							close();
							return false;
						}
					}
					continue;
				}
				else if (equalStartNocase("!NACC", keyword))
				{
					for (int counter = 1; counter <= numrecords; ++counter)
					{
						PSF::Acceptor temp_acceptor;
						temp_acceptor.number = counter;
						myFile >> temp_acceptor.atom1;
						myFile >> temp_acceptor.atom2;
						psf.acceptors.push_back(temp_acceptor);
						if (myFile.fail() || (myFile.eof() && counter < numrecords))
						{
							myFile.setstate(std::ios::failbit);
							close();
							return false;
						}
					}
					continue;
				}
				else if (equalStartNocase("!NNB", keyword))
				{
					for (int counter = 1; counter <= numrecords; ++counter)
					{
						PSF::Nonbonded temp_nonbonded;
						temp_nonbonded.number = counter;
						myFile >> temp_nonbonded.atom1;
						psf.nonbondeds.push_back(temp_nonbonded);
						if (myFile.fail() || (myFile.eof() && counter < numrecords))
						{
							report << recoverable << "[PSF::read] Expecting " << numrecords << " NNB, read " << counter - 1 << ", reached end of file." << endr;
							close();
							return false;
						}
					}
					continue;
				}
			}

			if (index == 2)
			{
				if (equalStartNocase("!NGRP", keyword))
				{
					for (int counter = 1; counter <= numrecords; ++counter)
					{
						PSF::Ngrp temp_ngrp;
						temp_ngrp.number = counter;
						myFile >> temp_ngrp.atom1;
						myFile >> temp_ngrp.atom2;
						myFile >> temp_ngrp.atom3;
						psf.ngrp.push_back(temp_ngrp);
						if (myFile.fail() || (myFile.eof() && counter < numrecords))
						{
							report << recoverable << "[PSF::read] Expecting " << numrecords << " NGRP, read " << counter - 1 << ", reached end of file." << endr;
							close();
							return false;
						}
					}
					continue;
				}
			}

			report << recoverable << "[PSF::read] Record " << keyword << " with " << numrecords << " entries not recognized." << endr;
		}

		close();
		return !myFile.fail();
	}

	PSF* PSFReader::orphanPSF()
	{
		PSF* tmp = myPSF;
		myPSF = NULL;
		return tmp;
	}

	PSFReader& operator>>(PSFReader& psfReader, PSF& psf)
	{
		psfReader.read(psf);
		return psfReader;
	}
}
