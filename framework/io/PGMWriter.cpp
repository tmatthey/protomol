#include "PGMWriter.h"

#include "systemutilities.h"
#include "Report.h"

using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_________________________________________________________________PGMWriter

	PGMWriter::PGMWriter(): Writer()
	{
	}

	PGMWriter::PGMWriter(const std::string& filename): Writer(filename)
	{
	}

	bool PGMWriter::write(const PGM& pgm)
	{
		if (!open())
			return false;

		myFile << "P5\n# " << "!ProtoMol (built on " << __DATE__ << " at " << __TIME__ << ") generated this PGM file by " << getUserName() << ". " << myComment << "\n" << pgm.width() << " " << pgm.height() << "\n255\n";
		myFile.write(reinterpret_cast<char*>(pgm.begin()), pgm.size());

		close();
		return !myFile.fail();
	}

	bool PGMWriter::write(const PPM& ppm)
	{
		PGM pgm;
		pgm = ppm;
		return write(pgm);
	}

	PGMWriter& operator<<(PGMWriter& pgmWriter, const PGM& pgm)
	{
		pgmWriter.write(pgm);
		return pgmWriter;
	}

	PGMWriter& operator<<(PGMWriter& pgmWriter, const PPM& ppm)
	{
		pgmWriter.write(ppm);
		return pgmWriter;
	}
}
