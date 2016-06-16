#include "PPMWriter.h"

#include "systemutilities.h"
#include "Report.h"

using namespace ProtoMol::Report;

namespace ProtoMol
{
	//_________________________________________________________________PPMWriter

	PPMWriter::PPMWriter(): Writer()
	{
	}

	PPMWriter::PPMWriter(const std::string& filename): Writer(filename)
	{
	}

	bool PPMWriter::write(const PPM& ppm)
	{
		if (!open())
			return false;

		myFile << "P6\n# " << "!ProtoMol (built on " << __DATE__ << " at " << __TIME__ << ") generated this PPM file by " << getUserName() << ". " << myComment << "\n" << ppm.width() << " " << ppm.height() << "\n255\n";
		myFile.write(reinterpret_cast<char*>(ppm.begin()), ppm.size());

		close();
		return !myFile.fail();
	}

	bool PPMWriter::write(const PGM& pgm)
	{
		PPM ppm;
		ppm = pgm;
		return write(ppm);
	}

	PPMWriter& operator<<(PPMWriter& ppmWriter, const PPM& ppm)
	{
		ppmWriter.write(ppm);
		return ppmWriter;
	}

	PPMWriter& operator<<(PPMWriter& ppmWriter, const PGM& pgm)
	{
		ppmWriter.write(pgm);
		return ppmWriter;
	}
}
