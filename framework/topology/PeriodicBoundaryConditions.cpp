#include "PeriodicBoundaryConditions.h"
#include "mathutilities.h"
#include "Report.h"
#include <algorithm>
using namespace ProtoMol::Report;
using std::string;
using std::vector;
using std::sort;

namespace ProtoMol
{
	//_____________________________________________________________________ normVector3DOp
	static bool normVector3DOp(const Vector3D& v1, const Vector3D& v2)
	{
		return (v1.normSquared() < v2.normSquared());
	}

	//_________________________________________________________________ PeriodicBoundaryConditions
	const string PeriodicBoundaryConditions::keyword("Periodic");

	PeriodicBoundaryConditions::PeriodicBoundaryConditions(): myE1(0, 0, 0), myE2(0, 0, 0), myE3(0, 0, 0),
	                                                          myE1r(0, 0, 0), myE2r(0, 0, 0), myE3r(0, 0, 0),
	                                                          myOrigin(0, 0, 0), myMin(0, 0, 0), myMax(0, 0, 0),
	                                                          myV(0), myOrthogonal(true)
	{
	}

	PeriodicBoundaryConditions::PeriodicBoundaryConditions(const Vector3D& e1, const Vector3D& e2,
	                                                       const Vector3D& e3, const Vector3D& origin)
	{
		set(e1, e2, e3, origin);
	}

	void PeriodicBoundaryConditions::set(const Vector3D& e1, const Vector3D& e2,
	                                     const Vector3D& e3, const Vector3D& origin)
	{
		myE1 = e1;
		myE2 = e2;
		myE3 = e3;
		myOrigin = origin;

		myV = fabs((e1.cross(e2)).dot(e3));
		if (myV < Constant::EPSILON)
		{
			report << error << "[PeriodicBoundaryConditions::set] No volume, aborting." << endr;
		}
		if (myV >= Constant::REAL_INFINITY)
		{
			report << error << "[PeriodicBoundaryConditions::set] Infinite volume, aborting." << endr;
		}

		myOrthogonal = !(e1.y != 0.0 || e1.z != 0.0 || e2.x != 0.0 || e2.z != 0.0 || e3.x != 0.0 || e3.y != 0.0);

		Vector3D a1(e2.cross(e3));
		myE1r = a1 / e1.dot(a1);
		Vector3D a2(e3.cross(e1));
		myE2r = a2 / e2.dot(a2);
		Vector3D a3(e1.cross(e2));
		myE3r = a3 / e3.dot(a3);

		Vector3D a(origin - (e1 + e2 + e3) * 0.5);
		Vector3D b(origin + (e1 + e2 + e3) * 0.5);
		myMin.x = std::min(a.x, b.x);
		myMin.y = std::min(a.y, b.y);
		myMin.z = std::min(a.z, b.z);
		myMax.x = std::max(a.x, b.x);
		myMax.y = std::max(a.y, b.y);
		myMax.z = std::max(a.z, b.z);

		myDX = power<2>(e1.x * 0.5);
		myDY = power<2>(e2.y * 0.5);
		myDZ = power<2>(e3.z * 0.5);
		myD = std::min(myDX, std::min(myDY, myDZ));
		myH = myMax - myMin;
		myH2 = myH * 0.5;
		report << debug(2) << "[PeriodicBoundaryConditions] maximal safe distance=" << myD << endr;
	}


	vector<Vector3D> PeriodicBoundaryConditions::buildLatticeVectors(Real cutoff) const
	{
		vector<Vector3D> lattice;
		if (myV >= Constant::REAL_INFINITY)
			return lattice;
		Vector3D a = myE1 + myE2 + myE3;
		Real lx = fabs(a.x);
		Real ly = fabs(a.y);
		Real lz = fabs(a.z);
		int boundK = (int)floor(cutoff / lx + 1.0);
		int boundL = (int)floor(cutoff / ly + 1.0);
		int boundM = (int)floor(cutoff / lz + 1.0);

		for (int k = -boundK; k <= boundK; k++)
		{
			for (int l = -boundL; l <= boundL; l++)
			{
				for (int m = -boundM; m <= boundM; m++)
				{
					if (k != 0 || l != 0 || m != 0)
					{
						Vector3D coord(myE1 * k + myE2 * l + myE3 * m);
						if (coord.norm() <= 2 * cutoff)
							lattice.push_back(coord);
					}
				}
			}
		}
		sort(lattice.begin(), lattice.end(), normVector3DOp);

		return lattice;
	}

	void PeriodicBoundaryConditions::getParameters(vector<Parameter>& parameters) const
	{
		parameters.push_back(Parameter("cellBasisVector1", Value(myE1, ConstraintValueType::NotZero())));
		parameters.push_back(Parameter("cellBasisVector2", Value(myE2, ConstraintValueType::NotZero())));
		parameters.push_back(Parameter("cellBasisVector3", Value(myE3, ConstraintValueType::NotZero())));
		parameters.push_back(Parameter("cellorigin", Value(myOrigin, Value::undefined)));
	}

	PeriodicBoundaryConditions PeriodicBoundaryConditions::make(string& errMsg, vector<Value> values)
	{
		Vector3D e1, e2, e3, o;
		values[0].get(e1);
		values[1].get(e2);
		values[2].get(e3);
		values[3].get(o);
		if (!values[0].valid() || !values[1].valid() || !values[2].valid() || !values[3].valid())
		{
			errMsg += keyword + " boundary conditions not valid: " + values[0].getString() + "," + values[1].getString() + "," + values[2].getString() + "," + values[3].getString() + ".";
			return PeriodicBoundaryConditions();
		}
		return PeriodicBoundaryConditions(e1, e2, e3, o);
	}

	vector<Parameter> PeriodicBoundaryConditions::getDefaults(const Vector3DBlock& positions) const
	{
		Vector3D a, b;
		Vector3D d(Constant::PERIODIC_BOUNDARY_TOLERANCE / 2.0, Constant::PERIODIC_BOUNDARY_TOLERANCE / 2.0, Constant::PERIODIC_BOUNDARY_TOLERANCE / 2.0);
		positions.boundingbox(a, b);
		a -= d;
		b += d;
		Vector3D c(b - a);
		vector<Parameter> parameters;
		parameters.push_back(Parameter("cellBasisVector1", Value(Vector3D(c.x, 0.0, 0.0), ConstraintValueType::NotZero())));
		parameters.push_back(Parameter("cellBasisVector2", Value(Vector3D(0.0, c.y, 0.0), ConstraintValueType::NotZero())));
		parameters.push_back(Parameter("cellBasisVector3", Value(Vector3D(0.0, 0.0, c.z), ConstraintValueType::NotZero())));
		parameters.push_back(Parameter("cellorigin", Value(Vector3D(a + c * 0.5))));
		return parameters;
	}
}
