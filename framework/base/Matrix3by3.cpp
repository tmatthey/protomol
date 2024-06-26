#include "Matrix3by3.h"
#include "mathutilities.h"
#include "Report.h"
using namespace ProtoMol::Report;

namespace ProtoMol
{
	//________________________________________________________________ Matrix3by3


	const Matrix3by3 Matrix3by3::ourIdentity = Matrix3by3(1, 0, 0, 0, 1, 0, 0, 0, 1);

	Matrix3by3::Matrix3by3()
	{
		// zeroMatrix(); 
		// for performance reasons, an empty Matrix3by3 object will contain
		// anything left in memory (not zeroed any more). user should avoid
		// using x += y type operations until the matrices (x) have been
		// properly initialized. A direct assian operation, x = y, will
		// replace whatever left in x with the new values.
	}

	Matrix3by3::Matrix3by3(Real x00, Real x01, Real x02,
	                       Real x10, Real x11, Real x12,
	                       Real x20, Real x21, Real x22)
	{
		m00 = x00;
		m01 = x01;
		m02 = x02;
		m10 = x10;
		m11 = x11;
		m12 = x12;
		m20 = x20;
		m21 = x21;
		m22 = x22;
	}

	Matrix3by3::Matrix3by3(float mat[9])
	{
		m00 = *mat++;
		m01 = *mat++;
		m02 = *mat++;
		m10 = *mat++;
		m11 = *mat++;
		m12 = *mat++;
		m20 = *mat++;
		m21 = *mat++;
		m22 = *mat++;
	}

	Matrix3by3::Matrix3by3(double mat[9])
	{
		m00 = *mat++;
		m01 = *mat++;
		m02 = *mat++;
		m10 = *mat++;
		m11 = *mat++;
		m12 = *mat++;
		m20 = *mat++;
		m21 = *mat++;
		m22 = *mat++;
	}

	//construct a matrix that is the outer products of two coordinates
	Matrix3by3::Matrix3by3(const Vector3D& a, const Vector3D& b)
	{
		(*this)(a.x * b.x, a.x * b.y, a.x * b.z,
		        a.y * b.x, a.y * b.y, a.y * b.z,
		        a.z * b.x, a.z * b.y, a.z * b.z);
	}

	Matrix3by3::Matrix3by3(const Vector3D& v1, const Vector3D& v2,
	                       const Vector3D& v3)
	{
		m00 = v1.x;
		m01 = v1.y;
		m02 = v1.z;
		m10 = v2.x;
		m11 = v2.y;
		m12 = v2.z;
		m20 = v3.x;
		m21 = v3.y;
		m22 = v3.z;
	}

	void Matrix3by3::identity()
	{
		m00 = 1;
		m01 = 0;
		m02 = 0;
		m10 = 0;
		m11 = 1;
		m12 = 0;
		m20 = 0;
		m21 = 0;
		m22 = 1;
	}

	void Matrix3by3::zeroMatrix()
	{
		m00 = 0;
		m01 = 0;
		m02 = 0;
		m10 = 0;
		m11 = 0;
		m12 = 0;
		m20 = 0;
		m21 = 0;
		m22 = 0;
	}

	bool Matrix3by3::zero() const
	{
		return *this == Matrix3by3(0, 0, 0, 0, 0, 0, 0, 0, 0);
	}

	Real Matrix3by3::operator()(int i, int j) const
	{
		if (i == 0)
		{
			if (j == 0) return m00;
			else if (j == 1) return m01;
			else if (j == 2) return m02;
		}
		else if (i == 1)
		{
			if (j == 0) return m10;
			else if (j == 1) return m11;
			else if (j == 2) return m12;
		}
		else if (i == 2)
		{
			if (j == 0) return m20;
			else if (j == 1) return m21;
			else if (j == 2) return m22;
		}
		report << recoverable << "[Matrix3by3::operator()] index out of range." << endr;
		return 0;
	}

	void Matrix3by3::operator()(int i, int j, Real x)
	{
		if (i == 0)
		{
			if (j == 0)
			{
				m00 = x;
				return;
			}
			else if (j == 1)
			{
				m01 = x;
				return;
			}
			else if (j == 2)
			{
				m02 = x;
				return;
			}
		}
		else if (i == 1)
		{
			if (j == 0)
			{
				m10 = x;
				return;
			}
			else if (j == 1)
			{
				m11 = x;
				return;
			}
			else if (j == 2)
			{
				m12 = x;
				return;
			}
		}
		else if (i == 2)
		{
			if (j == 0)
			{
				m20 = x;
				return;
			}
			else if (j == 1)
			{
				m21 = x;
				return;
			}
			else if (j == 2)
			{
				m22 = x;
				return;
			}
		}
		report << recoverable << "[Matrix3by3::operator()] index out of range." << endr;
	}


	void Matrix3by3::operator()(Real x00, Real x01, Real x02,
	                            Real x10, Real x11, Real x12,
	                            Real x20, Real x21, Real x22)
	{
		m00 = x00;
		m01 = x01;
		m02 = x02;
		m10 = x10;
		m11 = x11;
		m12 = x12;
		m20 = x20;
		m21 = x21;
		m22 = x22;
	}

	bool Matrix3by3::operator==(const Matrix3by3& tm) const
	{
		return (equal(m00, tm.m00) && equal(m01, tm.m01) && equal(m02, tm.m02) &&
			equal(m10, tm.m10) && equal(m11, tm.m11) && equal(m12, tm.m12) &&
			equal(m20, tm.m20) && equal(m21, tm.m21) && equal(m22, tm.m22));
	}

	bool Matrix3by3::operator!=(const Matrix3by3& tm) const
	{
		return !(*this == tm);
	}

	Matrix3by3& Matrix3by3::operator*=(const Real tm)
	{
		m00 *= tm;
		m01 *= tm;
		m02 *= tm;
		m10 *= tm;
		m11 *= tm;
		m12 *= tm;
		m20 *= tm;
		m21 *= tm;
		m22 *= tm;
		return *this;
	}

	Matrix3by3& Matrix3by3::operator/=(const Real tm)
	{
		m00 /= tm;
		m01 /= tm;
		m02 /= tm;
		m10 /= tm;
		m11 /= tm;
		m12 /= tm;
		m20 /= tm;
		m21 /= tm;
		m22 /= tm;
		return *this;
	}


	Matrix3by3& Matrix3by3::operator*=(const Matrix3by3& tm)
	{
		Matrix3by3 temp(*this);

		m00 = temp.m00 * tm.m00 + temp.m01 * tm.m10 + temp.m02 * tm.m20;
		m01 = temp.m00 * tm.m01 + temp.m01 * tm.m11 + temp.m02 * tm.m21;
		m02 = temp.m00 * tm.m02 + temp.m01 * tm.m12 + temp.m02 * tm.m22;

		m10 = temp.m10 * tm.m00 + temp.m11 * tm.m10 + temp.m12 * tm.m20;
		m11 = temp.m10 * tm.m01 + temp.m11 * tm.m11 + temp.m12 * tm.m21;
		m12 = temp.m10 * tm.m02 + temp.m11 * tm.m12 + temp.m12 * tm.m22;

		m20 = temp.m20 * tm.m00 + temp.m21 * tm.m10 + temp.m22 * tm.m20;
		m21 = temp.m20 * tm.m01 + temp.m21 * tm.m11 + temp.m22 * tm.m21;
		m22 = temp.m20 * tm.m02 + temp.m21 * tm.m12 + temp.m22 * tm.m22;

		return *this;
	}


	Matrix3by3 Matrix3by3::operator*(const Real tm) const
	{
		Matrix3by3 res;
		res.m00 = m00 * tm;
		res.m01 = m01 * tm;
		res.m02 = m02 * tm;

		res.m10 = m10 * tm;
		res.m11 = m11 * tm;
		res.m12 = m12 * tm;

		res.m20 = m20 * tm;
		res.m21 = m21 * tm;
		res.m22 = m22 * tm;

		return res;
	}


	Matrix3by3 Matrix3by3::operator/(const Real tm) const
	{
		Matrix3by3 res;
		res.m00 = m00 / tm;
		res.m01 = m01 / tm;
		res.m02 = m02 / tm;

		res.m10 = m10 / tm;
		res.m11 = m11 / tm;
		res.m12 = m12 / tm;

		res.m20 = m20 / tm;
		res.m21 = m21 / tm;
		res.m22 = m22 / tm;
		return res;
	}


	Matrix3by3 Matrix3by3::operator*(const Matrix3by3& tm) const
	{
		Matrix3by3 res;
		res.m00 = m00 * tm.m00 + m01 * tm.m10 + m02 * tm.m20;
		res.m01 = m00 * tm.m01 + m01 * tm.m11 + m02 * tm.m21;
		res.m02 = m00 * tm.m02 + m01 * tm.m12 + m02 * tm.m22;

		res.m10 = m10 * tm.m00 + m11 * tm.m10 + m12 * tm.m20;
		res.m11 = m10 * tm.m01 + m11 * tm.m11 + m12 * tm.m21;
		res.m12 = m10 * tm.m02 + m11 * tm.m12 + m12 * tm.m22;

		res.m20 = m20 * tm.m00 + m21 * tm.m10 + m22 * tm.m20;
		res.m21 = m20 * tm.m01 + m21 * tm.m11 + m22 * tm.m21;
		res.m22 = m20 * tm.m02 + m21 * tm.m12 + m22 * tm.m22;
		return res;
	}


	Vector3D Matrix3by3::operator*(const Vector3D& tm) const
	{
		Vector3D res;
		res.x = m00 * tm.x + m01 * tm.y + m02 * tm.z;
		res.y = m10 * tm.x + m11 * tm.y + m12 * tm.z;
		res.z = m20 * tm.x + m21 * tm.y + m22 * tm.z;
		return res;
	}


	Matrix3by3& Matrix3by3::operator-=(const Matrix3by3& tm)
	{
		m00 -= tm.m00;
		m01 -= tm.m01;
		m02 -= tm.m02;

		m10 -= tm.m10;
		m11 -= tm.m11;
		m12 -= tm.m12;

		m20 -= tm.m20;
		m21 -= tm.m21;
		m22 -= tm.m22;
		return *this;
	}

	Matrix3by3 Matrix3by3::operator-(const Matrix3by3& tm) const
	{
		Matrix3by3 res;

		res.m00 = m00 - tm.m00;
		res.m01 = m01 - tm.m01;
		res.m02 = m02 - tm.m02;

		res.m10 = m10 - tm.m10;
		res.m11 = m11 - tm.m11;
		res.m12 = m12 - tm.m12;

		res.m20 = m20 - tm.m20;
		res.m21 = m21 - tm.m21;
		res.m22 = m22 - tm.m22;
		return res;
	}


	Matrix3by3 Matrix3by3::operator-(void) const
	{
		Matrix3by3 res;

		res.m00 = -m00;
		res.m01 = -m01;
		res.m02 = -m02;
		res.m10 = -m10;
		res.m11 = -m11;
		res.m12 = -m12;
		res.m20 = -m20;
		res.m21 = -m21;
		res.m22 = -m22;
		return res;
	}

	Matrix3by3 Matrix3by3::transposed() const
	{
		Matrix3by3 res;

		res.m00 = m00;
		res.m01 = m10;
		res.m02 = m20;

		res.m10 = m01;
		res.m11 = m11;
		res.m12 = m21;

		res.m20 = m02;
		res.m21 = m12;
		res.m22 = m22;
		return res;
	}

	void Matrix3by3::transpose(const Matrix3by3& tm)
	{
		// NB: make sure that *this == tm also works
		Real a = tm.m10;
		m10 = tm.m01;
		m01 = a;

		Real b = tm.m20;
		m20 = tm.m02;
		m02 = b;

		Real c = tm.m21;
		m21 = tm.m12;
		m12 = c;

		m00 = tm.m00;
		m11 = tm.m11;
		m22 = tm.m22;
	}

	void Matrix3by3::transpose()
	{
		Real a = m10;
		m10 = m01;
		m01 = a;

		Real b = m20;
		m20 = m02;
		m02 = b;

		Real c = m21;
		m21 = m12;
		m12 = c;
	}


	Matrix3by3& Matrix3by3::operator+=(const Matrix3by3& tm)
	{
		m00 += tm.m00;
		m01 += tm.m01;
		m02 += tm.m02;

		m10 += tm.m10;
		m11 += tm.m11;
		m12 += tm.m12;

		m20 += tm.m20;
		m21 += tm.m21;
		m22 += tm.m22;
		return *this;
	}

	Matrix3by3 Matrix3by3::operator+(const Matrix3by3& tm) const
	{
		Matrix3by3 res;

		res.m00 = m00 + tm.m00;
		res.m01 = m01 + tm.m01;
		res.m02 = m02 + tm.m02;

		res.m10 = m10 + tm.m10;
		res.m11 = m11 + tm.m11;
		res.m12 = m12 + tm.m12;

		res.m20 = m20 + tm.m20;
		res.m21 = m21 + tm.m21;
		res.m22 = m22 + tm.m22;
		return res;
	}


	Real Matrix3by3::det()
	{
		Real a00 = m11 * m22 - m12 * m21;
		Real a10 = m10 * m22 - m12 * m20;
		Real a20 = m10 * m21 - m11 * m20;
		return m00 * a00 - m01 * a10 + m02 * a20;
	}

	bool Matrix3by3::invert()
	{
		Matrix3by3 tmp;
		Real det;

		//
		// Calculate the determinant of submatrix A (optimized version:
		// don,t just compute the determinant of A)
		//
		tmp.m00 = m11 * m22 - m12 * m21;
		tmp.m10 = m10 * m22 - m12 * m20;
		tmp.m20 = m10 * m21 - m11 * m20;

		tmp.m01 = m01 * m22 - m02 * m21;
		tmp.m11 = m00 * m22 - m02 * m20;
		tmp.m21 = m00 * m21 - m01 * m20;

		tmp.m02 = m01 * m12 - m02 * m11;
		tmp.m12 = m00 * m12 - m02 * m10;
		tmp.m22 = m00 * m11 - m01 * m10;

		det = m00 * tmp.m00 - m01 * tmp.m10 + m02 * tmp.m20;

		//
		// singular matrix ?
		//
		if (fabs(det) < Constant::EPSILON * Constant::EPSILON)
			return false;

		det = 1 / det;

		//
		// inverse(A) = adj(A)/det(A)
		//
		tmp.m00 *= det;
		tmp.m02 *= det;
		tmp.m11 *= det;
		tmp.m20 *= det;
		tmp.m22 *= det;

		det = -det;

		tmp.m01 *= det;
		tmp.m10 *= det;
		tmp.m12 *= det;
		tmp.m21 *= det;

		*this = tmp;
		return true;
	}

	/*
	 * Concat scale matrix to the current transformation.
	 */

	void Matrix3by3::scale(Real s)
	{
		m00 *= s;
		m01 *= s;
		m02 *= s;
		m10 *= s;
		m11 *= s;
		m12 *= s;
		m20 *= s;
		m21 *= s;
		m22 *= s;
	}


	void Matrix3by3::scale(Real sx, Real sy, Real sz)
	{
		m00 *= sx;
		m01 *= sy;
		m02 *= sz;
		m10 *= sx;
		m11 *= sy;
		m12 *= sz;
		m20 *= sx;
		m21 *= sy;
		m22 *= sz;
	}

	void Matrix3by3::scale(const Vector3D& scaleFactor)
	{
		m00 *= scaleFactor.x;
		m01 *= scaleFactor.y;
		m02 *= scaleFactor.z;
		m10 *= scaleFactor.x;
		m11 *= scaleFactor.y;
		m12 *= scaleFactor.z;
		m20 *= scaleFactor.x;
		m21 *= scaleFactor.y;
		m22 *= scaleFactor.z;
	}

	void Matrix3by3::rotate(const Vector3D& axis, Real sinAlpha, Real cosAlpha)
	{
		if (axis.normSquared() > 0.0)
		{
			Real n1, n2, n3;

			Vector3D n(axis.normalized());

			n1 = n.x;
			n2 = n.y;
			n3 = n.z;

			m00 = n1 * n1 + (1. - n1 * n1) * cosAlpha;
			m01 = n1 * n2 * (1. - cosAlpha) + n3 * sinAlpha;
			m02 = n1 * n3 * (1. - cosAlpha) - n2 * sinAlpha;

			m10 = n1 * n2 * (1. - cosAlpha) - n3 * sinAlpha;
			m11 = n2 * n2 + (1. - n2 * n2) * cosAlpha;
			m12 = n2 * n3 * (1. - cosAlpha) + n1 * sinAlpha;

			m20 = n1 * n3 * (1. - cosAlpha) + n2 * sinAlpha;
			m21 = n2 * n3 * (1. - cosAlpha) - n1 * sinAlpha;
			m22 = n3 * n3 + (1. - n3 * n3) * cosAlpha;
		}
		else
		{
			identity();
		}
	}

	void Matrix3by3::rotate(const Vector3D& axis, Real alpha)
	{
		Real sinAlpha, cosAlpha;
		sincos(alpha, sinAlpha, cosAlpha);
		rotate(axis, sinAlpha, cosAlpha);
	}

	void Matrix3by3::rotate(const Vector3D& from, const Vector3D& to)
	{
		if (from.normSquared() > 0.0 && to.normSquared() > 0.0)
		{
			Vector3D a(from.normalized());
			Vector3D b(to.normalized());

			// linear 
			if ((a - b).normSquared() < Constant::EPSILON)
			{
				// same direction
				identity();
			}
			else if ((a + b).normSquared() < Constant::EPSILON)
			{
				// opposite direction
				Real n1 = fabs(a.x);
				Real n2 = fabs(a.y);
				Real n3 = fabs(a.z);
				if (n1 >= max(n2, n3))
				{
					m00 = -1.0;
					m01 = 0.0;
					m02 = 0.0;

					m10 = 0.0;
					m11 = -1.0;
					m12 = 0.0;

					m20 = 0.0;
					m21 = 0.0;
					m22 = 1.0;
				}
				else if (n2 >= max(n1, n3))
				{
					m00 = 1.0;
					m01 = 0.0;
					m02 = 0.0;

					m10 = 0.0;
					m11 = -1.0;
					m12 = 0.0;

					m20 = 0.0;
					m21 = 0.0;
					m22 = -1.0;
				}
				else
				{
					m00 = -1.0;
					m01 = 0.0;
					m02 = 0.0;

					m10 = 0.0;
					m11 = 1.0;
					m12 = 0.0;

					m20 = 0.0;
					m21 = 0.0;
					m22 = -1.0;
				}
			}
			else
			{
				Vector3D axis(b.cross(a));
				axis.normalize();
				rotate(axis, a.cross(b).norm(), a.dot(b));
			}
		}
		else
		{
			identity();
		}
	}


	//______________________________________________________________________ friends

	std::ostream& operator<<(std::ostream& os, const Matrix3by3& tm)
	{
		os << "[[" << tm.m00 << ", " << tm.m01 << ", " << tm.m02 << "], ";
		os << "[" << tm.m10 << ", " << tm.m11 << ", " << tm.m12 << "], ";
		os << "[" << tm.m20 << ", " << tm.m21 << ", " << tm.m22 << "]]";
		return os;
	}

	Vector3D operator*(const Vector3D& point, const Matrix3by3& tm)
	{
		return Vector3D(tm.m00 * point.x + tm.m10 * point.y + tm.m20 * point.z,
		                tm.m01 * point.x + tm.m11 * point.y + tm.m21 * point.z,
		                tm.m02 * point.x + tm.m12 * point.y + tm.m22 * point.z);
	}

	void convert(const Matrix3by3& from, float to[9])
	{
		to[0] = static_cast<float>(from.m00);
		to[1] = static_cast<float>(from.m01);
		to[2] = static_cast<float>(from.m02);
		to[3] = static_cast<float>(from.m10);
		to[4] = static_cast<float>(from.m11);
		to[5] = static_cast<float>(from.m12);
		to[6] = static_cast<float>(from.m20);
		to[7] = static_cast<float>(from.m21);
		to[8] = static_cast<float>(from.m22);
	}

	void convert(const Matrix3by3& from, double to[16])
	{
		to[0] = static_cast<double>(from.m00);
		to[1] = static_cast<double>(from.m01);
		to[2] = static_cast<double>(from.m02);
		to[3] = static_cast<double>(from.m10);
		to[4] = static_cast<double>(from.m11);
		to[5] = static_cast<double>(from.m12);
		to[6] = static_cast<double>(from.m20);
		to[7] = static_cast<double>(from.m21);
		to[8] = static_cast<double>(from.m22);
	}
}
