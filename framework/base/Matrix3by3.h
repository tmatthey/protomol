/*  -*- c++ -*-  */
#ifndef MATRIX3by3_H
#define MATRIX3by3_H

#include <iostream>
#include "Vector3D.h"

namespace ProtoMol
{
	/*________________________________________________________________ Matrix3by3
	 *
	 *  [ m00 m01 m02 ]   
	 *  [ m10 m11 m12 ]   representation of a 3 by 3 matrix
	 *  [ m20 m21 m22 ]   
	 */

	class Matrix3by3
	{
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Constructors, destructors, assignment
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		Matrix3by3();
		// Matrix3by3(const Matrix3by3&);  // Use default verions.
		Matrix3by3(Real x00, Real x01, Real x02,
		           Real x10, Real x11, Real x12,
		           Real x20, Real x21, Real x22);
		Matrix3by3(float mat[9]);
		Matrix3by3(double mat[9]);
		Matrix3by3(const Vector3D& a, const Vector3D& b);
		// construct the outer-product matrix = a .outer. b;
		Matrix3by3(const Vector3D& v1, const Vector3D& v2, const Vector3D& v3);
		// construct the matrix with three Vector3D as rows.
	public:
		// Matrix3by3& operator=(const Matrix3by3&);  // Use default version.

	public:
		// ~Matrix3by3();  // Use default version.

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New methods of class Matrix3by3
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		const Matrix3by3& getIdentity(); // Return the identity matrix.
		void identity(); // Set 'this' to the identity matrix.
		void zeroMatrix();// Set 'this' to the zero matrix.
		bool zero() const; // Check if all elements equals zero.
		Real det();

		Real operator()(int i, int j) const;
		// Extract an element.
		void operator()(int i, int j, Real x);
		// Set an element.
		void operator()(Real x00, Real x01, Real x02,
		                Real x10, Real x11, Real x12,
		                Real x20, Real x21, Real x22);
		// set a whole matrix

		bool operator==(const Matrix3by3& tm) const;
		bool operator!=(const Matrix3by3& tm) const;

		Matrix3by3& operator*=(const Matrix3by3& tm);
		Matrix3by3 operator*(const Matrix3by3& tm) const;
		// Matrix multiplication.

		Vector3D operator*(const Vector3D& tm) const;
		// matrix mult a coordinate, return a coordinate

		Matrix3by3& operator*=(const Real tm);
		Matrix3by3 operator*(const Real tm) const;
		// Matrix mult constant.

		Matrix3by3& operator/=(const Real tm);
		Matrix3by3 operator/(const Real tm) const;
		// Matrix div constant.

		Matrix3by3& operator+=(const Matrix3by3& tm);
		Matrix3by3 operator+(const Matrix3by3& tm) const;
		// Matrix addition.

		Matrix3by3& operator-=(const Matrix3by3& tm);
		Matrix3by3 operator-(const Matrix3by3& tm) const;
		// Matrix subtraction.

		Matrix3by3 operator-(void) const;
		// negate operator that returns (- (*this)).

		void transpose(); // transpose 'this'. 
		void transpose(const Matrix3by3& tm); // Set 'this' to tm transposed, leave tm unchanged.
		Matrix3by3 transposed() const;
		// Return the transposed matrix, leave the original unchanged.
		// Matrix transpose

		bool invert();
		// Compute the inverse of the given matrix. 
		// Returns: false -> singular matrix
		//          true  -> success.

		void scale(Real sx, Real sy, Real sz);
		void scale(Real s);
		void scale(const Vector3D& scaleFactor);
		// Concat scale matrix to the current transformation.


		/// Rotation matrix along axis by alpha (rad)
		void rotate(const Vector3D& axis, Real alpha);
		void rotate(const Vector3D& axis, Real sinAlpha, Real cosAlpha);
		/// Rotation matrix such that from -> to
		/// and rotate along z-axes by beta (rad)
		void rotate(const Vector3D& from, const Vector3D& to);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Friends of class Matrix3by3
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
		//friend MyStreamer& operator<< (MyStreamer& OS, const Matrix3by3& tm);
		friend std::ostream& operator<<(std::ostream& os, const Matrix3by3& tm);

		friend Vector3D operator*(const Vector3D& point, const Matrix3by3& tm);

		// computes and retuns the outer product of two coordinates.

		friend void convert(const Matrix3by3& from, double to[9]);
		friend void convert(const Matrix3by3& from, float to[9]);
		// Conversion to built-in types.

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// My data members
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	private:
		Real m00, m01, m02;
		Real m10, m11, m12;
		Real m20, m21, m22;

	private:
		static const Matrix3by3 ourIdentity;
	};

	//____________________________________________________________________ friend
	//MyStreamer& operator<< (MyStreamer& OS, const Matrix3by3& tm);
	std::ostream& operator<<(std::ostream& os, const Matrix3by3& tm);
	Vector3D operator*(const Vector3D& point, const Matrix3by3& tm);
	void convert(const Matrix3by3& from, double to[9]);
	void convert(const Matrix3by3& from, float to[9]);

	//______________________________________________________________________ INLINES


	inline const Matrix3by3& Matrix3by3::getIdentity()
	{
		return ourIdentity;
	}
}
#endif // MATRIX3D_H
