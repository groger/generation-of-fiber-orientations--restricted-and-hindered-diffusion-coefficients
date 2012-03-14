#ifndef MATRIX_H
#define MATRIX_H

/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <iostream>

namespace geometry{

class Point ;
class Vector ;

class Matrix{
protected:
	double data[4][4] ;
public:
	struct Det3{
		double data[3][3] ;
		double operator()() const ;
	} ;
	
	Matrix() ;
	explicit Matrix(
		const double d00, const double d01, const double d02, const double d03,
		const double d10, const double d11, const double d12, const double d13,
		const double d20, const double d21, const double d22, const double d23,
		const double d30, const double d31, const double d32, const double d33
	) ;
	
	void set(
		const double d00, const double d01, const double d02, const double d03,
		const double d10, const double d11, const double d12, const double d13,
		const double d20, const double d21, const double d22, const double d23,
		const double d30, const double d31, const double d32, const double d33
	) ;
	
	double *operator[](const int rowNum) ;
	const double *operator[](const int rowNum) const ;
	
	bool isZero() const ;
	bool isIdentity() const ;
	
	double determinant() const ;		//Compute the determinant of the matrix
	Det3 getDet3(const int r, const int c) const ;
} ;

std::ostream &operator<<(std::ostream &os, const Matrix &m) ;

bool operator==(const Matrix &m1, const Matrix &m2) ;
Matrix operator*(const Matrix &m1, const Matrix &m2) ;

Matrix inverse(const Matrix& m) ;
Matrix identity() ;

Matrix lookAt(const Point &cameraPos ,const Point &centerPos, const Vector &up) ;
Matrix orthoProjection(const double left, const double right, const double bottom, const double top, const double near, const double far) ;

Matrix rotateAroundX(const double angle) ;
Matrix rotateAroundY(const double angle) ;
Matrix rotateAroundZ(const double angle) ;

Matrix translate(const double tx, const double ty, const double tz) ;
Matrix scale(const double fx, const double fy, const double fz) ;

}
#endif
