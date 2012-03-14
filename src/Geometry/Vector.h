#ifndef VECTOR_H
#define VECTOR_H

/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <iostream>

namespace geometry{

class Point ;
class Matrix ;

class Vector{
protected:
	double x, y, z ;
public:
	explicit Vector(const double x = 0, const double y = 0, const double z = 0) ;
	Vector(const Point &tail, const Point &tip) ;
	double magnitude() const ;
	double magnitudeSquare() const ;
	void normalize() ;
	bool isZero() const ;
	Vector operator-() const ;
	double getX() const {return x ;}
	double getY() const {return y ;}
	double getZ() const {return z ;}
	double operator[](const int index) const ;
	double &operator[](const int index) ;
	double &getRef(const int index) ;
} ;

std::ostream &operator<<(std::ostream &os, const Vector &v) ;

Vector operator-(const Point &tip, const Point &tail) ;

Vector operator*(const Vector &v, const double factor) ;
Vector operator*(const double factor, const Vector &v) ;

bool operator==(const Vector &v1, const Vector &v2) ;
Vector operator+(const Vector &v1, const Vector &v2) ;
Vector operator-(const Vector &v1, const Vector &v2) ;
double operator*(const Vector &v1, const Vector &v2) ;
Vector crossProduct(const Vector &v1, const Vector &v2) ;
Vector symmetric(const Vector &a, const Vector &bisector) ;	//To get the symmetric vector of a relative to the bisector

Vector operator*(const Matrix &m, const Vector &v) ;

}
#endif
