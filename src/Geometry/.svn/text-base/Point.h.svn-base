#ifndef POINT_H
#define POINT_H

/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <iostream>

namespace geometry{

class Vector ;
class Matrix ;

class Point{
protected:
	double x, y, z ;
public:
	explicit Point(const double x = 0, const double y = 0, const double z = 0) ;
	double getX() const {return x ;}
	double getY() const {return y ;}
	double getZ() const {return z ;}
	double operator[](const int index) const ;
	double &operator[](const int index) ;
	double &getRef(const int index) ;
} ;

std::ostream &operator<<(std::ostream &os, const Point &p) ;

Point operator+(const Point &p, const Vector &v) ;
Point operator-(const Point &p, const Vector &v) ;

bool operator==(const Point &p1, const Point &p2) ;
double distance(const Point &p1, const Point &p2) ;

Point operator*(const Matrix &m, const Point &p) ;

}
#endif
