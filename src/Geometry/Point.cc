/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <Point.h>
#include <Vector.h>
#include <Matrix.h>
#include <Geometry_Common.h>
#include <iomanip>

namespace geometry{

Point::Point(const double x, const double y, const double z){
	this->x = x ;
	this->y = y ;
	this->z = z ;
}

double Point::operator[](const int index) const{
	switch (index){
	case AXIS_X:
		return x ;
		break ;
	case AXIS_Y:
		return y ;
		break ;
	case AXIS_Z:
		return z ;
		break ;
	default:
		throw IndexOutOfRange("Point::operator[]") ;
		break ;
	} ;
}

double &Point::operator[](const int index){
	switch (index){
	case AXIS_X:
		return x ;
		break ;
	case AXIS_Y:
		return y ;
		break ;
	case AXIS_Z:
		return z ;
		break ;
	default:
		throw IndexOutOfRange("Point::operator[]") ;
		break ;
	} ;
}

double &Point::getRef(const int index){
	switch (index){
	case AXIS_X:
		return x ;
		break ;
	case AXIS_Y:
		return y ;
		break ;
	case AXIS_Z:
		return z ;
		break ;
	default:
		throw IndexOutOfRange("Point::getRef") ;
		break ;
	} ; 
}

std::ostream &operator<<(std::ostream &os, const Point &p){
	os << "Point: (" << std::fixed << p.getX() << ", " << p.getY() << ", " << p.getZ() << ")" ;
	return os ;
}

Point operator+(const Point &p, const Vector &v){
	return Point(p.getX() + v.getX(), p.getY() + v.getY(), p.getZ() + v.getZ()) ;
}

Point operator-(const Point &p, const Vector &v){
	return Point(p.getX() - v.getX(), p.getY() - v.getY(), p.getZ() - v.getZ()) ;
}

bool operator==(const Point &p1, const Point &p2){
	return p1.getX() == p2.getX() && p1.getY() == p2.getY() && p1.getZ() == p2.getZ() ;
}

double distance(const Point &p1, const Point &p2){
	return Vector(p1, p2).magnitude() ;
}

Point operator*(const Matrix &m, const Point &p){	//NOTICE:Homogeneous component of a point is 1
	double r[4] ;
	for (int i = 0; i < 4; i++){
		r[i] = 0 ;
		for (int j = 0; j < 3; j++)
			r[i] += m[i][j] * p[j] ;
		r[i] += m[i][3] * 1 ;
	}
	if (r[3] != 0){
		for (int i = 0; i < 3; i++)
			r[i] /= r[3] ;
	}
	return Point(r[0], r[1], r[2]) ;
}

}
