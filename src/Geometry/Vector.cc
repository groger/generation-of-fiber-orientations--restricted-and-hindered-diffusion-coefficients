/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <Vector.h>
#include <Point.h>
#include <Matrix.h>
#include <Geometry_Common.h>
#include <cmath>
#include <iomanip>

namespace geometry{

Vector::Vector(const double x, const double y, const double z){
	this->x = x ;
	this->y = y ;
	this->z = z ;
}

Vector::Vector(const Point &tail, const Point &tip){
	x = tip.getX() - tail.getX() ;
	y = tip.getY() - tail.getY() ;
	z = tip.getZ() - tail.getZ() ;
}

double Vector::magnitude() const{
	double squareMagnitude = x * x + y * y + z * z ;
	return sqrt(squareMagnitude) ;
}

double Vector::magnitudeSquare() const{
	return x * x + y * y + z * z ;
}

bool Vector::isZero() const{
	return (0 == magnitude()) ;
}

void Vector::normalize(){
	double mag = magnitude() ;
	if (mag > 0){
		x = x / mag ;
		y = y / mag ;
		z = z / mag ;
	}
	else{
		throw ZeroVector("Vector::normalize") ;
	}
}

Vector Vector::operator -() const{
	return Vector(-x, -y, -z) ;
}

double Vector::operator[](const int index) const{
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
		throw IndexOutOfRange("Vector::operator[]") ;
		break ;
	} ;
}

double &Vector::operator[](const int index){
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
		throw IndexOutOfRange("Vector::operator[]") ;
		break ;
	} ;
}

double &Vector::getRef(const int index){
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
		throw IndexOutOfRange("Vector::getRef") ;
		break ;
	} ;
}

std::ostream &operator<<(std::ostream &os, const Vector &v){
	os << "Vector: (" << std::fixed << v.getX() << ", " << v.getY() << ", " << v.getZ() << ")" ;
	return os ;
}

Vector operator-(const Point &tip, const Point &tail){
	return Vector(tail, tip) ;
}

Vector operator*(const Vector &v, const double factor){
	return Vector(v.getX() * factor, v.getY() * factor, v.getZ() * factor) ;
}

Vector operator*(const double factor, const Vector &v){
	return v * factor ;
}

bool operator==(const Vector &v1, const Vector &v2){
	return v1.getX() == v2.getX() && v1.getY() == v2.getY() && v1.getZ() == v2.getZ() ;
}

Vector operator+(const Vector &v1, const Vector &v2){
	return Vector(v1.getX() + v2.getX(), v1.getY() + v2.getY(), v1.getZ() + v2.getZ()) ;
}

Vector operator-(const Vector &v1, const Vector &v2){
	return Vector(v1.getX() - v2.getX(), v1.getY() - v2.getY(), v1.getZ() - v2.getZ()) ;
}

double operator*(const Vector &v1, const Vector &v2){
	return v1.getX() * v2.getX() + v1.getY() * v2.getY() + v1.getZ() * v2.getZ() ;
}

Vector crossProduct(const Vector &v1, const Vector &v2){
	return Vector(	v1.getY() * v2.getZ() - v1.getZ() * v2.getY(),
			v1.getZ() * v2.getX() - v1.getX() * v2.getZ(),
			v1.getX() * v2.getY() - v1.getY() * v2.getX()) ;
}

Vector symmetric(const Vector &a, const Vector &bisector){
	//   a    b    c    a,c are tails of r,v. Connect ac, ac intersects n at b, and |ob| = cos(theta)|v|
	//   ^\   ^n  /^
	//  r  \  |  / v
	//      \ | /
	//       \|/
	//        o
	// |n| = 1, namely a unit vector. r = v + (r - v) = v + 2 * (ob - v) = v + 2 * ((n * v)n - v)
	Vector n = bisector ;
	n.normalize() ;
	return 2 * ((n * a) * n) - a ;
}

Vector operator*(const Matrix &m, const Vector &v){	//NOTICE:Homogeneous component of a vector is 0!
	double r[3] ;
	for (int i = 0; i < 3; i++){
		r[i] = 0 ;
		for (int j = 0; j < 3; j++)
			r[i] += m[i][j] * v[j] ;
	}
	return Vector(r[0], r[1], r[2]) ;
}

}
