#ifndef SPHERE_H
#define SPHERE_H

/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <Point.h>
#include <Vector.h>

namespace geometry{

class Sphere{
protected:
	Point center ;
	double radius ;
public:
	Sphere(const Point &center, const double radius) ;
	Vector normal(const Point &pointOnSurface) const ;
	bool intersect(const Point &raySource, const Point &pointOnRay, Point &intersectPoint) const ;
	bool intersect(const Point &raySource, const Vector &direction, Point &intersectPoint) const ;
} ;

}
#endif
