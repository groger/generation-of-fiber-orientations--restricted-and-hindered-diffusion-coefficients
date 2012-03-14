#ifndef GEOMETRY_COMMON_H
#define GEOMETRY_COMMON_H

/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <Exception.h>

namespace geometry{

enum {AXIS_X = 0, AXIS_Y = 1, AXIS_Z = 2} ;

#define MIN(a, b) (((a) < (b))?(a) : (b))
#define MAX(a, b) (((a) > (b))?(a) : (b))

extern const double EPSILON ;

extern const double PI ;

class Geometry_Exception:public common::Exception{
} ;

class ZeroVector:public Geometry_Exception{
protected:
	std::string message ;
public:
	ZeroVector(const std::string &msg) ;
	std::string getMessage() const ;
} ;

class IndexOutOfRange:public Geometry_Exception{
protected:
	std::string message ;
public:
	IndexOutOfRange(const std::string &msg) ;
	std::string getMessage() const ;
} ;

}
#endif
