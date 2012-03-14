/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <Geometry_Common.h>

namespace geometry{

const double EPSILON = 0.000001 ;

const double PI = 3.141592653 ;

ZeroVector::ZeroVector(const std::string& msg): message(msg){
}

std::string ZeroVector::getMessage() const{
	return "Zero vector : " + message ;
}

IndexOutOfRange::IndexOutOfRange(const std::string& msg): message(msg){
}

std::string IndexOutOfRange::getMessage() const{
	return "Index out of range : " + message ;
}

}
