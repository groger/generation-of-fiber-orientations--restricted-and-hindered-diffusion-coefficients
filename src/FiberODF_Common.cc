/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <FiberODF_Common.h>

IndexOutOfRange::IndexOutOfRange(const std::string& msg): message(msg){
}

std::string IndexOutOfRange::getMessage() const{
	return "Index out of range : " + message ;
}

DimensionsDoNotMatch::DimensionsDoNotMatch(const std::string& msg): message(msg){
}

std::string DimensionsDoNotMatch::getMessage() const{
	return "Dimensions do not match : " + message ;
}

EndOfTraversal::EndOfTraversal(const std::string& msg): message(msg){
}

std::string EndOfTraversal::getMessage() const{
	return "End of traversal : " + message ;
}

InvalidDataRange::InvalidDataRange(const std::string& msg): message(msg){
}

std::string InvalidDataRange::getMessage() const{
	return "Data range is invalid : " + message ;
}
