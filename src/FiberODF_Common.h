#ifndef FIBERODF_COMMON_H
#define FIBERODF_COMMON_H

/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <Exception.h>

class IndexOutOfRange:public common::Exception{
protected:
	std::string message ;
public:
	IndexOutOfRange(const std::string &msg) ;
	std::string getMessage() const ;
} ;

class DimensionsDoNotMatch:public common::Exception{
protected:
	std::string message ;
public:
	DimensionsDoNotMatch(const std::string &msg) ;
	std::string getMessage() const ;
} ;

class EndOfTraversal:public common::Exception{
protected:
	std::string message ;
public:
	EndOfTraversal(const std::string &msg) ;
	std::string getMessage() const ;
} ;

class InvalidDataRange:public common::Exception{
protected:
	std::string message ;
public:
	InvalidDataRange(const std::string &msg) ;
	std::string getMessage() const ;
} ;

#endif
