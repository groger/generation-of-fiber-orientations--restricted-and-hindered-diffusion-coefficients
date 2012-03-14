#ifndef EXCEPTION_H
#define EXCEPTION_H

/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <string>
#include <iostream>

namespace common{

class Exception{
public:
	virtual std::string getMessage() const = 0 ;
	virtual ~Exception(){}
} ;

std::ostream &operator<<(std::ostream &os, const Exception &e) ;

}
#endif
