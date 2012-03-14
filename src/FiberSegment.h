#ifndef FIBERSEGMENT_H
#define FIBERSEGMENT_H

/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <Point.h>

struct FiberSegment{
	geometry::Point p1, p2 ;
    FiberSegment();
	FiberSegment(const geometry::Point &a, const geometry::Point &b) ;
	FiberSegment& operator=(const FiberSegment &fs) ;
} ;

std::ostream &operator<<(std::ostream &os, const FiberSegment &fs) ;

#endif
