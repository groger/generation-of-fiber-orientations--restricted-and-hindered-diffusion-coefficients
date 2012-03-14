/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <FiberSegment.h>

FiberSegment::FiberSegment(){
}
FiberSegment::FiberSegment(const geometry::Point& a, const geometry::Point& b):p1(a), p2(b){
}

FiberSegment& FiberSegment::operator=(const FiberSegment& fs){
    p1 = fs.p1;
    p2 = fs.p2;
    return *this;
}

std::ostream &operator<<(std::ostream &os, const FiberSegment &fs){
	return os << "Segment: [" << fs.p1 << ", " << fs.p2 << "]" ;
}
