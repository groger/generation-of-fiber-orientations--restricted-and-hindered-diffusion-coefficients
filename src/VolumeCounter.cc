/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <VolumeCounter.h>

PositionInVolume::PositionInVolume(const int x, const int y, const int z){
	val[0] = x ;
	val[1] = y ;
	val[2] = z ;
}

PositionInVolume& PositionInVolume::operator=(const PositionInVolume &p){
	for (int i = 0; i < VOLUME_DIMENSION; i++){
        val[i] = p[i];
    }
    return *this;
}

int &PositionInVolume::operator[](const int idx){
	if (idx < 0 || idx >= VOLUME_DIMENSION){
		throw IndexOutOfRange("PositionInVolume::operator[]") ;
	}
	
	return val[idx] ;
}

int PositionInVolume::operator[](const int idx) const{
	if (idx < 0 || idx >= VOLUME_DIMENSION){
		throw IndexOutOfRange("PositionInVolume::operator[]") ;
	}
	
	return val[idx] ;
}

bool PositionInVolume::isValid() const{
	return val[0] >= 0 && val[1] >= 0 && val[2] >= 0 ;
}

std::ostream &operator<<(std::ostream &os, const PositionInVolume &p){
	os << "PositionInVolume: (" << p[0] << ", " << p[1] << ", " << p[2] << ")" ;
	return os ;
}

bool operator==(const PositionInVolume &piv1, const PositionInVolume &piv2){
	for (int i = 0; i < VOLUME_DIMENSION; i++){
		if (piv1[i] != piv2[i])
			return false ;
	}
	return true ;
}
