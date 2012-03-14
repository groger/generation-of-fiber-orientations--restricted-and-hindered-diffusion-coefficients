#ifndef INDEXVALUEPAIR_H
#define INDEXVALUEPAIR_H

/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

namespace common{

template<typename T>
class IndexValuePair{
protected:
	int index ;
	T value ;
public:
	IndexValuePair(const int &idx, const T&v) ;
	int getIndex() const ;
	T getValue() const ;
} ;

template<typename T>
IndexValuePair<T>::IndexValuePair(const int& idx, const T& v):index(idx), value(v){
}

template<typename T>
int IndexValuePair<T>::getIndex() const{
	return index ;
}

template<typename T>
T IndexValuePair<T>::getValue() const{
	return value ;
}

template<typename T>
bool operator<(const IndexValuePair<T> &ivp1, const IndexValuePair<T> &ivp2){
	return ivp1.getValue() < ivp2.getValue() ;
}

}
#endif
