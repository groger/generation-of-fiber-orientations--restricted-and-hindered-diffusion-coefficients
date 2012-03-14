#ifndef COUNTERSTATISTICS_H
#define COUNTERSTATISTICS_H

/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <Counter.h>
#include <vnl/vnl_least_squares_function.h>

//Following classes are for Spherical K-means clustering

//In this SKmeans implementation, each point on the sphere is associated with a weight
struct PointOnSphere{
	geometry::Vector point ;
	double weight ;
	PointOnSphere(const double x, const double y, const double z, const double weight) ;
} ;

struct SKCluster{
	geometry::Vector mean ;
	double sigma ;	//sigma is currently the arclenth on a unit sphere, thus its value equals its angle in radians
} ;

std::ostream &operator<<(std::ostream &os, const SKCluster &skcluster) ;
std::ostream &operator<<(std::ostream &os, const std::vector<SKCluster> &skclusters) ;

class SKMeans{	//Perform the Spherical K-means
protected:
	std::vector<PointOnSphere> points ;
public:
	void AddPoint(const PointOnSphere& p) ;
	PointOnSphere GetPoint(const int pNum) const ;
	int GetNumberOfPoints() const ;
	
	std::vector<SKCluster> skmeans(const int k) ;
	
	~SKMeans(){}
} ;

//Following classes are for Levenberg Marquardt

class OppositeGaussianPairFunctor:public vnl_least_squares_function{
protected:
	Counter_SymmetricStatistics &counter ;

public:
	explicit OppositeGaussianPairFunctor(Counter_SymmetricStatistics &c) ;
	
	void f(vnl_vector<double> const& x, vnl_vector<double>& fx) ;
} ;

class TwoOppositeGaussianPairsFunctor:public vnl_least_squares_function{
protected:
	Counter_SymmetricStatistics &counter ;

public:
	explicit TwoOppositeGaussianPairsFunctor(Counter_SymmetricStatistics &c) ;
	
	void f(vnl_vector<double> const& x, vnl_vector<double>& fx) ;
} ;

#endif
