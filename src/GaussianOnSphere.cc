/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <GaussianOnSphere.h>
#include <cmath>

double gaussian(const double coefficient, const double disFromCenter, const double sigma){
	return coefficient * std::exp(- std::pow(disFromCenter, 2) / 2 / std::pow(sigma, 2)) ;
}

GaussianOnSphere::GaussianOnSphere(){
	coefficient = -1 ;
	phi = 0 ;
	theta = 0 ;
	sigma = -1 ;
}

bool GaussianOnSphere::IsValid() const{
	return (coefficient > 0 && sigma > 0) ;
}

std::ostream &operator<<(std::ostream &os, const GaussianOnSphere &g){
	os << "GaussianOnSphere: " << std::endl ;
	os << "coefficient: " << g.coefficient << std::endl ;
	os << "phi: " << g.phi << std::endl ;
	os << "theta: " << g.theta << std::endl ;
	os << "sigma: " << g.sigma ;
	return os ;
}

bool GaussianOnSphere_Two::IsValid() const{
	return (g1.IsValid() || g2.IsValid()) ;		//The two-gaussian group is deemed valid if one of them is valid
}
