#ifndef GAUSSIANONSPHERE_H
#define GAUSSIANONSPHERE_H

/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <iostream>

double gaussian(const double coefficient, const double disFromCenter, const double sigma) ;

struct GaussianOnSphere{
	double coefficient, phi, theta, sigma ;	//phi and theta indicate the location of the center of Gaussian
		//sigma is currently the arclenth on a unit sphere, thus its value equals its angle in radians
	GaussianOnSphere() ;
	bool IsValid() const ;
} ;

std::ostream &operator<<(std::ostream &os, const GaussianOnSphere &g) ;

struct GaussianOnSphere_Two{
	GaussianOnSphere g1, g2 ;
	bool IsValid() const ;
} ;

#endif
