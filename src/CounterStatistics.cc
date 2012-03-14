/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <CounterStatistics.h>
#include <FiberODF_Common.h>
#include <cmath>
#include <set>

PointOnSphere::PointOnSphere(const double x, const double y, const double z, const double weight): point(x, y, z){
	if (weight <= 0){
		throw InvalidDataRange("PointOnSphere::PointOnSphere: the weight must be positive!") ;
	}
	this->weight = weight ;
	point.normalize() ;
}

std::ostream &operator<<(std::ostream &os, const SKCluster &skcluster){
	os << "SKCluster: " << std::endl ;
	os << '\t' << "mean: " << skcluster.mean << std::endl ;
	os << '\t' << "sigma: " << skcluster.sigma ;
	return os ;
}

std::ostream &operator<<(std::ostream &os, const std::vector<SKCluster> &skclusters){
	for (unsigned int i = 0; i < skclusters.size(); i++){
		os << "Cluster " << i << " : " << skclusters[i] << std::endl ;
	}
	return os ;
}

void SKMeans::AddPoint(const PointOnSphere& p){
	points.push_back(p) ;
}

PointOnSphere SKMeans::GetPoint(const int pNum) const{
	if (pNum < 0 || pNum >= static_cast<int>(points.size()))
		throw IndexOutOfRange("SKMeans::GetPoint") ;
	return points[pNum] ;
}

int SKMeans::GetNumberOfPoints() const{
	return static_cast<int>(points.size()) ;
}

std::vector<SKCluster> SKMeans::skmeans(const int k){
  
	if (k <= 0){
		throw InvalidDataRange("SKMeans::skmeans: the number of clusters must be positive!") ;
	}
	
	if (k > GetNumberOfPoints()){
		throw InvalidDataRange("SKMeans::skmeans: the number of clusters must be less than or equal to the number of points!") ;
	}
  
	std::vector<SKCluster> clusters(k) ;
	
	std::vector<std::set<int> > clusterToPoint(k) ;		//Mapping from cluster to point
	std::vector<int> pointToCluster(GetNumberOfPoints()) ;	//Mapping from point to cluster
	
	//Initialize the partitions
	//NOTE: Might need more improvement
	for (int i = 0; i < GetNumberOfPoints(); i++){
		const int clusterID = i % k ;
		clusterToPoint[clusterID].insert(i) ;
		pointToCluster[i] = clusterID ;
	}
	
	bool somePointMoved = true ;
	
	while (somePointMoved){
		
		somePointMoved = false ;
		
		//Compute the new centroids
		//Note that the underlying probabilistic model is NOT gaussian
		for (int i = 0; i < k; i++){
		  
			geometry::Vector mean ;
		  
			for (std::set<int>::iterator j = clusterToPoint[i].begin(); j != clusterToPoint[i].end(); j++){
				mean = mean + points[*j].point * points[*j].weight ;
			}
			
			mean.normalize() ;
			clusters[i].mean = mean ;
		}

		//Compute the new assignment for each point on sphere
		for (int i = 0; i < GetNumberOfPoints(); i++){
		  
			//Unit vectors, so the larger the cos value, the closer the point to mean
			double maxCos = points[i].point * clusters[pointToCluster[i]].mean ;
			
			int moveToCluster = -1 ;
			
			for (int j = 0; j < k; j++){
			  
				if (points[i].point * clusters[j].mean > maxCos){
				  
					maxCos = points[i].point * clusters[j].mean ;
					moveToCluster = j ;
					somePointMoved = true ;
					clusterToPoint[pointToCluster[i]].erase(i) ;
				}
			}
			
			if (moveToCluster > -1){	//Move point to new cluster
		
				clusterToPoint[moveToCluster].insert(i) ;
				pointToCluster[i] = moveToCluster ;
			}
			
		}

	}	//No assignment changes, clustering finished
	
	for (int i = 0; i < k; i++){
	  
		double weightSum = 0 ;
		double stdDeviation = 0 ;
		
		for (std::set<int>::iterator j = clusterToPoint[i].begin(); j != clusterToPoint[i].end(); j++){
			double disOnSphere = std::acos(clusters[i].mean * points[*j].point) ;
			stdDeviation += disOnSphere * disOnSphere * points[*j].weight ;
			weightSum += points[*j].weight ;
		}
		
		stdDeviation = std::sqrt(stdDeviation / weightSum) ;
		
		clusters[i].sigma = stdDeviation ;
	}
	
	return clusters ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

OppositeGaussianPairFunctor::OppositeGaussianPairFunctor(Counter_SymmetricStatistics &c):
		vnl_least_squares_function(4,	//Four parameters: phi, theta, coefficient, sigma
					   4,	//Required by the vnl levenberg marquardt to be at least of the same dimension with the parameters
					   no_gradient),
		counter(c){
}

void OppositeGaussianPairFunctor::f(vnl_vector<double> const& x, vnl_vector<double>& fx){
	const VectorType center = counter.GetCartesianCoordinates(x[0], x[1]) ;
	const geometry::Vector vecToCenter(center[0], center[1], center[2]) ;
	
	double res = 0 ;
	
	for (int i = 0; i < counter.GetSize(); i++){
		const VectorType vertex = counter.GetBinCartesianCoordinates(i) ;
		const geometry::Vector vecToVertex(vertex[0], vertex[1], vertex[2]) ;	//This is guaranteed to be a unit vector
		
		//NOTICE that since the counter is symmetric and the two gaussians are opposite to each other, we assume that the distribution of the gaussian is only on the half sphere where the center of gaussian is on
		//NOTE: "<" is used here, namely the vertices on the equator which is perpendicular to vecToCenter are also taken into account
		if (vecToCenter * vecToVertex < 0)	
			continue ;
		
		const double arcLength = std::acos(vecToCenter * vecToVertex) ;
		res += std::pow(counter.GetBinValue(i) - gaussian(x[2], arcLength, x[3]), 2) ;
	}
	
	fx[0] = res ;
	fx[1] = fx[2] = fx[3] = 0 ;
}

TwoOppositeGaussianPairsFunctor::TwoOppositeGaussianPairsFunctor(Counter_SymmetricStatistics &c):
		vnl_least_squares_function(8,	//Parameters for two gaussians: phi1, theta1, coefficient1, sigma1, ...
					   8,
					   no_gradient),
		counter(c){
}

void TwoOppositeGaussianPairsFunctor::f(vnl_vector<double> const& x, vnl_vector<double>& fx){
	const VectorType center1 = counter.GetCartesianCoordinates(x[0], x[1]) ;
	const VectorType center2 = counter.GetCartesianCoordinates(x[4], x[5]) ;
	const geometry::Vector vecToCenters[] = {
		geometry::Vector(center1[0], center1[1], center1[2]),
		geometry::Vector(-center1[0], -center1[1], -center1[2]),
		geometry::Vector(center2[0], center2[1], center2[2]),
		geometry::Vector(-center2[0], -center2[1], -center2[2])
	} ;
  
	double res = 0 ;
	
	for (int i = 0; i < counter.GetSize(); i++){
		const VectorType vertex = counter.GetBinCartesianCoordinates(i) ;
		const geometry::Vector vecToVertex(vertex[0], vertex[1], vertex[2]) ;
		
		double modelSum = 0 ;	//Sum of models at this vertex
		
		for (int j = 0; j < 4; j++){
			if (vecToCenters[j] * vecToVertex >= 0){	//Each gaussian model spans a half sphere
				const double arcLength = std::acos(vecToCenters[j] * vecToVertex) ;
				modelSum += gaussian(x[j / 2 * 4 + 2], arcLength, x[j / 2 * 4 + 3]) ;
			}
		}
		
		res += std::pow(counter.GetBinValue(i) - modelSum, 2) ;
	}
	
	fx[0] = res ;
	for (int i = 1; i < 8; i++)
		fx[i] = 0 ;
}
