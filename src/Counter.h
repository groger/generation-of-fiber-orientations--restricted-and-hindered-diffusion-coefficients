#ifndef COUNTER_H
#define COUNTER_H

/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <SphereIkosahedron.h>
#include <Vector.h>
#include <Matrix.h>
#include <GaussianOnSphere.h>

typedef itk::SphereIkosahedron<double> SphereIkosahedronType ;

typedef double AccumulateType ;

//*******************************************************************************

class Counter{
	/*
	 * 
	 * Base class for icosahedron-based 3d histogram
	 * 
	 */
	public:
		static const int NO_ID = -1 ;

		static void Initialize(const int subdivisionLevel) ;		//Initialize the icosahedron model shared by all the Counters
										//Note that this must be invoked before any concrete Counter is instantiated
		static int GetIcosahedronSubdivisionLevel() ;
		static int GetSize() ;
		
		explicit Counter() ;
		
		virtual void Add(geometry::Vector direction, const double weight = 1) = 0 ;
		
		virtual vtkSmartPointer<vtkPolyData> GetVTKPolyData(const geometry::Matrix &m = geometry::identity(), const double weight = 1) const ;
		
		virtual vtkSmartPointer<vtkPolyData> GetVTKPolyDataWithLocalMaximaMarked(const geometry::Matrix &m = geometry::identity(), const int neighborLevel = 2) const ;
		
		void SetId(const int id) ;
		int GetId() const ;
		
		bool IsEmpty() const ;					//Check if the bins are empty
		
		int GetNumTotalVotes() const ;				//Get the total number of votes
		
		const std::vector<AccumulateType>& GetBins() const ;
		void SetBins(const std::vector<AccumulateType>& b) ;
		
		AccumulateType GetBinValue(const int binNum) const ;
		
		VectorType GetBinCartesianCoordinates(const int binNum) ;
		VectorType GetBinSphericalCoordinates(const int binNum) ;
		
		VectorType GetCartesianCoordinates(const double phi, const double theta) const ;		//Convert spherical coordinates to cartesian coordinates, where phi and theta are in radians
		VectorType RegularizeSphericalCoordinates(const double phi, const double theta) const ;		//Convert phi and theta to right range
		VectorType GetSphericalCoordinates(const double x, const double y, const double z) const ;	//Convert cartesian coordinates to spherical coordinates
		
		std::vector<double> GetFrequency() const ;
		
		void WriteCounterToVTKFile(const char *fname = "Counter.vtk") const ;
		
		void WriteCounterWithLocalMaximaMarkedToVTKFile(const char *fname = "CounterWithLocalMaxima.vtk") const ;
		
		virtual ~Counter(){}
		
	protected:
		static SphereIkosahedronType::Pointer icosahedron ;
		
		std::vector<AccumulateType> bins ;
		
		int id ;
	
	private:
		Counter(const Counter &) ;
		void operator=(const Counter &) ;
} ;

std::ostream &operator<<(std::ostream &os, Counter &c) ;

//*******************************************************************************

class Counter_NearestNeighborVertex:public Counter{
	public:
		explicit Counter_NearestNeighborVertex() ;
		
		void Add(geometry::Vector direction, const double weight = 1) ;
		
		~Counter_NearestNeighborVertex(){}
	private:
		Counter_NearestNeighborVertex(const Counter_NearestNeighborVertex &) ;
		void operator=(const Counter_NearestNeighborVertex &) ;
} ;

//*******************************************************************************

class Counter_WeightedVertices:public Counter{
	public:
		explicit Counter_WeightedVertices() ;
		
		void Add(geometry::Vector direction, const double weight = 1) ;

		~Counter_WeightedVertices(){}
	private:
		Counter_WeightedVertices(const Counter_WeightedVertices &) ;
		void operator=(const Counter_WeightedVertices &) ;
} ;

//*******************************************************************************

//NOTICE that this counter is for statistics use only.
//This counter is constructed with the bin values (possibily normalized) of another counter. It assumes that the bin values are SYMMETRIC with respect to the center of the icosahedron
class Counter_SymmetricStatistics:public Counter{
	public:
		explicit Counter_SymmetricStatistics(const std::vector<AccumulateType> &bins) ;
		
		void Add(geometry::Vector direction, const double weight = 1){}	//This does nothing
	  
		GaussianOnSphere FitOppositeGaussianPair() ;		//Fit a pair of Gaussian models that are opposite to each other on the sphere
									//The returned values are parameters for one of them
		GaussianOnSphere_Two FitTwoOppositeGaussianPairs() ;	//Fit two pairs of opposite Gaussian models on the sphere
									//The returned values are parameters for one model in each pair
		void SubtractOppositeGaussianPair(const GaussianOnSphere &g) ;	//Subtract an opposite Gaussian pair from bins
		
		virtual vtkSmartPointer<vtkPolyData> GetVTKStatResult_LevenbergMarquardt(const GaussianOnSphere_Two &g, const geometry::Matrix &m = geometry::identity()) const ;
		
		virtual vtkSmartPointer<vtkPolyData> GetVTKStatResult_SKmeans(const geometry::Matrix &m = geometry::identity()) const ;
		
		virtual ~Counter_SymmetricStatistics(){}
	private:
		Counter_SymmetricStatistics(const Counter_SymmetricStatistics &) ;
		void operator=(const Counter_SymmetricStatistics &) ;
} ;

#endif
