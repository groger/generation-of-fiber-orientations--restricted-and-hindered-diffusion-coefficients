/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <CounterStatistics.h>
#include <Triangle.h>
#include <FiberODF_Common.h>
#include <IndexValuePair.h>
#include <iomanip>
#include <cmath>
#include <set>
#include <algorithm>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkAppendPolyData.h>
#include <vnl/algo/vnl_levenberg_marquardt.h>
#include <Cylinder.h>

//*******************************************************************************

SphereIkosahedronType::Pointer Counter::icosahedron = SphereIkosahedronType::New() ;

void Counter::Initialize(const int subdivisionLevel){

	if (subdivisionLevel < 0){
		std::cerr << "The subdivision level must be at least 0!" << std::endl ;
		exit(1) ;
	}
  
	Counter::icosahedron->SetSubdivisionLevel(subdivisionLevel);
	Counter::icosahedron->Initialize();
}

int Counter::GetIcosahedronSubdivisionLevel(){
	return icosahedron->GetSubdivisionLevel() ;
}

int Counter::GetSize(){
	return icosahedron->GetNumberOfVertices() ;
}

Counter::Counter(){
	id = NO_ID ;		//Initialize the id
	
	if (GetSize() == 0){
		std::cerr << "Counter::Initialize must be invoked before any Counter is instantiated!" << std::endl ;
		exit(1) ;
	}
	
	bins.resize(GetSize()) ;
	for (int i = 0; i < GetSize(); i++)
		bins[i] = 0 ;
}

void Counter::SetId(const int id){
	this->id = id ;
}

int Counter::GetId() const{
	return id ;
}

bool Counter::IsEmpty() const{
	for (int i = 0; i < GetSize(); i++)
		if (bins[i] > 0)
			return false ;
	return true ;
}

int Counter::GetNumTotalVotes() const{
	double temp = 0 ;
	for (int i = 0; i < GetSize(); i++){
		temp += bins[i] ;
	}
	return static_cast<int>(std::floor(temp)) ;
}

const std::vector<AccumulateType>& Counter::GetBins() const{
	return bins ;
}

void Counter::SetBins(const std::vector<AccumulateType>& b){
	if (static_cast<int>(b.size()) != GetSize())
		throw DimensionsDoNotMatch("Counter::SetBins") ;
	bins = b ;
}

AccumulateType Counter::GetBinValue(const int binNum) const{
	if (binNum < 0 || binNum >= GetSize())
		throw IndexOutOfRange("Counter::GetBinValue") ;
	return bins[binNum] ;
}

VectorType Counter::GetBinCartesianCoordinates(const int binNum){
	return icosahedron->GetCoordinateTableatIndex(binNum) ;
}

VectorType Counter::GetBinSphericalCoordinates(const int binNum){
	return icosahedron->GetPhiThetaTableatIndex(binNum) ;
}

VectorType Counter::GetCartesianCoordinates(const double phi, const double theta) const{
	VectorType cartesian(3) ;
	cartesian[2] = std::cos(theta) ;
	
	double temp = theta ;
	if (std::sin(theta) < 0)
		temp = -theta ;	//Guarantee that the function works when theta is out of [0, PI]
		
	cartesian[0] = std::cos(phi) * std::sin(temp) ;
	cartesian[1] = std::sin(phi) * std::sin(temp) ;
	
	return cartesian ;
}

VectorType Counter::RegularizeSphericalCoordinates(const double phi, const double theta) const{
	VectorType spherical(2) ;

	double tempPhi = phi ;
	//Regularize the phi value
	while(tempPhi >= M_PI)
		tempPhi -= 2 * M_PI ;
	while(tempPhi < -M_PI)
		tempPhi += 2 * M_PI ;
	spherical[0] = tempPhi ;
	
	double tempTheta = theta ;
	if (std::sin(theta) < 0)
		tempTheta = -theta ;
	while(tempTheta > M_PI)
		tempTheta -= 2 * M_PI ;
	while(tempTheta < 0)
		tempTheta += 2 * M_PI ;
	spherical[1] = tempTheta ;
	
	return spherical ;
}

VectorType Counter::GetSphericalCoordinates(const double x, const double y, const double z) const{
	VectorType spherical(2) ;
	
	geometry::Vector v(x, y, z) ;
	v.normalize() ;		//Exception will be thrown when the vector is zero
	
	double tempPhi = std::atan2(v[1], v[0]) ;
	//Regularize the phi value
	while(tempPhi >= M_PI)
		tempPhi -= 2 * M_PI ;
	while(tempPhi < -M_PI)
		tempPhi += 2 * M_PI ;
	spherical[0] = tempPhi ;
	
	spherical[1] = std::acos(v[2]) ;	//Theta
	
	return spherical ;
}

std::vector<double> Counter::GetFrequency() const{
	std::vector<double> frequency(bins.size()) ;
	double total = 0 ;
	for (unsigned int i = 0; i < frequency.size(); i++){
		total += bins[i] ;
		frequency[i] = 0 ;
	}
	if (total > 0){
		for (unsigned int i = 0; i < frequency.size(); i++)
			frequency[i] = bins[i] / total ;
	}
	return frequency ;
}

vtkSmartPointer<vtkPolyData> Counter::GetVTKPolyData(const geometry::Matrix &m, const double weight) const{
	vtkSmartPointer<vtkPolyData> poly = icosahedron->CreateVTKPolyData() ;
	
	vtkSmartPointer<vtkPoints> pts = poly->GetPoints() ;
	
	double temp[3] ;
	for (vtkIdType i = 0; i < pts->GetNumberOfPoints(); i++){
		pts->GetPoint(i, temp) ;
		geometry::Point p1(temp[0], temp[1], temp[2]) ;
		geometry::Point p2 = m * p1 ;
		temp[0] = p2[0] ;
		temp[1] = p2[1] ;
		temp[2] = p2[2] ;
		pts->SetPoint(i, temp) ;
	}
	
	vtkSmartPointer<vtkPointData> pointData = poly->GetPointData() ;
	
	vtkSmartPointer<vtkDoubleArray> scalars = vtkSmartPointer<vtkDoubleArray>::New() ;
	scalars->SetName("Frequencies") ;
	scalars->SetNumberOfComponents(1) ;
	std::vector<double> frequencies = GetFrequency() ;
	scalars->SetNumberOfValues(static_cast<vtkIdType>(frequencies.size())) ;
	for (unsigned int i = 0; i < frequencies.size(); i++){
		scalars->SetValue(static_cast<vtkIdType>(i), frequencies[i] * weight) ;
	}
	
	pointData->SetScalars(scalars) ;
	
	return poly ;
}

vtkSmartPointer<vtkPolyData> Counter::GetVTKPolyDataWithLocalMaximaMarked(const geometry::Matrix &m, const int neighborLevel) const{
	vtkSmartPointer<vtkPolyData> poly = icosahedron->CreateVTKPolyData() ;
	
	vtkSmartPointer<vtkPoints> pts = poly->GetPoints() ;
	
	double temp[3] ;
	for (vtkIdType i = 0; i < pts->GetNumberOfPoints(); i++){
		pts->GetPoint(i, temp) ;
		geometry::Point p1(temp[0], temp[1], temp[2]) ;
		geometry::Point p2 = m * p1 ;
		temp[0] = p2[0] ;
		temp[1] = p2[1] ;
		temp[2] = p2[2] ;
		pts->SetPoint(i, temp) ;
	}
	
	vtkSmartPointer<vtkPointData> pointData = poly->GetPointData() ;
	
	std::vector<double> mark(GetSize(), 0) ;	//The mark data for local maxima
	
	for (int i = 0; i < GetSize(); i++){
		std::set<int> coveredVertices ;
		std::set<int> newlyInsertedVertices ;
		
		coveredVertices.insert(i) ;
		newlyInsertedVertices.insert(i) ;
		
		for (int lev = 0; lev < neighborLevel; lev++){
			std::set<int> temp ;
		  
			for (std::set<int>::iterator iter = newlyInsertedVertices.begin(); iter != newlyInsertedVertices.end(); iter++){
				IndexList triangleList = icosahedron->GetSurroundingTriangles(*iter) ;
				for (unsigned int j = 0; j < triangleList.size(); j++){
					IndexList vertices = icosahedron->GetTriangleVertices(triangleList[j]) ;
					for (unsigned int k = 0; k < vertices.size(); k++){
						if (coveredVertices.find(vertices[k]) == coveredVertices.end()){	//The vertex is new
							temp.insert(vertices[k]) ;
							coveredVertices.insert(vertices[k]) ;
						}
					}
				}
			}
			
			newlyInsertedVertices = temp ;
		}
		
		double maximum = -1 ;
		for (std::set<int>::iterator iter = coveredVertices.begin(); iter != coveredVertices.end(); iter++){
			if (bins[*iter] > maximum)
				maximum = bins[*iter] ;
		}
		if (maximum == bins[i])		//This is local maximum
			mark[i] = bins[i] ;	//At local maximum, the mark value equals the bin value. Otherwise the mark stays 0
	}
	
	//Find the top 4 local maxima, namely two pairs of local maxima
	const int NUM_LOCALMAXIMA = 4 ;
	std::vector<common::IndexValuePair<double> > indexMarkPairs ;
	for (int i = 0; i < GetSize(); i++){
		indexMarkPairs.push_back(common::IndexValuePair<double>(i, mark[i])) ;
	}
	std::sort(indexMarkPairs.begin(), indexMarkPairs.end()) ;
	for (int i = 0; i < GetSize() - NUM_LOCALMAXIMA; i++){
		mark[indexMarkPairs[i].getIndex()] = 0 ;	//Mask the other marks to 0
	}
	//Normalize the local maxima marks
	if (mark[indexMarkPairs[GetSize() - 1].getIndex()] > 0){
		for (int i = GetSize() - NUM_LOCALMAXIMA; i < GetSize(); i++){
			mark[indexMarkPairs[i].getIndex()] = mark[indexMarkPairs[i].getIndex()] / mark[indexMarkPairs[GetSize() - 1].getIndex()] ;
		}
	}

	vtkSmartPointer<vtkDoubleArray> scalars = vtkSmartPointer<vtkDoubleArray>::New() ;
	scalars->SetName("LocalMaxima") ;
	scalars->SetNumberOfComponents(1) ;
	scalars->SetNumberOfValues(static_cast<vtkIdType>(mark.size())) ;
	for (unsigned int i = 0; i < mark.size(); i++){
		scalars->SetValue(static_cast<vtkIdType>(i), mark[i]) ;
	}
	
	pointData->SetScalars(scalars) ;
	
	return poly ;

}

void Counter::WriteCounterToVTKFile(const char *fname) const{
	vtkSmartPointer<vtkPolyData> data = GetVTKPolyData() ;

	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New() ;
	writer->SetInput(data) ;
	writer->SetFileName(fname) ;
	std::stringstream ss ;
	ss << "Counter : id " << id ;
	writer->SetHeader(ss.str().c_str()) ;
	writer->Write() ;
}

void Counter::WriteCounterWithLocalMaximaMarkedToVTKFile(const char *fname) const{
	vtkSmartPointer<vtkPolyData> data = GetVTKPolyDataWithLocalMaximaMarked() ;

	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New() ;
	writer->SetInput(data) ;
	writer->SetFileName(fname) ;
	std::stringstream ss ;
	ss << "Counter with local maxima : id " << id ;
	writer->SetHeader(ss.str().c_str()) ;
	writer->Write() ;
}

std::ostream &operator<<(std::ostream &os, Counter &c){
	os << "Counter: " << c.GetId() << std::endl ;
	os << "Bin values:" << std::endl ;
	const std::vector<AccumulateType> &bins = c.GetBins() ;
	os << std::fixed ;
	for (unsigned int i = 0; i < bins.size(); i++){
		os << bins[i] ;
		if (i % 4 == 3)
			os << std::endl ;
		else
			os << '\t' ;
	}
	if (bins.size() % 4 != 0)
		os << std::endl ;
	os << "Bin Cartesian coordinates and spherical coordinates:" << std::endl ;
	for (unsigned int i = 0; i < bins.size(); i++){
		VectorType cartesian = c.GetBinCartesianCoordinates(static_cast<int>(i)) ;
		//VectorType spherical = c.GetBinSphericalCoordinates(static_cast<int>(i)) ;
		os << cartesian[0] << "  " << cartesian[1] << "  " << cartesian[2] << /*"  " << spherical[0] << "  " << spherical[1] <<*/ std::endl ;
	}
	return os ;
}

//**************************************************************************

Counter_NearestNeighborVertex::Counter_NearestNeighborVertex(){
}

void Counter_NearestNeighborVertex::Add(geometry::Vector direction, const double weight){
	if (direction.isZero()){
		return ;
	}
	else{
		direction.normalize() ;
	}
	
	double temp = -1 ;
	int nearest = -1 ;
	for (int i = 0; i < GetSize(); i++){
		VectorType vertex = icosahedron->GetCoordinateTableatIndex(i) ;
		assert(vertex.size() == 3) ;
		geometry::Vector v_vertex(vertex[0], vertex[1], vertex[2]) ;
		if (v_vertex * direction > temp){
			nearest = i ;
			temp = v_vertex * direction ;
		}
	}
	
	bins[nearest] += 1 * weight ;
}

//**************************************************************************

Counter_WeightedVertices::Counter_WeightedVertices(){
}

void Counter_WeightedVertices::Add(geometry::Vector direction, const double weight){
	if (direction.isZero())
		return ;
	else{
		direction.normalize() ;
	}
	
	double temp = -1 ;
	int nearest = -1 ;
	for (int i = 0; i < GetSize(); i++){
		VectorType vertex = icosahedron->GetCoordinateTableatIndex(i) ;
		assert(vertex.size() == 3) ;
		geometry::Vector v_vertex(vertex[0], vertex[1], vertex[2]) ;
		if (v_vertex * direction > temp){
			nearest = i ;
			temp = v_vertex * direction ;
		}
	}
	
	IndexList surroundingTriangles = icosahedron->GetSurroundingTriangles(nearest) ;
	
	assert(surroundingTriangles.size() == 5 || surroundingTriangles.size() == 6) ;
	
	for (unsigned int i = 0; i < surroundingTriangles.size(); i++){
		const std::vector<VectorType> &triangle = icosahedron->GetTriangle(surroundingTriangles[i]) ;
		assert(triangle.size() == 3) ;
		geometry::Triangle t_triangle(
			geometry::Point(triangle[0][0], triangle[0][1], triangle[0][2]),
			geometry::Point(triangle[1][0], triangle[1][1], triangle[1][2]),
			geometry::Point(triangle[2][0], triangle[2][1], triangle[2][2])
		) ;
		
		geometry::Point intersectionPoint ;
		
		if (t_triangle.intersect(geometry::Point(0, 0, 0), direction, intersectionPoint)){
			IndexList triangleVertices = icosahedron->GetTriangleVertices(surroundingTriangles[i]) ;
			assert(triangleVertices.size() == 3) ;
			geometry::TriangleBarycentricCoords baryCentric = t_triangle.barycentric(intersectionPoint) ;
			bins[triangleVertices[0]] += baryCentric.w1 * weight ;
			bins[triangleVertices[1]] += baryCentric.w2 * weight ;
			bins[triangleVertices[2]] += baryCentric.w3 * weight ;
			break ;		//Stop the iteration
		}
	}
}

//**************************************************************************

Counter_SymmetricStatistics::Counter_SymmetricStatistics(const std::vector<AccumulateType> &bins){
	SetBins(bins) ;
}

GaussianOnSphere Counter_SymmetricStatistics::FitOppositeGaussianPair(){
	GaussianOnSphere res ;
	
	if (IsEmpty())
		return res ;
	
	int maxBin = -1 ;
	double saveBinMax = -1 ;
	for (int i = 0; i < GetSize(); i++){
		if (bins[i] > saveBinMax){
			maxBin = i ;
			saveBinMax = bins[i] ;
		}
	}

//Debug//////////////
	std::cout << "Opposite gaussian pair fitting:" << std::endl ;
	VectorType cart = GetBinCartesianCoordinates(maxBin) ;
	geometry::Vector cartv(cart[0], cart[1], cart[2]) ;
	std::cout << "max bin value:" << bins[maxBin] << std::endl ;
	std::cout << "max bin vec:" << cartv << std::endl ;
////////////////////


	vnl_vector<double> param(4) ;	//Initialization for levenberg marquardt
	VectorType spherical = GetBinSphericalCoordinates(maxBin) ;
	param[0] = spherical[0] ;	//phi
	param[1] = spherical[1] ;	//theta
	param[2] = saveBinMax ;		//coefficient
		//Note that the maximum bin value is used to initialize the coefficient
	param[3] = 1 ;			//sigma

	OppositeGaussianPairFunctor func(*this) ;
	vnl_levenberg_marquardt lm(func) ;
	
	lm.minimize(param) ;

	res.coefficient = param[2] ;
	res.phi = param[0] ;
	res.theta = param[1] ;
	res.sigma = std::abs(param[3]) ;	//Guarantee that sigma is positive
	
	return res ;
}

GaussianOnSphere_Two Counter_SymmetricStatistics::FitTwoOppositeGaussianPairs(){
	GaussianOnSphere_Two res ;
	
	if (IsEmpty())
		return res ;

//Debug////////////
	std::cout << "Original counter:" << std::endl ;
	std::cout << *this << std::endl ;
///////////////////
	
	const std::vector<AccumulateType> binsBackup = GetBins() ;
	
	res.g1 = FitOppositeGaussianPair() ;
	
	if (!res.g1.IsValid())		//Failed while fitting the first gaussian pair
		return res ;
	
	SubtractOppositeGaussianPair(res.g1) ;

//Debug/////////////
	std::cout << "Gaussian 1:" << std::endl ;
	std::cout << res.g1 << std::endl ;
	std::cout << "center vec:" << std::endl ;
	{
	    VectorType cenv1 = GetCartesianCoordinates(res.g1.phi, res.g1.theta) ;
	    std::cout << cenv1[0] << "," << cenv1[1] << "," << cenv1[2] << std::endl ;
	}
	std::cout << *this << std::endl ;
////////////////////

	res.g2 = FitOppositeGaussianPair() ;
	
	if (!res.g2.IsValid())		//Failed while fitting the second gaussian pair
		return res ;

//Debug/////////////////
	std::cout << "Gaussian 2:" << std::endl ;
	std::cout << res.g2 << std::endl ;
	SubtractOppositeGaussianPair(res.g2);
	std::cout << "center vec:" << std::endl ;
	{
	    VectorType cenv1 = GetCartesianCoordinates(res.g2.phi, res.g2.theta) ;
	    std::cout << cenv1[0] << "," << cenv1[1] << "," << cenv1[2] << std::endl ;
	}
	std::cout << *this << std::endl ;
////////////////////////

	SetBins(binsBackup) ;		//Restore the bins to original values

	//Fit the two pairs together. The initialization is done with the independent fitting results of the two pairs
	vnl_vector<double> param(8) ;	//Initialization for levenberg marquardt
	param[0] = res.g1.phi ;
	param[1] = res.g1.theta ;
	param[2] = res.g1.coefficient ;
	param[3] = res.g1.sigma ;
	param[4] = res.g2.phi ;
	param[5] = res.g2.theta ;
	param[6] = res.g2.coefficient ;
	param[7] = res.g2.sigma ;
	
	TwoOppositeGaussianPairsFunctor func(*this) ;
	vnl_levenberg_marquardt lm(func) ;
	
	lm.minimize(param) ;
	
	res.g1.coefficient = param[2] ;
	res.g1.phi = param[0] ;
	res.g1.theta = param[1] ;
	res.g1.sigma = std::abs(param[3]) ;
	res.g2.coefficient = param[6] ;
	res.g2.phi = param[4] ;
	res.g2.theta = param[5] ;
	res.g2.sigma = std::abs(param[7]) ;
	
	return res ;
}

void Counter_SymmetricStatistics::SubtractOppositeGaussianPair(const GaussianOnSphere &g){
	if (!g.IsValid()){
		std::cerr << "Counter_SymmetricStatistics::SubtractOppositeGaussianPair: the gaussian model is invalid!" << std::endl ;
		exit(1) ;
	}
	
	const VectorType center1 = GetCartesianCoordinates(g.phi, g.theta) ;
	const geometry::Vector vecToCenter1(center1[0], center1[1], center1[2]) ;
	const geometry::Vector vecToCenter2 = -vecToCenter1 ;
	
	for (int i = 0; i < GetSize(); i++){
		const VectorType vertex = GetBinCartesianCoordinates(i) ;
		const geometry::Vector vecToVertex(vertex[0], vertex[1], vertex[2]) ;
		
		//NOTICE that ">=" is used here, namely vertices on the equator that is perpendicular to vecToCenter1 are subtracted by both gaussian models
		if (vecToVertex * vecToCenter1 >= 0){
			const double arcLength = std::acos(vecToVertex * vecToCenter1) ;
			bins[i] -= gaussian(g.coefficient, arcLength, g.sigma) ;
		}
		if (vecToVertex * vecToCenter2 >= 0){
			const double arcLength = std::acos(vecToVertex * vecToCenter2) ;
			bins[i] -= gaussian(g.coefficient, arcLength, g.sigma) ;
		}
		
		//Regularize the bin values so that they are >= 0
		if (bins[i] < 0)
			bins[i] = 0 ;
	}
	
}

vtkSmartPointer<vtkPolyData> Counter_SymmetricStatistics::GetVTKStatResult_LevenbergMarquardt(const GaussianOnSphere_Two &g, const geometry::Matrix &m) const{

	if (!g.g1.IsValid() && !g.g2.IsValid())
		throw InvalidDataRange("Counter_SymmetricStatistics::GetVTKStatResult_LevenbergMarquardt : both gaussian models are invalid!") ;
  
	vtkSmartPointer<vtkAppendPolyData> ap = vtkSmartPointer<vtkAppendPolyData>::New() ;
	vtkSmartPointer<vtkPolyData> res = ap->GetOutput() ;
	
	if (g.g1.IsValid()){
		VectorType c1 = GetCartesianCoordinates(g.g1.phi, g.g1.theta) ;
		geometry::Vector center1(c1[0], c1[1], c1[2]) ;
		Cylinder cylin1(2, std::sin(g.g1.sigma)) ;
		ap->AddInput(cylin1.GetVTKPolyData(center1, m)) ;
	}
	
	if (g.g2.IsValid()){
		VectorType c2 = GetCartesianCoordinates(g.g2.phi, g.g2.theta) ;
		geometry::Vector center2(c2[0], c2[1], c2[2]) ;
		Cylinder cylin2(2, std::sin(g.g2.sigma)) ;
		ap->AddInput(cylin2.GetVTKPolyData(center2, m)) ;
	}
	
	ap->Update() ;
	
	return res ;
}

vtkSmartPointer<vtkPolyData> Counter_SymmetricStatistics::GetVTKStatResult_SKmeans(const geometry::Matrix &m) const{

	vtkSmartPointer<vtkAppendPolyData> ap = vtkSmartPointer<vtkAppendPolyData>::New() ;
	vtkSmartPointer<vtkPolyData> res = ap->GetOutput() ;
	
	ap->Update() ;

	return res ;
}
