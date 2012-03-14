#ifndef VOLUMECOUNTER_H
#define VOLUMECOUNTER_H

/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <FiberODF_Common.h>
#include <Counter.h>
#include <FiberSegment.h>
#include <vtkAppendPolyData.h>
#include <teem/nrrd.h>
#include <algorithm>
#include <cmath>
#include <sstream>

const int VOLUME_DIMENSION = 3 ;

class PositionInVolume{
public:
	PositionInVolume(const int x = -1, const int y = -1, const int z = -1) ;
	PositionInVolume& operator=(const PositionInVolume &p) ;

	int &operator[](const int idx) ;
	int operator[](const int idx) const ;
	bool isValid() const ;
	
protected:
	
	int val[VOLUME_DIMENSION] ;
	
	template<typename CounterType>
	friend class VolumeCounter ;
} ;

bool operator==(const PositionInVolume &piv1, const PositionInVolume &piv2) ;
std::ostream &operator<<(std::ostream &os, const PositionInVolume &p) ;

//*******************************************************************************
  
template<typename CounterType>
class VolumeCounter{
	/*
	 * 
	 * NOTICE that the VolumeCounter employs RAS coordinate system
	 * One Counter is assigned to each voxel
	 * 
	 */
public:
	explicit VolumeCounter(const short subdivisionLevel, const char *volumeNrrdFileName) ;
	
	explicit VolumeCounter(const char *volumeCounterNrrdFileName) ;							//Load from the volume counter rasterized to disk
	
	int GetSize() const ;												//Get the number of voxels
	
	int GetNumVoxelsAlongAxis(const int axis) const ;								//Get the number of voxels along a certain axis
    //Get Origin Along Axis
    int GetOriginAlongAxis(const int axis) const;
    //Get Spacing Along Axis
    int GetSpacingAlongAxis(const int axis) const;
	PositionInVolume GetVoxelPositionInVolume(const geometry::Point &p) const ;					//Get location of the voxel where the geometry point resides in
	
	PositionInVolume GetVoxelPositionInVolume(const int x, const int y, const int z) const ;			//Indicate the location of the voxel. The three coordinates are listed in the same order with the space(domain) axes in the volume nrrd file
	
	geometry::Matrix GetTranslationMatrixToVoxel(const PositionInVolume &piv) const ;				//Get the translation matrix from the origin to the indicated voxel
	geometry::Matrix GetTranslationMatrixToVoxel(const int x, const int y, const int z) const ;
	
	CounterType &GetCounterForVoxel(const PositionInVolume &piv) ;							//Get the counter associated with the specified voxel
	CounterType &GetCounterForVoxel(const int x, const int y, const int z) ;
	
	CounterType &GetCounterForResidingVoxel(const geometry::Point &p) ;
	
	vtkSmartPointer<vtkPolyData> GetVTKPolyData(const bool omitEmptyCounters = true) ;
	vtkSmartPointer<vtkPolyData> GetVTKPolyData(const int axis, const int sliceNumber,
						    const bool omitEmptyCounters = true) ;				//Get the indicated slice
	vtkSmartPointer<vtkPolyData> GetVTKPolyData(const PositionInVolume &corner1, const PositionInVolume &corner2,
						    const bool omitEmptyCounters = true) ;				//Get the indicated region
	vtkSmartPointer<vtkPolyData> GetVTKPolyData(const int c1x, const int c1y, const int c1z,
						    const int c2x, const int c2y, const int c2z,
						    const bool omitEmptyCounters = true) ;
	
	void WriteVolumeCounterToVTKFile(const char *fname = "VolumeCounter.vtk", const bool omitEmptyCounters = true) ;
	void WriteVolumeCounterSliceToVTKFile(const int axis, const int sliceNumber,
					      const char *fname = "VolumeCounterSlice.vtk",
					      const bool omitEmptyCounters = true) ;
	void WriteVolumeCounterRegionToVTKFile(const int c1x, const int c1y, const int c1z,
					       const int c2x, const int c2y, const int c2z,
					       const char *fname = "VolumeCounterRegion.vtk",
					       const bool omitEmptyCounters = true) ;

	vtkSmartPointer<vtkPolyData> GetVTKPolyDataWithLocalMaximaMarked(const int axis, const int sliceNumber,
									 const bool omitEmptyCounters = true) ;

	void WriteVolumeCounterSliceWithLocalMaximaMarkedToVTKFile(const int axis, const int sliceNumber,
								   const char *fname = "VolumeCounterSliceWithLocalMaximaMarked.vtk",
								   const bool omitEmptyCounters = true) ;
	
	void WriteVolumeCounterToNrrdFile(const char *fname = "VolumeCounter.nrrd") const  ;
	
	void WriteNumTotalVotesOfEachVoxelToNrrdFile(const char *fname = "NumTotalVotesOfEachVoxel.nrrd") const ;
	
	virtual ~VolumeCounter() ;
	
protected:
	CounterType *counters ;
	geometry::Point origin ;
	geometry::Vector axisDirections[VOLUME_DIMENSION] ;	//Unit vectors here
	double spacingsAlongAxis[VOLUME_DIMENSION] ;
	int numElementsAlongAxis[VOLUME_DIMENSION] ;
	double b_value ;
	
	template<typename T>
	friend std::ostream &operator<<(std::ostream &os, const VolumeCounter<T> &c) ;

private:
	VolumeCounter(const VolumeCounter &) ;
	void operator=(const VolumeCounter &) ;
} ;

template<typename CounterType>
VolumeCounter<CounterType>::VolumeCounter(const char *volumeCounterNrrdFileName){

	const int DATA_DIMENSION = VOLUME_DIMENSION + 1 ;
  
	Nrrd *volumeCounter = nrrdNew() ;
	if (nrrdLoad(volumeCounter, volumeCounterNrrdFileName, NULL)){
		const char *err = biffGetDone(NRRD) ;
		std::cerr << "Failed to load serialized Volume Counter " << volumeCounterNrrdFileName << " : " << err << std::endl ;
		delete err ;
		exit(1) ;
	}
	
	if (static_cast<int>(volumeCounter->dim) != DATA_DIMENSION){
		std::cerr << "The dimension of the serialized data must be " << DATA_DIMENSION << std::endl ;
		exit(1) ;
	}
	
	if (volumeCounter->type != nrrdTypeDouble){
		std::cerr << "The type of rasterized data must be double!" << std::endl ;
		exit(1) ;
	}
	
	if (volumeCounter->space != nrrdSpaceRightAnteriorSuperior){
		std::cerr << "The word coordinate frame of the serialized data must be RAS!" << std::endl ;
		exit(1) ;
	}
	
	//Find the list type axis, namely the gradient axis
	int listAxis = -1 ;
	for (int i = 0; i < DATA_DIMENSION; i++){
		if (volumeCounter->axis[i].kind == nrrdKindList || volumeCounter->axis[i].kind == nrrdKindVector || volumeCounter->axis[i].kind == nrrdKindPoint){
			if (listAxis != -1){
				std::cerr << "Too many list axes in the data!" << std::endl ;
				exit(1) ;
			}
			listAxis = i ;
		}
		else if (volumeCounter->axis[i].kind != nrrdKindDomain && volumeCounter->axis[i].kind != nrrdKindSpace){
			std::cerr << "Unrecognizable axis kind: axis " << i << " is of kind " << volumeCounter->axis[i].kind << std::endl ;
			exit(1) ;
		}
	}
	if (listAxis == -1){
		std::cerr << "Can not find the list axis!" << std::endl ;
		exit(1) ;
	}

	//Compute the permutation
	std::vector<unsigned int> permutation(DATA_DIMENSION) ;
	unsigned int permuteCounter = 0 ;
	permutation[0] = static_cast<unsigned int>(listAxis) ;		//Shift the list axis to the fastest axis
	for (int i = 1; i < DATA_DIMENSION; i++){
		if (permuteCounter == static_cast<unsigned int>(listAxis))
			permuteCounter++ ;
		permutation[i] = permuteCounter++ ;
	}

	Nrrd *temp = nrrdNew() ;
	if (nrrdAxesPermute(temp, volumeCounter, &permutation[0])){		//Permute the data
		const char *txt = biffGetDone(NRRD) ;
		std::cerr << "Failed while permuting the data! : " << txt << std::endl ;
		delete txt ;
		exit(1) ;
	}
	nrrdNuke(volumeCounter) ;
	volumeCounter = temp ;
	temp = NULL ;

	//Extract the key/value info
	int subdivisionLevel ;
	if (nrrdKeyValueSize(volumeCounter) < 3){
		std::cerr << "Wrong number of nrrd key/value pairs!" << std::endl ;
		exit(1) ;
	}
	else{
		char *key, *value ;
		
		nrrdKeyValueIndex(volumeCounter, &key, &value, 0) ;
		if (strcmp(key, "Icosahedron Subdivision Level") != 0){
			std::cerr << "No subdivision level found in file!" << std::endl ;
			exit(1) ;
		}
		else{
			sscanf(value, "%d", &subdivisionLevel) ;
		}
		delete key ;
		delete value ;

		nrrdKeyValueIndex(volumeCounter, &key, &value, 2) ;
		if (strcmp(key, "DWMRI_b-value") != 0){
			std::cerr << "No b value found in file!" << std::endl ;
			exit(1) ;
		}
		else{
			sscanf(value, "%lf", &b_value) ;
		}
		delete key ;
		delete value ;
	}
	
	CounterType::Initialize(subdivisionLevel) ;
	
	//Check the validity of the serialized data
	if (static_cast<int>(volumeCounter->axis[0].size) != CounterType::GetSize()){
		std::cerr << "Number of elements in the rasterized data's first axis and voxel Counter size do not match!" << std::endl ;
		exit(1) ;
	}
	
	//Get geometric info
	for (int i = 0; i < DATA_DIMENSION - 1; i++){
		numElementsAlongAxis[i] = static_cast<int>(volumeCounter->axis[i + 1].size) ;
		origin.getRef(i) = volumeCounter->spaceOrigin[i] ;
		for (int j = 0; j < DATA_DIMENSION - 1; j++){
			axisDirections[i].getRef(j) = volumeCounter->axis[i + 1].spaceDirection[j] ;
		}
		spacingsAlongAxis[i] = axisDirections[i].magnitude() ;
		axisDirections[i].normalize() ;		//To make the axis directions unit vectors
	}	

	counters = new CounterType[GetSize()] ;		//Allocate memory for the counters

	//Get the data
	const double * const data = reinterpret_cast<double *>(volumeCounter->data) ;
	std::vector<double> bins(CounterType::GetSize()) ;
	for (int i = 0; i < GetSize(); i++){
		memcpy(&bins[0], data + CounterType::GetSize() * i, sizeof(double) * CounterType::GetSize()) ;
		counters[i].SetBins(bins) ;
	}
	
	nrrdNuke(volumeCounter) ;
}

template<typename CounterType>
VolumeCounter<CounterType>::VolumeCounter(const short subdivisionLevel, const char *volumeNrrdFileName){
  
	const int DATA_DIMENSION = VOLUME_DIMENSION + 1 ;
  
	CounterType::Initialize(subdivisionLevel) ;		//Initialize the counters
	
	Nrrd *volume = nrrdNew() ;
	if (nrrdLoad(volume, volumeNrrdFileName, NULL)){
		const char *err = biffGetDone(NRRD) ;
		std::cerr << "Failed to open volume " << volumeNrrdFileName << " : " << err << std::endl ;
		delete err ;
		exit(1) ;
	}
	
	if (static_cast<int>(volume->dim) != DATA_DIMENSION){
		std::cerr << "The dimension of data must be " << DATA_DIMENSION << std::endl ;
		exit(1) ;
	}
	
	if (volume->space != nrrdSpaceRightAnteriorSuperior &&
	  volume->space != nrrdSpaceLeftAnteriorSuperior &&
	  volume->space != nrrdSpaceLeftPosteriorSuperior){
		std::cerr << "Can only handle RAS, LAS and LPS world coordinate frames!" << std::endl ;
		exit(1) ;
	}
	
	//Find the list type axis, namely the gradient axis
	int listAxis = -1 ;
	for (int i = 0; i < DATA_DIMENSION; i++){
		if (volume->axis[i].kind == nrrdKindList || volume->axis[i].kind == nrrdKindVector || volume->axis[i].kind == nrrdKindPoint){
			if (listAxis != -1){
				std::cerr << "Too many list axes in the data!" << std::endl ;
				exit(1) ;
			}
			listAxis = i ;
		}
		else if (volume->axis[i].kind != nrrdKindDomain && volume->axis[i].kind != nrrdKindSpace){
			std::cerr << "Unrecognizable axis kind: axis " << i << " is of kind " << volume->axis[i].kind << std::endl ;
			exit(1) ;
		}
	}
	if (listAxis == -1){
		std::cerr << "Can not find the list axis!" << std::endl ;
		exit(1) ;
	}

	//Compute the permutation
	std::vector<unsigned int> permutation(DATA_DIMENSION) ;
	unsigned int permuteCounter = 0 ;
	permutation[0] = static_cast<unsigned int>(listAxis) ;		//Shift the list axis to the fastest axis
	for (int i = 1; i < DATA_DIMENSION; i++){
		if (permuteCounter == static_cast<unsigned int>(listAxis))
			permuteCounter++ ;
		permutation[i] = permuteCounter++ ;
	}

	Nrrd *temp = nrrdNew() ;
	if (nrrdAxesPermute(temp, volume, &permutation[0])){		//Permute the data
		const char *txt = biffGetDone(NRRD) ;
		std::cerr << "Failed while permuting the data! : " << txt << std::endl ;
		delete txt ;
		exit(1) ;
	}
	nrrdNuke(volume) ;
	volume = temp ;
	temp = NULL ;

	//Get b-value from the volume
	bool foundBValue = false ;
	for (unsigned int i = 0; i < nrrdKeyValueSize(volume); i++){
		char *key, *value ;
		nrrdKeyValueIndex(volume, &key, &value, i) ;
		if (strcmp(key, "DWMRI_b-value") == 0){
			foundBValue = true ;
			sscanf(value, "%lf", &b_value) ;
			break ;
		}
		delete key ;
		delete value ;
	}
	if (!foundBValue){
		std::cerr << "Can not find the b value!" << std::endl ;
		exit(1) ;
	}
	
	//Convert to RAS coordinates
	if (volume->space != nrrdSpaceRightAnteriorSuperior){
		const bool usesLAS = (volume->space == nrrdSpaceLeftAnteriorSuperior) ;
		const bool usesLPS = (volume->space == nrrdSpaceLeftPosteriorSuperior) ;
		for (int i = 1; i < DATA_DIMENSION; i++){
			if (usesLAS || usesLPS){
				volume->axis[i].spaceDirection[0] = -volume->axis[i].spaceDirection[0] ;
				volume->measurementFrame[i - 1][0] = -volume->measurementFrame[i - 1][0] ;
			}
			if (usesLPS){
				volume->axis[i].spaceDirection[1] = -volume->axis[i].spaceDirection[1] ;
				volume->measurementFrame[i - 1][1] = -volume->measurementFrame[i - 1][1] ;
			}
		}
    
		if (usesLAS || usesLPS){
			volume->spaceOrigin[0] = -volume->spaceOrigin[0] ;
		}
		if (usesLPS){
			volume->spaceOrigin[1] = -volume->spaceOrigin[1] ;
		}

		volume->space = nrrdSpaceRightAnteriorSuperior ;
	}
	
	//Guarantee that the space(domain) axes are cell centered
	for (int i = 0; i < VOLUME_DIMENSION; i++){
		if (volume->axis[i + 1].center != nrrdCenterCell){
			std::cerr << "Data along the space(domain) axes must be cell centered!" << std::endl ;
			exit(1) ;
		}
	}
	
	//Extract geometric info of the volume
	for (int i = 0; i < DATA_DIMENSION - 1; i++){
		numElementsAlongAxis[i] = static_cast<int>(volume->axis[i + 1].size) ;
		origin.getRef(i) = volume->spaceOrigin[i] ;
		for (int j = 0; j < DATA_DIMENSION - 1; j++){
			axisDirections[i].getRef(j) = volume->axis[i + 1].spaceDirection[j] ;
		}
		spacingsAlongAxis[i] = axisDirections[i].magnitude() ;
		axisDirections[i].normalize() ;		//To make the axis directions unit vectors
	}
	
	nrrdNuke(volume) ;
	
	if (GetSize() <= 0){
		std::cerr << "The volume is empty!" << std::endl ;
		exit(1) ;
	}
	
	counters = new CounterType[GetSize()] ;		//Allocate for the individual Counters
}

template<typename CounterType>
int VolumeCounter<CounterType>::GetSize() const{
	int temp = 1 ;
	for (int i = 0; i < VOLUME_DIMENSION; i++)
		temp *= numElementsAlongAxis[i] ;
	return temp ;
}

template<typename CounterType>
int VolumeCounter<CounterType>::GetNumVoxelsAlongAxis(const int axis) const{
	if (axis < 0 || axis >= VOLUME_DIMENSION)
		throw IndexOutOfRange("VolumeCounter::GetNumVoxelsAlongAxis") ;
	
	return numElementsAlongAxis[axis] ;
}

template<typename CounterType>
int VolumeCounter<CounterType>::GetOriginAlongAxis(const int axis) const{
	if (axis < 0 || axis >= VOLUME_DIMENSION)
		throw IndexOutOfRange("VolumeCounter::GetOriginAlongAxis") ;
	
	return origin[axis] ;
}

template<typename CounterType>
int VolumeCounter<CounterType>::GetSpacingAlongAxis(const int axis) const{
	if (axis < 0 || axis >= VOLUME_DIMENSION)
		throw IndexOutOfRange("VolumeCounter::GetSpacingAlongAxis") ;
	
	return spacingsAlongAxis[axis] ;
}


template<typename CounterType>
PositionInVolume VolumeCounter<CounterType>::GetVoxelPositionInVolume(const geometry::Point &p) const{
	PositionInVolume pos ;
	geometry::Vector v(origin, p) ;
	for (int i = 0; i < VOLUME_DIMENSION; i++){
		const double distanceAlongAxis = v * axisDirections[i] ;
		const int index = static_cast<int>(round(distanceAlongAxis / spacingsAlongAxis[i])) ;
			//Note that round() is used here, which uses the assumption that data along the axis is cell centered
		if (index < 0 || index >= numElementsAlongAxis[i])
			return pos ;
		else
			pos[i] = index ;
	}
	return pos ;
}

template<typename CounterType>
PositionInVolume VolumeCounter<CounterType>::GetVoxelPositionInVolume(const int x, const int y, const int z) const{
	PositionInVolume pos ;
	if (x < 0 || x >= numElementsAlongAxis[0] ||
	    y < 0 || y >= numElementsAlongAxis[1] ||
	    z < 0 || z >= numElementsAlongAxis[2])
		return pos ;
	pos[0] = x ;
	pos[1] = y ;
	pos[2] = z ;
	return pos ;
}

template<typename CounterType>
geometry::Matrix VolumeCounter<CounterType>::GetTranslationMatrixToVoxel(const PositionInVolume &piv) const{
	if (!piv.isValid())
		throw IndexOutOfRange("VolumeCounter::GetTranslationMatrixToVoxel") ;
	
	geometry::Vector transVector(origin[0], origin[1], origin[2]) ;
	for (int i = 0; i < VOLUME_DIMENSION; i++){
		transVector = transVector + axisDirections[i] * spacingsAlongAxis[i] * piv[i] ;
	}
	
	//This factor is to scale the geometry entity for the voxel to proper size
	const double scalingFactor = std::min(std::min(spacingsAlongAxis[0], spacingsAlongAxis[1]), spacingsAlongAxis[2]) / 2 ;
	
	return geometry::translate(transVector[0], transVector[1], transVector[2]) * geometry::scale(scalingFactor, scalingFactor, scalingFactor) ;
}

template<typename CounterType>
geometry::Matrix VolumeCounter<CounterType>::GetTranslationMatrixToVoxel(const int x, const int y, const int z) const{
	return GetTranslationMatrixToVoxel(GetVoxelPositionInVolume(x, y, z)) ;
}

template<typename CounterType>
CounterType &VolumeCounter<CounterType>::GetCounterForVoxel(const PositionInVolume &piv){
	if (!piv.isValid())
		throw IndexOutOfRange("VolumeCounter::GetCounterForVoxel") ;
	/*
	 * NOTICE that the Counters are stored in similar order with the NRRD format, namely the axes are from
	 * fastest to slowest
	 */
	return counters[piv[2] * numElementsAlongAxis[0] * numElementsAlongAxis[1] + piv[1] * numElementsAlongAxis[0] + piv[0]] ;
}

template<typename CounterType>
CounterType &VolumeCounter<CounterType>::GetCounterForVoxel(const int x, const int y, const int z){
	return GetCounterForVoxel(GetVoxelPositionInVolume(x, y, z)) ;
}

template<typename CounterType>
CounterType &VolumeCounter<CounterType>::GetCounterForResidingVoxel(const geometry::Point &p){
	return GetCounterForVoxel(GetVoxelPositionInVolume(p)) ;
}

template<typename CounterType>
vtkSmartPointer<vtkPolyData> VolumeCounter<CounterType>::GetVTKPolyData(const bool omitEmptyCounters){
	vtkSmartPointer<vtkAppendPolyData> ap = vtkSmartPointer<vtkAppendPolyData>::New() ;
	for (int i = 0; i < numElementsAlongAxis[0]; i++)
		for (int j = 0; j < numElementsAlongAxis[1]; j++)
			for (int k = 0; k < numElementsAlongAxis[2]; k++){
				if (!omitEmptyCounters || !GetCounterForVoxel(i, j, k).IsEmpty())
					ap->AddInput(GetCounterForVoxel(i, j, k).GetVTKPolyData(GetTranslationMatrixToVoxel(i, j, k))) ;
			}
	vtkSmartPointer<vtkPolyData> res = ap->GetOutput() ;
	ap->Update() ;
	return res ;
}

template<typename CounterType>
vtkSmartPointer<vtkPolyData> VolumeCounter<CounterType>::GetVTKPolyData(const int axis, const int sliceNumber,
									const bool omitEmptyCounters){
	if (axis < 0 || axis >= VOLUME_DIMENSION || sliceNumber < 0 || sliceNumber >= numElementsAlongAxis[axis])
		throw IndexOutOfRange("VolumeCounter::GetVTKPolyData") ;
	
	vtkSmartPointer<vtkAppendPolyData> ap = vtkSmartPointer<vtkAppendPolyData>::New() ;
	vtkSmartPointer<vtkPolyData> res = ap->GetOutput() ;
	
	if (axis == 0){
		for (int i = 0; i < numElementsAlongAxis[1]; i++)
			for (int j = 0; j < numElementsAlongAxis[2]; j++)
				if (!omitEmptyCounters || !GetCounterForVoxel(sliceNumber, i, j).IsEmpty())
					ap->AddInput(GetCounterForVoxel(sliceNumber, i, j).GetVTKPolyData(GetTranslationMatrixToVoxel(sliceNumber, i, j))) ;
	}
	else if (axis == 1){
		for (int i = 0; i < numElementsAlongAxis[0]; i++)
			for (int j = 0; j < numElementsAlongAxis[2]; j++)
				if (!omitEmptyCounters || !GetCounterForVoxel(i, sliceNumber, j).IsEmpty())
					ap->AddInput(GetCounterForVoxel(i, sliceNumber, j).GetVTKPolyData(GetTranslationMatrixToVoxel(i, sliceNumber, j))) ;
	}
	else{
		for (int i = 0; i < numElementsAlongAxis[0]; i++)
			for (int j = 0; j < numElementsAlongAxis[1]; j++)
				if (!omitEmptyCounters || !GetCounterForVoxel(i, j, sliceNumber).IsEmpty())
					ap->AddInput(GetCounterForVoxel(i, j, sliceNumber).GetVTKPolyData(GetTranslationMatrixToVoxel(i, j, sliceNumber))) ;
	}
	ap->Update() ;
	return res ;
}

template<typename CounterType>
vtkSmartPointer<vtkPolyData> VolumeCounter<CounterType>::GetVTKPolyData(const PositionInVolume &corner1, const PositionInVolume &corner2,
									const bool omitEmptyCounters){
	if (!corner1.isValid() || !corner2.isValid() || corner1[0] == corner2[0] || corner1[1] == corner2[1] || corner1[2] == corner2[2])
		throw IndexOutOfRange("VolumeCounter::GetVTKPolyData") ;
	
	vtkSmartPointer<vtkAppendPolyData> ap = vtkSmartPointer<vtkAppendPolyData>::New() ;
	vtkSmartPointer<vtkPolyData> res = ap->GetOutput() ;
	
	for (int i = std::min(corner1[0], corner2[0]); i <= std::max(corner1[0], corner2[0]); i++)
		for (int j = std::min(corner1[1], corner2[1]); j <= std::max(corner1[1], corner2[1]); j++)
			for (int k = std::min(corner1[2], corner2[2]); k <= std::max(corner1[2], corner2[2]); k++){
				if (!omitEmptyCounters || !GetCounterForVoxel(i, j, k).IsEmpty())
					ap->AddInput(GetCounterForVoxel(i, j, k).GetVTKPolyData(GetTranslationMatrixToVoxel(i, j, k))) ;
			}
	
	ap->Update() ;
	return res ;
}

template<typename CounterType>
vtkSmartPointer<vtkPolyData> VolumeCounter<CounterType>::GetVTKPolyData(const int c1x, const int c1y, const int c1z,
									const int c2x, const int c2y, const int c2z,
									const bool omitEmptyCounters){
	return GetVTKPolyData(GetVoxelPositionInVolume(c1x, c1y, c1z), GetVoxelPositionInVolume(c2x, c2y, c2z), omitEmptyCounters) ;
}

template<typename CounterType>
void VolumeCounter<CounterType>::WriteVolumeCounterToVTKFile(const char *fname, const bool omitEmptyCounters){
	vtkSmartPointer<vtkPolyData> data = GetVTKPolyData(omitEmptyCounters) ;

	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New() ;
	writer->SetInput(data) ;
	writer->SetFileName(fname) ;
	std::stringstream ss ;
	ss << "Volume Counter" ;
	writer->SetHeader(ss.str().c_str()) ;
	writer->Write() ;
}

template<typename CounterType>
void VolumeCounter<CounterType>::WriteVolumeCounterSliceToVTKFile(const int axis, const int sliceNumber, const char *fname,
								  const bool omitEmptyCounters){
	vtkSmartPointer<vtkPolyData> data = GetVTKPolyData(axis, sliceNumber, omitEmptyCounters) ;

	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New() ;
	writer->SetInput(data) ;
	writer->SetFileName(fname) ;
	std::stringstream ss ;
	ss << "Volume Counter Slice : axis " << axis << ", slice " << sliceNumber ;
	writer->SetHeader(ss.str().c_str()) ;
	writer->Write() ;
}

template<typename CounterType>
void VolumeCounter<CounterType>::WriteVolumeCounterRegionToVTKFile(const int c1x, const int c1y, const int c1z,
								   const int c2x, const int c2y, const int c2z,
								   const char *fname, const bool omitEmptyCounters){
	vtkSmartPointer<vtkPolyData> data = GetVTKPolyData(c1x, c1y, c1z, c2x, c2y, c2z, omitEmptyCounters) ;

	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New() ;
	writer->SetInput(data) ;
	writer->SetFileName(fname) ;
	std::stringstream ss ;
	ss << "Volume Counter Region : (" << c1x << ", " << c1y << ", " << c1z << "), (" << c2x << ", " << c2y << ", " << c2z << ")" ;
	writer->SetHeader(ss.str().c_str()) ;
	writer->Write() ;
}

template<typename CounterType>
vtkSmartPointer<vtkPolyData> VolumeCounter<CounterType>::GetVTKPolyDataWithLocalMaximaMarked(const int axis, const int sliceNumber,
											     const bool omitEmptyCounters){
	if (axis < 0 || axis >= VOLUME_DIMENSION || sliceNumber < 0 || sliceNumber >= numElementsAlongAxis[axis])
		throw IndexOutOfRange("VolumeCounter::GetVTKPolyDataWithLocalMaximaMarked") ;
	
	vtkSmartPointer<vtkAppendPolyData> ap = vtkSmartPointer<vtkAppendPolyData>::New() ;
	vtkSmartPointer<vtkPolyData> res = ap->GetOutput() ;
	
	if (axis == 0){
		for (int i = 0; i < numElementsAlongAxis[1]; i++)
			for (int j = 0; j < numElementsAlongAxis[2]; j++)
				if (!omitEmptyCounters || !GetCounterForVoxel(sliceNumber, i, j).IsEmpty())
					ap->AddInput(GetCounterForVoxel(sliceNumber, i, j).GetVTKPolyDataWithLocalMaximaMarked(GetTranslationMatrixToVoxel(sliceNumber, i, j))) ;
	}
	else if (axis == 1){
		for (int i = 0; i < numElementsAlongAxis[0]; i++)
			for (int j = 0; j < numElementsAlongAxis[2]; j++)
				if (!omitEmptyCounters || !GetCounterForVoxel(i, sliceNumber, j).IsEmpty())
					ap->AddInput(GetCounterForVoxel(i, sliceNumber, j).GetVTKPolyDataWithLocalMaximaMarked(GetTranslationMatrixToVoxel(i, sliceNumber, j))) ;
	}
	else{
		for (int i = 0; i < numElementsAlongAxis[0]; i++)
			for (int j = 0; j < numElementsAlongAxis[1]; j++)
				if (!omitEmptyCounters || !GetCounterForVoxel(i, j, sliceNumber).IsEmpty())
					ap->AddInput(GetCounterForVoxel(i, j, sliceNumber).GetVTKPolyDataWithLocalMaximaMarked(GetTranslationMatrixToVoxel(i, j, sliceNumber))) ;
	}
	ap->Update() ;
	return res ;
}

template<typename CounterType>
void VolumeCounter<CounterType>::WriteVolumeCounterSliceWithLocalMaximaMarkedToVTKFile(const int axis, const int sliceNumber, const char *fname,
											const bool omitEmptyCounters){
	vtkSmartPointer<vtkPolyData> data = GetVTKPolyDataWithLocalMaximaMarked(axis, sliceNumber, omitEmptyCounters) ;

	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New() ;
	writer->SetInput(data) ;
	writer->SetFileName(fname) ;
	std::stringstream ss ;
	ss << "Volume Counter Slice with local maxima marked : axis " << axis << ", slice " << sliceNumber ;
	writer->SetHeader(ss.str().c_str()) ;
	writer->Write() ;
}

template<typename CounterType>
void VolumeCounter<CounterType>::WriteVolumeCounterToNrrdFile(const char *fname) const{
	Nrrd *volumeCounter = nrrdNew() ;
	
	if (nrrdAlloc_va(volumeCounter, nrrdTypeDouble, static_cast<unsigned int>(VOLUME_DIMENSION + 1),	//Allocate memory for the rasterized data
		static_cast<size_t>(CounterType::GetSize()),
		static_cast<size_t>(numElementsAlongAxis[0]),
		static_cast<size_t>(numElementsAlongAxis[1]),
		static_cast<size_t>(numElementsAlongAxis[2])
	)){
		const char *err = biffGetDone(NRRD) ;
		std::cerr << "Failed to allocate data for rasterizing VolumeCounter : " << err << std::endl ;
		delete err ;
		exit(1) ;
	}
	
	const char *content = "Volume Counter" ;
	volumeCounter->content = new char[strlen(content) + 1] ;
	strcpy(volumeCounter->content, content) ;
	
	nrrdSpaceSet(volumeCounter, nrrdSpaceRightAnteriorSuperior) ;
	
	for (int i = 0; i < VOLUME_DIMENSION; i++)
		volumeCounter->spaceOrigin[i] = origin[i] ;
	
	volumeCounter->axis[0].center = nrrdCenterUnknown ;
	volumeCounter->axis[0].kind = nrrdKindList ;
	
	for (int i = 0; i < VOLUME_DIMENSION; i++){
		volumeCounter->axis[i + 1].center = nrrdCenterCell ;
		volumeCounter->axis[i + 1].kind = nrrdKindSpace ;
		geometry::Vector dir = axisDirections[i] * spacingsAlongAxis[i] ;
		for (int j = 0; j < VOLUME_DIMENSION; j++){
			volumeCounter->axis[i + 1].spaceDirection[j] = dir[j] ;
		}
	}
	
	//Put the key/value pairs in nrrd
	char valuestr[2000] ;
	char keystr[2000] ;
	sprintf(valuestr, "%d", CounterType::GetIcosahedronSubdivisionLevel()) ;
	nrrdKeyValueAdd(volumeCounter, "Icosahedron Subdivision Level", valuestr) ;
	nrrdKeyValueAdd(volumeCounter, "modality", "DWMRI") ;
	sprintf(valuestr, "%f", b_value) ;
	nrrdKeyValueAdd(volumeCounter, "DWMRI_b-value", valuestr) ;
	vtkSmartPointer<vtkPolyData> counterVTK = counters[0].GetVTKPolyData() ;
	vtkSmartPointer<vtkPoints> vertices = counterVTK->GetPoints() ;
	for (int i = 0; i < CounterType::GetSize(); i++){
		sprintf(keystr, "DWMRI_gradient_%04d", i) ;
		double vertex[3] ;
		vertices->GetPoint(static_cast<vtkIdType>(i), vertex) ;
		sprintf(valuestr, "% .8f    % .8f    % .8f", vertex[0], vertex[1], vertex[2]) ;
		nrrdKeyValueAdd(volumeCounter, keystr, valuestr) ;
	}
	
	//Put the data in nrrd
	double * const data = reinterpret_cast<double *>(volumeCounter->data) ;
	for (int i = 0; i < GetSize(); i++){
		const std::vector<AccumulateType> &counterBin = counters[i].GetBins() ;
		memcpy(data + CounterType::GetSize() * i, &counterBin[0], sizeof(double) * CounterType::GetSize()) ;
	}	
	
	if (nrrdSave(fname, volumeCounter, NULL)){
		const char *err = biffGetDone(NRRD) ;
		std::cerr << "Failed to save rasterized VolumeCounter to " << fname << " : " << err << std::endl ;
		delete err ;
		exit(1) ;
	}
	
	nrrdNuke(volumeCounter) ;
}

template<typename CounterType>
void VolumeCounter<CounterType>::WriteNumTotalVotesOfEachVoxelToNrrdFile(const char *fname) const{
	Nrrd *volumeCounter = nrrdNew() ;
	
	if (nrrdAlloc_va(volumeCounter, nrrdTypeInt, static_cast<unsigned int>(VOLUME_DIMENSION),
		static_cast<size_t>(numElementsAlongAxis[0]),
		static_cast<size_t>(numElementsAlongAxis[1]),
		static_cast<size_t>(numElementsAlongAxis[2])
	)){
		const char *err = biffGetDone(NRRD) ;
		std::cerr << "Failed to allocate data for number of total votes of each voxel : " << err << std::endl ;
		delete err ;
		exit(1) ;
	}
	
	const char *content = "Number of total votes of each voxel" ;
	volumeCounter->content = new char[strlen(content) + 1] ;
	strcpy(volumeCounter->content, content) ;
	
	nrrdSpaceSet(volumeCounter, nrrdSpaceRightAnteriorSuperior) ;
	
	for (int i = 0; i < VOLUME_DIMENSION; i++)
		volumeCounter->spaceOrigin[i] = origin[i] ;
	
	for (int i = 0; i < VOLUME_DIMENSION; i++){
		volumeCounter->axis[i].center = nrrdCenterCell ;
		volumeCounter->axis[i].kind = nrrdKindSpace ;
		geometry::Vector dir = axisDirections[i] * spacingsAlongAxis[i] ;
		for (int j = 0; j < VOLUME_DIMENSION; j++){
			volumeCounter->axis[i].spaceDirection[j] = dir[j] ;
		}
	}
	
	//Put the data in nrrd
	int * const data = reinterpret_cast<int *>(volumeCounter->data) ;
	for (int i = 0; i < GetSize(); i++){
		data[i] = counters[i].GetNumTotalVotes() ;
	}

	if (nrrdSave(fname, volumeCounter, NULL)){
		const char *err = biffGetDone(NRRD) ;
		std::cerr << "Failed to save number of total votes of each voxel to " << fname << " : " << err << std::endl ;
		delete err ;
		exit(1) ;
	}
	
	nrrdNuke(volumeCounter) ;
}

template<typename CounterType>
VolumeCounter<CounterType>::~VolumeCounter(){
	delete [] counters ;
}

template<typename CounterType>
std::ostream &operator<<(std::ostream &os, const VolumeCounter<CounterType> &c){
	os << std::fixed ;
	os << "Origin:" << std::endl ;
	os << c.origin << std::endl ;
	os << "Sizes:" << std::endl ;
	os << c.numElementsAlongAxis[0] << ", " << c.numElementsAlongAxis[1] << ", " << c.numElementsAlongAxis[2] << std::endl ;
	os << "Spacings:" << std::endl ;
	os << c.spacingsAlongAxis[0] << ", " << c.spacingsAlongAxis[1] << ", " << c.spacingsAlongAxis[2] << std::endl ;
	os << "Axis directions:" << std::endl ;
	os << c.axisDirections[0] << std::endl ;
	os << c.axisDirections[1] << std::endl ;
	os << c.axisDirections[2] << std::endl ;
	return os ;
}

#endif
