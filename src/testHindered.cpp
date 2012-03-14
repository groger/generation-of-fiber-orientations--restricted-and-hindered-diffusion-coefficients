//CMAKELIST/////
/*cmake_minimum_required(VERSION 2.6)

project(FiberSig)

option(BUILD_SHARED_LIBS "Build with shared libraries." OFF)

option(BUILD_FIBERODF_AS_LIBRARY "Build fiberodf as a library rather than an executable." OFF)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(Teem REQUIRED)
include(${Teem_USE_FILE})

find_package(VXL REQUIRED)
include(${VXL_CMAKE_DIR}/UseVXL.cmake)

add_definitions(-Wall)

include_directories(${PROJECT_SOURCE_DIR}/Common)
include_directories(${PROJECT_SOURCE_DIR}/Geometry)
include_directories(${PROJECT_SOURCE_DIR})
add_subdirectory(Common)
add_subdirectory(Geometry)

set(FIBERODF_SOURCE SphereIkosahedron.h Counter.h Counter.cc CounterStatistics.h CounterStatistics.cc Cylinder.h Cylinder.cc FiberODF_Common.h FiberODF_Common.cc FiberSegment.h FiberSegment.cc FiberSegmentExtractor.h FiberSegmentExtractor.cc GaussianOnSphere.h GaussianOnSphere.cc VolumeCounter.h VolumeCounter.cc DisplayFunc.h Stats.h FiberOrientation.h)

if(BUILD_FIBERODF_AS_LIBRARY)
  add_library(fiberodf ${FIBERODF_SOURCE})
else(BUILD_FIBERODF_AS_LIBRARY)
  add_executable(fibersig ${FIBERODF_SOURCE} main.cc)
endif(BUILD_FIBERODF_AS_LIBRARY)

target_link_libraries(fibersig common geometry vtkIO vtkGraphics ${ITK_LIBRARIES} teem vnl)*/



#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <string.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
//#include "itkMaskImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkDataArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkMath.h>
#include <vtkTensorGlyph.h>
#include <vtkProperty.h>
#include <vtkDataSet.h>

#include <DisplayFunc.h>
#include <VolumeCounter.h>
#include <FiberSegmentExtractor.h>
//#include <vtkPolyDataReader.h>
#include <Counter.h>

#define DEBUG false
#define COUNTDISP false 

int main(int argc, char **argv){

const unsigned int Dimension = 3;
    typedef std::vector< double > PixelType;
    typedef itk::Image< PixelType, Dimension > ImageType;
    typedef itk::ImageFileWriter< ImageType > WriterType;
    typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
    typedef itk::VectorImage< PixelType , 3 > VectorImageType ;

    VectorImageType::PointType EigenOrigin;
    VectorImageType::SpacingType EigenSpacing;
    VectorImageType::SizeType EigenSize;
    
	for (int i = 0; i < VOLUME_DIMENSION; i++){

	EigenOrigin[i]  = faOrigin[i];
        EigenSpacing[i] = faSpacing[i];
        EigenSize[i] = faSize[i];
    }
	
    
    VectorImageType::Pointer EigenImage = VectorImageType::New();
    //ImageType::Pointer EigenImage = ImageType::New();
    VectorImageType::RegionType EigenRegion;
    EigenRegion.SetSize(EigenSize);
    EigenImage->SetSpacing(EigenSpacing);
    EigenImage->SetOrigin(EigenOrigin);
    EigenImage->SetRegions(EigenRegion);
    EigenImage->Allocate();

    VectorImageType::IndexType pixelIndex;
    VectorImageType::PixelType pixelValue;
 
    
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    std::cout << "Reading " << fiberImgFileName << std::endl;
    reader->SetFileName(fiberImgFileName);
    reader->Update();
 
    // Extract the polydata
    vtkSmartPointer<vtkPolyData> polydata =
    reader->GetOutput();
    double init[9];
    for (i=0;i<9;i++){	init[i]=0;}
    int numbpoints=polydata->GetNumberOfPoints();
    double ResultAndsize[10];
    for (i=0;i<10;i++){	ResultAndsize[i]=0;}
    EigenImage->FillBuffer( init );

//Get the coordinates of points in the polydata 
for(vtkIdType i = 0; i < polydata->GetNumberOfPoints(); i++)
    {
    double p[3];
    polydata->GetPoint(i,p);
    vtkIdType IndexPoly =polydata->FindPoint(p);//index in the polydata of the point
    pixelIndex[0] = p[0]/EigenSpacing[0];//index in EigenImage of the point				 
    pixelIndex[1] = p[1]/EigenSpacing[1];
    pixelIndex[2] = p[2]/EigenSpacing[2];
    double tuple1[9];
    double tuple2[9];
    double temp[10];
    double Result[9];
    polydata->GetPointData()->GetTensors()->GetTuple(IndexPoly,tuple1);//get the tuple of the point
    temp=EigenImage->GetPixel(pixelIndex);//get the value store in the pixel containing the current point
    for (i=0;i<9;i++){	tuple2[i]=temp[i];}
    if (tuple2==init){Result=tuple1;}
    else{ vtkMath::Multiply3x3(tuple1,tuple2,Result);}
    for (i=0;i<9;i++){	ResultAndsize[i]=Result[i];}
    ResultAndsize[9]=temp[9]+1;//we store the number of tuples multiplied in each voxel
    EigenImage->SetPixel(pixelIndex,ResultAndsize);
    // This is identical to:
    // polydata->GetPoints()->GetPoint(i,p);
    //std::cout << "Point " << i << " : (" << p[0] << " " << p[1] << " " << p[2] << ")" << std::endl;
    }
 // Write the EigenImage image for Hindered part to a text file
    output("Writing the image to text file...please be patient");
 	std::ofstream Eigen(EigenFileName);
    Eigen<<"begin header"<<endl;
    Eigen<<"dimension="<<VOLUME_DIMENSION<<endl;
    Eigen<<"type=double"<<endl;
    Eigen<<"size="<<EigenSize[0]<<" "<<EigenSize[1]<<" "<<EigenSize[2]<<endl;
    Eigen<<"origin="<<foOrigin[0]<<" "<<EigenOrigin[1]<<" "<<EigenOrigin[2]<<endl;
    Eigen<<"spacing="<<foSpacing[0]<<" "<<EigenSpacing[1]<<" "<<EigenSpacing[2]<<endl;
    Eigen<<"end header"<<endl;
    Eigen<<endl;
    Eigen<<"begin data"<<endl;

//To have the index of a point in the polydata
    long int size = 0;
    long int size2 = 0;
    ConstIteratorType Eigen_it( EigenImage, EigenImage->GetLargestPossibleRegion());
    Eigen_it.GoToBegin();
    
    while(!Eigen_it.IsAtEnd()){
	double *w = new double[3];
	double temp[9];
	pixelIndex = Eigen_it.GetIndex();
        pixelValue = Eigen_it.Get();//on recupere le tuple average stocke precedemment
	Eigen<<"voxel="<<pixelIndex[0]<<" "<<pixelIndex[1]<<" "<<pixelIndex[2]<<endl;
        size = pixelValue[9];//le nombre de tuple multiplies precedemment est stocke dans le dernier element du vecteur pixelvalue
	for(i=0;i<9;i++){temp[i]=pow(pixelValue[i],(1/size));}//on calcule l'average definitif par voxel}
	w=EstimateEigenValue(temp);//we estimate eigenvalues
	EigenImage->SetPixel(pixelIndex,w);
	pixelValue = Eigen_it.Get();//on check si la nouvelle valeur du pixel est bien un vector de taille 3 contenant les eigenvalues de notre average de tuples
	size2 = pixelValue.size();
	Eigen<<"size="<<size2<<endl;
	if(size2==3){std::cout<<"on a la bonne taille de voxel"<<std::endl;}
	else{std::cout<<"mauvaise taille de voxel"<<std::endl;}
	for(unsigned int i = 0; i < size2; i++){
            Eigen<<pixelValue[i]<<endl;
        }
        ++Eigen_it;
    }
    Eigen<<"end data"<<endl;
    return 0;
}	
