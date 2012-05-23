#ifndef FIBERORIENTATION_H
#define FIBERORIENTATION_H


#include <itkImageIOBase.h>
#include "itkArray2D.h"
#include "itkImageDuplicator.h"
#include <itkImageFileWriter.h>
#include "itkVariableSizeMatrix.h"
#include "itkAffineTransform.h"
#include <itksys/SystemTools.hxx>
#include <itkImage.h>
#include <itkVectorImage.h> 
#include <itkVariableLengthVector.h>
#include <itkMetaDataObject.h> 
#include <itkVector.h>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include <itkVectorImage.h> 
#include <itkVariableLengthVector.h>
#include <itkImageAdaptor.h>
#include <itkPoint.h>

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
#include <Counter.h>
#include <sstream>
#include <stdio.h> 
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include <cmath>
#define _USE_MATH_DEFINES
#include "math.h"
#include <string.h>
#define DEBUG false
#define COUNTDISP false 

//Fiber Orientation Image meta info struct 
struct FOMeta{
    double origin[3];
    double spacing[3];
    unsigned int size[3];
};

// FA Image
const unsigned int FADimension = 3;
typedef double FAPixelType;
typedef itk::Image< FAPixelType, FADimension > FAImageType;
typedef itk::ImageFileReader<FAImageType> FAReaderType;

class FiberOrientation{
	public:
	explicit FiberOrientation(const char* fiberFile, const char* dwiFile, const char* faFile, const char* maskFile);
		
	double *EstimateEigenValue(double *tuple);
        int writeFiberOrientationsPerVoxelTextFile(const char* foImgFileName, const char* EigenFileName,const char* EigenHinderedFileName);
        int readFOMeta(const char* foImgFileName, FOMeta *m);

        template <class ImagePointer>
        int getFOImageHandler(ImagePointer &foImage, const char* foImageFileName);

        template <class ImagePointer>
        int updatePixelValue(ImagePointer &foImage, PositionInVolume p, geometry::Vector direction);
	
	template <class ImagePointer>
        int updatePixelValue3(ImagePointer &foImage, PositionInVolume p, geometry::Vector direction);
  	
	template <class ImagePointer>
	int updatePixelValue2(ImagePointer &foImage, PositionInVolume p, double tuple[9]);

		virtual ~FiberOrientation(){}
		
	protected:
        //Filenames
        const char* fiberImgFileName;
        const char* dwiImgFileName;
        const char* faImgFileName;
        const char* maskImgFileName;
        const char* fiberSegmentsPerVoxelImgFileName;
	const char* EigenFileName;
	const char* EigenHinderedFileName;
        //FA Image
        FAReaderType::Pointer faReader;
        FAImageType::Pointer faImage;
        FAImageType::PointType faOrigin;
        FAImageType::SpacingType faSpacing;
        FAImageType::SizeType faSize;

        //Mask Image
        FAReaderType::Pointer maskReader;
        FAImageType::Pointer maskImage;

        //Max segments in fiber file
        long int maxFiberSegments;

        //Subdivision Level
        int subdivisionLevel;
} ;

FiberOrientation::FiberOrientation(const char *fiberFile, const char* dwiFile, const char* faFile, const char* maskFile){
    //Assign default values for variables
    maxFiberSegments = 0;
    //Assign the subdivisionLevel
    subdivisionLevel = 6;

    //Assign filenames for future use
    fiberImgFileName = fiberFile;
    dwiImgFileName = dwiFile;
    faImgFileName = faFile;
    maskImgFileName = maskFile;

    faReader = ReaderType::New();
    faReader->SetFileName(faImgFileName);
    faReader->Update();
    faImage = faReader->GetOutput();
    
    maskReader = ReaderType::New();
    maskReader->SetFileName(maskImgFileName);
    maskReader->Update();
    maskImage = maskReader->GetOutput();

    faOrigin = faImage->GetOrigin();
    faSpacing = faImage->GetSpacing();
    faSize = faImage->GetLargestPossibleRegion().GetSize();
}



double * FiberOrientation::EstimateEigenValue(double *tuple){

       double **a;//the matrix which will contain tensors
       a = new double*[3];
       double *w = new double[3];//create a vector for the eigenvalues to be stored in
       double **v;//create a matrix for the eigenvectors to be stored in
       v = new double*[3];
 
       for( unsigned int j=0; j<3; j++)
       {
         a[j] = new double[3];
         v[j] = new double[3];
       }
 
       for( int nl=0; nl<3; nl++)
       {
         for( int nc=0; nc<3; nc++)
         {
           a[nl][nc] = tuple[nl*3+nc];//we "reshape" tuple to make a 3*3matrix
         }
       }

        vtkMath::Jacobi (a, w, v); 
	return w;
}



int FiberOrientation::writeFiberOrientationsPerVoxelTextFile(const char* foImgFileName, const char* EigenFileName,const char* EigenHinderedFileName){
    const unsigned int Dimension = 3;
    typedef std::vector< double > PixelType;
    typedef itk::Image< PixelType, Dimension > ImageType;
    typedef itk::ImageFileWriter< ImageType > WriterType;
    typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;

    ImageType::PointType foOrigin;
    ImageType::SpacingType foSpacing;
    ImageType::SizeType foSize;
    ImageType::PointType EigenOrigin;
    ImageType::SpacingType EigenSpacing;
    ImageType::SizeType EigenSize;
    ImageType::PointType EigenHinderedOrigin;
    ImageType::SpacingType EigenHinderedSpacing;
    ImageType::SizeType EigenHinderedSize;

	for (int i = 0; i < VOLUME_DIMENSION; i++){
	foOrigin[i]  = faOrigin[i];
        foSpacing[i] = faSpacing[i];
        foSize[i] = faSize[i];
	EigenOrigin[i]  = faOrigin[i];
        EigenSpacing[i] = faSpacing[i];
        EigenSize[i] = faSize[i];
	EigenHinderedOrigin[i]  = faOrigin[i];
        EigenHinderedSpacing[i] = faSpacing[i];
        EigenHinderedSize[i] = faSize[i];

    }
	
    ImageType::Pointer foImage = ImageType::New();
    ImageType::RegionType foRegion;
    foRegion.SetSize(foSize);
    foImage->SetSpacing(foSpacing);
    foImage->SetOrigin(foOrigin);
    foImage->SetRegions(foRegion);
    foImage->Allocate();

    ImageType::Pointer EigenImage = ImageType::New();
    ImageType::RegionType EigenRegion;
    EigenRegion.SetSize(EigenSize);
    EigenImage->SetSpacing(EigenSpacing);
    EigenImage->SetOrigin(EigenOrigin);
    EigenImage->SetRegions(EigenRegion);
    EigenImage->Allocate();

    
    ImageType::Pointer EigenHinderedImage = ImageType::New();
    std::vector< double >  initialValue(10,0);
    EigenHinderedImage->FillBuffer( initialValue );
    ImageType::RegionType EigenHinderedRegion;
    EigenHinderedRegion.SetSize(EigenHinderedSize);
    EigenHinderedImage->SetSpacing(EigenHinderedSpacing);
    EigenHinderedImage->SetOrigin(EigenHinderedOrigin);
    EigenHinderedImage->SetRegions(EigenHinderedRegion);
    EigenHinderedImage->Allocate(); 

    ImageType::IndexType pixelIndex;
    ImageType::PixelType pixelValue;
    ImageType::PixelType pixelValue2;
    
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    std::cout << "Reading " << fiberImgFileName << std::endl;
    reader->SetFileName(fiberImgFileName);
    reader->Update();
 
    // Extract the polydata
    vtkSmartPointer<vtkPolyData> polydata =
    reader->GetOutput();
    

    // Pass the vtk polydata fiber file which is generated from ukf tractography 
    output("Reading fiber file started");
    FiberSegmentExtractor fse(fiberImgFileName) ;
    output("Reading fiber file finished");
    // Pass the original DWI file 
    output("Reading DWI file started");
    VolumeCounter<Counter_WeightedVertices> vcc(6, dwiImgFileName) ;
    output("Reading DWI file finished");
    // Compute the total number of fibers 
    int totalfibers = (fse.GetNumberOfFibers());
    output2("Number of fibers = ", totalfibers);
    // Compute the total number of segments 
    int totalSegments = fse.GetNumberOfSegments();
    output2("Number of segments = ", totalSegments);

    output("Computing Fiber Orientations per voxel");
    if(COUNTDISP){
        cout<<"No. of Fiber Segments read = ";
    }
    long int segmentCount = 1, wildCount = 0;
    bool fiberStart = true, fiberEnd = false;
    FiberSegment currSegment, prevSegment, nextSegment;
    PositionInVolume point1, point2;

    while(fse){
	/*we define the current segment(=a line of a fiber) which is between two points: point1 and point2*/
        currSegment = fse.Next(fiberEnd);
        point1 = vcc.GetVoxelPositionInVolume(currSegment.p1);/*we define point1 and point2*/
        point2 = vcc.GetVoxelPositionInVolume(currSegment.p2);
	double tuple1[9];double tuple2[9];/*we define two tuples (tensors) at point1 and point2 that we take from the vtk file*/
	double x[3];double y[3];double *w = new double[3];
    	x[0]=currSegment.p1[0];x[1]=currSegment.p1[1];x[2]=currSegment.p1[2];
	y[0]=currSegment.p2[0];y[1]=currSegment.p2[1];y[2]=currSegment.p2[2];
    	vtkIdType Index1 =polydata->FindPoint(x); 
	vtkIdType Index2 =polydata->FindPoint(y);	

        polydata->GetPointData()->GetTensors()->GetTuple(Index1,tuple1);//Get the data tuple at ith location by filling in a user-provided array, we store it in 'tuple'
	polydata->GetPointData()->GetTensors()->GetTuple(Index2,tuple2);

        if (!point1.isValid() || !point2.isValid()){
            wildCount++ ;
            continue ;
        }
        //Print debug
        if(DEBUG){
            displayPoint(point1, "p1 = ");
            displayPoint(point2, "p2 = ");
            displayFiberSegment(currSegment);
        }
        
        //Handle special case of fiberStart
        if(fiberStart){	/*if we are on the first segment of a fiber*/ 
            geometry::Vector direction = (currSegment.p2 - currSegment.p1);/*we define fiber orientation vector*/
	    geometry::Vector eigenvalue;
	    w=EstimateEigenValue(tuple1);/*we estimate eigenvalues of the tuple1(=tensor of the point1)*/
	    eigenvalue[0]=w[0];eigenvalue[1]=w[1];eigenvalue[2]=w[2];/*we store these eigenvalues in eigenvalue vector*/
            updatePixelValue(foImage, point1, direction);/*we update our foImage= file containing fiber orientation vectors*/
	    updatePixelValue3(EigenImage, point1, eigenvalue);/*we update our EigenImage= file containing eigenvalues of each fiber*/
	    updatePixelValue2(EigenHinderedImage, point1, tuple1);/*we update our EigenHinderedImage= file containing eigenvalues and eigenvectors of tensors average*/
            fiberStart = false;
            prevSegment = currSegment;
            segmentCount++;
            continue;
        }

        // Normal case of computing tangential vector at current fiber point using neighbor points
        geometry::Vector direction = (currSegment.p2 - prevSegment.p1);/*we define fiber orientation vector*/
        geometry::Vector eigenvalue;
        w=EstimateEigenValue(tuple1);/*we estimate eigenvalues of the tuple1(=tensor of the point1)*/
	eigenvalue[0]=w[0];eigenvalue[1]=w[1];eigenvalue[2]=w[2];/*we store these eigenvalues in eigenvalue vector*/
        updatePixelValue(foImage, point1, direction);/*we update our foImage= file containing fiber orientation vectors*/
	updatePixelValue3(EigenImage, point1, eigenvalue);/*we update our EigenImage= file containing eigenvalues of each fiber*/
	updatePixelValue2(EigenHinderedImage, point1, tuple1);/*we update our EigenHinderedImage= file containing eigenvalues and eigenvectors of tensors average*/
        prevSegment = currSegment;
        
        //Handle special case of fiberEnd
        if(fiberEnd){/*if we are at the end of a fiber*/
            geometry::Vector direction = (currSegment.p2 - currSegment.p1);/*we define fiber orientation vector*/
	    geometry::Vector eigenvalue;
	    w=EstimateEigenValue(tuple2);/*we estimate eigenvalues of the tuple2(=tensor of the point2)*/
	    eigenvalue[0]=w[0];eigenvalue[1]=w[1];eigenvalue[2]=w[2];/*we store these eigenvalues in eigenvalue vector*/
            updatePixelValue(foImage, point2, direction);/*we update our foImage= file containing fiber orientation vectors*/
	    updatePixelValue3(EigenImage, point2, eigenvalue);/*we update our EigenImage= file containing eigenvalues of each fiber*/
	    updatePixelValue2(EigenHinderedImage, point2, tuple2);/*we update our EigenHinderedImage= file containing eigenvalues and eigenvectors of tensors average*/
            fiberStart = true;
            fiberEnd = false;
        }
        segmentCount++ ;
        if(COUNTDISP){
            display(segmentCount);
        }
    }

    if(COUNTDISP){
        cout<<endl;
    }
    // Format of Fiber Orientation File
    // ================================
    // begin header
    // dimension=3
    // type=double
    // size=maxx maxy maxz
    // origin=ox oy oz
    // spacing=sx sy sz
    // end header
    //
    // begin data
    // voxel=0 0 0
    // size=3
    // 0.1
    // 0.8
    // 0.1
    // voxel=1 0 0
    // ...
    // ...
    // end data
    
     ///////////// Write the Fiber orientation image to a text file//////////////
    output("Writing the image to text file...please be patient");
 	std::ofstream fo(foImgFileName);
    fo<<"begin header"<<endl;
    fo<<"dimension="<<VOLUME_DIMENSION<<endl;
    fo<<"type=double"<<endl;
    fo<<"size="<<foSize[0]<<" "<<foSize[1]<<" "<<foSize[2]<<endl;
    fo<<"origin="<<foOrigin[0]<<" "<<foOrigin[1]<<" "<<foOrigin[2]<<endl;
    fo<<"spacing="<<foSpacing[0]<<" "<<foSpacing[1]<<" "<<foSpacing[2]<<endl;
    fo<<"end header"<<endl;
    fo<<endl;
    fo<<"begin data"<<endl;

    long int size = 0;
    ConstIteratorType fo_it( foImage, foImage->GetLargestPossibleRegion());
    fo_it.GoToBegin();
    while(!fo_it.IsAtEnd()){
        pixelIndex = fo_it.GetIndex();
        pixelValue = fo_it.Get();
        size = pixelValue.size();
        fo<<"voxel="<<pixelIndex[0]<<" "<<pixelIndex[1]<<" "<<pixelIndex[2]<<endl;
        fo<<"size="<<size<<endl;
        for(unsigned int i = 0; i < size; i++){
            fo<<pixelValue[i]<<endl;
        }
        ++fo_it;
    }
    fo<<"end data"<<endl;
    
     /////////////// Write the EigenImage image to a text file/////////////////////
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

    long int size2 = 0;
    ConstIteratorType Eigen_it( EigenImage, EigenImage->GetLargestPossibleRegion());
    Eigen_it.GoToBegin();
    while(!Eigen_it.IsAtEnd()){
        pixelIndex = Eigen_it.GetIndex();
        pixelValue = Eigen_it.Get();
        size2 = pixelValue.size();
        Eigen<<"voxel="<<pixelIndex[0]<<" "<<pixelIndex[1]<<" "<<pixelIndex[2]<<endl;
        Eigen<<"size="<<size2<<endl;
        if(size2!=0){std::cout<<"new pixel "<<std::endl;}
        for(unsigned int i = 0; i < size2; i++){
            Eigen<<pixelValue[i]<<endl;
            //std::cout<<"pixelvalue : "<<pixelValue[i]<<std::endl;
        }
        ++Eigen_it;
    }
    Eigen<<"end data"<<endl; 

    ///////////// Write the EigenImage Hindered part image to a text file///////////////
    output("Writing the image to text file...please be patient");
    std::ofstream EigenHindered(EigenHinderedFileName);
    EigenHindered<<"begin header"<<endl;
    EigenHindered<<"dimension="<<VOLUME_DIMENSION<<endl;
    EigenHindered<<"type=double"<<endl;
    EigenHindered<<"size="<<EigenHinderedSize[0]<<" "<<EigenHinderedSize[1]<<" "<<EigenHinderedSize[2]<<endl;
    EigenHindered<<"origin="<<EigenHinderedOrigin[0]<<" "<<EigenHinderedOrigin[1]<<" "<<EigenHinderedOrigin[2]<<endl;
    EigenHindered<<"spacing="<<EigenHinderedSpacing[0]<<" "<<EigenHinderedSpacing[1]<<" "<<EigenHinderedSpacing[2]<<endl;
    EigenHindered<<"end header"<<endl;
    EigenHindered<<endl;
    EigenHindered<<"begin data"<<endl;

//To have the index of a point in the polydata
    double size3 = 0;
    double size4 = 0; 
    typedef itk::ImageRegionIterator< ImageType > IteratorType;
    IteratorType EigenHindered_it( EigenHinderedImage, EigenHinderedImage->GetLargestPossibleRegion().GetSize() );
    EigenHindered_it.GoToBegin();

    while(!EigenHindered_it.IsAtEnd()){//through the EigenImage
	std::vector<double> pixelValue(10);
	pixelIndex = EigenHindered_it.GetIndex();
        pixelValue = EigenHindered_it.Get();//we take the value of the voxel(a tuple) stored before 

	EigenHindered<<"voxel="<<pixelIndex[0]<<" "<<pixelIndex[1]<<" "<<pixelIndex[2]<<endl;
        size3 = pixelValue[9];//size is stored in the last component of the tuple(pixelValue) stored in the voxel
	std::vector< double > eigeninfo(6); //#Dxx # Dxy # Dxz# Dyy # Dyz # Dzz 
	if(size3!=0){//if there is a tensor for this voxel(voxel not empty)
		
		for (int i=0; i < 3; i++){
			eigeninfo[i] = pixelValue[i];
		}
		for (int i=4; i < 6; i++){
			eigeninfo[i-1] = pixelValue[i];
		}
		eigeninfo[5] = pixelValue[8];
		/*
		double matrix[3][3];double V[3][3];double w[3];
		for( int nl=0; nl<3; nl++)
		{
			for( int nc=0; nc<3; nc++)
			{
				matrix[nl][nc] = pixelValue[nl*3+nc];//we "reshape" tuple(tensor) to make a 3*3matrix
				
			}
		}
		*/
		//vtkMath::Diagonalize3x3(matrix,w,V);/*we estimate eigenvalues(w) and eigenvectors (V) of our current tensor(average tensor)*/

		/*for(int i=0;i<3;i++){
			if(w[i]<0){ w[i]=fabs(w[i]);}
		}
		
		std::vector< double > wv(3);
		for(int i=0;i<3;i++){
			wv[i]=pow(w[i],(1/size3));//we apply ^(1/numbertensormultipliedinthisvoxel) on singular values
			}

		//we will store first 3eigenvalues then first eigenvector then second eigenvector then third eigenvector
		for(int i=0;i<3;i++){eigeninfo[i]=wv[i];}//first we store three eigenvalues in eigeninfo vector*/
		
		/*for(int j=0;j<3;j++){
			for(int i=0;i<3;i++){
				eigeninfo[3+3*j+i]=V[i][j];//then we store the three eigenvectors in eigeninfo
			}
		}*/
		
		
		
		EigenHindered_it.Set(eigeninfo);//we store eigenvalues and eigenvectors results in EigenImage
		//size4=12;/*this is the size of information containing in a voxel in our EigenImage(3coordinates of ;eigenvalues, three eigenvectors=12)*/
		size4=6;//#Dxx # Dxy # Dxz# Dyy # Dyz # Dzz 
	}

	else{ //if there is no tensor in the voxel (empty)
		size4=0;
              }

	EigenHindered<<"size="<<size4<<endl;
	pixelValue = EigenHindered_it.Get();	
	for(unsigned int i = 0; i < size4; i++){
            EigenHindered<<pixelValue[i]<<endl;
        }
        ++EigenHindered_it;
    }
    EigenHindered<<"end data"<<endl;
    return 0;
}



 /*specific function to get information from our txt files containing fiber orientation, eigenvalues and eigenvectors*/ 
template <class ImagePointer>
int FiberOrientation::getFOImageHandler(ImagePointer &foImage, const char* foImgFileName){
    //Read the Fiber Orientation Image
    const unsigned int FODimension = 3;
    typedef std::vector< double > FOPixelType;
    typedef itk::Image< FOPixelType, FODimension > FOImageType;

    typedef itk::ImageFileWriter< FOImageType > FOWriterType;
    typedef itk::ImageRegionConstIterator< FOImageType > FOConstIteratorType;

    FOImageType::PointType foOrigin;
    FOImageType::SpacingType foSpacing;
    FOImageType::SizeType foSize;

    FOMeta *meta = new FOMeta();
    readFOMeta(foImgFileName, meta);

    for(int i = 0; i < VOLUME_DIMENSION; i++){
        foOrigin[i] = meta->origin[i];
        foSpacing[i] = meta->spacing[i];
        foSize[i] = meta->size[i];
    }
    
   // FOImageType::Pointer foImage = FOImageType::New();
    FOImageType::RegionType foRegion;
    foRegion.SetSize(foSize);
    foImage->SetSpacing(foSpacing);
    foImage->SetOrigin(foOrigin);
    foImage->SetRegions(foRegion);
    foImage->Allocate();

    FOImageType::IndexType pixelIndex;
    FOImageType::PixelType pixelValue;

    pixelIndex[0] = 1, pixelIndex[1] = 1; pixelIndex[2] = 1;
    pixelValue = foImage->GetPixel(pixelIndex);
    pixelValue.push_back(0.234);
    foImage->SetPixel(pixelIndex, pixelValue);

    ifstream inFile(foImgFileName);
    std::string line;
    bool headerFinish = false;
    bool keyFinish = false;
    unsigned int i, count;
    int j, k;
    char value[100];

    while(getline(inFile, line)){
        if(line.find("end header") != std::string::npos){
            headerFinish = true; 
        }
        if(!headerFinish){
            continue;
        }
        
        //Process voxels
        if(line.find("voxel") != std::string::npos){
            keyFinish = false;
            j = 0;
            k = 0;
            for(i = 0; i < line.length(); i++){
               if(line[i] == '='){
                   keyFinish = true;
                   continue;
               }
               if(!keyFinish){
                   continue;
               }
               if(line[i] == ' '){
                   value[j] = '\0';
                   j = 0;
                   pixelIndex[k++] = (unsigned int)atoi(value);
               }else{
                   value[j++] = line[i];
               }
            }
            value[j] = '\0';
            pixelIndex[k] = (unsigned int)atoi(value);

            pixelValue = foImage->GetPixel(pixelIndex);

            //Compute size 
            getline(inFile, line);
            if(line.find("size") != std::string::npos){
                keyFinish = false;
                j = 0;
                for(i = 0; i < line.length(); i++){
                   if(line[i] == '='){
                       keyFinish = true;
                       continue;
                   }
                   if(!keyFinish){
                       continue;
                   }
                   value[j++] = line[i];
                }
                value[j] = '\0';
                count = (unsigned int)atoi(value);
            }
            while(count > 0){
                getline(inFile, line);
                pixelValue.push_back(strtod(line.c_str(), NULL));
                count--;
            }
            foImage->SetPixel(pixelIndex, pixelValue);
        }
    }
    return 0;
}
/*specific function to read our txt files containing fiber orientation, eigenvalues and eigenvectors*/
int FiberOrientation::readFOMeta(const char* foImgFileName, FOMeta* m){
    ifstream inFile(foImgFileName);
    std::string line;
    bool keyFinish = false;
    unsigned int i;
    int j, k;
    char value[100];

    while(getline(inFile, line)){
        if(line.find("size") != std::string::npos){
            keyFinish = false; 
            j = 0;
            k = 0;
            for(i = 0; i < line.length(); i++){
               if(line[i] == '='){
                   keyFinish = true;
                   continue;
               }
               if(!keyFinish){
                   continue;
               }
               if(line[i] == ' '){
                   value[j] = '\0';
                   j = 0;
                   m->size[k++] = (unsigned int)atoi(value);
               }else{
                   value[j++] = line[i];
               }
            }
            value[j] = '\0';
            m->size[k] = (unsigned int)atoi(value);
        }

        if(line.find("origin") != std::string::npos){
            keyFinish = false; 
            j = 0;
            k = 0;
            for(i = 0; i<line.length(); i++){
               if(line[i] == '='){
                   keyFinish = true;
                   continue;
               }
               if(!keyFinish){
                   continue;
               }
               if(line[i] == ' '){
                   value[j] = '\0';
                   j = 0;
                   m->origin[k++] = strtod(value, NULL);
               }else{
                   value[j++] = line[i];
               }
            }
            value[j] = '\0';
            m->origin[k] = strtod(value, NULL);
        }
        if(line.find("spacing") != std::string::npos){
            keyFinish = false; 
            j = 0;
            k = 0;
            for(i = 0; i < line.length(); i++){
               if(line[i] == '='){
                   keyFinish = true;
                   continue;
               }
               if(!keyFinish){
                   continue;
               }
               if(line[i] == ' '){
                   value[j] = '\0';
                   j = 0;
                   m->spacing[k++] = strtod(value, NULL);
               }else{
                   value[j++] = line[i];
               }
            }
            value[j] = '\0';
            m->spacing[k] = strtod(value, NULL);
        }
        if(line.find("end header") != std::string::npos){
            break;
        }
    }
    return 0;
}

template <class ImagePointer>/*function to update pixel(to write) for the fiber orientation image*/
int FiberOrientation::updatePixelValue(ImagePointer &foImage, PositionInVolume p, geometry::Vector direction){
    const unsigned int Dimension = 3;
    typedef std::vector< double > PixelType;
    typedef itk::Image< PixelType, Dimension > ImageType;

    ImageType::IndexType pixelIndex;
    ImageType::PixelType pixelValue;

    if(!direction.isZero()){/*we normalize our vector direction*/
        direction.normalize();
        if(DEBUG){
            displayVect(direction, "Direction Vector = ");
        }
        //Compute the pixel position, retrieve the vector, update it and set the pixel value again
        for(int i = 0; i < VOLUME_DIMENSION; i++){
            pixelIndex[i] = (int)p[i];
        }
        pixelValue = foImage->GetPixel(pixelIndex);

        for(int i = 0; i < VOLUME_DIMENSION; i++){
            pixelValue.push_back(direction[i]);
        }
        foImage->SetPixel(pixelIndex, pixelValue);
	
	pixelValue = foImage->GetPixel(pixelIndex);

    }
    return 0;
}
template <class ImagePointer>/*function to update pixel(to write) for the restricted diffusion*/
int FiberOrientation::updatePixelValue3(ImagePointer &foImage, PositionInVolume p, geometry::Vector direction){
    const unsigned int Dimension = 3;
    typedef std::vector< double > PixelType;
    typedef itk::Image< PixelType, Dimension > ImageType;

    ImageType::IndexType pixelIndex;
    ImageType::PixelType pixelValue;

    if(!direction.isZero()){
       
        //Compute the pixel position, retrieve the vector, update it and set the pixel value again
        for(int i = 0; i < VOLUME_DIMENSION; i++){
            pixelIndex[i] = (int)p[i];
        }
        pixelValue = foImage->GetPixel(pixelIndex);
       
        for(int i = 0; i < VOLUME_DIMENSION; i++){
            pixelValue.push_back(direction[i]);
        }
        foImage->SetPixel(pixelIndex, pixelValue);
	
	pixelValue = foImage->GetPixel(pixelIndex);	

    }
    return 0;
}

template <class ImagePointer>/*function to update pixel for the hindered diffusion: we estimate tensor average and store the number of tensor that we average in each voxel*/
int FiberOrientation::updatePixelValue2(ImagePointer &foImage, PositionInVolume p, double tuple[9]){
    
    int verif=0;
    const unsigned int Dimension = 3;
    typedef std::vector< double > PixelType;
    typedef itk::Image< PixelType, Dimension > ImageType;
    std::vector< double > ResultAndsize(10);
    std::vector< double > Result(9);
    ImageType::IndexType pixelIndex;
    ImageType::PixelType pixelValue1;
        //Compute the pixel position, retrieve the vector, update it and set the pixel value again
        for(int i = 0; i < VOLUME_DIMENSION; i++){
            pixelIndex[i] = (int)p[i];
        }
        std::vector< double >  pixelValue(10);
	pixelValue = foImage->GetPixel(pixelIndex);

	for (int i=0;i<9;i++)
	{
		if (pixelValue[i]==0){ ++verif;}
	}
	

       if (verif==9){//if tuple2=0, the current voxel is empty in EigenImage, we are on the first tensor of the voxel, we keep it and store it in the voxel

		for (int i=0;i<9;i++){
			ResultAndsize[i]=tuple[i];
			}
		ResultAndsize[9]=1;//the size stored is one, we have only one tensor
		foImage->SetPixel(pixelIndex,ResultAndsize);//we store the tuple in EigenImage
	}

	else{
	double matrix1[3][3];double matrix2[3][3];double matrixResult[3][3];
	//reshape tuple in matrix
	 for( int nl=0; nl<3; nl++)
       {
         for( int nc=0; nc<3; nc++)
         {
           matrix1[nl][nc] = tuple[nl*3+nc];//we "reshape" tuple1 to make a 3*3matrix
         }
       }
	for( int nl=0; nl<3; nl++)
       {
         for( int nc=0; nc<3; nc++)
         {
           matrix2[nl][nc] = pixelValue[nl*3+nc];//we "reshape" tuple2 to make a 3*3matrix
         }
       }
	vtkMath::Multiply3x3(matrix1,matrix2,matrixResult);//we multiply tensors

	for( int nl=0; nl<3; nl++)
       {
         for( int nc=0; nc<3; nc++)
         {
            Result[nl*3+nc]= matrixResult[nl][nc];//we "reshape" Result matrix to make a tuple
         }
       } 
       for (int i=0;i<9;i++){	ResultAndsize[i]=Result[i];}
       ResultAndsize[9]=pixelValue[9]+1;/*we store the number of tuples multiplied in each voxel in the last component of ResultAndsize vector*/

       foImage->SetPixel(pixelIndex, ResultAndsize);/*we set the vector in our image*/

    }
  
    return 0;
}

#endif
