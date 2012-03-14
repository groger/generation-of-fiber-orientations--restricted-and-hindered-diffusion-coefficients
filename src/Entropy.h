#ifndef ENTROPY_H
#define ENTROPY_H

/*
 *	Author: Ravikiran J
 *	Email : ravikirn@cs.unc.edu
 */
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include <VolumeCounter.h>
#include <FiberSegmentExtractor.h>
#include <iostream>
#include <fstream>
#include <vtkPolyDataReader.h>
#include <Counter.h>
#include <DisplayFunc.h>

class Entropy{
	public:
		explicit Entropy(const char* volumeCounterFileName);
		virtual ~Entropy(){}
        void writeToTextFile(const char* fname);
        int writeToGrayscaleImage(const char* fname);
		
	protected:
        std::vector< std::vector < std:: vector < double > > > entropy;
        int dims[3];
        int spacing[3];
        int origin[3];
} ;

Entropy::Entropy(const char* volumeCounterFileName){
    output("Reading the Volume Counter");
    VolumeCounter<Counter_WeightedVertices> vc(volumeCounterFileName);
    // Initialize the dimensions, origin and axes
    for(int i = 0; i < VOLUME_DIMENSION; i++){
        dims[i] = vc.GetNumVoxelsAlongAxis(i);
        origin[i] = vc.GetOriginAlongAxis(i);
        spacing[i] = vc.GetSpacingAlongAxis(i);
    }

    // Setup sizes of the 3D entropy vector 
    entropy.resize(dims[0]);
    for (int i = 0; i < dims[0]; i++) {
       entropy[i].resize(dims[1]);
       for (int j = 0; j < dims[1]; j++){
          entropy[i][j].resize(dims[2]);
       }
    }

    double e = 0.0;
    double total = 0.0;
    double p = 0.0;

    output("Computing Entropy");
    for (int i = 0; i < dims[0]; i++){
        for (int j = 0; j < dims[1]; j++){
            for (int k = 0; k < dims[2]; k++){
                Counter &c = vc.GetCounterForVoxel(i, j, k) ;
                if (c.IsEmpty()){
                    e = 0.0;
                    total = 0.0;
                }else{
                    std::vector<double> res = c.GetBins() ;
                    total = 0.0 ;
                    for (unsigned int l = 0; l < res.size(); l++){
                        total += (double) res[l] ;
                    }
                    e = 0.0;
                    for (unsigned int l = 0; l < res.size(); l++){
                        p = (double)res[l] / total;
                        // if p < 0, you need to add 0 to entropy, hence skip it
                        if(p > 0){
                            e += p * log(p); 
                        }
                    }
                    e *= -1;
                }
                entropy[i][j][k] = e;
            }
        }
    }
}
void Entropy::writeToTextFile(const char* fname){
    output("Writing out the entropy result to text file");
 	std::ofstream entropyOutput(fname);
    for (int i = 0; i < dims[0]; i++){
        for (int j = 0; j < dims[1]; j++){
            for (int k = 0; k < dims[2]; k++){
                entropyOutput << i << "," << j << "," << k << ","<< entropy[i][j][k] << std::endl ;
            }
        }
    }
}

int Entropy::writeToGrayscaleImage(const char* fname){
    typedef double PixelType;
    const unsigned int Dimension = 3;
    typedef itk::Image< PixelType, Dimension > gImageType;
    typedef itk::ImageFileWriter< gImageType > gWriterType;

    gImageType::PointType gOrigin;
    gImageType::SpacingType gSpacing;
    gImageType::SizeType gSize;

	for (int i = 0; i < VOLUME_DIMENSION; i++){
		gOrigin[i]  = origin[i];
        gSpacing[i] = spacing[i];
        gSize[i] = dims[i];
    }

    gImageType::Pointer gImage = gImageType::New();
    
    gImageType::RegionType gRegion;
    gRegion.SetSize(gSize);
    gImage->SetSpacing(gSpacing);
    gImage->SetOrigin(gOrigin);
    gImage->SetRegions(gRegion);
    gImage->Allocate();
     
    gImageType::IndexType gPixelIndex;
    gImageType::PixelType gPixelValue;

    // Write the entropy to the image
    output("Writing out the entropy image");
    for (int i = 0; i < dims[0]; i++){
        for (int j = 0; j < dims[1]; j++){
            for (int k = 0; k < dims[2]; k++){
                gPixelIndex[0] = i; gPixelIndex[1] = j; gPixelIndex[2] = k;
                gPixelValue = entropy[i][j][k];
                gImage->SetPixel(gPixelIndex, gPixelValue);
            }
        }
    }

    // Update the writer 
    gWriterType::Pointer gWriter = gWriterType::New();
    gWriter->SetFileName(fname);
    gWriter->SetInput(gImage);
    try{
        gWriter->Update();
    }
    catch( itk::ExceptionObject & err ){
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }
    return 0;
}

#endif
