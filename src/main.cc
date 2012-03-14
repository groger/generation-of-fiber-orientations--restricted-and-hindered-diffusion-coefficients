#include "mainCLP.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVariableLengthVector.h"
#include "itkRGBPixel.h"
#include <VolumeCounter.h>
#include <FiberSegmentExtractor.h>
#include <iostream>
#include <fstream>
#include <vtkPolyDataReader.h>
#include <Counter.h>
//#include <Entropy.h>
#include <Stats.h>
#include <FiberOrientation.h>
#include<vtkPolynomialSolversUnivariate.h>
#include <string>
int computeVolumeCounter(const char *, const char *);
void computeFAImageParams(const char *);
void readFOImage(const char *);

int main(int argc, char **argv){
    PARSE_ARGS;
    //Input files 
    const char *dwiImgFileName=dwiImgFilename.c_str();
    const char *maskImgFileName =maskImgFilename.c_str();
    const char *EigenFileName=EigenFilename.c_str();
    const char *EigenHinderedImgFileName=EigenHinderedImgFilename.c_str();
    const char *foImgFileName=foImgFilename.c_str();
    const char *faImgFileName =faImgFilename.c_str();
    const char *fiberFileName=fiberFilename.c_str();
   
    FiberOrientation fo(fiberFileName, dwiImgFileName, faImgFileName, maskImgFileName);

    fo.writeFiberOrientationsPerVoxelTextFile(foImgFileName,EigenFileName,EigenHinderedImgFileName);
	
    output("Done");
    return 0;   
} 
