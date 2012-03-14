#ifndef STATS_H
#define STATS_H

/*
 *	Author: Ravikiran J
 *	Email : ravikirn@cs.unc.edu
 */
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkMaskImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <DisplayFunc.h>

const unsigned int Dimension = 3;
typedef double PixelType;
typedef itk::Image< PixelType, Dimension > ImageType;
typedef itk::ImageFileReader<ImageType> ReaderType;
typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
typedef itk::ImageFileWriter< ImageType > WriterType;

class Stats{
	public:
		explicit Stats(const char* faImgFileName, const char* enImgFileName, const char* maskImgFileName);
        void computeMeanAndSD();
        double computePearsonCorrelation();
        int writeNormalizedAndDiffImages(const char* faNormImgFileName, const char* enNormImgFileName);
        int writeNormalizedCrossCorrImage();
        double computeNormalizedCrossCorrelation();
        void displayStats();
		virtual ~Stats(){}
		
	protected:
        ReaderType::Pointer faReader;
        ImageType::Pointer faImage;
        ImageType::PointType faOrigin;
        ImageType::SpacingType faSpacing;
        ImageType::SizeType faSize;

        ReaderType::Pointer enReader;
        ImageType::Pointer enImage;
        ImageType::PointType enOrigin;
        ImageType::SpacingType enSpacing;
        ImageType::SizeType enSize;

        ReaderType::Pointer maskReader;
        ImageType::Pointer maskImage;
        ImageType::PointType maskOrigin;
        ImageType::SpacingType maskSpacing;
        ImageType::SizeType maskSize;

        double meanFA;
        double sdFA;

        double meanEn;
        double sdEn;

        double correlation;
        double normCrossCorr;

        long int N;
} ;

Stats::Stats(const char *faImgFileName, const char* enImgFileName, const char* maskImgFileName){
    meanFA = meanEn = 0.0;
    sdFA = sdEn = 0.0;
    correlation = 0.0;
    normCrossCorr = 0.0;
    N = 0;

    faReader = ReaderType::New();
    faReader->SetFileName(faImgFileName);
    faReader->Update();
    enReader = ReaderType::New();
    enReader->SetFileName(enImgFileName);
    enReader->Update();
    maskReader = ReaderType::New();
    maskReader->SetFileName(maskImgFileName);
    maskReader->Update();

    
    faImage = faReader->GetOutput();
    enImage = enReader->GetOutput();
    maskImage = maskReader->GetOutput();

    faOrigin = faImage->GetOrigin();
    faSpacing = faImage->GetSpacing();
    faSize = faImage->GetLargestPossibleRegion().GetSize();

    enOrigin = enImage->GetOrigin();
    enSpacing = enImage->GetSpacing();
    enSize = enImage->GetLargestPossibleRegion().GetSize();

    maskOrigin = maskImage->GetOrigin();
    maskSpacing = maskImage->GetSpacing();
    maskSize = maskImage->GetLargestPossibleRegion().GetSize();

    computeMeanAndSD();
    displayStats();
}

void Stats::computeMeanAndSD(){
    //Initialize iterators
    ConstIteratorType fa_it( faImage, faImage->GetLargestPossibleRegion());
    ConstIteratorType en_it( enImage, enImage->GetLargestPossibleRegion());
    ConstIteratorType mask_it( maskImage, maskImage->GetLargestPossibleRegion());

    //Compute Mean
    fa_it.GoToBegin();
    en_it.GoToBegin();
    mask_it.GoToBegin();

    ImageType::IndexType pixelIndex;
    ImageType::PixelType enValue;
    ImageType::PixelType faValue;
    ImageType::PixelType maskValue;
    N = 0;
    double sumFA = 0.0, sumEn = 0.0;
    while(!fa_it.IsAtEnd() && !en_it.IsAtEnd() && !mask_it.IsAtEnd()){
        pixelIndex = fa_it.GetIndex();
        maskValue = mask_it.Get();
        if(maskValue > 0){
            faValue = fa_it.Get();
            enValue = en_it.Get();
            sumFA += faValue;
            sumEn += enValue;
            N++;
        }
        ++fa_it;
        ++en_it;
        ++mask_it;
    }
    meanFA = sumFA/N;
    meanEn = sumEn/N;

    //Compute SD
    fa_it.GoToBegin();
    en_it.GoToBegin();
    mask_it.GoToBegin();

    sumFA = 0.0, sumEn = 0.0;
    while(!fa_it.IsAtEnd() && !en_it.IsAtEnd() && !mask_it.IsAtEnd()){ 
        pixelIndex = fa_it.GetIndex(); 
        maskValue = mask_it.Get();
        if(maskValue > 0){
            faValue = fa_it.Get();
            enValue = en_it.Get();
            sumFA += pow(faValue-meanFA,2);
            sumEn += pow(enValue-meanEn,2);
        }
        ++fa_it;
        ++en_it;
        ++mask_it;
    }
    sumFA = sumFA / N;
    sumEn = sumEn / N;
    sdFA = sqrt(sumFA);
    sdEn = sqrt(sumEn);
}

double Stats::computePearsonCorrelation(){
    // Pearson Correlation is computed as per the definition found at http://en.wikipedia.org/wiki/Correlation_and_dependence
    //Initialize iterators
    ConstIteratorType fa_it( faImage, faImage->GetLargestPossibleRegion());
    ConstIteratorType en_it( enImage, enImage->GetLargestPossibleRegion());
    ConstIteratorType mask_it( maskImage, maskImage->GetLargestPossibleRegion());

    fa_it.GoToBegin();
    en_it.GoToBegin();
    mask_it.GoToBegin();

    ImageType::IndexType pixelIndex;
    ImageType::PixelType enValue;
    ImageType::PixelType faValue;
    ImageType::PixelType maskValue;
    ImageType::PixelType xy = 0.0, x = 0.0, y = 0.0, xpow2 = 0.0, ypow2 = 0.0;

    while(!fa_it.IsAtEnd() && !en_it.IsAtEnd() && !mask_it.IsAtEnd()){
        pixelIndex = fa_it.GetIndex();
        maskValue = mask_it.Get();
        if(maskValue > 0){
            faValue = fa_it.Get();
            enValue = en_it.Get();
            x += faValue;
            y += enValue;
            xy += faValue * enValue;
            xpow2 += pow(faValue,2);
            ypow2 += pow(enValue,2);
        }
        ++fa_it;
        ++en_it;
        ++mask_it;
    }
    correlation = ((N * xy) - x * y) / (sqrt(N * xpow2 - pow(x,2)) * sqrt(N * ypow2 - pow(y,2)));
    return correlation;
}


int Stats::writeNormalizedCrossCorrImage(){
    //Initialize iterators
    ConstIteratorType fa_it( faImage, faImage->GetLargestPossibleRegion());
    ConstIteratorType en_it( enImage, enImage->GetLargestPossibleRegion());
    ConstIteratorType mask_it( maskImage, maskImage->GetLargestPossibleRegion());

    // Set same as FA image
    ImageType::Pointer ncc = ImageType::New();
    ImageType::RegionType nccRegion;
    nccRegion.SetSize(faSize);
    ncc->SetSpacing(faSpacing);
    ncc->SetOrigin(faOrigin);
    ncc->SetRegions(nccRegion);
    ncc->Allocate();

    fa_it.GoToBegin();
    en_it.GoToBegin();
    mask_it.GoToBegin();

    ImageType::IndexType pixelIndex, newIndex;
    ImageType::PixelType enValue;
    ImageType::PixelType faValue;
    ImageType::PixelType maskValue;
    ImageType::PixelType corrValue;
    ImageType::PixelType xy = 0.0, x = 0.0, y = 0.0, xpow2 = 0.0, ypow2 = 0.0;
    int count = 0, neighborLevel = 3;
    output2("Computing Normalized Cross Correlation with neighborhood level of ", neighborLevel);

    while(!mask_it.IsAtEnd()){
        pixelIndex = mask_it.GetIndex();
        maskValue = mask_it.Get();
        if(maskValue > 0){
            xy = 0.0, x = 0.0, y = 0.0, xpow2 = 0.0, ypow2 = 0.0, count = 0;
            for(int i = -neighborLevel/2; i <= neighborLevel/2 ; i++){
                for(int j = -neighborLevel/2; j <= neighborLevel/2; j++){
                    for(int k = -neighborLevel/2; k <= neighborLevel/2; k++){
                        newIndex[0] = pixelIndex[0]+i;
                        newIndex[1] = pixelIndex[1]+j;
                        newIndex[2] = pixelIndex[2]+k;
                        if(newIndex[0] < 0 || newIndex[0] > (unsigned int)faSize[0] || newIndex[1] < 0 || newIndex[1] > (unsigned int)faSize[1] || newIndex[2] < 0 || newIndex[2] > (unsigned int)faSize[2]){
                            faValue = 0.0;
                            enValue = 0.0;
                        }else{
                            faValue = faImage->GetPixel(newIndex);
                            enValue = enImage->GetPixel(newIndex);
                        }
                        x += faValue;
                        y += enValue;
                        xy += faValue * enValue;
                        xpow2 += pow(faValue,2);
                        ypow2 += pow(enValue,2);
                        count++;
                    }
                }
            }
            if( (xpow2 == 0.0 && x == 0.0) || (ypow2 == 0.0 && y == 0.0) ){
                corrValue = 0.0;
            }else{
                corrValue = fabs(((count * xy) - x * y) / (sqrt(count * xpow2 - pow(x,2)) * sqrt(count * ypow2 - pow(y,2))));
            }
        }else{
            corrValue = 0.0;
        }
        ncc->SetPixel(pixelIndex, corrValue);
        ++mask_it;
    }
    output("Writing Normalized Cross Correlation Image");
    WriterType::Pointer nccWriter = WriterType::New();
    nccWriter->SetFileName("../normalizedCrossCorrImage.nrrd");
    nccWriter->SetInput(ncc);
    try{
        nccWriter->Update();
    }
    catch( itk::ExceptionObject & err ){
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }
    return 0;
}

int Stats::writeNormalizedAndDiffImages(const char* faNormImgFileName, const char* enNormImgFileName){
    //Initialize iterators
    ConstIteratorType fa_it( faImage, faImage->GetLargestPossibleRegion());
    ConstIteratorType en_it( enImage, enImage->GetLargestPossibleRegion());
    ConstIteratorType mask_it( maskImage, maskImage->GetLargestPossibleRegion());

    ImageType::Pointer faNormal = ImageType::New();
    ImageType::RegionType faRegion;
    faRegion.SetSize(faSize);
    faNormal->SetSpacing(faSpacing);
    faNormal->SetOrigin(faOrigin);
    faNormal->SetRegions(faRegion);
    faNormal->Allocate();

    //Set origin same as FA - hack
    ImageType::Pointer enNormal = ImageType::New();
    ImageType::RegionType enRegion;
    enRegion.SetSize(enSize);
    enNormal->SetSpacing(enSpacing);
    enNormal->SetOrigin(faOrigin);
    enNormal->SetRegions(enRegion);
    enNormal->Allocate();
    
    // Set same as FA image
    ImageType::Pointer diffNormal = ImageType::New();
    ImageType::RegionType diffRegion;
    diffRegion.SetSize(faSize);
    diffNormal->SetSpacing(faSpacing);
    diffNormal->SetOrigin(faOrigin);
    diffNormal->SetRegions(diffRegion);
    diffNormal->Allocate();

    fa_it.GoToBegin();
    en_it.GoToBegin();
    mask_it.GoToBegin();

    ImageType::IndexType pixelIndex;
    ImageType::PixelType enValue;
    ImageType::PixelType faValue;
    ImageType::PixelType maskValue;
    ImageType::PixelType diffValue;
    output("Computing Normalized FA, Entropy and Diff Images");
    while(!fa_it.IsAtEnd() && !en_it.IsAtEnd() && !mask_it.IsAtEnd()){
        pixelIndex = fa_it.GetIndex();
        maskValue = mask_it.Get();
        if(maskValue > 0){
            faValue = (fa_it.Get()-meanFA) / sdFA;
            enValue = (en_it.Get()-meanEn) / sdEn;
            diffValue = fabs(faValue-enValue);
        }else {
            faValue = 0.0;
            enValue = 0.0;
            diffValue = 0.0;
        }
        faNormal->SetPixel(pixelIndex, faValue);
        enNormal->SetPixel(pixelIndex, enValue);
        diffNormal->SetPixel(pixelIndex, diffValue);
        ++fa_it;
        ++en_it;
        ++mask_it;
    }
    output("Writing out Normalized FA, Entropy and Diff Images");
    WriterType::Pointer faWriter = WriterType::New();
    WriterType::Pointer enWriter = WriterType::New();
    WriterType::Pointer diffWriter = WriterType::New();

    faWriter->SetFileName(faNormImgFileName);
    enWriter->SetFileName(enNormImgFileName);
    diffWriter->SetFileName("../diff_fa_en.nrrd");

    faWriter->SetInput(faNormal);
    enWriter->SetInput(enNormal);
    diffWriter->SetInput(diffNormal);

    try{
        faWriter->Update();
        enWriter->Update();
        diffWriter->Update();
    }
    catch( itk::ExceptionObject & err ){
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }
    return 0;
}

double Stats::computeNormalizedCrossCorrelation(){
    //Initialize iterators
    ConstIteratorType fa_it( faImage, faImage->GetLargestPossibleRegion());
    ConstIteratorType en_it( enImage, enImage->GetLargestPossibleRegion());
    ConstIteratorType mask_it( maskImage, maskImage->GetLargestPossibleRegion());

    fa_it.GoToBegin();
    en_it.GoToBegin();
    mask_it.GoToBegin();

    ImageType::IndexType pixelIndex;
    ImageType::PixelType enValue;
    ImageType::PixelType faValue;
    ImageType::PixelType maskValue;
    normCrossCorr = 0.0;
    while(!fa_it.IsAtEnd() && !en_it.IsAtEnd() && !mask_it.IsAtEnd()){
        pixelIndex = fa_it.GetIndex();
        maskValue = mask_it.Get();
        if(maskValue > 0){
            faValue = (fa_it.Get()-meanFA) / sdFA;
            enValue = (en_it.Get()-meanEn) / sdEn;
            normCrossCorr += faValue * enValue;  
        }
        ++fa_it;
        ++en_it;
        ++mask_it;
    }
    normCrossCorr = normCrossCorr / (N-1);
    return normCrossCorr;
}

void Stats::displayStats(){
    output2("Number of Non-zero values in the mask = ", N);
    output2("Mean of FA Image = ", meanFA);
    output2("Mean of Entropy Image = ", meanEn);
    output2("Standard Deviation of FA Image = ", sdFA);
    output2("Standard Deviation of Entropy Image = ", sdEn);
}
#endif
