// Copyright (c) General Electric Company, 2017.  All rights reserved.

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h" 

#include "itkWaveletNucleiSegmentationFilter.h"

int main(int argc, char* argv[])
{
	typedef itk::Image< unsigned short, 2 > ShortImageType;
	typedef itk::Image< unsigned short, 2 > LongImageType;

	typedef itk::Image< float, 2 > DoubleImageType;
	typedef itk::ImageFileReader< ShortImageType > ReaderType;
	typedef itk::ImageFileWriter< ShortImageType > WriterType;
	
	if (argc < 3)
	{
		std::cerr << "Usage: " << argv[0]<< " InputImageFileName OutputImageFileName [lowerDetailLevel upperDetailLevel numberOfLevels smoothingSigma]"
<< std::endl;
		return 0;
	}
	// cast

	std::cout<<argv[0];

	typedef itk::CastImageFilter <ShortImageType, DoubleImageType> ShortToDoubleType;
	ShortToDoubleType::Pointer toDoubleImage = ShortToDoubleType::New();

	typedef itk::WaveletNucleiSegmentationFilter<DoubleImageType, ShortImageType>
WaveletNucleiSegmentationFilterType;
	//typedef itk::WaveletBlobDetectionSegmentationFilter< DoubleImageType, DoubleImageType > WaveletBlobDetectionFilterType;
	WaveletNucleiSegmentationFilterType::Pointer  nucleiSegFilter = 
WaveletNucleiSegmentationFilterType::New();

  ReaderType::Pointer reader = ReaderType::New();

	std::cout<<argv[1]<<std::endl;

	reader->SetFileName(argv[1]);

	reader->Update();
	
	toDoubleImage->SetInput(reader->GetOutput());
	toDoubleImage->Update();


	//segment
	int lowerDetailLevel = 3;
	int upperDetailLevel = 3;
	int numberOfLevels   = 3;
	float sigmaMultiplier = 2;
	if( argc > 3 )
	{
      lowerDetailLevel = atoi(argv[3]);
	}
	if( argc > 4 )
	{
	  upperDetailLevel = atoi(argv[4]);
	}
	if( argc > 5 )
	{
	  numberOfLevels = atoi(argv[5]);
	}
	if( argc > 6 )
	{
	  sigmaMultiplier = atof(argv[6]);
	}

	nucleiSegFilter->SetDetailImagesToUse_X(lowerDetailLevel); // was 4
	nucleiSegFilter->SetDetailImagesToUse_Y(upperDetailLevel); // was 4
	nucleiSegFilter->SetNumberOfLevels(numberOfLevels); // was 5
	nucleiSegFilter->SetSigmaMultiplier(sigmaMultiplier);
	nucleiSegFilter->SetInput(toDoubleImage->GetOutput());
	nucleiSegFilter->Update();

	// cast back
	//toLongImage->SetInput(nucleiSegFilter->GetOutput());
	//toLongImage->Update();
	//write
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput(nucleiSegFilter->GetOutput());
	writer->SetFileName( argv[2] );
	try{writer->Update();}catch (itk::ExceptionObject &e){std::cerr << "Can't write image: "<<  e << std::endl;}    
	std::cout << "File Output: " << argv[2] << std::endl;
	std::cout<<"Wavelet Nuclei Segmentation... DONE"<<std::endl;


	return 0;
}
