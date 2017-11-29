// Copyright (c) General Electric Company, 2017.  All rights reserved.

#ifndef _itkShapeWaterShedImageFilter_txx
#define _itkShapeWaterShedImageFilter_txx

#include <iostream>
 
#include "itkShapeWaterShedImageFilter.h"
#include "Utilities/itkMinimumMaximumImageCalculatorWrapper.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkReconstructionByDilationImageFilter.h"
#include "Utilities/itkInvertIntensityImageFilterWrapper.h"
#include "itkBinaryThresholdImageFilter.h"

// Handle filters that were deprecated in ITKv4
#if ITK_VERSION_MAJOR < 4
  #include "itkSubtractConstantFromImageFilter.h"
  #include "itkAddConstantToImageFilter.h"
#else
  #include "itkSubtractImageFilter.h"
  #include "itkAddImageFilter.h"
#endif

#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
//#include "itkImageFileWriterWrapper.h"
#include "itkCastImageFilter.h"

#include "Utilities/itkRegionalExtremaImageFilter.h"
#include "Utilities/itkImposeMinimaImageFilter.h"
//#include "Utilities/WatershedSequenceImageFilter.h" //commented by Yousef to avoid linking errors
#include "Utilities/itkImageFileWriterWrapper.h"
#include "itkMorphologicalWatershedImageFilter.h" //added by Mirabela
/**
Original Code by: Dirk Paddfield, 
Addapted by Xiaofeng Liu, August 2010

*/ 

namespace itk
{
//#define DEBUG

/**
 *    Constructor
 */
template <class TInputImage, class TOutputImage>
ShapeWaterShedImageFilter<TInputImage, TOutputImage>
::ShapeWaterShedImageFilter()
{
	m_ImageSpacing = 0.5; 

	this->Superclass::SetNumberOfRequiredInputs( 1 );
	this->Superclass::SetNumberOfRequiredOutputs( 1 );
	this->Superclass::SetNthOutput( 0, OutputImageType::New() );
}


 
template <class TInputImage,class TOutputImage>
void ShapeWaterShedImageFilter<TInputImage,TOutputImage>
::GenerateData()
{
	// The input image should be a binary image.  
	// The function operates on distance maps generated from the binary image.

	typedef itk::Image<float, 2> FloatImageType; 
	typedef itk::Image< unsigned short, 2 > UShortImageType;
	typedef itk::Image< unsigned long, 2 > ULongImageType;
	typedef itk::Image< unsigned char, 2 > BinaryImageType;

	float radius = (2.0 / 3.0) / this->m_ImageSpacing; // This is specified in physical spacing.

	typename InputImageType::ConstPointer inputImage = this->GetInput();
	
	// step 1. distance map of the inverted mask
	// 1.1 Invert the mask image.
	float minimum;
	float maximum;
	
	// calculate min and max value
	typedef itk::MinimumMaximumImageCalculator<InputImageType> CalculatorType;
	typename CalculatorType::Pointer calculator = CalculatorType::New();
	calculator->SetImage( inputImage );
	calculator->Compute();
	minimum = calculator->GetMinimum();
	maximum = calculator->GetMaximum();

	// invert image
	typedef itk::InvertIntensityImageFilter<TInputImage> InvertIntensityType;
	typename InvertIntensityType::Pointer inverter = InvertIntensityType::New();
	inverter->SetInput( inputImage );
	inverter->SetMaximum( maximum + minimum);
	inverter->Update();  

	InputImagePointer invertedImage = inverter->GetOutput();


	// 1.2 compute distance map
	typedef itk::DanielssonDistanceMapImageFilter< InputImageType, FloatImageType >  DistanceFilterType;
	typename DistanceFilterType::Pointer distanceFilter = DistanceFilterType::New();
	distanceFilter->SetInput( invertedImage );
	distanceFilter->Update();
	invertedImage = NULL;

#if DEBUG
	{
		ImageFileWriterWrapper<FloatImageType>( distanceFilter->GetOutput(), "D:/TEMP/DistanceMap.mhd" );
	}
#endif

	// step 2. Saturate the distance image.
	// 2.1 subtract the distance map by radius
#if ITK_VERSION_MAJOR < 4
	typedef itk::SubtractConstantFromImageFilter<FloatImageType,float,FloatImageType> SubtractFilter;
#else
  typedef itk::SubtractImageFilter<FloatImageType,FloatImageType,FloatImageType> SubtractFilter;
#endif
	typename SubtractFilter::Pointer subtracter = SubtractFilter::New();
	subtracter->SetInput( distanceFilter->GetOutput() );
	subtracter->SetConstant( vnl_math_rnd(radius) );
	subtracter->Update();

	// 2.2 reconstruct. marker image: subtracter->GetOutput(); mask image: distanceImage
	typedef itk::ReconstructionImageFilter<FloatImageType, FloatImageType, std::greater<float> > ReconstructImageFilterType; 
	typename ReconstructImageFilterType::Pointer reconFilter = ReconstructImageFilterType::New(); 
	reconFilter->SetMarkerImage( subtracter->GetOutput());
	reconFilter->SetMaskImage(distanceFilter->GetOutput());
	reconFilter->SetFullyConnected(true); 
	reconFilter->Update(); 
	FloatImageType::Pointer saturatedIm = reconFilter->GetOutput(); 

#if DEBUG
	{
		ImageFileWriterWrapper<FloatImageType>( saturatedIm, "D:/TEMP/Saturated.mhd" );
	}
#endif


	// step 3. compute costImage: max(saturatedIm) - saturatedIm
	// 3.1 compute max
	typedef itk::MinimumMaximumImageCalculator<FloatImageType> MinMaxCalculator; 
	typename MinMaxCalculator::Pointer minMax = MinMaxCalculator::New();
	minMax->SetImage(saturatedIm);
	minMax->Compute();
	//float minVal = minMax->GetMinimum();
	float maxVal = minMax->GetMaximum(); 
	
	// 3.2 compute -saturatedIm
	typedef itk::InvertIntensityImageFilter<FloatImageType> InvertImageFilter;
	typename InvertImageFilter::Pointer invertImage = InvertImageFilter::New();
	invertImage->SetInput(saturatedIm);
	invertImage->SetMaximum(0.0);
	invertImage->Update();

	// 3.3 compute max-saturatedIm
#if ITK_VERSION_MAJOR < 4
	typedef itk::AddConstantToImageFilter<FloatImageType, float, FloatImageType> AddToImageFilterType;
#else
  typedef itk::AddImageFilter<FloatImageType,FloatImageType,FloatImageType> AddToImageFilterType;
#endif
	typename AddToImageFilterType::Pointer adder = AddToImageFilterType::New();
	adder->SetInput(invertImage->GetOutput());
	adder->SetConstant(maxVal);
	adder->Update(); 
	FloatImageType::Pointer costIm = adder->GetOutput(); 

#if DEBUG
	{
		ImageFileWriterWrapper<FloatImageType>( costIm, "D:/TEMP/cost_image.mhd" );
	}
#endif

	// step 4. watershed on costIm
	typedef itk::MorphologicalWatershedImageFilter<FloatImageType, ULongImageType> WatershedFilterType;
	typename WatershedFilterType::Pointer watershedFilter = WatershedFilterType::New();
	watershedFilter->SetInput(costIm); 
	watershedFilter->SetFullyConnected(true); 
	watershedFilter->Update();
	ULongImageType::Pointer watershedIm = watershedFilter->GetOutput(); 

#if DEBUG
	{
		ImageFileWriterWrapper<ULongImageType>( watershedIm, "D:/TEMP/shape_watersheded.mhd" );
	}
#endif


	// step 5. multiply watershed result with the mask image
	typedef itk::ImageRegionIterator<ULongImageType> ULongImageIterator; 
	typedef itk::ImageRegionConstIterator<InputImageType> InputImageConstIterator; 
	ULongImageIterator outImItr (watershedIm, watershedIm->GetBufferedRegion());
	InputImageConstIterator inputImItr (inputImage, inputImage->GetBufferedRegion() );
	for (outImItr.GoToBegin(), inputImItr.GoToBegin();
		!outImItr.IsAtEnd();
		++outImItr, ++inputImItr)
	{
		if (inputImItr.Get() == 0)
			outImItr.Set(0);
	}
#if DEBUG
	{
		ImageFileWriterWrapper<ULongImageType>( watershedIm, "D:/TEMP/shape_watershed_masked.mhd" );
	}
#endif

	
	// step 6. produce output
	typename  OutputImageType::Pointer outputImage       =  this->GetOutput(0);
	typedef itk::CastImageFilter<ULongImageType, OutputImageType> ImageCasterType;
	typename ImageCasterType::Pointer caster = ImageCasterType::New();
	caster->SetInput(watershedIm);
	caster->Update();
	outputImage->Graft(caster->GetOutput()); 
	
	return;
}

template <class TInputImage, class TOutputImage>
void 
ShapeWaterShedImageFilter<TInputImage,TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  
  os << " FIXME - needs to be implemented  " << std::endl;

}

	

} // end namespace itk

#endif
