// Copyright (c) General Electric Company, 2017.  All rights reserved.

#ifndef _itkWaveletNucleiSegmentationFilter_txx
#define _itkWaveletNucleiSegmentationFilter_txx

#include <iostream>

#include "itkWaveletNucleiSegmentationFilter.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkShapeWaterShedImageFilter.h"
#include "itkGrayscaleFillholeImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkFastChamferDistanceImageFilter.h"
#include "itkRegionalMaximaImageFilter.h"
#include "itkRegionalMinimaImageFilter.h"
#include "itkValuedRegionalMinimaImageFilter.h"
#include "itkMorphologicalWatershedImageFilter.h"
#include "itkMorphologicalWatershedFromMarkersImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkReconstructionImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkMultiplyImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "vnl/vnl_sample.h"
#include "itkSimpleFilterWatcher.h"

// From MIA
#include "Utilities/itkOpeningClosingFilterWrapper.h"

// From ITK
#include "Utilities/itkATrousWaveletBlobSegmentationImageFilter.h"

#define DEBUG_WAVELET 0

#if DEBUG_WAVELET
#include "itkImageFileWriterWrapper.h"
#endif

#include <vnl/vnl_double_3x3.h>
#include <vector>
#include <math.h>

#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_real_eigensystem.h>
/**
Original Code by: Dirk Paddfield, 
Adapted by Alberto Santamaria/Xiaofeng Liu, August 2010
*/

namespace itk
{

/**
 *    Constructor
 */
template <class TInputImage, class TOutputImage, class TLabelImageType>
WaveletNucleiSegmentationFilter<TInputImage, TOutputImage, TLabelImageType>::WaveletNucleiSegmentationFilter()
{

  m_NumberOfLevels  = 5;
  m_SigmaMultiplier = 3;
  m_DetailImagesToUse_X = 4;
  m_DetailImagesToUse_Y = 4;

  this->Superclass::SetNumberOfRequiredInputs( 1 );
  this->Superclass::SetNumberOfRequiredOutputs( 1 );
  this->Superclass::SetNthOutput( 0, OutputImageType::New() );

}


/**
 *    GenerateData 
 */
template <class TInputImage, class TOutputImage, class TLabelImageType>
void WaveletNucleiSegmentationFilter<TInputImage,TOutputImage,TLabelImageType>::GenerateData()
{

  this->GetOutput(0)->SetBufferedRegion( this->GetOutput(0)->GetRequestedRegion() );
  this->GetOutput(0)->Allocate();

  typename  InputImageType::ConstPointer inputImage =  this->GetInput(0);
  typename  OutputImageType::Pointer outImage       =  this->GetOutput(0);

  // step 1, do  morphological operations on the input image
  // 1.1 erode
  typedef itk::BinaryBallStructuringElement<double, 2> StructuringElementType;
  StructuringElementType   se;
  se.SetRadius ( 3);
  se.CreateStructuringElement();
  typedef itk::GrayscaleErodeImageFilter<InputImageType, InputImageType, StructuringElementType> ErodeImageFilter;
  typename ErodeImageFilter::Pointer erodeFilter = ErodeImageFilter::New();
  erodeFilter->SetInput(inputImage);
  erodeFilter->SetKernel(se); 
  erodeFilter->Update();
  
  // 1.2 reconstruction
  
  typedef itk::ReconstructionImageFilter<InputImageType, InputImageType, std::greater<typename InputImageType::PixelType> > ReconstructImageFilterType; 
  typename ReconstructImageFilterType::Pointer reconFilter = ReconstructImageFilterType::New(); 
  reconFilter->SetMarkerImage(erodeFilter->GetOutput());
  reconFilter->SetMaskImage(inputImage);
  reconFilter->SetFullyConnected(true); 
  reconFilter->Update(); 

  //ImageFileWriterWrapper< InputImageType >( reconFilter->GetOutput(), "ReconstructionImageFilter.mhd" );

  // step 2 detect blobs using wavelet
  typedef itk::ATrousWaveletBlobSegmentationImageFilter<InputImageType, InputImageType> WaveletBlobDetector;
  typename WaveletBlobDetector::Pointer blobDetector = WaveletBlobDetector::New();
  // This is a temporary fix: subtract 1 from the lower and upper detail levels to conform with the
  // bug that was in the MIA implementation of the blob detector.
  // Eventually, we would like to not subtract 1, but this is used so that the experiments already running
  // will not break.
  int lowerLevel = this->GetDetailImagesToUse_X() - 1;
  int upperLevel = this->GetDetailImagesToUse_Y() - 1;
  blobDetector->SetLowerLevel( lowerLevel );
  blobDetector->SetUpperLevel( upperLevel );
  blobDetector->SetSigmaMultiplier(this->GetSigmaMultiplier()); 
  blobDetector->SetInput(reconFilter->GetOutput());
  blobDetector->Update();

  typename InputImageType::Pointer blobIm = blobDetector->GetOutput(); 

  // 2.1 mask: by thresholding the blobimage
  // no need because the output of blobdetection is binary image
  typedef itk::BinaryThresholdImageFilter<InputImageType, BinaryImageType> ThresholdInputImageFilter;
  typename ThresholdInputImageFilter::Pointer inputThresholder = ThresholdInputImageFilter::New(); 
  inputThresholder->SetLowerThreshold( 0.0001 );
  inputThresholder->SetInsideValue( 1 );
  inputThresholder->SetOutsideValue( 0 );
  inputThresholder->SetInput(blobIm); 
  inputThresholder->Update(); 
  BinaryImageType::Pointer maskIm = inputThresholder->GetOutput(); 

#if DEBUG_WAVELET
	{
		ImageFileWriterWrapper<InputImageType>( blobIm, "image_After_Wavelet.mhd" );
		ImageFileWriterWrapper<BinaryImageType>(maskIm, "mask_After_Wavelet.tif" );
	}
#endif
	// step 3 apply watershed
	typename ULongImageType::Pointer watershedIm = WatershedSegmentation(reconFilter->GetOutput(), maskIm); 

	// step 4. multiply the watershed image with the mask, and threshold it
	// 4.1 watershedIm(blobIm == 0) = 0; 
	typedef itk::ImageRegionIterator<ULongImageType> ULongImageIterator; 
	typedef itk::ImageRegionIterator<InputImageType> InputImageIterator; 
	ULongImageIterator watershedImItr (watershedIm, watershedIm->GetBufferedRegion() );
	InputImageIterator blobItr (blobIm, blobIm->GetBufferedRegion() );
	for (watershedImItr.GoToBegin(), blobItr.GoToBegin();
		!watershedImItr.IsAtEnd();
		++watershedImItr, ++blobItr)
	{
		if (blobItr.Get() == 0)
			watershedImItr.Set(0); 
	}
	
#if DEBUG_WAVELET
	{
		ImageFileWriterWrapper<ULongImageType>( watershedIm, "watershed_masked.mhd" );
	}
#endif
	
	// 4.2 threshold watershedIm (watershedIm>0)
	typedef itk::BinaryThresholdImageFilter<ULongImageType,BinaryImageType> ThresholderType;
	ThresholderType::Pointer thresholder = ThresholderType::New();
	thresholder->SetLowerThreshold( 1 );
	thresholder->SetInsideValue( 1 );
	thresholder->SetOutsideValue( 0 );
	thresholder->SetInput(watershedIm); 
	thresholder->Update(); 

#if DEBUG_WAVELET
	{
		ImageFileWriterWrapper<BinaryImageType>( thresholder->GetOutput(), "watershed_thresholded.tif" );
	}
#endif
  
  // step 5. post-processing
  typename OutputImageType::Pointer processedSegmentedImage = OutputImageType::New();
  processedSegmentedImage = this->Segmentation_PostProcessing( reconFilter->GetOutput(), thresholder->GetOutput() ); 
  
  // output
  outImage->Graft(processedSegmentedImage); 

	typedef itk::ImageRegionIterator<InputImageType> InputImageIterator; 
	typedef itk::ImageRegionIterator<BinaryImageType> BinaryImageIterator; 

}

template <class TInputImage, class TOutputImage, class TLabelImageType>
void 
WaveletNucleiSegmentationFilter<TInputImage,TOutputImage,TLabelImageType>::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  
  os << " FIXME - needs to be implemented  " << std::endl;

}

// apply watershed on the distance map computed on the image after wavelet, 
template <class TInputImage, class TOutputImage, class TLabelImageType>
itk::Image<unsigned long,2>::Pointer 
WaveletNucleiSegmentationFilter<TInputImage,TOutputImage,TLabelImageType>::WatershedSegmentation(InputImagePointer image,
                                                                                 BinaryImagePointer maskIm)
{
	// step 1. compute distance map after inverting the image, and apply an "imimposemin" function. 
	// distImageInverted = imimposemin (-double(inImage), outImageSmall>0);
	
	// 1.1 invert the input image
	//InputImageType::ConstPointer image = this->GetInput();
	typedef itk::InvertIntensityImageFilter<InputImageType, InputImageType> InvertInputImageFilterType;
	typename InvertInputImageFilterType::Pointer invertInputFilter = InvertInputImageFilterType::New();
	invertInputFilter->SetInput(image);
	invertInputFilter->SetMaximum(0.0);
	invertInputFilter->Update();

#if DEBUG_WAVELET
	{
	ImageFileWriterWrapper<InputImageType>( invertInputFilter->GetOutput(), "ws_invertInput.mhd" );
	}
#endif

	// 1.2 Implement the matlab function: imimposemin

	// 1.2.1 no matter what, cast the image into float type
	typedef itk::CastImageFilter<InputImageType, FloatImageType> CastFloatImageFilter; 
	typename CastFloatImageFilter::Pointer floatCaster = CastFloatImageFilter::New();
	floatCaster->SetInput(invertInputFilter->GetOutput());
	floatCaster->Update();

	// 1.2.2 convert the image, so that: for pixels in the mask, set to a big negative number;
	// for pixels outside of the image, set to the image value. 
	float minInfVal = -1.0e5;
	float maxInfVal = 1.0e5; // later: need to change this to relate to the min and max value of the image
	typedef itk::ImageRegionIterator<FloatImageType> FloatImageIterator; 
	typedef itk::ImageRegionIterator<InputImageType> InputImageIterator; 
	typedef itk::ImageRegionIterator<BinaryImageType> BinaryImageIterator; 
	FloatImageIterator invertImItr (floatCaster->GetOutput(), floatCaster->GetOutput()->GetBufferedRegion()); 
	BinaryImageIterator maskImItr (maskIm, maskIm->GetBufferedRegion()); 
	
	FloatImageType::Pointer reconMask = FloatImageType::New();
	reconMask->SetRegions(maskIm->GetBufferedRegion());
	reconMask->Allocate();
	reconMask->FillBuffer(minInfVal);
	FloatImageType::Pointer reconMarker = FloatImageType::New();
	reconMarker->SetRegions(maskIm->GetBufferedRegion());
	reconMarker->Allocate();
	reconMarker->FillBuffer(minInfVal);

	FloatImageIterator reconMaskItr(reconMask, reconMask->GetBufferedRegion());
	FloatImageIterator reconMarkerItr(reconMarker, reconMarker->GetBufferedRegion());

	for (reconMaskItr.GoToBegin(), reconMarkerItr.GoToBegin(), maskImItr.GoToBegin(), invertImItr.GoToBegin();
		!reconMaskItr.IsAtEnd();
		++reconMaskItr, ++reconMarkerItr, ++maskImItr, ++invertImItr)
	{
		if (maskImItr.Get() == 0)
		{
			reconMaskItr.Set(invertImItr.Get());
			reconMarkerItr.Set(maxInfVal); 
		}
	}
#if DEBUG_WAVELET
	{
	ImageFileWriterWrapper<FloatImageType>( reconMarker, "ws_reconMarker.mhd" );
	ImageFileWriterWrapper<FloatImageType>( reconMask, "ws_reconMask.mhd" );
	}
#endif
	
	// 1.2.2 compute the complement of recon mask and marker images
	typedef itk::InvertIntensityImageFilter<FloatImageType, FloatImageType> InvertFloatImageFilter;
	InvertFloatImageFilter::Pointer floatImageInverter = InvertFloatImageFilter::New();
	floatImageInverter->SetInput(reconMask);
	floatImageInverter->SetMaximum(0.0);
	floatImageInverter->Update();
	FloatImageType::Pointer compReconMask = floatImageInverter->GetOutput();

	floatImageInverter = InvertFloatImageFilter::New();
	floatImageInverter->SetInput(reconMarker);
	floatImageInverter->SetMaximum(0.0);
	floatImageInverter->Update();
	FloatImageType::Pointer compReconMarker = floatImageInverter->GetOutput();

#if DEBUG_WAVELET
	{
	ImageFileWriterWrapper<FloatImageType>( compReconMarker, "ws_CompreconMarker.mhd" );
	ImageFileWriterWrapper<FloatImageType>( compReconMask, "ws_CompreconMask.mhd" );
	}
#endif

	// 1.2.3 apply image reconstruction
	typedef itk::ReconstructionImageFilter<FloatImageType, FloatImageType, std::greater<float> > ReconstructionFilter; 
	ReconstructionFilter::Pointer reconFilter = ReconstructionFilter::New();
	reconFilter->SetMarkerImage ( compReconMarker);
	reconFilter->SetMaskImage ( compReconMask);
	reconFilter->SetFullyConnected (true);
	reconFilter->Update();

#if DEBUG_WAVELET
	{
	ImageFileWriterWrapper<FloatImageType>( reconFilter->GetOutput(), "ws_reconImage.mhd" );
	}
#endif
	// 1.2.4 invert image
	floatImageInverter = InvertFloatImageFilter::New();
	floatImageInverter->SetInput(reconFilter->GetOutput());
	floatImageInverter->SetMaximum(0.0);
	floatImageInverter->Update();


	// 2. apply watershed on the image after 
	// first way: use filter	MorphologicalWatershedImageFilter
	typedef itk::MorphologicalWatershedImageFilter<FloatImageType, ULongImageType> WatershedFilterType;
	WatershedFilterType::Pointer watershedFilter = WatershedFilterType::New();
	watershedFilter->SetInput(floatImageInverter->GetOutput()); 
	watershedFilter->SetFullyConnected(true); 
	watershedFilter->Update();

#if DEBUG_WAVELET
	{
	//ImageFileWriterWrapper<itk::Image<unsigned long,2>>( watershedFilter->GetOutput(), "ws_result.mhd" );
	}
#endif

	//3. generate output; 
	return watershedFilter->GetOutput(); 
}


// This function is supposed to replace the counter part in the matlab code
// However, maybe it need not do more than hole filling.
// The watershed method may be should not be mixed with the wavelet part
template <class TInputImage, class TOutputImage, class TLabelImageType>
typename TOutputImage::Pointer
WaveletNucleiSegmentationFilter<TInputImage,TOutputImage,TLabelImageType>::Segmentation_PostProcessing(typename InputImageType::Pointer image,
                                                                                       typename BinaryImageType::Pointer mask)
{
	// step 1. fill holes
	typedef itk::GrayscaleFillholeImageFilter<BinaryImageType, BinaryImageType> HoleFillingFilter;
	HoleFillingFilter::Pointer holeFiller = HoleFillingFilter::New();
	holeFiller->SetInput(mask);
	holeFiller->Update();

#if DEBUG_WAVELET
	{
		ImageFileWriterWrapper<BinaryImageType>( mask, "mask.tif" );
		ImageFileWriterWrapper( holeFiller->GetOutput(), "holeFilling.tif" );
	}
#endif

	//std::cout<<"step 2. run shape based watershed..."<<std::endl;
	// step 2. run shape based watershed 
	float imageSpacing = 0.5; 
	float meanCellRadiusInPixels = 5.0 / imageSpacing;

	typedef itk::ShapeWaterShedImageFilter<BinaryImageType, ULongImageType> ShapeWatershedFilter;
	ShapeWatershedFilter::Pointer shapeWatershed =  ShapeWatershedFilter::New();
	
	shapeWatershed->SetInput(holeFiller->GetOutput());
	shapeWatershed->SetImageSpacing(imageSpacing );
	shapeWatershed->Update();
	ULongImageType::Pointer shapeWSIm = shapeWatershed->GetOutput(); 


	//std::cout<<"step 3. remove small objects..."<<std::endl;
	// step 3. remove small objects
	// 3.1 binarize image
	typedef itk::BinaryThresholdImageFilter<ULongImageType, BinaryImageType> ThresholdULongImageFilter;
	ThresholdULongImageFilter::Pointer ulongThresholder = ThresholdULongImageFilter::New(); 
	ulongThresholder->SetLowerThreshold( 1 );
	ulongThresholder->SetInsideValue( 1 );
	ulongThresholder->SetOutsideValue( 0 );
	ulongThresholder->SetInput(shapeWSIm); 
	ulongThresholder->Update(); 

	//ImageFileWriterWrapper< BinaryImageType >(ulongThresholder->GetOutput(), "ulongThresholder.tif" );

	int Radius = 1;
	BinaryImageType::Pointer TempImage = itk::OpeningClosingFilterWrapper<BinaryImageType>( ulongThresholder->GetOutput(),  Radius);

	//ImageFileWriterWrapper< BinaryImageType >(TempImage, "ulongThresholder_2.tif" );

	//std::cout<<" 3.2 connected component analysis..."<<std::endl;
	// 3.2 connected component analysis. 
	typedef itk::ConnectedComponentImageFilter<BinaryImageType, TLabelImageType> ConnectedComponentType;
	typename ConnectedComponentType::Pointer ccl = ConnectedComponentType::New(); 

	//ccl->SetInput(ulongThresholder->GetOutput());
	ccl->SetInput( TempImage );
  ccl->Update(); 

  // Suppress small objects
  typedef itk::RelabelComponentImageFilter<TLabelImageType, OutputImageType> RelabelType;
  typename RelabelType::Pointer relabel = RelabelType::New();
  relabel->SetInput( ccl->GetOutput());
	float areaThreshold = 0.2 * (vnl_math::pi *  vcl_pow(meanCellRadiusInPixels,(float)2.0));
  relabel->SetMinimumObjectSize( areaThreshold );
  try
    {
    relabel->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Relabel: exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }
  return relabel->GetOutput();

#if 0
	UShortImageType::Pointer labeledIm = ccl->GetOutput(); 
	UShortImageType::ConstPointer constlabeledIm = ccl->GetOutput(); 
	// 3.3 compute component area
	vnl_vector<unsigned int>  area = ConnectedComponentArea<UShortImageType> (constlabeledIm); 

	// 3.4 threshold based on area. 
	float areaThreshold = 0.2 * (vnl_math::pi *  vcl_pow(meanCellRadiusInPixels,(float)2.0));
	typedef itk::ImageRegionIterator<UShortImageType> UShortImageIterator; 
	UShortImageIterator labeledItr (labeledIm, labeledIm->GetBufferedRegion());
	unsigned short t; 
	for (labeledItr.GoToBegin(); !labeledItr.IsAtEnd(); ++labeledItr)
	{	
		if ( (t = labeledItr.Get()) == 0) continue;
		if (area[t-1] <= areaThreshold )
			labeledItr.Set(0);
	}
	// 3.4 threshold labeledIm to get the mask
	typedef itk::BinaryThresholdImageFilter<UShortImageType, OutputImageType> ThresholdILabelImageFilter;
	ThresholdILabelImageFilter::Pointer thresholder = ThresholdILabelImageFilter::New();
	thresholder->SetLowerThreshold (1);
	thresholder->SetInsideValue(1);
	thresholder->SetOutsideValue(0);
	thresholder->SetInput (labeledIm);
	thresholder->Update();

	return thresholder->GetOutput(); 
#endif
}

} // end namespace itk

#endif
