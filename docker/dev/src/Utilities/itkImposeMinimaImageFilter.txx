// Copyright (c) General Electric Company, 2017.  All rights reserved.
/*=========================================================================

  Author: Dirk Padfield
  date: 03/03/2007

  =========================================================================*/

#ifndef _itkImposeMinimaImageFilter_txx
#define _itkImposeMinimaImageFilter_txx

#include "itkImposeMinimaImageFilter.h"
#include <time.h>

namespace itk
{

template <class TInputImage, class TMarkerImage, class TOutputImage>
ImposeMinimaImageFilter< TInputImage, TMarkerImage, TOutputImage>
::ImposeMinimaImageFilter()
{
  m_FullyConnected = false;
}

template <class TInputImage, class TMarkerImage, class TOutputImage>
void 
ImposeMinimaImageFilter< TInputImage, TMarkerImage, TOutputImage>
::SetMarkerImage(const MarkerImageType* markerImage)
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(1, const_cast<MarkerImageType *>( markerImage ));
}


template <class TInputImage, class TMarkerImage, class TOutputImage>
const typename ImposeMinimaImageFilter< TInputImage, TMarkerImage, TOutputImage>::MarkerImageType *
ImposeMinimaImageFilter< TInputImage, TMarkerImage, TOutputImage>
::GetMarkerImage()
{
  return static_cast<MarkerImageType*>(const_cast<DataObject *>(this->ProcessObject::GetInput(1)));
  //return this->GetInput(1);
}


template <class TInputImage, class TMarkerImage, class TOutputImage>
void 
ImposeMinimaImageFilter< TInputImage, TMarkerImage, TOutputImage>
::SetMaskImage(const MaskImageType* maskImage)
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(0, const_cast<MaskImageType *>( maskImage ));
}


template <class TInputImage, class TMarkerImage, class TOutputImage>
const typename ImposeMinimaImageFilter< TInputImage, TMarkerImage, TOutputImage>::MaskImageType *
ImposeMinimaImageFilter< TInputImage, TMarkerImage, TOutputImage>
::GetMaskImage()
{
  return static_cast<MaskImageType*>(const_cast<DataObject *>(this->ProcessObject::GetInput(0)));
  //return this->GetInput(0);
}

template <class TInputImage, class TMarkerImage, class TOutputImage>
void
ImposeMinimaImageFilter< TInputImage, TMarkerImage, TOutputImage>
::GenerateData()
{
//   // Start a timer.
//   time_t start,end;
//   double diff;
//   time (&start);

  // Allocate the output
  this->AllocateOutputs();


  // Invert the marker image and set the values to the extremes of the
  // pixel type.
  // Any positive values in the maker image will be set to the
  // smallest value of the PixelType.
  // Zero values will be set to the largest value of the PixelType.
  /*
  typedef itk::InvertIntensityImageFilter< TInputImage, TInputImage > InvertType;
  InvertType::Pointer inverter = InvertType::New();
  inverter->SetInput( this->GetMarkerImage() );
  inverter->Print(std::cout);

  typedef itk::ShiftScaleImageFilter< TInputImage, TInputImage > ShiftScaleType;
  ShiftScaleType::Pointer inverter = ShiftScaleType::New();
  inverter->SetShift( itk::NumericTraits<TInputImage::PixelType>::max() );
  inverter->SetScale( -1 );
  */


  typename TInputImage::Pointer invertedImage = TInputImage::New();
  invertedImage->SetRegions( this->GetMarkerImage()->GetLargestPossibleRegion() );
  invertedImage->CopyInformation( this->GetMarkerImage() );
  invertedImage->Allocate();
  itk::ImageRegionConstIterator< MarkerImageType > markerIt( this->GetMarkerImage(), this->GetMarkerImage()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< TInputImage > inverterIt( invertedImage, invertedImage->GetLargestPossibleRegion() );
  for( markerIt.Begin(), inverterIt.Begin(); !markerIt.IsAtEnd(); ++markerIt, ++inverterIt )
    {
    if( markerIt.Get() )      
      inverterIt.Set( itk::NumericTraits<typename TInputImage::PixelType>::NonpositiveMin() );
    else
      inverterIt.Set( itk::NumericTraits<typename TInputImage::PixelType>::max() );
    }



  // Add 1 to the mask image.  This ensures that none of the pixels in the image are the smallest possible value.
  typename TInputImage::Pointer constantImage = TInputImage::New();
  constantImage->SetRegions( this->GetMaskImage()->GetLargestPossibleRegion() );
  constantImage->CopyInformation( this->GetMaskImage() );
  constantImage->Allocate();
  constantImage->FillBuffer( 1 );

  typedef itk::AddImageFilter< TInputImage, TInputImage, TInputImage > AddType;
  typename AddType::Pointer adder = AddType::New();
  adder->SetInput1( this->GetMaskImage() );
  adder->SetInput2( constantImage );

  //Find the min of these two images.  The background pixels of the marker image have no effect because they have been set to the maximum of the pixel type.  The marker pixels are always kept because they correspond with the smallest values of the image type.  So, the result is the same as the mask image with min values placed at all of the marker locations.  This is the new mask image.
  typedef itk::MinimumImageFilter< TInputImage, TInputImage, TInputImage > MinType;
  typename MinType::Pointer minFilter = MinType::New();
  minFilter->SetInput1( invertedImage );
  minFilter->SetInput2( adder->GetOutput() );

  typedef itk::ReconstructionByErosionImageFilter< TInputImage, TOutputImage > ReconstructionType;
  typename ReconstructionType::Pointer reconstructionFilter = ReconstructionType::New();
  reconstructionFilter->SetMarkerImage( invertedImage );
  reconstructionFilter->SetMaskImage( minFilter->GetOutput() );
  reconstructionFilter->Update();


  /** graft the minipipeline output back into this filter's output */
  this->GraftOutput( reconstructionFilter->GetOutput() );

//   // Stop the timer
//   time (&end);
//   diff = difftime (end,start);
//   std::cout << "Time taken: " << diff << std::endl;
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TMarkerImage, class TOutputImage>
void
ImposeMinimaImageFilter< TInputImage, TMarkerImage, TOutputImage >
::PrintSelf(
  std::ostream& os,
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "FullyConnected: "  << m_FullyConnected << std::endl;

}

} // end namespace itk

#endif
