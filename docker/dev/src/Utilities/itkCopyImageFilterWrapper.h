// Copyright (c) General Electric Company, 2017.  All rights reserved.
#ifndef __itkCopyImageFilterWrapper_h
#define __itkCopyImageFilterWrapper_h

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk 
{
/** \class itkCopyImageFilterWrapper.h
 *
 * Author: Dirk Padfield, GE Global Research, padfield@research.ge.com
 *         Alberto Santamaria, GE Global Research, santamar@research.ge.com
 * Date: August 2011
 *
 */

template < class TInputImage, class TOutputImage >
typename TOutputImage::Pointer CopyImageFilterWrapper(TInputImage const * inputImage )
{
  typename TOutputImage::Pointer outputImage = TOutputImage::New();
  outputImage->SetRegions( inputImage->GetLargestPossibleRegion() );
  //outputImage->SetLargestPossibleRegion( inputImage->GetLargestPossibleRegion() );
  //outputImage->SetBufferedRegion( inputImage->GetBufferedRegion() );
  //outputImage->SetRequestedRegion( inputImage->GetRequestedRegion() );
  outputImage->Allocate();
  outputImage->CopyInformation( inputImage );
  // Graft copies the information, sets up the regions, and copies the pixel container.
  // However, Graft does not work in this case for copying the data because it simply sets the pointer to the memory location of the data.
  // So we need to do a deep copy of the image using iterators.
  itk::ImageRegionConstIterator<TInputImage> inputIt( inputImage, inputImage->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<TOutputImage> outputIt( outputImage, outputImage->GetLargestPossibleRegion() );
  for( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++outputIt )
  {
    outputIt.Set( inputIt.Get() );
  }

  return outputImage;	
}
} // end namespace itk

#endif
