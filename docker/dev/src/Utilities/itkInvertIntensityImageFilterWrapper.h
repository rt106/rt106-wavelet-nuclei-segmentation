// Copyright (c) General Electric Company, 2017.  All rights reserved.
#ifndef __itkInvertIntensityImageFilterWrapper_h
#define __itkInvertIntensityImageFilterWrapper_h

#include "itkInvertIntensityImageFilter.h"
#include "itkMinimumMaximumImageCalculatorWrapper.h"

namespace itk
{
/** \class itkInvertIntensityImageFilterWrapper.h
 *
 * Author: Dirk Padfield, GE Global Research, padfield@research.ge.com
 *
 */

// Version that takes the maximum intensity as an input rather than computing it.
template <class TInputImage >
typename TInputImage::Pointer InvertIntensityImageFilterWrapper( TInputImage const * inputImage, typename TInputImage::PixelType maximum )
{
  typedef itk::InvertIntensityImageFilter<TInputImage> InvertIntensityType;
  typename InvertIntensityType::Pointer inverter = InvertIntensityType::New();
  inverter->SetInput( inputImage );
  inverter->SetMaximum( maximum );

  try
  {
    inverter->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
  }

  typename TInputImage::Pointer invertedImage = inverter->GetOutput();
  invertedImage->DisconnectPipeline();
  
  return invertedImage;
}

// Version that computes the maximum intensity directly from the image.
template <class TInputImage >
typename TInputImage::Pointer InvertIntensityImageFilterWrapper( TInputImage const * inputImage )
{
  typename TInputImage::PixelType minimum;
  typename TInputImage::PixelType maximum;
  MinimumMaximumImageCalculatorWrapper<TInputImage>( inputImage, minimum, maximum );

  return InvertIntensityImageFilterWrapper<TInputImage>(inputImage, maximum);
}


} // end namespace itk

#endif

