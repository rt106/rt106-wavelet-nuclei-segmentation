// Copyright (c) General Electric Company, 2017.  All rights reserved.
#ifndef __itkIntensityWindowingImageFilterWrapper_h
#define __itkIntensityWindowingImageFilterWrapper_h

#include "itkIntensityWindowingImageFilter.h"
#include "itkMinimumMaximumImageCalculatorWrapper.h"

namespace itk
{
/** \class itkIntensityWindowingImageFilterWrapper.h
 *
 * Author: Dirk Padfield, GE Global Research, padfield@research.ge.com
 *
 */

template <class TInputImage, class TOutputImage>
typename TOutputImage::Pointer IntensityWindowingImageFilterWrapper( TInputImage const * inputImage, typename TInputImage::PixelType inputMinimum, typename TInputImage::PixelType inputMaximum, typename TOutputImage::PixelType outputMinimum, typename TOutputImage::PixelType outputMaximum )
{
  typedef typename itk::IntensityWindowingImageFilter< TInputImage, TOutputImage > WindowingFilter;
  typename WindowingFilter::Pointer windower = WindowingFilter::New();
  windower->SetInput( inputImage );
  windower->SetWindowMinimum( inputMinimum  );
  windower->SetWindowMaximum( inputMaximum  );
  windower->SetOutputMinimum( outputMinimum );
  windower->SetOutputMaximum( outputMaximum );
  windower->Update();

  return windower->GetOutput();
}

// Overloaded function.
// If only two parameters are passed in, we automatically compute the input minimum and maximum and use those for the
// window minimum and maximum.
template <class TInputImage, class TOutputImage>
typename TOutputImage::Pointer IntensityWindowingImageFilterWrapper( TInputImage const * inputImage, typename TOutputImage::PixelType outputMinimum, typename TOutputImage::PixelType outputMaximum )
{
  typedef typename itk::IntensityWindowingImageFilter< TInputImage, TOutputImage > WindowingFilter;
  typename WindowingFilter::Pointer windower = WindowingFilter::New();
  windower->SetInput( inputImage );

  // Compute the input minimum and maximum.
  typename TInputImage::PixelType inputMinimum;
  typename TInputImage::PixelType inputMaximum;
  MinimumMaximumImageCalculatorWrapper<TInputImage>(inputImage,inputMinimum,inputMaximum);

  windower->SetWindowMinimum( inputMinimum  );
  windower->SetWindowMaximum( inputMaximum  );
  windower->SetOutputMinimum( outputMinimum );
  windower->SetOutputMaximum( outputMaximum );
  windower->Update();

  return windower->GetOutput();
}


} // end namespace itk

#endif


