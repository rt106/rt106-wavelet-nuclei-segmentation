// Copyright (c) General Electric Company, 2017.  All rights reserved.
#ifndef __itkBinaryThresholdImageFilterWrapper_h
#define __itkBinaryThresholdImageFilterWrapper_h

#include "itkBinaryThresholdImageFilter.h"

namespace itk
{
/** \class itkBinaryThresholdImageFilterWrapper.h
 * 
 * This is a simple function that enables the 
 * itkBinaryThresholdImageFilter to be called as a function and using
 * just one line of code rather than 8.
 *
 * It should be called as follows:
 * outputImage = itk::BinaryThresholdImageFilterWrapper<TInputImage,TOutputImage>( inputImage, lowerThreshold, upperThreshold, insideValue, outsideValue )
 *
 * It would be better and more intuitive to the user if we could use
 * operator overloading in this filter.
 *
 * Author: Dirk Padfield, GE Global Research, padfield@research.ge.com
 *
 */

template <class TInputImage, class TOutputImage>
typename TOutputImage::Pointer BinaryThresholdImageFilterWrapper( TInputImage const * inputImage, typename TInputImage::PixelType lowerThreshold, typename TInputImage::PixelType upperThreshold, typename TOutputImage::PixelType insideValue = 255, typename TOutputImage::PixelType outsideValue = 0)
{
    typedef itk::BinaryThresholdImageFilter< TInputImage, TOutputImage > ThresholdType;
    typename ThresholdType::Pointer thresholder = ThresholdType::New();
    thresholder->SetInput( inputImage );
    thresholder->SetLowerThreshold( lowerThreshold );
    thresholder->SetUpperThreshold( upperThreshold );
    thresholder->SetInsideValue( insideValue );
    thresholder->SetOutsideValue( outsideValue );
    thresholder->Update();

    return thresholder->GetOutput();
}

} // end namespace itk

#endif
