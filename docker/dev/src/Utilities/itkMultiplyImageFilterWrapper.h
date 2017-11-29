// Copyright (c) General Electric Company, 2017.  All rights reserved.
#ifndef __itkMultiplyImageFilterWrapper_h
#define __itkMultiplyImageFilterWrapper_h

#include "itkMultiplyImageFilter.h"

namespace itk
{
/** \class itkMultiplyImageFilterWrapper.h
 *
 * This is a simple function that enables the
 * itkMultiplyImageFilter to be called as a function and using
 * just one line of code.
 * 
 * It should be called as follows:
 * outputImage = itk::MultiplyImageFilterWrapper<InputImageType1,InputImageType2,OutputImageType>(inputImage1,inputImage2);
 *
 * Author: Dirk Padfield, GE Global Research, padfield@research.ge.com
 *
 */

template <typename TInputImage1, typename TInputImage2, typename TOutputImage>
typename TOutputImage::Pointer MultiplyImageFilterWrapper( TInputImage1 const * inputImage1, TInputImage2 const * inputImage2 )
{
  typedef typename itk::MultiplyImageFilter< TInputImage1, TInputImage2, TOutputImage > MultiplyImageType;
  typename MultiplyImageType::Pointer multiplyFilter = MultiplyImageType::New();
  multiplyFilter->SetInput1( inputImage1 );
  multiplyFilter->SetInput2( inputImage2 );
  multiplyFilter->Update();

  return multiplyFilter->GetOutput();
}


// Overloaded method taking a constant as the second argument.
template <typename TInputImage1, typename TInputImage2, typename TOutputImage>
typename TOutputImage::Pointer MultiplyImageFilterWrapper( TInputImage1 const * inputImage1, typename TInputImage2::PixelType inputConstant2 )
{
  typedef typename itk::MultiplyImageFilter< TInputImage1, TInputImage2, TOutputImage > MultiplyImageType;
  typename MultiplyImageType::Pointer multiplyFilter = MultiplyImageType::New();
  multiplyFilter->SetInput1( inputImage1 );
  multiplyFilter->SetConstant2( inputConstant2 );
  multiplyFilter->Update();

  return multiplyFilter->GetOutput();
}
} // end namespace itk


#endif
