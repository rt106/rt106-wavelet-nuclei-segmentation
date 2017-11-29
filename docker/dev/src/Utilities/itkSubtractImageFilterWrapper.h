// Copyright (c) General Electric Company, 2017.  All rights reserved.

#ifndef __itkSubtractImageFilterWrapper_h
#define __itkSubtractImageFilterWrapper_h

#include "itkSubtractImageFilter.h"

namespace itk
{
/** \class itkSubtractImageFilterWrapper.h
 *
 * Author: Dirk Padfield, GE Global Research, padfield@research.ge.com
 *
 */

template <class TInputImage1, class TInputImage2, class TOutputImage >
typename TOutputImage::Pointer SubtractImageFilterWrapper( TInputImage1 const * inputImage1, TInputImage2 const * inputImage2 )
{
  typedef itk::SubtractImageFilter< TInputImage1, TInputImage2, TOutputImage > FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1( inputImage1 ); 
  filter->SetInput2( inputImage2 );

  try
    {
    filter->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

  return filter->GetOutput();
}

} // end namespace itk

#endif






