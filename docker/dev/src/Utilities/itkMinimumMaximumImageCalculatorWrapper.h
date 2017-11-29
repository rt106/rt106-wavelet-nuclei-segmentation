// Copyright (c) General Electric Company, 2017.  All rights reserved.
#ifndef __itkMinimumMaximumImageCalculatorWrapper_h
#define __itkMinimumMaximumImageCalculatorWrapper_h

#include "itkMinimumMaximumImageCalculator.h"

namespace itk
{
/** \class itkMinimumMaximumImageCalculatorWrapper.h
 *
 * Author: Dirk Padfield, GE Global Research, padfield@research.ge.com
 *
 */

template <class TInputImage>
bool MinimumMaximumImageCalculatorWrapper( TInputImage const * inputImage, typename TInputImage::PixelType & minimum, typename TInputImage::PixelType & maximum )
{
  typedef itk::MinimumMaximumImageCalculator<TInputImage> CalculatorType;
  typename CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetImage( inputImage );
  calculator->Compute();

  minimum = calculator->GetMinimum();
  maximum = calculator->GetMaximum();

  return true;
}

} // end namespace itk

#endif
