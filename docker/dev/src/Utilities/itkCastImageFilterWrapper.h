// Copyright (c) General Electric Company, 2017.  All rights reserved.
#ifndef __itkCastImageFilterWrapper_h
#define __itkCastImageFilterWrapper_h

#include "itkCastImageFilter.h"

namespace itk
{
/** \class itkCastImageFilterWrapper.h
 *
 * Author: Dirk Padfield, GE Global Research, padfield@research.ge.com
 *
 */

template <class TInputImage, class TOutputImage >
typename TOutputImage::Pointer CastImageFilterWrapper( TInputImage * inputImage )
{
  typedef itk::CastImageFilter<TInputImage,TOutputImage> FilterType;
  typename FilterType::Pointer caster = FilterType::New();
  caster->SetInput( inputImage );
  caster->Update();
  return caster->GetOutput();
}

} // end namespace itk

#endif


