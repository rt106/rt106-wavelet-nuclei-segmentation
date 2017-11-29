// Copyright (c) General Electric Company, 2017.  All rights reserved.

#ifndef __itkRecursiveGaussianImageFilterWrapper_h
#define __itkRecursiveGaussianImageFilterWrapper_h

#include "itkImage.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkCastImageFilter.h"

namespace itk
{
/** \class itkRecursiveGaussianImageFilterWrapper.h
 *
 * Author: Dirk Padfield, GE Global Research, padfield@research.ge.com
 *
 */

template <class TInputImage, class TOutputImage>
typename TOutputImage::Pointer RecursiveGaussianImageFilterWrapper( TInputImage const * inputImage, double sigma )
{
  const unsigned int ImageDimension = TInputImage::ImageDimension;
  typedef itk::Image<float, ImageDimension> RealImageType;
  typedef TInputImage InputImageType;
  typedef TOutputImage OutputImageType;

  typedef itk::CastImageFilter< RealImageType, OutputImageType > CastType;
  typename CastType::Pointer caster = CastType::New();

  typedef itk::RecursiveGaussianImageFilter<InputImageType,RealImageType> GaussianType1;
  typedef itk::RecursiveGaussianImageFilter<RealImageType,RealImageType> GaussianType2;
  typename GaussianType1::Pointer gaussianX = GaussianType1::New();
  gaussianX->SetInput( inputImage );
  gaussianX->SetSigma( sigma );
  gaussianX->SetDirection( 0 );   // 0 --> X direction
  gaussianX->ReleaseDataFlagOn();
  gaussianX->Update();

  typename GaussianType2::Pointer gaussianY = GaussianType2::New();
  gaussianY->SetInput( gaussianX->GetOutput() );
  gaussianY->SetSigma( sigma );
  gaussianY->SetDirection( 1 );   // 1 --> Y direction
  gaussianY->ReleaseDataFlagOn();
  gaussianY->Update();
  
  if( ImageDimension == 2 )
    {
    // The input is a 2D image
    caster->SetInput( gaussianY->GetOutput() );
    caster->Update();
    }
  else
    {
    // The input is a 3D image.
    typename GaussianType2::Pointer gaussianZ = GaussianType2::New();
    gaussianZ->SetInput( gaussianY->GetOutput() );
    gaussianZ->SetSigma( sigma );
    gaussianZ->SetDirection( 2 );   // 2 --> Z direction
    gaussianZ->ReleaseDataFlagOn();
    gaussianZ->Update();
    
    caster->SetInput( gaussianZ->GetOutput() );
    caster->Update();
    }
  
  return caster->GetOutput();
}

} // end namespace itk

#endif


