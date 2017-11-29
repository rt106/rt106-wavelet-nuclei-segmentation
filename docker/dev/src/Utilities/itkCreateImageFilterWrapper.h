// Copyright (c) General Electric Company, 2017.  All rights reserved.
#ifndef __itkCreateImageFilterWrapper_h
#define __itkCreateImageFilterWrapper_h

#include "itkImage.h"

namespace itk
{
/** \class itkCreateImageFilterWrapper.h
 *
 * This is a simple function to create an image of a specified size
 * filled by a specified value.
 * 
 * Author: Dirk Padfield, GE Global Research, padfield@research.ge.com
 *
 */

template <class TOutputImage >
typename TOutputImage::Pointer CreateImageFilterWrapper( typename TOutputImage::SizeType imageSize, typename TOutputImage::PixelType pixelValue = 0 )
{
  typedef TOutputImage OutputImageType;

  // Create an image.
  typename OutputImageType::IndexType imageStart;
  imageStart.Fill( 0 );
  typename OutputImageType::RegionType imageRegion;
  imageRegion.SetSize( imageSize );
  imageRegion.SetIndex( imageStart );
    
  typename OutputImageType::Pointer outputImage = OutputImageType::New();
  outputImage->SetRegions( imageRegion );
  outputImage->Allocate();
  outputImage->FillBuffer( pixelValue );

  return outputImage;
}

} // end namespace itk

#endif
