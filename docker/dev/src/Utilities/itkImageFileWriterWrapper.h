// Copyright (c) General Electric Company, 2017.  All rights reserved.
#ifndef __itkImageFileWriterWrapper_h
#define __itkImageFileWriterWrapper_h

#include "itkImageFileWriter.h"
#include "itkIntensityWindowingImageFilterWrapper.h"

namespace itk 
{
/** \class itkImageFileWriterWrapper.h
 *
 * Author: Dirk Padfield, GE Global Research, padfield@research.ge.com
 *
 */

template < class TInputImage >
bool ImageFileWriterWrapper(TInputImage const * inputImage, std::string filename, bool useCompression = true)
{
  typedef ImageFileWriter<TInputImage> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput( inputImage );
  writer->SetUseCompression( useCompression );

  std::cout << "Writing out image " << filename << std::endl;
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

  return true;
}


template < class TInputImage >
bool SimpleIntensityRescaleImageFileWriter(TInputImage const * inputImage, std::string filename, bool useCompression = true)
{
  typedef itk::Image<unsigned char, TInputImage::ImageDimension> TOutputImage;
  typename TOutputImage::PixelType outputMinimum = 0;
  typename TOutputImage::PixelType outputMaximum = 255;
  typename TOutputImage::Pointer rescaledImage = IntensityWindowingImageFilterWrapper<TInputImage,TOutputImage>( inputImage, outputMinimum, outputMaximum );

  ImageFileWriterWrapper<TOutputImage>(rescaledImage,filename,useCompression);

  return true;
}


} // end namespace itk

#endif
