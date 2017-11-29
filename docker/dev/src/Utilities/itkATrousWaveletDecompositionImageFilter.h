// Copyright (c) General Electric Company, 2017.  All rights reserved.
/*=========================================================================

  Author: Dirk Padfield
  Date: 03/19/2012

=========================================================================*/
#ifndef __itkATrousWaveletDecompositionImageFilter_h
#define __itkATrousWaveletDecompositionImageFilter_h

#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include <vcl_vector.h>

namespace itk
{
/** \class ATrousWaveletDecompositionImageFilter

Decomposition: 
approximation image:   	I_i(x,y) = \sum_{m,n} h(m,n) I_{i-1}(x-2^{i-1}m, y-2^{i-1}n)
detail images: 		W_i(x,y) = I_i(x,y)- I_{i+1}(x,y)
			h(m,n): scaling function (use B3 spline) 
			 I_0(x,y) is the original image
Reconstruction
  I_0(x,y)=I_s(x,y)+\sum_{i=1}^s W_i(x,y)

Advantages over standard wavelet:
no need to downsampling and upsampling

The output is the final approximation image.  It can also be obtained using the GetApproximationImage method.
The detail images is a vector of images corresponding to the different levels of the decomposition.
They are obtained using the GetDetailImages() method.

Reference: 
 */

template <class TInputImage, class TOutputImage>
class ATrousWaveletDecompositionImageFilter :
public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Run-time type information (and related methods). */
  itkTypeMacro(ATrousWaveletDecompositionImageFilter, ImageToImageFilter);

  /** Extract dimension from input and output image. */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

  /** Standard class typedefs. */
  typedef ATrousWaveletDecompositionImageFilter Self;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;

  /** Some typedefs associated with the input images. */
  typedef TInputImage InputImageType;
  typedef typename InputImageType::Pointer          InputImagePointer;
  typedef typename InputImageType::ConstPointer     InputImageConstPointer;
  typedef typename InputImageType::RegionType       InputImageRegionType;

  typedef TOutputImage OutputImageType;
  typedef typename OutputImageType::Pointer         OutputImagePointer;
  typedef typename OutputImageType::RegionType      OutputImageRegionType; 
  typedef typename OutputImageType::PixelType       OutputImagePixelType;
  typedef typename vcl_vector<OutputImagePointer>	VectorOutputImagePointer;

  typedef itk::Image< float, ImageDimension > InternalImageType;
  typedef typename InternalImageType::Pointer   InternalImagePointer;
  typedef typename InternalImageType::SizeType  InternalSizeType;
  typedef typename InternalImageType::PixelType  InternalPixelType;
  typedef typename InternalImageType::IndexType  InternalIndexType;
  typedef typename InternalImageType::RegionType  InternalRegionType;
  typedef typename InternalImageType::SpacingType InternalSpacingType;

  typedef typename InternalImageType::PixelType     RealPixelType;
  typedef Image< RealPixelType, ImageDimension>     RealImageType;
  typedef typename RealImageType::Pointer           RealImagePointer;
  typedef typename RealImageType::IndexType         RealIndexType;
  typedef typename RealImageType::SizeType          RealSizeType;
  typedef typename RealImageType::RegionType        RealRegionType;

  typedef Image< std::complex<RealPixelType>, ImageDimension >  FFTImageType;
  typedef typename FFTImageType::Pointer                        FFTImagePointer;


  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  VectorOutputImagePointer GetDetailImages()
  {
    return this->m_DetailImages;
    //return static_cast<vcl_vector<OutputImageType*>>(const_cast<vcl_vector<OutputImageType*>>(m_DetailImages));
  };

  OutputImageType* GetApproximationImage()
  {
    return static_cast<OutputImageType*>(const_cast<typename OutputImageType::Pointer>(this->GetOutput(0)));
  };

  const InputImageType* GetInputImage()
  {
    return static_cast<InputImageType*>(const_cast<DataObject *>(this->ProcessObject::GetInput(0)));
  };

  itkSetMacro(LowerLevel, unsigned short);
  itkGetConstReferenceMacro(LowerLevel, unsigned short);

  itkSetMacro(UpperLevel, unsigned short);
  itkGetConstReferenceMacro(UpperLevel, unsigned short);


  typedef vnl_matrix<typename OutputImageType::PixelType> MatrixType;

protected:
  ATrousWaveletDecompositionImageFilter();
  virtual ~ATrousWaveletDecompositionImageFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;
  void GenerateData();

  InternalImagePointer ConvolveImageVNL( InternalImageType* inputImage, int level );
  InternalImagePointer ConvolveImageITK( InternalImageType * inputImage, int level );
  InternalImagePointer ConvolveImageFast( InputImageType const * inputImage, int level );

  vnl_vector<double> GenerateGaussianKernel(float sigma, unsigned int kernelRadius, float dimSpacing);


private:
  ATrousWaveletDecompositionImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implementaed

  unsigned short m_LowerLevel;
  unsigned short m_UpperLevel;
  VectorOutputImagePointer m_DetailImages;

  InternalSpacingType m_ImageSpacing;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkATrousWaveletDecompositionImageFilter.hxx"
#endif

#endif
