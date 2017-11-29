// Copyright (c) General Electric Company, 2017.  All rights reserved.
#ifndef __itkATrousWaveletBlobSegmentationImageFilter_h
#define __itkATrousWaveletBlobSegmentationImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkConstantBoundaryCondition.h"

#include <vcl_vector.h>

// ITK
#include "itkATrousWaveletDecompositionImageFilter.h"

namespace itk
{
/** \class itkATrousWaveletBlobSegmentationImageFilter.h
  This function decomposes the input image using a trous wavelets.  Then,
  using a noise estimate calculated from a representative background patch,
  it thresholds the wavelet coefficients.  Finally, it combines several
  levels to reduce the noise.

  LowerLevel and UpperLevel define the lower and upper levels of the decomposition to combine.
  The numbering of the levels starts with 1, not zero.

  Author: Dirk Padfield, GE Global Research, padfield@research.ge.com

  References:
  [1] Dirk Padfield. "Segmentation and Tracking Algorithms for Monitoring Cellular Motion and Function."  PhD thesis, Rensselaer Polytechnic Institute, 2009.
  [2] Dirk Padfield, Jens Rittscher, and Badrinath Roysam. "Coupled minimum-cost flow cell tracking."  In Information Processing in Medical Imaging (IPMI), 2009.
*/

template <class TInputImage,class TOutputImage>
class ATrousWaveletBlobSegmentationImageFilter :
public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef ATrousWaveletBlobSegmentationImageFilter    Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  itkNewMacro(Self);
  itkTypeMacro( ATrousWaveletBlobSegmentationImageFilter, ImageToImageFilter );

  /** Type for input image. */
  typedef   TInputImage                        InputImageType;
  itkStaticConstMacro(ImageDimension, unsigned int, InputImageType::ImageDimension);
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::PointType      InputPointType;
  typedef typename InputPointType::CoordRepType   InputCoordType;
  typedef typename InputImageType::RegionType     InputRegionType;
  typedef typename InputImageType::IndexType      InputIndexType;
  typedef typename InputIndexType::IndexValueType  InputIndexValueType;
  typedef typename InputImageType::SizeType       InputSizeType;
  typedef typename InputImageType::PixelType      InputPixelType;

  typedef   TOutputImage      OutputImageType;
  typedef typename OutputImageType::Pointer        OutputImagePointer;
  typedef typename OutputImageType::ConstPointer   OutputConstImagePointer;
  typedef typename OutputImageType::PointType      OutputPointType;
  typedef typename OutputPointType::CoordRepType   OutputCoordType;
  typedef typename OutputImageType::RegionType     OutputRegionType;
  typedef typename OutputImageType::IndexType      OutputIndexType;
  typedef typename OutputIndexType::IndexValueType OutputIndexValueType;
  typedef typename OutputImageType::SizeType       OutputSizeType;
  typedef typename OutputImageType::PixelType      OutputPixelType;

  typedef itk::Image<float,ImageDimension>         RealImageType;
  typedef typename RealImageType::Pointer      RealImagePointer;

  //typedef InputImageType    BinaryImageType;
  typedef typename itk::Image<unsigned char, ImageDimension>    BinaryImageType;
  typedef typename BinaryImageType::Pointer                     BinaryImagePointer;


  typedef itk::ATrousWaveletDecompositionImageFilter<InputImageType,RealImageType> DecompositionType;
  
  /** Boundary condition type for the neighborhood iterator */
  typedef ConstantBoundaryCondition< InputImageType > ConstBoundaryConditionType;

  /** Neighborhood iterator type */
  typedef NeighborhoodIterator<InputImageType, ConstBoundaryConditionType> NeighborhoodIteratorType;

  /** Neighborhood type */
  typedef typename NeighborhoodIteratorType::NeighborhoodType NeighborhoodType;

  typedef ImageRegionConstIterator< InputImageType	   >  ConstIteratorType;

  itkSetMacro(SigmaMultiplier, float);
  itkGetMacro(SigmaMultiplier, float);

  itkSetMacro(LowerLevel, unsigned short);
  itkGetMacro(LowerLevel, unsigned short);

  itkSetMacro(UpperLevel, unsigned short);
  itkGetMacro(UpperLevel, unsigned short);

  /** Set/Get whether image spacing should be used in computing distances. */
  /** Set On/Off whether spacing is used. */
  itkSetMacro( UseImageSpacing, bool );
  itkGetConstReferenceMacro( UseImageSpacing, bool );
  itkBooleanMacro( UseImageSpacing );
  
  typename DecompositionType::VectorOutputImagePointer GetDetailImagesThresholded()
  {
    return this->m_DetailImagesThresholded;
  };

  itkSetMacro( DecomposeByDifferenceOfGaussians, bool );
  itkGetConstReferenceMacro( DecomposeByDifferenceOfGaussians, bool );
  itkBooleanMacro( DecomposeByDifferenceOfGaussians );

  void SetSmoothingVariances( std::vector<float> variances )
  {
  if( variances.size() < 2 )
    {
    std::cerr << "The SmoothingVariances vector must have at least two entries!" << std::endl;
    return;
    }
  m_SmoothingVariances.clear();
  for( std::vector<float>::iterator it = variances.begin(); it != variances.end(); ++it )
    {
    if( *it <= 0.0 )
      {
      m_SmoothingVariances.push_back( vnl_math::eps );
      }
    else
      {
      m_SmoothingVariances.push_back( *it );
      }
    }
  // Ensure that values are sorted by increasing value.
  std::sort(m_SmoothingVariances.begin(),m_SmoothingVariances.end());
  }

protected:
  ATrousWaveletBlobSegmentationImageFilter();
  ~ATrousWaveletBlobSegmentationImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateData(void);

  BinaryImagePointer CalculateBackgroundMask( const InputImageType * inputImage, size_t & numberOfBackgroundPixels );
  float ABE_Thresholder(float inputPixelValue, float noiseSigma, float sigmaMultiplier = 3);

private:   
  ATrousWaveletBlobSegmentationImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  float				m_SigmaMultiplier;
  unsigned short	m_LowerLevel;
  unsigned short	m_UpperLevel;

  typename DecompositionType::VectorOutputImagePointer m_DetailImagesThresholded;

  bool m_UseImageSpacing;

  // For decomposition using difference of Gaussians.
  bool m_DecomposeByDifferenceOfGaussians;
  std::vector<float> m_SmoothingVariances;
  typename DecompositionType::VectorOutputImagePointer DifferenceOfGaussians(const InputImageType * inputImage);
  

}; // end of ATrousWaveletBlobSegmentationImageFilter class

} //end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkATrousWaveletBlobSegmentationImageFilter.hxx"
#endif


#endif
