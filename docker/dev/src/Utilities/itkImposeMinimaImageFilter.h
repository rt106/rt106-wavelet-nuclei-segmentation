// Copyright (c) General Electric Company, 2017.  All rights reserved.
/*=========================================================================

  Author: Dirk Padfield
  date: 03/03/2007

=========================================================================*/
#ifndef __itkImposeMinimaImageFilter_h
#define __itkImposeMinimaImageFilter_h

#include "itkImage.h"
//#include "itkReconstructionImageFilter.h"
#include "itkReconstructionByErosionImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkMinimumImageFilter.h"
#include "itkShiftScaleImageFilter.h"

//#include "FileDumper.h"

namespace itk
{
/** \class ImposeMinimaImageFilter

This filter imposes minima defined in a marker image on the grayscale
mask image.

 */

template < class TInputImage, class TMarkerImage, class TOutputImage >
class ITK_EXPORT ImposeMinimaImageFilter :
    public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Run-time type information (and related methods). */
  itkTypeMacro(ImposeMinimaImageFilter, ImageToImageFilter);

  /** Extract dimension from input and output image. */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

  /** Some typedefs associated with the input images. */
  typedef TInputImage InputImageType;
  typedef typename InputImageType::Pointer          InputImagePointer;
  typedef typename InputImageType::ConstPointer     InputImageConstPointer;
  typedef typename InputImageType::RegionType       InputImageRegionType;

  typedef TOutputImage OutputImageType;
  typedef typename OutputImageType::Pointer         OutputImagePointer;
  typedef typename OutputImageType::RegionType      OutputImageRegionType; 


  /** Image typedef support. */
  typedef typename InputImageType::PixelType InputImagePixelType;
  typedef typename OutputImageType::PixelType OutputImagePixelType;


  /** Standard class typedefs. */
  typedef ImposeMinimaImageFilter Self;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  
  /** Some convenient typedefs. */
  typedef TMarkerImage MarkerImageType;
  typedef TInputImage MaskImageType;


  typedef typename InputImageType::SizeType        ISizeType;
  typedef typename MarkerImageType::Pointer        MarkerImagePointer;
  typedef typename MarkerImageType::ConstPointer   MarkerImageConstPointer;
  typedef typename MarkerImageType::RegionType     MarkerImageRegionType;
  typedef typename MarkerImageType::PixelType      MarkerImagePixelType;
  typedef typename InputImageType::IndexType       InputImageIndexType;
  typedef typename MaskImageType::Pointer          MaskImagePointer;
  typedef typename MaskImageType::ConstPointer     MaskImageConstPointer;
  typedef typename MaskImageType::RegionType       MaskImageRegionType;
  typedef typename MaskImageType::PixelType        MaskImagePixelType;
  typedef typename OutputImageType::ConstPointer   OutputImageConstPointer;
  typedef typename OutputImageType::IndexType      OutputImageIndexType;



  /** Method for creation through the object factory. */
  itkNewMacro(Self);


  /** Set/Get the marker image. Traditionally, the marker image must
   * be pixelwise less than or equal to the mask image (for dilation),
   * however this filter implicitly applies a mask to force the
   * constraint to hold. */
  void SetMarkerImage(const MarkerImageType *);
  const MarkerImageType* GetMarkerImage();

  /** Set/Get the mask image. The mask image is used to "mask" the
   * dilated marker image. The mask operation is a pixelwise
   * minimum. */
  void SetMaskImage(const MaskImageType *);
  const MaskImageType* GetMaskImage();

  /**
   * Set/Get whether the connected components are defined strictly by
   * face connectivity or by face+edge+vertex connectivity.  Default is
   * FullyConnectedOff.  For objects that are 1 pixel wide, use
   * FullyConnectedOn.
   */
  itkSetMacro(FullyConnected, bool);
  itkGetConstReferenceMacro(FullyConnected, bool);
  itkBooleanMacro(FullyConnected);


protected:
  ImposeMinimaImageFilter();
  virtual ~ImposeMinimaImageFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;
  void  GenerateData ();

private:
  ImposeMinimaImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  bool                m_FullyConnected;



};

} // end namespace itk

#ifndef GEITK_MANUAL_INSTANTIATION
#include "itkImposeMinimaImageFilter.txx"
#endif

#endif
