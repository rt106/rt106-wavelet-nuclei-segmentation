/*=========================================================================
  Copyright (c) General Electric Company, 2017.  All rights reserved.
  Author: Dirk Padfield
  date: 06/02/2006
=========================================================================*/
#ifndef __itkRegionalExtremaImageFilter_h
#define __itkRegionalExtremaImageFilter_h

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryThresholdImageFunction.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkFloodFilledImageFunctionConditionalIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkShiftScaleImageFilter.h"

#include "itkReconstructionByErosionImageFilter.h"
#include "itkReconstructionByDilationImageFilter.h"
#include "itkSubtractImageFilter.h"

#include <stack>

namespace itk
{
/** \class RegionalExtremaImageFilter

Algorithm reference: E. Breen and R. Jones, "Attribute openings, thinnings, and
granulometries," Computer Vision and Image Understanding, vol 64 n 3,
pp. 337-389.  Implementation here follows the description in P. Soille,
Morphological Image Analysis: Principles and Applications, Springer, 1999,
p. 169.


Give the description of what "regional" means.

Regional maxima are connected components of pixels with the same
    intensity value, t, whose external boundary pixels all have a value
    less than t.

Describe the default settings.

To change whether finding maxima or minima, use FindMinimaOff() and
FindMinimaOn().  The default is FindMinimaOn().

In addition to being
faster, the fast method allows the setting of an arbitrary sized
neighborhood.  
The slow method uses reconstruction after shifting
the image down by 1, so it is not as flexible.  Also, if the image
spacing is not 1, then shifting by 1 will not reconstruct by 1 pixel.
If the spacing is less than 1, then it will reconstruct by several
pixels, and if it is greater than 1, it will reconstruct by less than
1 pixel.


 */
template < class TInputImage, class TOutputImage >
class ITK_EXPORT RegionalExtremaImageFilter :
    public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
   /** Standard class typedefs. */
  typedef RegionalExtremaImageFilter Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(RegionalExtremaImageFilter, ImageToImageFilter);

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

  typedef typename InputImageType::SizeType RadiusType;
  typedef typename InputImageType::IndexType InputIndexType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Set/Get the radius of the neighborhood. */
  itkSetMacro(Radius, RadiusType);
  itkGetConstReferenceMacro(Radius, RadiusType);

  /** Set/Get the foreground value */
  itkSetMacro(ForegroundValue, OutputImagePixelType);
  itkGetConstReferenceMacro(ForegroundValue, OutputImagePixelType);

  /** Set/Get the background value*/
  itkSetMacro(BackgroundValue, OutputImagePixelType);
  itkGetConstReferenceMacro(BackgroundValue, OutputImagePixelType);

  /** Set/Get whether to find the minima */
  itkSetMacro(FindMinima,bool);
  itkBooleanMacro(FindMinima);
  itkGetConstReferenceMacro(FindMinima,bool);

  /** Set/Get whether to use the fast method  */
  itkSetMacro(UseFastMethod,bool);
  itkBooleanMacro(UseFastMethod);
  itkGetConstReferenceMacro(UseFastMethod,bool);

  /** Set/Get whether the output is binary or the extrema have the
   * intensity of the input image.  */
  itkSetMacro(OutputIsBinary,bool);
  itkBooleanMacro(OutputIsBinary);
  itkGetConstReferenceMacro(OutputIsBinary,bool);


protected:
  RegionalExtremaImageFilter();
  virtual ~RegionalExtremaImageFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;
  void  GenerateData ();

private:
  RegionalExtremaImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  RadiusType m_Radius;
  OutputImagePixelType m_ForegroundValue;
  OutputImagePixelType m_BackgroundValue;

  bool m_FindMinima;
  bool m_UseFastMethod;

  bool m_OutputIsBinary;

};

} // end namespace itk

#ifndef GEITK_MANUAL_INSTANTIATION
#include "itkRegionalExtremaImageFilter.txx"
#endif

#endif
