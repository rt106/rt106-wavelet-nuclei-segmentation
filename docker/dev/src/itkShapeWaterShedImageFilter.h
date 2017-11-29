// Copyright (c) General Electric Company, 2017.  All rights reserved.

#ifndef _itkShapeWaterShedImageFilter_H
#define _itkShapeWaterShedImageFilter_H

#include <itkImageToImageFilter.h>
// This is implemented by Dirk and re-organized by Xiaofeng

// The input should be a binary image. 

namespace itk
{
/** \class itkShapeWaterShedImageFilter.h
*/

template <class TInputImage,class TOutputImage>
class ShapeWaterShedImageFilter :
    public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef ShapeWaterShedImageFilter    Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  itkNewMacro(Self);
  itkTypeMacro( ShapeWaterShedImageFilter, ImageToImageFilter );

  /** Type for input image. */
  typedef   TInputImage                        InputImageType;
  itkStaticConstMacro(Dimension, unsigned int, InputImageType::ImageDimension);
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputConstImagePointer;
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
  typedef typename OutputIndexType::IndexValueType  OutputIndexValueType;
  typedef typename OutputImageType::SizeType       OutputSizeType;
  typedef typename OutputImageType::PixelType      OutputPixelType;
  

 

  
  /** Boundary condition type for the neighborhood iterator */
  typedef ConstantBoundaryCondition< InputImageType > ConstBoundaryConditionType;
  
  /** Neighborhood iterator type */
  typedef NeighborhoodIterator<InputImageType, ConstBoundaryConditionType> NeighborhoodIteratorType;
  
  /** Neighborhood type */
  typedef typename NeighborhoodIteratorType::NeighborhoodType NeighborhoodType;
  
  typedef   ImageRegionConstIterator< InputImageType	   >  ConstIteratorType;

  typedef vnl_matrix< InputPixelType > MatrixType;


  itkSetMacro(ImageSpacing, float);
  itkGetMacro(ImageSpacing, float);
  
  itkSetMacro(AreaThreshold, unsigned int);
  itkGetMacro(AreaThreshold, unsigned int);


protected:
  ShapeWaterShedImageFilter();
  ~ShapeWaterShedImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateData(void);

private:   
	float m_ImageSpacing; 
	unsigned short m_AreaThreshold;
  




}; // end of WaveletNucleiSegmentationFilter class

} //end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkShapeWaterShedImageFilter.txx"
#endif


#endif
