// Copyright (c) General Electric Company, 2017.  All rights reserved.

#ifndef __itkWaveletNucleiSegmentationFilter_h
#define __itkWaveletNucleiSegmentationFilter_h

#include <itkNeighborhoodIterator.h>
#include <itkImageToImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkConstantBoundaryCondition.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <vnl/vnl_matrix.h>
#include <vcl_vector.h>

namespace itk
{
/** \class itkWaveletNucleiSegmentationFilter.h
*/
//#define DEBUG

template <class TInputImage, class TOutputImage, class TLabelImageType = itk::Image<unsigned short,2> >
class WaveletNucleiSegmentationFilter :
    public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef WaveletNucleiSegmentationFilter    Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  itkNewMacro(Self);
  itkTypeMacro( WaveletNucleiSegmentationFilter, ImageToImageFilter );

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
  

  typedef itk::Image<float, 2> FloatImageType; 
  typedef itk::Image< unsigned short, 2 > UShortImageType;
  typedef itk::Image<unsigned char, 2> BinaryImageType;
    typedef itk::Image<unsigned long, 2> ULongImageType; 

  typedef typename BinaryImageType::Pointer  BinaryImagePointer;
  
  /** Boundary condition type for the neighborhood iterator */
  typedef ConstantBoundaryCondition< InputImageType > ConstBoundaryConditionType;
  
  /** Neighborhood iterator type */
  typedef NeighborhoodIterator<InputImageType, ConstBoundaryConditionType> NeighborhoodIteratorType;
  
  /** Neighborhood type */
  typedef typename NeighborhoodIteratorType::NeighborhoodType NeighborhoodType;
  
  typedef   ImageRegionConstIterator< InputImageType	   >  ConstIteratorType;

  typedef vnl_matrix< InputPixelType > MatrixType;


  itkSetMacro(NumberOfLevels, unsigned short);
  itkGetMacro(NumberOfLevels, unsigned short);

  itkSetMacro(SigmaMultiplier, float);
  itkGetMacro(SigmaMultiplier, float);

  itkSetMacro(DetailImagesToUse_X, unsigned short);
  itkGetMacro(DetailImagesToUse_X, unsigned short);

  itkSetMacro(DetailImagesToUse_Y, unsigned short);
  itkGetMacro(DetailImagesToUse_Y, unsigned short);


protected:
  WaveletNucleiSegmentationFilter();
  ~WaveletNucleiSegmentationFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateData(void);

  //itk::Image<unsigned long,2>::Pointer WatershedSegmentation(typename InputImageType::Pointer, typename BinaryImageType::Pointer);
  itk::Image<unsigned long,2>::Pointer WatershedSegmentation(InputImagePointer, BinaryImagePointer);
  
  //typename TOutputImage::Pointer  Segmentation_PostProcessing(typename InputImageType::Pointer, typename BinaryImageType::Pointer );
  typename TOutputImage::Pointer Segmentation_PostProcessing(typename InputImageType::Pointer, typename BinaryImageType::Pointer );




private:   
  WaveletNucleiSegmentationFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  unsigned short	m_NumberOfLevels;
  float				m_SigmaMultiplier;
  unsigned short	m_DetailImagesToUse_X;
  unsigned short	m_DetailImagesToUse_Y;






}; // end of WaveletNucleiSegmentationFilter class

template <class TImageType>
vnl_vector<unsigned int> 
ConnectedComponentArea( typename TImageType::ConstPointer image)
{
	// compute the number of labels
	typedef itk::MinimumMaximumImageCalculator<TImageType>  MinMaxCalculatorType;
	typename MinMaxCalculatorType::Pointer minMaxCalculator = MinMaxCalculatorType::New();
	minMaxCalculator->SetImage(image);
	minMaxCalculator->Compute();
	double numberLabel = minMaxCalculator->GetMaximum();

	vnl_vector<unsigned int> area(numberLabel); 
	vnl_vector< unsigned int >::iterator areaItr;
	typedef itk::ImageRegionConstIterator<TImageType> ConstTImageIterator; 
	ConstTImageIterator imageItr(image, image->GetLargestPossibleRegion());
	for (areaItr = area.begin(); areaItr != area.end(); ++areaItr)
		*areaItr = 0; 
	unsigned int t; 
	for (imageItr.GoToBegin(); !imageItr.IsAtEnd(); ++imageItr)
	{
		if ((t = imageItr.Get())>0)
			area[t-1] ++; 
	}

	return area; 
}

} //end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkWaveletNucleiSegmentationFilter.txx"
#endif


#endif
