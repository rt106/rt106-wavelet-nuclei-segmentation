// Copyright (c) General Electric Company, 2017.  All rights reserved.
/*=========================================================================

 Author: Dirk Padfield
 Date: 03/19/2012

 =========================================================================*/

#ifndef _itkATrousWaveletDecompositionImageFilter_txx
#define _itkATrousWaveletDecompositionImageFilter_txx

#ifndef CONVOLVETYPE
// 0 = ConvolveImageFast, 1 = ConvolveImageITK, 2 = ConvolveImageVNL
#define CONVOLVETYPE 0
#endif

#include "itkATrousWaveletDecompositionImageFilter.h"

// ITK
#include "itkCastImageFilter.h"
#include "itkPeriodicBoundaryCondition.h"
#include "itkTimeProbe.h"
#include "itkMirrorPadImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"

// ITK
#include "itkCreateImageFilterWrapper.h"
#include "itkCopyImageFilterWrapper.h"
#include "itkSubtractImageFilterWrapper.h"
#include "itkImageFileWriterWrapper.h"
#include "ITKImageToVNLMatrixToITKImage.h"
#include "itkRecursiveGaussianImageFilterWrapper.h"

#ifdef __linux__
// Memory usage
#include "sys/types.h"
#include "sys/sysinfo.h"
#endif

// Constructor
namespace itk {
template<class TInputImage, class TOutputImage>
ATrousWaveletDecompositionImageFilter<TInputImage, TOutputImage>
::ATrousWaveletDecompositionImageFilter() 
{
  m_LowerLevel = 3;
  m_UpperLevel = 3;
}

template<class TInputImage, class TOutputImage>
void ATrousWaveletDecompositionImageFilter<TInputImage, TOutputImage>
::GenerateData()
 {
#if DEBUG
  itk::TimeProbe clock;
  itk::TimeProbe::TimeStampType prevTotal = 0;
  std::ofstream timingOutput;
#endif

  InputImageConstPointer inputImage = this->GetInput();
  m_ImageSpacing = inputImage->GetSpacing();

  // Ensure that m_UpperLevel is greater than or equal to m_LowerLevel.
  // Otherwise, the program will crash.
  if( m_UpperLevel < m_LowerLevel )
  {
    m_UpperLevel = m_LowerLevel;
  }

  // We need to work with an image of real type.
  typedef itk::CastImageFilter<InputImageType,InternalImageType> CasterType;
  typename CasterType::Pointer caster = CasterType::New();
  caster->SetInput( inputImage );
  caster->Update();
  InternalImagePointer approximationImage0 = caster->GetOutput();
  InternalImagePointer approximationImage = approximationImage0;
  OutputImagePointer detail_img;

  m_DetailImages.resize(m_UpperLevel - m_LowerLevel + 1);

#if CONVOLVETYPE==0
#if DEBUG
    clock.Start();
#endif
  // In this algorithm, the first approximation image is the approximation image at level = m_LowerLevel, not the inputImage.
    approximationImage0 = this->ConvolveImageFast( inputImage, m_LowerLevel );
    approximationImage = approximationImage0;
#if DEBUG
    clock.Stop();
    std::cout << "Difference time: " << clock.GetTotal() - prevTotal << std::endl;
    timingOutput.open("Timing.csv", std::ios::out | std::ios::app);
    timingOutput << clock.GetTotal() - prevTotal << ",";
    timingOutput.close();
    //std::cout << "Total time ITK Fast: " << clock.GetTotal() << std::endl;
    prevTotal = clock.GetTotal();
#endif
#endif

  // The number of detail levels is one more than the upper level because the levels are numbered from 0, and the 
  // upper level will be used.  For example, if the upper level is 3, we need four detail images.
  unsigned int numberOfLevels = m_UpperLevel + 1;
  for( unsigned int level = 0; level < numberOfLevels; level++ )
  {
#if DEBUG
    clock.Start();
    std::cout << "Level: " << level << std::endl;
#endif

#if WRITEIMAGES
    std::string numberString = static_cast<std::ostringstream*>( &(std::ostringstream() << level) )->str();
    std::cout << "Number string: " << numberString << std::endl;
    std::string approximationImageFileName = "ApproximationImage" + numberString + ".png";
    itk::IntensityRescaleImageFileWriterWrapper<InternalImageType>(approximationImage,approximationImageFileName);
#endif

    // Compute the approximation image
#if CONVOLVETYPE==0
    // In this algorithm, there is no need to do anything when the level is less than the m_LowerLevel.
    // So we can simply continue the loop.
    if( level < m_LowerLevel )
    {
      continue;
    }
    approximationImage = this->ConvolveImageFast( inputImage, level+1 );
#elif CONVOLVETYPE==1
    approximationImage = this->ConvolveImageITK(approximationImage0, level + 1);
#elif CONVOLVETYPE==2
    approximationImage = this->ConvolveImageVNL(approximationImage0, level + 1);
#endif

#if DEBUG
    clock.Stop();
    std::cout << "Difference time: " << clock.GetTotal() - prevTotal << std::endl;
    timingOutput.open("Timing.csv", std::ios::out | std::ios::app);
    timingOutput << clock.GetTotal() - prevTotal << ",";
    timingOutput.close();
    prevTotal = clock.GetTotal();
#endif


    if( level >= m_LowerLevel )
    {
      //compute detail image
      detail_img = SubtractImageFilterWrapper<InternalImageType,InternalImageType,OutputImageType>(approximationImage0,approximationImage);
      // Add to the image vector if this is one of the desired levels.
      m_DetailImages[level-m_LowerLevel] = detail_img;
    }
    //set approximation image for the next level
    approximationImage0 = approximationImage;
  }

  // Cast the approximation image to the desired output type.
  typedef itk::CastImageFilter<InternalImageType,OutputImageType> CastBackType;
  typename CastBackType::Pointer castBack = CastBackType::New();
  castBack->SetInput( approximationImage );
  castBack->Update();

  this->GetOutput(0)->Graft( castBack->GetOutput() );

#if DEBUG
  timingOutput.open("Timing.csv", std::ios::out | std::ios::app);
  timingOutput << std::endl;
  timingOutput.close();
  std::ofstream memoryOutput("Memory.csv", std::ios::out | std::ios::app);
  memoryOutput << std::endl;
  memoryOutput.close();
#endif

 }

template<class TInputImage, class TOutputImage>
typename ATrousWaveletDecompositionImageFilter<TInputImage, TOutputImage>::InternalImagePointer
ATrousWaveletDecompositionImageFilter<TInputImage, TOutputImage>
::ConvolveImageVNL(InternalImageType * inputImage,int level)
 {
  // Setup scaling function.
  //the scaling function is a B3 spline
  const static int scalingFunctionLength = 5;
  const static float scalingFunction[scalingFunctionLength*scalingFunctionLength] = {
      0.00390625, 0.015625, 0.0234375, 0.015625, 0.00390625,
      0.01562500, 0.062500, 0.0937500, 0.062500, 0.01562500,
      0.02343750, 0.093750, 0.1406250, 0.093750, 0.02343750,
      0.01562500, 0.062500, 0.0937500, 0.062500, 0.01562500,
      0.00390625, 0.015625, 0.0234375, 0.015625, 0.00390625
  };

  vnl_matrix<float> scalingFunctionVector;
  scalingFunctionVector.set_size(scalingFunctionLength,scalingFunctionLength);
  scalingFunctionVector.set(scalingFunction);

  MatrixType inputMatrix = ConvertITKImageToVNLMatrix<InternalImageType, MatrixType>( inputImage );

  int nrows = inputMatrix.rows();
  int ncols = inputMatrix.columns();
  MatrixType outputMatrix;
  outputMatrix.set_size(nrows, ncols);

  int step = 1 << (level - 1); // vcl_pow(2.0, level-1.0);

  // Note: the matlab version did an imreconstruct here, but it is not done here.

  //pad image borders (by mirroring) for filtering
  int filtersz = (scalingFunctionLength - 1) * step + 1;
  int border = (int) (filtersz / 2);
  MatrixType img_pad;
  imagePaddingMirror(inputMatrix, img_pad, border);

  itk::TimeProbe clock;
  clock.Start();
  //2-d filtering with filter using approximation equation
  // I_i(x,y) = \sum_{m,n} h(m,n) I_{i-1}(x-2^{i-1}m, y-2^{i-1}n)
  int y, x, m, n, y0, x0;
  //the indexing here may look a bit strange, since we don't center the filter (m=-2 to 2, n=-2 to 2)
  //for filtering. However, because we pad the image border by half of the size of the
  //filter, the center of the filter falls right into the correct indexing (y0,x0) of the output image
  for (y = border * 2, y0 = 0; y < border * 2 + nrows; y++, y0++)
    for (x = border * 2, x0 = 0; x < border * 2 + ncols; x++, x0++) {
      outputMatrix(y0, x0) = 0;
      for (m = 0; m < scalingFunctionLength; m++)
        for (n = 0; n < scalingFunctionLength; n++)
          outputMatrix(y0, x0) += img_pad(y - step * m, x - step * n)
          * scalingFunctionVector[m][n];
    }
  clock.Stop();
  std::cout << "ConvolveImage time: " << clock.GetTotal() << std::endl;

  InternalImagePointer outputImage = ConvertVNLMatrixToITKImage<InternalImageType, MatrixType>( outputMatrix );
  outputImage->CopyInformation( inputImage );
  return outputImage;
 }


template<class TInputImage, class TOutputImage>
typename ATrousWaveletDecompositionImageFilter<TInputImage, TOutputImage>::InternalImagePointer
ATrousWaveletDecompositionImageFilter<TInputImage, TOutputImage>
::ConvolveImageITK(InternalImageType * inputImage,int level)
 {
#if DEBUG
#ifdef __linux__
  const double gigabyte = 1024 * 1024 * 1024;
  struct sysinfo memInfo;
  long long initialFreeMem;
  long long currentFreeMem;
  sysinfo (&memInfo);
  initialFreeMem = memInfo.freeram;
  //std::cout << "Initial free memory: " << initialFreeMem / gigabyte << std::endl;
#endif
#endif

  // Level 0 is the original image.  For processing, the level must be greater than 0.  Otherwise, just return the original image.
  if( level <= 0 )
  {
    return itk::CopyImageFilterWrapper<InternalImageType,InternalImageType>(inputImage);
  }

  // When step is 1, we use a 5*5 neighborhood.  As it increases, we skip more values, but we still use only 25 values.
  unsigned short step = 1 << (level - 1); // vcl_pow(2.0, level-1.0);
  // We will iterate over negative radius values, so elementRadius cannot be unsigned.
  int elementRadius = 2*step;

  // Pad image borders (by mirroring) for filtering
  InternalSizeType padBorder;
  padBorder.Fill( elementRadius );
  typedef itk::MirrorPadImageFilter<InternalImageType,InternalImageType> PadType;
  typename PadType::Pointer padder = PadType::New();
  padder->SetInput( inputImage );
  padder->SetPadUpperBound( padBorder);
  padder->SetPadLowerBound( padBorder );
  padder->Update();

#if DEBUG
#ifdef __linux__
  // For a fair comparison to the fast convolution method, we should only look at the additional memory required by the padder.
  // Everything else in this function is a copy of the padder.
  // Multiple copies are also used by the fast method as part of the RecursiveGaussian filter, but those copies are not tracked either.
  // So the correct comparison is with the memory used by the padder.
  std::cout << "In itkATrousWaveletDecompositionImageFilter: padder size: " << padder->GetOutput()->GetLargestPossibleRegion().GetSize() << std::endl;
  sysinfo (&memInfo);
  currentFreeMem = memInfo.freeram;
  std::cout << "Used memory: " << (initialFreeMem - currentFreeMem)/gigabyte << std::endl;
  std::ofstream memoryOutput("Memory.csv", std::ios::out | std::ios::app);
  memoryOutput << (initialFreeMem - currentFreeMem)/gigabyte << ",";
  memoryOutput.close();
#endif
#endif


  InternalImagePointer convolvedImage = itk::CopyImageFilterWrapper<InternalImageType,InternalImageType>( padder->GetOutput() );
  // Save some memory by pointing the convolvedImage to the padder instead of copying it.
  //InternalImagePointer convolvedImage = padder->GetOutput();
  //convolvedImage->DisconnectPipeline();
  // We cannot also point the tempImage to the padder output because the convolvedImage and tempImage must have separate memory locations since one is an
  // input and one is an output in the processing of the separable filters.
  InternalImagePointer tempImage = itk::CopyImageFilterWrapper<InternalImageType,InternalImageType>( padder->GetOutput() );

  // We cannot use a shaped neighborhood iterator because that iterator creates iterators for all pixels in the neighborhood even if many of them are 
  // not activated.  This becomes very innefficient in our case because we always have the same number of iterators that are turned on even though the
  // neighborhood increases by about double in each dimension for each iterator.
  const typename InternalImageType::OffsetValueType * offsetTable = padder->GetOutput()->GetOffsetTable();

  // Implementation of convolution.
  // Loop over each dimension separately.  
  // This improves computation time by computing in each dimension separately.
  // This is possible since the scaling function is separable.
  
  //itk::TimeProbe clock;
  //clock.Start();

  // Sigma calculated based on the B3-spline approximations to the Gaussian.
  // A B3 spline has a kernel consisting of the following values: 1/16 * [1, 4, 6, 4, 1], 
  // and this sigma best approximates that kernel by considering the Gaussian fit to that kernel.
  // See "The Magic Sigma", CVPR 2011 for more details.
  //float sigma = sqrt(4.0/(2*log(6.0/1.0))); // this is approximately 1.05
  //float sigma = 1/(vcl_sqrt(vcl_log(4.0)));
  // For simplicity and comparison with the fast version, we set sigma = 1
  float sigma = 1;

#if DEBUG
  std::cout << "In itkATrousWaveletDecompositionImageFilter: variance: " << sigma*sigma << std::endl;
#endif

  // The sigma should depend on the image spacing.
  // As a reference, we choose the spacing in x.
  // Sigma should be consistent across dimensions.  However, the distance in physical spacing
  // is calculated in GenerateGaussianKernel so that less weight is given for pixels that are
  // farther away in physical space.
  sigma = sigma * m_ImageSpacing[0];

  //float magicSigma = sqrt(1/(std::log(4.0)));
  unsigned int kernelRadius = 2;

  // The input to the iterator is the padded image.  But the region over which to iterate is the region of the input image before padding.  
  // Thus, the iterator is only over the valid region (not including the padded region).
  InternalSizeType neighborhoodRadius;
  neighborhoodRadius.Fill( 0 );
  // We use a neighborhood iterator only because it gives us access to the GetCenterValue method, which enables faster computation by doing iterator arithmetic
  // with the help of the offset table.
  itk::ConstNeighborhoodIterator<InternalImageType> inputIt;
  // The output iterator also iterates only over the valid region of the image, which is the same as the region of the input image.
  itk::ImageRegionIterator<InternalImageType> outputIt;
  // Set up an iterator for the scaling function.
  std::vector<typename InternalImageType::OffsetValueType> neighborhoodOffsets( kernelRadius*2 + 1 );
  typename std::vector<typename InternalImageType::OffsetValueType>::iterator neighborhoodOffsetsIt;
  vnl_vector<double>::iterator gaussianKernelIt;

  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    // Precompute the offsets in the image given the indices.
    // This should work for ND not just 2D!
    // Since the origin of the padded image is already shifted, we need to shift twice as far to get the offsets.
    int j = -2;
    for( neighborhoodOffsetsIt = neighborhoodOffsets.begin(); neighborhoodOffsetsIt != neighborhoodOffsets.end(); ++neighborhoodOffsetsIt )
    {
      // Only change the indices for the current dimension.
      *neighborhoodOffsetsIt = j * step * offsetTable[i];
      j++;
    }
 
    // The distance in physical spacing
    // is calculated in GenerateGaussianKernel so that less weight is given for pixels that are
    // farther away in physical space.
    vnl_vector<double> gaussianKernel = GenerateGaussianKernel(sigma, kernelRadius,(float)(m_ImageSpacing[i]));  

    // The region over which to iterate changes based on the dimension over which we are iterating.
    // All dimension should have the same region as the padded image except for the dimension over which
    // we are currently iterating, which should have a restricted region corresponding to the input size.
    InternalSizeType iterationSize = convolvedImage->GetLargestPossibleRegion().GetSize();
    iterationSize[i] = inputImage->GetLargestPossibleRegion().GetSize()[i];
    InternalIndexType iterationIndex = convolvedImage->GetLargestPossibleRegion().GetIndex();
    iterationIndex[i] = inputImage->GetLargestPossibleRegion().GetIndex()[i];
    InternalRegionType iterationRegion;
    iterationRegion.SetSize( iterationSize );
    iterationRegion.SetIndex( iterationIndex );

    inputIt = itk::ConstNeighborhoodIterator<InternalImageType>( neighborhoodRadius, convolvedImage, iterationRegion );
    outputIt = itk::ImageRegionIterator<InternalImageType>( tempImage, iterationRegion );

    for (inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++outputIt)
    {
      InternalPixelType * centerValueAddress = inputIt.GetCenterValue();
      float convolutionValue = 0;    
      for( neighborhoodOffsetsIt = neighborhoodOffsets.begin(), gaussianKernelIt = gaussianKernel.begin(); 
        neighborhoodOffsetsIt != neighborhoodOffsets.end(); ++neighborhoodOffsetsIt, ++gaussianKernelIt )
      {
        // Compute the values at the desired locations by finding the value at the offsets from the center value.
        // Multiply these by the corresponding values of the scaling function.
        convolutionValue += *(centerValueAddress+(*neighborhoodOffsetsIt)) * (*gaussianKernelIt);
      }
      outputIt.Set( convolutionValue );
    }

    // Copy the temporary image back to the convolvedImage for the next iteration since we are using separable filters.
    convolvedImage = itk::CopyImageFilterWrapper<InternalImageType,InternalImageType>( tempImage );
  }
  //clock.Stop();
  //std::cout << "ConvolveImageITK time: " << clock.GetTotal() << std::endl;

  // Remove the tempImage since it is no longer needed.
  tempImage = NULL;

  // Copy the processed image into the output image.  The processed image has a larger region than the output image.
  typedef itk::RegionOfInterestImageFilter<InternalImageType,InternalImageType> ROIType;
  typename ROIType::Pointer roiFilter = ROIType::New();
  roiFilter->SetInput( convolvedImage );
  roiFilter->SetRegionOfInterest( inputImage->GetLargestPossibleRegion() );
  roiFilter->Update();

  return roiFilter->GetOutput();
}

template<class TInputImage, class TOutputImage>
typename ATrousWaveletDecompositionImageFilter<TInputImage, TOutputImage>::InternalImagePointer
ATrousWaveletDecompositionImageFilter<TInputImage,TOutputImage>
::ConvolveImageFast(InputImageType const * inputImage, int level)
{
#if DEBUG
#ifdef __linux__
  const double gigabyte = 1024 * 1024 * 1024;
  struct sysinfo memInfo;
  long long initialFreeMem;
  long long currentFreeMem;
  sysinfo (&memInfo);
  initialFreeMem = memInfo.freeram;
  //std::cout << "Initial free memory: " << initialFreeMem / gigabyte << std::endl;
#endif
#endif

  // Level 0 is the original image.  For processing, the level must be greater than 0.  Otherwise, just return the original image.
  if( level <= 0 )
  {
    return itk::CopyImageFilterWrapper<InputImageType,InternalImageType>(inputImage);
    //return static_cast<InternalImageType*>(const_cast<DataObject *>(inputImage));
  }

  // Determine the Gaussian sigma for smoothing.  This is based on the relationship that a smoothing of 
  // G_\sigma = G_\sigma_1 * G_\sigma_2
  // will have sigma values determined by the following relationship
  // \sigma^2 = \sigma_1^2 + \sigma_2^2
  // For each step of the decomposition, the smoothing sigma doubles, which we represent here by the "scale" parameter.
  // We could use the magic sigma for the smoothing.
  // See "The Magic Sigma", CVPR 2011 for more details.
  // We use sigma = 1 for simplicity since there is no need for integer kernels here.
  // We work with variances instead of standard deviations in building up the equivalent standard deviation because
  // that is easier than continually taking square roots and then squares again.
  //float sigma = sqrt(4.0/(2*log(6.0/1.0)));
  //float sigma = sqrt(1/(2*log(6.0/4.0)));
  float sigma = 1;
  float initialVariance = sigma*sigma;
  float variance = initialVariance;
  float scale = 1;
  for( int i = 1; i < level; i++ )
  {
    scale = 2*scale;
    variance = variance + scale*scale*initialVariance;
  }
  float varianceInPixels = variance;

#if DEBUG
  //std::cout << "In itkATrousWaveletDecompositionImageFilter: base sigma: " << sigma << std::endl;
  //std::cout << "In itkATrousWaveletDecompositionImageFilter: variance: " << variance << std::endl;
#endif

  // The sigma should depend on the image spacing.
  // As a reference, we choose the spacing in x.
  // Sigma should be consistent across dimensions.
  // Multiply variance by the square of the image spacing since variance is the square of sigma.
  // We use x as the reference, so we only use index 0.
  variance = variance * m_ImageSpacing[0] * m_ImageSpacing[0];

#if DEBUG
  std::cout << "Level: " << level << std::endl;
  std::cout << "Variance in pixels: " << varianceInPixels << std::endl;
  std::cout << "Variance in physical spacing: " << variance << std::endl;
#endif

  float equivalentSigma = std::sqrt(variance);
  InternalImagePointer convolvedImage= itk::RecursiveGaussianImageFilterWrapper<InputImageType,InternalImageType>(inputImage,equivalentSigma);

#if DEBUG
#ifdef __linux__
  sysinfo (&memInfo);
  currentFreeMem = memInfo.freeram;
  std::cout << "Used memory: " << (initialFreeMem - currentFreeMem)/gigabyte << std::endl;
  std::ofstream memoryOutput("Memory.csv", std::ios::out | std::ios::app);
  memoryOutput << (initialFreeMem - currentFreeMem)/gigabyte << ",";
  memoryOutput.close();
#endif
#endif

  return convolvedImage;
}

template<class TInputImage, class TOutputImage>
vnl_vector<double>
ATrousWaveletDecompositionImageFilter<TInputImage, TOutputImage>
::GenerateGaussianKernel(float sigma, unsigned int kernelRadius, float dimSpacing)
{
  //const double factor = 1.0 / (vcl_sqrt( 2.0 * vnl_math::pi ) * sigma );

  vnl_vector<double> gaussianKernel(kernelRadius*2+1,0);
  vnl_vector<double>::iterator gaussianKernelIt;
  vnl_vector<double> x(kernelRadius*2+1,0);
  vnl_vector<double>::iterator xIt;
  unsigned int counter = 0;
  double sum = 0;
  for( xIt = x.begin(), gaussianKernelIt = gaussianKernel.begin(); xIt != x.end(); ++xIt, ++gaussianKernelIt )
    {
    *xIt = (double)counter - (double)kernelRadius;
    // The weight is computed based on physical distance instead of number of pixels.
    *xIt = *xIt * dimSpacing;
    *gaussianKernelIt = vcl_exp(-0.5 * vnl_math_sqr(*xIt) / vnl_math_sqr(sigma) );
    sum += *gaussianKernelIt;
    counter++;
    }

  // Normalize the kernel to have area 1.  It is better to do it this
  // way rather than using the theoretical factor because of numerical
  // errors.
  for( gaussianKernelIt = gaussianKernel.begin(); gaussianKernelIt != gaussianKernel.end(); ++gaussianKernelIt )
    {
    *gaussianKernelIt /= sum;
    }
  
  return gaussianKernel;
};



/**
 * Standard "PrintSelf" method
 */
template<class TInputImage, class TOutputImage>
void ATrousWaveletDecompositionImageFilter<TInputImage, TOutputImage>::PrintSelf(
    std::ostream& os, Indent indent) const {
  Superclass::PrintSelf(os, indent);
}

} // end namespace itk

#endif
