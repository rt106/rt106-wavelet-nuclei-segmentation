// Copyright (c) General Electric Company, 2017.  All rights reserved.
#ifndef _itkATrousWaveletBlobSegmentationImageFilter_txx
#define _itkATrousWaveletBlobSegmentationImageFilter_txx

#include <iostream>

#include "itkATrousWaveletBlobSegmentationImageFilter.h"

// ITK
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkFlatStructuringElement.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkReconstructionImageFilter.h"
#include "itkImageRandomConstIteratorWithIndex.h"

// ITK
#include "ITKImageToVNLMatrixToITKImage.h"
#include "atrous_waveletTemplate.h"
#include "itkImageFileWriterWrapper.h"
#include "itkCreateImageFilterWrapper.h"
#include "itkMultiplyImageFilterWrapper.h"
#include "itkBinaryThresholdImageFilterWrapper.h"
#include "itkMinimumMaximumImageCalculatorWrapper.h"
#include "itkCastImageFilterWrapper.h"

#include <vnl/vnl_double_3x3.h>

#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_real_eigensystem.h>

/*
Author: Dirk Padfield, GE Global Research, padfield@research.ge.com
*/

namespace itk
{

/**
 *    Constructor
 */
template <class TInputImage, class TOutputImage>
ATrousWaveletBlobSegmentationImageFilter<TInputImage, TOutputImage>
::ATrousWaveletBlobSegmentationImageFilter()
  {
  m_SigmaMultiplier = 3;
  m_LowerLevel = 3;
  m_UpperLevel = 3;

  m_SmoothingVariances.push_back(5);
  m_SmoothingVariances.push_back(21);
  m_DecomposeByDifferenceOfGaussians = false;

  this->Superclass::SetNumberOfRequiredInputs( 1 );
  this->Superclass::SetNumberOfRequiredOutputs( 1 );
  this->Superclass::SetNthOutput( 0, OutputImageType::New() );
  }


/**
 *    DifferenceOfGaussians
 */
  template< class TInputImage, class TOutputImage >
  typename ATrousWaveletBlobSegmentationImageFilter<TInputImage, TOutputImage>::DecompositionType::VectorOutputImagePointer ATrousWaveletBlobSegmentationImageFilter<TInputImage, TOutputImage>
  ::DifferenceOfGaussians(const InputImageType * inputImage)
  {
  typename InputImageType::SpacingType imageSpacing = inputImage->GetSpacing();
  typename DecompositionType::VectorOutputImagePointer differenceImages;

  RealImagePointer lowerConvolvedImage, upperConvolvedImage;

  // Loop up to one less than the maximum because we will be offsetting by 1 in the loop.
  for( int i = 0; i < m_SmoothingVariances.size()-1; ++i )
    {
    // The sigma should depend on the image spacing.
    // As a reference, we choose the spacing in x.
    // Sigma should be consistent across dimensions.
    // Multiply variance by the square of the image spacing since variance is the square of sigma.
    // We use x as the reference, so we only use index 0.
    float lowerSigma = std::sqrt(m_SmoothingVariances[i] * imageSpacing[0] * imageSpacing[0]);
    float upperSigma = std::sqrt(m_SmoothingVariances[i+1] * imageSpacing[0] * imageSpacing[0]);

#if DEBUG
    std::cout << "(lowerSigma,upperSigma): " << lowerSigma << ", " << upperSigma << std::endl;
#endif

    // If the variance is set to less than or equal to 0, we simply use the input image instead of smoothing the image.
    // This avoids exceptions in the smoothing filter, and it also speeds up the processing.
    if( i == 0 )
      {
        lowerConvolvedImage = itk::RecursiveGaussianImageFilterWrapper<InputImageType,RealImageType>(inputImage,lowerSigma);
      }
    upperConvolvedImage = itk::RecursiveGaussianImageFilterWrapper<InputImageType,RealImageType>(inputImage,upperSigma);
    RealImagePointer differenceImage = SubtractImageFilterWrapper<RealImageType,RealImageType,RealImageType>(lowerConvolvedImage,upperConvolvedImage);

    lowerConvolvedImage = upperConvolvedImage;
    differenceImages.push_back( differenceImage );
    }

  return differenceImages;
  }
  

/**
 *    GenerateData
 */
template <class TInputImage, class TOutputImage>
void ATrousWaveletBlobSegmentationImageFilter<TInputImage, TOutputImage>
::GenerateData()
  {
  InputImageConstPointer inputImage = this->GetInput();


#if DEBUG
  itk::TimeProbe clock;
  itk::TimeProbe::TimeStampType prevTotal = 0;
  clock.Start();
#endif

  // Calculate the mask image.
  size_t numberOfBackgroundPixels;
  BinaryImagePointer backgroundMask = this->CalculateBackgroundMask( inputImage, numberOfBackgroundPixels );

#if DEBUG
  clock.Stop();
  std::cout << "Compute background time: " << clock.GetTotal() - prevTotal << std::endl;
  prevTotal = clock.GetTotal();
#endif

#if WRITEIMAGES
  itk::IntensityRescaleImageFileWriterWrapper<BinaryImageType>(backgroundMask,"BackgroundMask.png");
#endif

#if DEBUG
  clock.Start();
#endif

  if( m_DecomposeByDifferenceOfGaussians )
    {
    // Decompose the image using difference of Gaussians
    m_DetailImagesThresholded = this->DifferenceOfGaussians(inputImage);
    }
  else
    {
    // Decompose the image using a trous wavelets.
    typename DecompositionType::Pointer waveletDecomposer = DecompositionType::New();
    waveletDecomposer->SetInput( inputImage );
    waveletDecomposer->SetLowerLevel( m_LowerLevel );
    waveletDecomposer->SetUpperLevel( m_UpperLevel );
    waveletDecomposer->Update();
    m_DetailImagesThresholded = waveletDecomposer->GetDetailImages();
    }

#if WRITEIMAGES
  for( unsigned int i = 0; i < m_DetailImagesThresholded.size(); i++ )
    {
    std::string numberString = static_cast<std::ostringstream*>( &(std::ostringstream() << i) )->str();
    std::cout << "Number string: " << numberString << std::endl;
    std::string detailImageFileName = "DetailImage" + numberString + ".tif";
    itk::ImageFileWriterWrapper<RealImageType>(m_DetailImagesThresholded[i],detailImageFileName);
    }
#endif

#if DEBUG
  clock.Stop();
  std::cout << "Wavelet decomposition time: " << clock.GetTotal() - prevTotal << std::endl;
  prevTotal = clock.GetTotal();
#endif

#if DEBUG
  clock.Start();
#endif

  // Denoise each level of the decomposition using an amplitude invariant Bayes estimator.
  itk::ImageRegionIterator<RealImageType> detailImageIt;
  itk::ImageRegionConstIterator<BinaryImageType> backgroundMaskIt(backgroundMask,backgroundMask->GetLargestPossibleRegion());

#if DEBUG
  std::cout << "Number of detail images: " << m_DetailImagesThresholded.size() << std::endl;
#endif
  
  for( unsigned int i = 0; i < m_DetailImagesThresholded.size(); i++ )
  {
    // The creation of the detailImageThresholded must be inside the loop so that a new image is created
    // on each iteration because these images are stored in the m_DetailImagesThresholded.
    // If this new image pointer is only created before the loop, all entries in the m_DetailImagesThresholded will
    // contain the same pointer pointing to the same data.
    detailImageIt = ImageRegionIterator<RealImageType>( m_DetailImagesThresholded[i],m_DetailImagesThresholded[i]->GetLargestPossibleRegion() );

    // Calculate the standard deviation of the noise using the background mask.
    vnl_vector< float > backgroundPixels(numberOfBackgroundPixels);
    unsigned int counter = 0;
    for( detailImageIt.GoToBegin(), backgroundMaskIt.GoToBegin(); !detailImageIt.IsAtEnd(); ++detailImageIt, ++backgroundMaskIt )
    {
      if( backgroundMaskIt.Get() == 1 )
      {
        backgroundPixels[counter] = detailImageIt.Get();
        counter++;
      }
    }

    float N = backgroundPixels.size();
    float meanX = backgroundPixels.mean();
    float rmsX = backgroundPixels.rms(); // RMS = sqrt(1/N * sum(X_i^2));
    float noiseSigma = vcl_sqrt(N/(N-1+vnl_math::eps) * (rmsX*rmsX - meanX*meanX));

    for( detailImageIt.GoToBegin(); !detailImageIt.IsAtEnd(); ++detailImageIt )
    {
      // detailImageThresholded image is overwritten on each iteration.
      detailImageIt.Set( ABE_Thresholder(detailImageIt.Get(), noiseSigma, m_SigmaMultiplier) );

      // Clamp the thresholded value to 0 if it is less than 0.
      if( detailImageIt.Get() < 0 )
      {
        detailImageIt.Set( 0 );
      }
    }

#if WRITEIMAGES
    std::string numberString = static_cast<std::ostringstream*>( &(std::ostringstream() << i) )->str();
    std::cout << "Number string: " << numberString << std::endl;
    std::string detailImageFileName = "DetailImageThresholded" + numberString + ".tif";
    itk::ImageFileWriterWrapper<RealImageType>(m_DetailImagesThresholded[i],detailImageFileName);
#endif
  }

#if DEBUG
  clock.Stop();
  std::cout << "Wavelet denoise time: " << clock.GetTotal() - prevTotal << std::endl;
  prevTotal = clock.GetTotal();
#endif

#if DEBUG
  clock.Start();
#endif

  // Combine the requested levels.
  // Original outputMatrix must be set to one.
  RealImagePointer combinedImage = m_DetailImagesThresholded[0];
  combinedImage->CopyInformation( inputImage );
  // Since the combinedImage is already pointing to the first element of the detailImages vector, we only need to start from index 1 to multiply the rest.
  for( unsigned int i = 1; i < m_DetailImagesThresholded.size(); i++ )
  {
    combinedImage = MultiplyImageFilterWrapper<RealImageType,RealImageType,RealImageType>(combinedImage,m_DetailImagesThresholded[i]);
  } 

#if DEBUG
  clock.Stop();
  std::cout << "Wavelet level combine time: " << clock.GetTotal() - prevTotal << std::endl;
  prevTotal = clock.GetTotal();
#endif

#if WRITEIMAGES
 itk::ImageFileWriterWrapper<RealImageType>(combinedImage,"CombinedImage.mhd");
#endif

#if DEBUG
  clock.Start();
#endif

  // Threshold the image to result in a binary image.
  // All values less than 0 are set to 0, and all values greater than 0 are set to 1.
  // outputImage should already be non-negative.
  typename RealImageType::PixelType minimum, maximum;
  MinimumMaximumImageCalculatorWrapper<RealImageType>(combinedImage,minimum,maximum);
  OutputImagePointer outputImage = BinaryThresholdImageFilterWrapper<RealImageType,OutputImageType>( combinedImage, minimum, 0, 0, 1 );
  combinedImage = NULL;
  // We need to copy the input information so that it can be passed to other filters.
  outputImage->CopyInformation( inputImage );

#if DEBUG
  clock.Stop();
  std::cout << "Wavelet threshold time: " << clock.GetTotal() - prevTotal << std::endl;
  prevTotal = clock.GetTotal();
#endif


#if WRITEIMAGES
  itk::ImageFileWriterWrapper<OutputImageType>(outputImage,"OutputImage.mhd");
#endif

  this->GetOutput(0)->Graft( outputImage );

  return;
 }

template <class TInputImage, class TOutputImage>
float 
ATrousWaveletBlobSegmentationImageFilter<TInputImage, TOutputImage>
::ABE_Thresholder(float inputPixelValue, float noiseSigma, float sigmaMultiplier)
 {
  // Amplitude-scale-invariant Bayes Estimator
  // From "Wavelet-Based Image Estimation: An Empirical Bayes Approach Using
  // Jeffrey's Noninformative Prior."

  // Avoid problems with divide by zero.
  inputPixelValue = inputPixelValue == 0 ? vnl_math::eps : inputPixelValue;

  float numerator;
  //numerator = vcl_pow((float)inputPixelValue,(float)2.0) - sigmaMultiplier * vcl_pow((float)noiseSigma,(float)2.0);
  numerator = inputPixelValue*inputPixelValue - sigmaMultiplier * noiseSigma*noiseSigma;

  // Clip all values lower than 0 to 0.
  numerator = numerator < 0 ? 0 : numerator;

  float outputPixelValue = numerator / inputPixelValue;

  return outputPixelValue;
 }


template <class TInputImage, class TOutputImage>
typename ATrousWaveletBlobSegmentationImageFilter<TInputImage, TOutputImage>::BinaryImagePointer
ATrousWaveletBlobSegmentationImageFilter<TInputImage, TOutputImage>
::CalculateBackgroundMask( const InputImageType * inputImage, size_t & numberOfBackgroundPixels )
 {

#if DEBUG
  itk::TimeProbe clock;
  itk::TimeProbe::TimeStampType prevTotal = 0;
  clock.Start();
#endif

  // Sort the inputImage pixels.
  // Sample the input image rather than selecting all pixels in order to improve speed (and memory) of the processing.
  unsigned int step = 10;
  unsigned int numberOfPixels = inputImage->GetLargestPossibleRegion().GetNumberOfPixels();
  // This will truncate the number of pixels to sample if the result is floating point.  Thus, numberOfPixelsToSample * step
  // will always be <= numberOfPixels.
  unsigned int numberOfPixelsToSample = numberOfPixels / step;

  // Create a map of the pixel value to the number of times it is used.  This creates a histogram that only has bins that are filled.
  typedef typename std::map<InputPixelType,unsigned long> MapType;
  MapType mapper;

  unsigned int counter = 0;
  ImageRegionConstIterator<InputImageType> inputImageIt( inputImage,inputImage->GetLargestPossibleRegion() );
  for( inputImageIt.GoToBegin(); !inputImageIt.IsAtEnd(); ++inputImageIt )
  {
    counter++;
    if( counter == step )
    {
      mapper[inputImageIt.Get()]++;
      counter = 0;
    }
  }

#if DEBUG
  std::cout << "numberOfPixels: " << numberOfPixels << std::endl;
  std::cout << "numberOfPixelsToSample: " << numberOfPixelsToSample << std::endl;
#endif


#if DEBUG
  clock.Stop();
  std::cout << "map pixels time: " << clock.GetTotal() - prevTotal << std::endl;
  prevTotal = clock.GetTotal();
#endif

  // Accumulate the bins of the mapper until the desired number of pixels have been read.
  //unsigned int approximateNumberOfBackgroundPixels = vnl_math_rnd(0.05*numberOfPixelsToSample);
  unsigned int approximateNumberOfBackgroundPixels = vnl_math_rnd(0.05*numberOfPixelsToSample);
  InputPixelType backgroundThreshold = 0;

  unsigned int accumulator = 0;
  typename MapType::iterator mapIt;
  for( mapIt = mapper.begin(); mapIt != mapper.end(); ++mapIt )
  {
    accumulator += mapIt->second;
    if( accumulator > approximateNumberOfBackgroundPixels )
    {
      backgroundThreshold = mapIt->first;
      break;
    }
  }

#if DEBUG
  std::cout << "Number of background pixels: " << approximateNumberOfBackgroundPixels << std::endl;
  std::cout << "Background threshold: " << backgroundThreshold << std::endl;
#endif



#if DEBUG
  clock.Start();
#endif

  // Create an output mask.
  BinaryImagePointer maskImage = CreateImageFilterWrapper<BinaryImageType>(inputImage->GetLargestPossibleRegion().GetSize());
  ImageRegionIterator<BinaryImageType> maskImageIt( maskImage,maskImage->GetLargestPossibleRegion() );

  // This value may be slightly different from approximateNumberOfBackgroundPixels because of rounding errors in calculating approximateNumberOfBackgroundPixels.
  numberOfBackgroundPixels = 0;

  for( inputImageIt.GoToBegin(), maskImageIt.GoToBegin(); !inputImageIt.IsAtEnd(); ++inputImageIt, ++maskImageIt )
  {
    if( inputImageIt.Get() <= backgroundThreshold )
    {
      maskImageIt.Set( 1 );
      ++numberOfBackgroundPixels;
    }
    else
    {
      maskImageIt.Set( 0 );
    }
  }

#if DEBUG
  clock.Stop();
  std::cout << "Mask time: " << clock.GetTotal() - prevTotal << std::endl;
  prevTotal = clock.GetTotal();
#endif

  return maskImage;
 }

template <class TInputImage, class TOutputImage>
void 
ATrousWaveletBlobSegmentationImageFilter<TInputImage,TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
 {
  Superclass::PrintSelf(os,indent);

  os << " FIXME - needs to be implemented  " << std::endl;
 }

} // end namespace itk


#endif



