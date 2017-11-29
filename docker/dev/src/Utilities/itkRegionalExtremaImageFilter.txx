/*=========================================================================
  Copyright (c) General Electric Company, 2017.  All rights reserved.
  Author: Dirk Padfield
  date: 06/02/2006
  =========================================================================*/

#ifndef _itkRegionalExtremaImageFilter_txx
#define _itkRegionalExtremaImageFilter_txx

#include "itkRegionalExtremaImageFilter.h"
#include <time.h>

namespace itk
{

template <class TInputImage, class TOutputImage>
RegionalExtremaImageFilter< TInputImage, TOutputImage>
::RegionalExtremaImageFilter()
{

  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs(1);

  // Take into account image spacing.
  m_Radius.Fill(1);
  m_ForegroundValue = 1;
  m_BackgroundValue = 0;


  m_FindMinima = true;
  m_UseFastMethod = true;

  m_OutputIsBinary = true;

}



template <class TInputImage, class TOutputImage>
void
RegionalExtremaImageFilter< TInputImage, TOutputImage>
::GenerateData()
{
  // Start a timer.
  time_t start,end;
  double diff;
  time (&start);

  // Allocate the output
  this->AllocateOutputs();

  if( !m_UseFastMethod )
    {
    std::cout << "Using the slow method." << std::endl;

    // construct a marker image to manipulate using reconstruction by
    // erosion. the marker image is the input image minus the height
    // parameter.
    typedef ShiftScaleImageFilter<TInputImage, TInputImage>
      ShiftFilterType;
    typename ShiftFilterType::Pointer shift = ShiftFilterType::New();
    shift->SetInput( this->GetInput() );
    shift->SetShift( 1 );

    typedef SubtractImageFilter< TInputImage, TInputImage, TOutputImage > SubtractType;
    typename SubtractType::Pointer subtracter = SubtractType::New();

    if( m_FindMinima )
      {
      std::cout << "Finding the minima. " << std::endl;

      // Delegate to a ReconstructionByErosion filter.
      typename ReconstructionByErosionImageFilter<TInputImage, TInputImage>::Pointer
        erode = ReconstructionByErosionImageFilter<TInputImage, TInputImage>::New();
      
//       // Create a process accumulator for tracking the progress of this minipipeline
//       ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
//       progress->SetMiniPipelineFilter(this);
//       progress->RegisterInternalFilter(erode,1.0f);
      
      // set up the erode filter
      erode->SetMarkerImage( shift->GetOutput() );
      erode->SetMaskImage( this->GetInput() );
//      erode->SetFullyConnected( m_FullyConnected );
      
      erode->Update();
      
      subtracter->SetInput1( erode->GetOutput() );
      subtracter->SetInput2( this->GetInput() );      
      subtracter->Update();

      // graft the output of the erode filter back onto this filter's
      // output. this is needed to get the appropriate regions passed
      // back.
      this->GraftOutput( subtracter->GetOutput() );
      }
    else
      {
      std::cout << "Finding the maxima. " << std::endl;

      // Delegate to a ReconstructionByDilation filter.
      typename ReconstructionByDilationImageFilter<TInputImage, TInputImage>::Pointer
        dilate = ReconstructionByDilationImageFilter<TInputImage, TInputImage>::New();
      
//       // Create a process accumulator for tracking the progress of this minipipeline
//       ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
//       progress->SetMiniPipelineFilter(this);
//       progress->RegisterInternalFilter(erode,1.0f);
      
      // set up the erode filter
      dilate->SetMarkerImage( this->GetInput() );
      dilate->SetMaskImage( shift->GetOutput() );
//      dilate->SetFullyConnected( m_FullyConnected );

      subtracter->SetInput1( shift->GetOutput() );
      subtracter->SetInput2( dilate->GetOutput() );
      subtracter->Update();

      // graft the output of the erode filter back onto this filter's
      // output. this is needed to get the appropriate regions passed
      // back.
      this->GraftOutput( subtracter->GetOutput() );
      }
    }
  else
    {

    InputImageConstPointer inputImage = this->GetInput();
    OutputImagePointer outputImage = this->GetOutput();

    // Set the output image to the ForegroundValue.
    OutputImageRegionType region = outputImage->GetRequestedRegion() ;
    outputImage->SetBufferedRegion( region );
    outputImage->Allocate();
    outputImage->FillBuffer ( m_ForegroundValue );
  
    // Set up a neighborhood iterator for the input image.
    typedef itk::ConstNeighborhoodIterator< InputImageType > ConstNeighborhoodIteratorType;
    ConstNeighborhoodIteratorType neighborhoodIt;

    // Set up a neighborhood iterator for the output image.
    typedef itk::NeighborhoodIterator< OutputImageType > NeighborhoodIteratorType;
    NeighborhoodIteratorType outputNeighborhoodIt;

    // Set up an iterator for the output image.
    typedef itk::ImageRegionIterator< OutputImageType > IteratorType;
    IteratorType outputIt;

    typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< InputImageType > FaceCalculatorType;
    FaceCalculatorType faceCalculator;
    typename FaceCalculatorType::FaceListType faceList;
    typename FaceCalculatorType::FaceListType::iterator fit;

    faceList = faceCalculator( inputImage, 
                               inputImage->GetRequestedRegion(), 
                               m_Radius );


//     // Set up the flood filled iterator.
//     typedef itk::BinaryThresholdImageFunction< InputImageType > FunctionType;
//     FunctionType::Pointer function = FunctionType::New();
//     function->SetInputImage( inputImage );

//     typedef itk::FloodFilledImageFunctionConditionalIterator< OutputImageType, FunctionType> FloodIteratorType;

    // Iterate over the image and check whether the neighbors of each
    // pixel have intensity less than the pixel.
    // Loop through each face
    bool foundExtremum;
    typename ConstNeighborhoodIteratorType::PixelType centerPixel;
    typename ConstNeighborhoodIteratorType::PixelType neighPixel;
    for ( fit=faceList.begin(); fit != faceList.end(); ++fit)
      {
      neighborhoodIt = ConstNeighborhoodIteratorType( m_Radius, inputImage , *fit );
      unsigned int neighborhoodSize = neighborhoodIt.Size();

      outputNeighborhoodIt = NeighborhoodIteratorType( m_Radius, outputImage, *fit );

      outputIt = IteratorType( outputImage, *fit);

      // Loop through each neighborhood
      for (neighborhoodIt.GoToBegin(), outputIt.GoToBegin(); !neighborhoodIt.IsAtEnd(); ++neighborhoodIt, ++outputIt)
        {
        // All of the pixels in the output image are initially set to be
        // the foreground value.  If the current
        // pixel is a background value, then it has been set by the
        // flood filled iterator and doesn't need to be visited again.
        if( outputIt.Get() == m_BackgroundValue )
          {
//        std::cout << "Found a background value" << std::endl;
          continue;
          }


//       if( outputIt.GetIndex()[0] == 1.0 )
//         {
//        std::cout << "Image index: " << outputIt.GetIndex() << std::endl;
//         }

        // Assume that an extremum (min or max depending on the flag)
        // has been found.
        foundExtremum = true;

        centerPixel = neighborhoodIt.GetCenterPixel();

        // Loop through each active element of each neighborhood.
        for (unsigned int indexValue = 0; indexValue < neighborhoodSize; indexValue++)
          {
          neighPixel = neighborhoodIt.GetPixel(indexValue);
          // Find the minima in the image.
          if( m_FindMinima )
            {
            if( neighPixel < centerPixel )
              {
              foundExtremum = false;
              break;
              }
            }
          // Find the maxima in the image.        
          else
            {
            if( neighPixel > centerPixel )
              {
              foundExtremum = false;
              break;
              }          
            }
          }


        // If the pixel we found was not an extremum, set it to zero along
        // will all connected pixels that have the same pixel value.
        if( !foundExtremum )
          {  
          // Store the initial neighborhood location in the image so
          // that we can go back to it.
          typename NeighborhoodIteratorType::IndexType currentImageIndex = neighborhoodIt.GetIndex();

        
          //std::vector< NeighborhoodIteratorType::IndexType > indexStack;
          std::stack < typename NeighborhoodIteratorType::IndexType > indexStack;
          indexStack.push( neighborhoodIt.GetIndex() );
          // centerPixel contains the intensity value we want to flood
          // with.
          outputIt.SetIndex( neighborhoodIt.GetIndex() ); 
          outputIt.Set(  m_BackgroundValue );
        
          while( !indexStack.empty() )
            {
            typename NeighborhoodIteratorType::IndexType neighborhoodIndex = indexStack.top();
            // Popping the vector does not return a value.
            indexStack.pop();
            neighborhoodIt.SetLocation( neighborhoodIndex );
          
            // Loop through each active element of the neighborhood.
            for (unsigned int indexValue = 0; indexValue < neighborhoodSize; indexValue++)
              {
              typename NeighborhoodIteratorType::IndexType imageIndex = neighborhoodIt.GetIndex(indexValue);

              // Check whether the index is out of bounds.
              // There should be a better way to do this.
              bool outOfBounds = false;
              for( unsigned int i = 0; i < ImageDimension; i++ )
                {
                if( imageIndex[i] < 0 || imageIndex[i] >= inputImage->GetBufferedRegion().GetSize()[i] )
                  {
                  outOfBounds = true;
                  }
                }
              if( outOfBounds )
                {
                continue;
                }

              typename IteratorType::PixelType outputPixel = outputIt.Get();
              outputIt.SetIndex( imageIndex );
              neighPixel = neighborhoodIt.GetPixel(indexValue);
              // If the output pixel is not already set not to be an
              // extremum and the input pixel has the same value as the
              // flood filling value.
              if( outputIt.Get() != m_BackgroundValue && neighPixel == centerPixel )
                {
                indexStack.push( imageIndex );
                outputIt.Set( m_BackgroundValue );
                }
              }
            }

          // Go back to the original neighborhood location.
          neighborhoodIt.SetLocation( currentImageIndex);
          outputIt.SetIndex( currentImageIndex );


//         while (stack.getSequenceLength() > 0)
//           {
//           pp = stack.pop();
//           nhSetWalkerLocation(walker, pp);
//           while (nhGetNextInboundsNeighbor(walker, &qq, NULL))
//             {
//             if ((BW[qq] != 0) && (F[qq] == val))
//               {
//               stack.push(qq);
//               BW[qq] = 0;
//               }
//             }
//           }


//See page 769 of 836 of the Insight Guide.



/*        
//         std::cout << "Found extremum at index: " << neighborhoodIt.GetIndex() << std::endl;
//         std::cout << "Intensity at extremum: " << centerPixel << std::endl;

// The function should be set up to retain only pixels with
// the same value as the current center pixel.
function->ThresholdBetween( centerPixel, centerPixel );

FloodIteratorType floodIt = FloodIteratorType( outputImage, function, neighborhoodIt.GetIndex() );  

floodIt.GoToBegin();
while( !floodIt.IsAtEnd())
{
//           std::cout << "Flood filling index: " << floodIt.GetIndex() << std::endl;
floodIt.Set( m_BackgroundValue );
++floodIt;
          }
*/
  }
        }
      
      }
    }

  // Stop the timer
  time (&end);
  diff = difftime (end,start);
  //std::cout << "Time taken: " << diff << std::endl;

  // If instead of a binary output, the user desires the regional
  // maxima to have their original intensities.
  if( !m_OutputIsBinary )
    {
    typedef itk::ImageRegionConstIterator< TInputImage > InputImageIteratorType;
    typedef itk::ImageRegionIterator< TOutputImage > OutputImageIteratorType;

    InputImageIteratorType inputIt( this->GetInput(), this->GetInput()->GetBufferedRegion() );
    OutputImageIteratorType outputIt( this->GetOutput(), this->GetOutput()->GetBufferedRegion() );
    
    for( inputIt.Begin(), outputIt.Begin(); !inputIt.IsAtEnd(); ++inputIt, ++outputIt )
      {
      if( outputIt.Get() )
        {
        outputIt.Set( inputIt.Get() );
        }
      }
    }

}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutputImage>
void
RegionalExtremaImageFilter< TInputImage, TOutputImage>
::PrintSelf(
  std::ostream& os,
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Radius: "
     << static_cast<typename NumericTraits<RadiusType>::PrintType>(m_Radius)
     << std::endl;
  os << indent << "ForegroundValue: "
     << static_cast<typename NumericTraits<OutputImagePixelType>::PrintType>(m_ForegroundValue)
     << std::endl;
  os << indent << "BackgroundValue: "
     << static_cast<typename NumericTraits<OutputImagePixelType>::PrintType>(m_BackgroundValue)
     << std::endl;
}

} // end namespace itk

#endif
