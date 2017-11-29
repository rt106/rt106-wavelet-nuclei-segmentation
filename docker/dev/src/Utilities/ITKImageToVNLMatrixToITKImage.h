// Copyright (c) General Electric Company, 2017.  All rights reserved.
#ifndef _ITKImageToVNLMatrixToITKImage_h
#define _ITKImageToVNLMatrixToITKImage_h

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

template < class TGenericImage, class MatrixType >
MatrixType ConvertITKImageToVNLMatrix( const TGenericImage * image )
{
  typedef TGenericImage ImageType;
  //typedef vnl_matrix< double > MatrixType;
  //typedef vnl_matrix< double >::iterator MatrixIteratorType;
  typedef typename MatrixType::iterator MatrixIteratorType;

  typename ImageType::RegionType region;
  region = image->GetBufferedRegion();
  typename ImageType::SizeType size = region.GetSize();

  MatrixType outputMatrix(size[1],size[0]);
  MatrixIteratorType outputMatrixIterator;

  typedef itk::ImageRegionConstIterator< ImageType> IteratorType;
  IteratorType it( image, region);
  it.GoToBegin();
  outputMatrixIterator = outputMatrix.begin();

  // Elements in VNL matrices are stored in row-major order, just like ITK images.
  // So we can simply iterate through the image and the matrix and copy the values directly.
  while( !it.IsAtEnd() )
    {
    *outputMatrixIterator = it.Get();
    ++it;
    ++outputMatrixIterator;
    }

  return outputMatrix;
}

template < class TGenericImage, class MatrixType >
typename TGenericImage::Pointer ConvertVNLMatrixToITKImage( MatrixType inputMatrix )
{
  typedef TGenericImage ImageType;
  //typedef vnl_matrix< double >::iterator MatrixIteratorType;
  typedef typename MatrixType::iterator MatrixIteratorType;

  MatrixIteratorType inputMatrixIterator;

  typename ImageType::Pointer image = ImageType::New();

  typename ImageType::SizeType size;
  size[0] = inputMatrix.columns();
  size[1] = inputMatrix.rows();

  typename ImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;

  typename ImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );

  double spacing[2];
  spacing[0] = 1;
  spacing[1] = 1;

  double origin[2];
  origin[0] = 0;
  origin[1] = 0;

  image->SetRegions( region);
  image->Allocate();
  image->SetSpacing( spacing );
  image->SetOrigin( origin );

  typedef itk::ImageRegionIterator< ImageType> IteratorType;
  IteratorType it( image, region);
  it.GoToBegin();
  inputMatrixIterator = inputMatrix.begin();

  // Elements in VNL matrices are stored in row-major order, just like ITK images.
  // So we can simply iterate through the image and the matrix and copy the values directly.
  while( !it.IsAtEnd() )
    {
    it.Set( *inputMatrixIterator );
    ++it;
    ++inputMatrixIterator;
    }

  return image;
}

#endif
