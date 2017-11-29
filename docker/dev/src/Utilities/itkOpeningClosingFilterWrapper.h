// Copyright (c) General Electric Company, 2017.  All rights reserved.
/*
 * itkOpeningClosingFilterWrapper.h
 *
 *  Created on: Jan 31, 2013
 *      Author: Alberto Santamaria
 *
 */

#ifndef itkOpeningClosingFilterWrapper_H_
#define itkOpeningClosingFilterWrapper_H_

#include "itkBinaryBallStructuringElement.h"
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"

namespace itk
{
	template < class TInputImage>
	typename TInputImage::Pointer OpeningClosingFilterWrapper( TInputImage * InputImage, int Radius)
	{

		typedef typename TInputImage::PixelType InputPixelType;
		/**
		 *  Removes individual pixels
		 */

		typedef  typename itk::BinaryBallStructuringElement< InputPixelType , 2 > SRType;
		SRType kernel;
		kernel.SetRadius(Radius); kernel.CreateStructuringElement();

	    typedef typename itk::GrayscaleMorphologicalOpeningImageFilter< TInputImage, TInputImage, SRType > ICloseType;
	    typename ICloseType::Pointer ICloseF = ICloseType::New();
		ICloseF->SetInput( InputImage ); ICloseF->SetKernel( kernel );
		ICloseF->Update();


		  try
			{
			  ICloseF->Update();
			}
		  catch( itk::ExceptionObject & excep )
			{
			std::cerr << "Exception caught !" << std::endl;
			std::cerr << excep << std::endl;
			}

		  return ICloseF->GetOutput();

	}

} // end namespace itk

#endif /* ITKOpeningClosingFilterWrapper_H_ */
