// Copyright (c) General Electric Company, 2017.  All rights reserved.
#ifndef img_atrous_wavelet_template_h_
#define img_atrous_wavelet_template_h_ 

// Copyright (C) 2009 General Electric Company
// 
// This software is intellectual property of General Electric Co.
// and may not be copied or redistributed without express written consent.

/// \file atrous_wavelet.h
/// \brief class for A Trous wavelet decomposition
/// \author Dashan Gao, gaoda@ge.com
/// \date 04/13/2009
/// 
/// 
/// 
/// \verbatim 
/// Modifications 
/// Apr. 2009 original version of the file 
/// \endverbatim 

// Put includes here.

#include <vnl/vnl_matrix.h>
#include <vcl_vector.h>
#include "ImagePadding.h"

//namespace gevxl {
//namespace img {
//namespace filter {

/********************************
Decomposition: 
approximation image:   	I_i(x,y) = \sum_{m,n} h(m,n) I_{i-1}(x-2^{i-1}m, y-2^{i-1}n)
detail images: 		W_i(x,y) = I_i(x,y)- I_{i+1}(x,y)
			h(m,n): scaling function (use B3 spline) 
			 I_0(x,y) is the original image
Reconstruction
  I_0(x,y)=I_s(x,y)+\sum_{i=1}^s W_i(x,y)

Advantages over standard wavelet:
no need to downsampling and upsampling

Reference: 

*******************************/

namespace atrous_wavelet_const{
  //the scaling function is a B3 spline
const static float scaling_func_[5][5] = { 
      {0.00390625,   0.015625,   0.0234375,   0.015625,   0.00390625},
      {0.01562500,   0.062500,   0.0937500,   0.062500,   0.01562500},
      {0.02343750,   0.093750,   0.1406250,   0.093750,   0.02343750},
      {0.01562500,   0.062500,   0.0937500,   0.062500,   0.01562500},
      {0.00390625,   0.015625,   0.0234375,   0.015625,   0.00390625}  };
const static int scaling_func_size_ = 5;
} //end namespace


template <class DataType>
class atrous_wavelet 
{
public:

  /// Default constructor.
  atrous_wavelet();
  
  atrous_wavelet(int nlevels);

  atrous_wavelet(int nlevels, bool enable_fft);

  //set defaults for the class
  void set_to_default();

  //set the number of levels, the finest level = nlevel-1 
  void set_levels(int nlevels) {nlevels_ = nlevels;}
  int get_levels() const {return nlevels_; }

  //enable/disable fft2 for image filtering 
  void set_enable_fft(bool flag) { enable_fft_ = flag; }
  bool get_enable_fft() const {return enable_fft_; }

  void decompose(const vnl_matrix<DataType> &img, vcl_vector< vnl_matrix<DataType> > &detail_images, 
                               vnl_matrix<DataType> &approx_image);

private: 
  
  int nlevels_; //number of pyramid level (finese level = nlevel-1)
  bool enable_fft_;

  //get the approximate image at some pyramid level
  // the input (level) is the approximate image on the previous level, level 0 is the original image
  void get_approx_image(const vnl_matrix<DataType> &img, int level_ind, vnl_matrix<DataType> &approx_image);
  void get_approx_image_filter2d(const vnl_matrix<DataType> &img, int level_ind, vnl_matrix<DataType> &approx_image);
  //void get_approx_image_fft2(const vnl_matrix<DataType> &img, int level_ind, vnl_matrix<DataType> &approx_image);

  void get_scaling_func(int level_ind, vnl_matrix<DataType> &scaling_func);



 
};
	
template <class DataType>
atrous_wavelet<DataType>
::atrous_wavelet()
{
  set_to_default();
}
  
template <class DataType>
atrous_wavelet<DataType>
::atrous_wavelet(int nlevels)
{
  set_levels(nlevels);
  set_enable_fft(false); //default do not use fft
}


template <class DataType>
atrous_wavelet<DataType>
::atrous_wavelet(int nlevels, bool enable_fft)
{
  set_levels(nlevels);
  set_enable_fft(enable_fft);
}


//set defaults for the class
template <class DataType>
void atrous_wavelet<DataType>
::set_to_default() 
{
  set_levels(3); //default 3 levels are used
  set_enable_fft(false); ////default do not use fft
}


/********************************
Decomposition: 
approximation image:   	I_i(x,y) = \sum_{m,n} h(m,n) I_{i-1}(x-2^{i-1}m, y-2^{i-1}n)
detail images: 		W_i(x,y) = I_i(x,y)- I_{i+1}(x,y)
			h(m,n): scaling function (use B3 spline) 
			 I_0(x,y) is the original image
Return: the detail images at level 0 to nlevels-1, and the approximation image at finest level, i.e. nlevels
***********************************/

template <class DataType>
void atrous_wavelet<DataType>
::decompose(const vnl_matrix<DataType> &img, vcl_vector<vnl_matrix<DataType> > &detail_images, 
                               vnl_matrix<DataType> &approx_image)
{
  int level;
  int nrows, ncols; 
  vnl_matrix<DataType> approx_img0, detail_img; 

  detail_images.clear();

  nrows = img.rows();
  ncols = img.columns();

  approx_image.set_size(nrows, ncols);
  detail_img.set_size(nrows,ncols);

  get_approx_image(img, 0, approx_img0);
/*  {
  itk::Image<double, 2>::Pointer tempImage = itk::Image<double, 2>::New();
  tempImage = ConvertVNLMatrixToITKImage<itk::Image<double, 2>,vnl_matrix<DataType>>( approx_img0 );
  itk::ImageFileWriterWrapper<itk::Image<double, 2>>( tempImage, "D:/temp/approximation0.mhd" );
  }
*/
  for (level=0;level<nlevels_;level++)
  {
    // get the approximation image
    get_approx_image(approx_img0,level+1,approx_image);
/*	{
	itk::Image<double, 2>::Pointer tempImage = itk::Image<double, 2>::New();
	tempImage = ConvertVNLMatrixToITKImage<itk::Image<double, 2>,vnl_matrix<DataType>>( approx_img0 );
	char str[30];
	sprintf(str, "approximation%d.mhd", level+1);
	itk::ImageFileWriterWrapper<itk::Image<double, 2>>( tempImage, str );
	}
  */
	//compute detail image
    detail_img = approx_img0 - approx_image;  
    //add to the image vector
    detail_images.push_back(detail_img);
    //set approximation image for the next level
    approx_img0 = approx_image;

////////////////// test //////////////////////////
/*    //save matrix
    vcl_string fname;
    vcl_ofstream outfile;

    fname = vcl_string("c:/tmp/img_")+ util::to_str(level);
    outfile.open(fname.c_str(), vcl_ios::out);
    outfile << img;
    outfile.close();

    fname = vcl_string("c:/tmp/detail_")+ util::to_str(level);
    outfile.open(fname.c_str(), vcl_ios::out);
    outfile << detail_img;
    outfile.close();
    
    vil_image_view<vxl_byte> tmp_byte;
     fname = vcl_string("c:/tmp/detail_")+ util::to_str(level) + vcl_string(".png");
     matrix_to_image(detail_img, tmp_byte, ncols,nrows);
     vil_save(tmp_byte,fname.c_str());

    fname = vcl_string("c:/tmp/approx_")+ util::to_str(level+1);
    outfile.open(fname.c_str(), vcl_ios::out);
    outfile <<  approx_image;
    outfile.close();

     fname = vcl_string("c:/tmp/approx_")+ util::to_str(level+1) + vcl_string(".png");
      matrix_to_image(approx_image, tmp_byte,  ncols, nrows);
      vil_save(tmp_byte, fname.c_str());
*/
////////////////////////////////////////////////////////
  }
}

template <class DataType>
void atrous_wavelet<DataType>
::get_approx_image(const vnl_matrix<DataType> &img, int level_ind, 
                                      vnl_matrix<DataType> &approx_image)
{
  if (enable_fft_) 
	  // For now, only enable calculation in the spatial domain.
    get_approx_image_filter2d(img, level_ind, approx_image);
    //get_approx_image_fft2(img, level_ind, approx_image);
  else
    get_approx_image_filter2d(img, level_ind, approx_image);
}


template <class DataType>
void atrous_wavelet<DataType>
::get_approx_image_filter2d(const vnl_matrix<DataType> &img, int level_ind, 
                                               vnl_matrix<DataType> &approx_image)
{

  int nrows = img.rows();
  int ncols = img.columns();
  approx_image.set_size(nrows,ncols);

  //level 0 is the original image
  if (level_ind<=0)      // can level_ind be less than 0?
  {
    approx_image = img;
    return;
  }

  int step = 1 << (level_ind-1);  // vcl_pow(2.0, level_ind-1.0);

  // Note: the matlab version did an imreconstruct here, but it is not done here.

  //pad image borders (by mirroring) for filtering
  int filtersz = (atrous_wavelet_const::scaling_func_size_-1)*step+1;
  int border = (int)(filtersz/2);
  vnl_matrix<DataType> img_pad;
  imagePaddingMirror(img, img_pad, border);

  //2-d filtering with filter using approximation equation
  // I_i(x,y) = \sum_{m,n} h(m,n) I_{i-1}(x-2^{i-1}m, y-2^{i-1}n)
  int y,x,m,n,y0,x0;
  //the indexing here may look a bit strang, since we don't center the filter (m=-2 to 2, n=-2 to 2) 
  //for filtering. However, because we pad the image border by half of the size of the 
  //filter, the center of the filter falls right into the correct indexing (y0,x0) of the output image
  for (y=border*2, y0=0; y<border*2+nrows; y++, y0++)
    for (x=border*2, x0=0; x<border*2+ncols; x++, x0++)
    {
      approx_image(y0,x0) = 0;
      for (m=0;m<atrous_wavelet_const::scaling_func_size_;m++)
		  for (n=0;n<atrous_wavelet_const::scaling_func_size_; n++)
		  	  approx_image(y0,x0) += img_pad(y-step*m, x-step*n)*atrous_wavelet_const::scaling_func_[m][n];
    }

}

/*
void atrous_wavelet::get_approx_image_fft2(const vnl_matrix<DataType> &img, int level_ind, 
                                           vnl_matrix<DataType> &approx_image)
{

  int nrows = img.rows();
  int ncols = img.columns();
  approx_image.set_size(nrows,ncols);

  //level 0 is the original image
  if (level_ind==0)
  {
    approx_image = img;
    return;
  }

  //get the scaling function at current level
  vnl_matrix<DataType> scaling_func;
  get_scaling_func(level_ind, scaling_func);

  filter_2d_fft(img, scaling_func, approx_image);
*/
////////////////////////test code ///////////////////
/*  //convert matrix to vil_image_view to try the fft
  vil_image_view<DataType> img_vil, scalfun_vil, approx_img_vil;
  img_vil.set_size(img.columns(),img.rows());
  for (int y=0;y<img_vil.nj();y++)
    for (int x=0;x<img_vil.ni();x++)
      img_vil(x,y) = img(y,x);

  scalfun_vil.set_size(scaling_func.columns(),scaling_func.rows());
  for (int y=0;y<scalfun_vil.nj();y++)
    for (int x=0;x<scalfun_vil.ni();x++)
      scalfun_vil(x,y) = scaling_func(y,x);

  vcl_cout << scaling_func.mean() << vcl_endl;

  approx_img_vil.set_size(approx_image.columns(),approx_image.rows());
  filter_2d_fft(img_vil, scalfun_vil, approx_img_vil);
  //copy result
  for (int y=0;y<approx_img_vil.nj();y++)
    for (int x=0;x<approx_img_vil.ni();x++)
      approx_image(y,x) = approx_img_vil(x,y);
*/
  //}


template <class DataType>
void atrous_wavelet<DataType>
::get_scaling_func(int level_ind, vnl_matrix<DataType> &scaling_func)
{
 
  int step = 1 << (level_ind-1);  // vcl_pow(2.0, level_ind-1.0);
  int filtersz = (atrous_wavelet_const::scaling_func_size_-1)*step+1;
  scaling_func.set_size(filtersz,filtersz);
  scaling_func.fill(0);
  
  for (int m=0;m<atrous_wavelet_const::scaling_func_size_;m++)
    for (int n=0;n<atrous_wavelet_const::scaling_func_size_;n++)
      scaling_func(step*m, step*n) = atrous_wavelet_const::scaling_func_[m][n];
}


#endif
