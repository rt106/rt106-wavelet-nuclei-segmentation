// Copyright (c) General Electric Company, 2017.  All rights reserved.
#ifndef img_image_padding_h_
#define img_image_padding_h_ 

// Copyright (C) 2009 General Electric Company
// 
// This software is intellectual property of General Electric Co.
// and may not be copied or redistributed without express written consent.

/// \file image_padding.h
/// \brief Pad image boundaries for filtering and other operations
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

//#include <vil/vil_image_view.h>
//#include <vnl/vnl_matrix.h>

// namespace gevxl {
// namespace img {
// namespace filter {

	
/*** functions to pad the boundaries of an image (or a matrix) for filtering (or convolution) 
**/

//declarations

/*
//pad the image by zeros
template<class T>
void image_padding_zero(const vil_image_view<T> &img_in, 
                        vil_image_view<T> &img_out, int border_width);

//pad the image by copying the border pixels
template<class T>
void image_padding_copy(const vil_image_view<T> &img_in, 
                        vil_image_view<T> &img_out, int border_width);
                        
//pad the image by mirroring the border pixels
template<class T>
void image_padding_mirror(const vil_image_view<T> &img_in, 
                        vil_image_view<T> &img_out, int border_width);

//pad the image by mirroring the border pixels
template<class T, class destT >
void image_padding_mirror(const vil_image_view<T> &img_in, 
                        vil_image_view<destT > &img_out, int border_width);
*/

//pad an image matrix by zeros
template<class T>
void imagePaddingZero(const vnl_matrix<T> &img_in, 
                        vnl_matrix<T> &img_out, int border_width); 

//pad an image matrix by copying the border pixels
template<class T>
void imagePaddingCopy(const vnl_matrix<T> &img_in, 
                        vnl_matrix<T> &img_out, int border_width);

//pad an image matrix by mirroring the border pixels
template<class T>
void imagePaddingMirror(const vnl_matrix<T> &img_in, 
                        vnl_matrix<T> &img_out, int border_width);

template<class T, class destT >
void imagePaddingMirror(const vnl_matrix<T> &img_in, 
                        vnl_matrix<destT > &img_out, int border_width);


////////////////////////////////////////////////////////////////

/*** functions to pad the boundaries of an image (or a matrix) for filtering (or convolution) */

//definitions:

/*
//pad the image by zeros
template<class T>
void image_padding_zero(const vil_image_view<T> &img_in, 
                        vil_image_view<T> &img_out, int border_width)
{
  int nx,ny, nx_new, ny_new;
  nx = img_in.ni();
  ny = img_in.nj();
  nx_new = nx + border_width*2;
  ny_new = ny + border_width*2;

  //adjust output image size only when it's smaller than padded size
  if (img_out.ni()<nx_new || img_out.nj()<ny_new)
    img_out.set_size(nx_new, ny_new);
  img_out.fill(0);

  unsigned x,x_new, y, y_new;
  for (y=0,y_new=border_width; y<ny; y++, y_new++)
    for (x=0,x_new=border_width; x<nx; x++, x_new++)
      img_out(x_new, y_new) = img_in(x,y);
  
}


//pad the image by mirroring the border pixels
template<class T>
void image_padding_mirror(const vil_image_view<T> &img_in, 
                        vil_image_view<T> &img_out, int border_width)
{
  int nx,ny, nx_new, ny_new;
  int h, w, r1, r2, r3, r4;

  nx = img_in.ni();
  ny = img_in.nj();
  nx_new = nx + border_width*2;
  ny_new = ny + border_width*2;

  //adjust output image size only when it's smaller than padded size
  if (img_out.ni()<nx_new || img_out.nj()<ny_new)
    img_out.set_size(nx_new, ny_new);
  else
    img_out.fill(0);

  //padding the boundaries
  //top rows
  r1 = nx+border_width; r2 = nx+2*border_width;
	for (h=0;h<border_width;h++) 
  {
	    for (w=0;w<border_width;w++) 
        img_out(w,h) = img_in(border_width-1-w, border_width-1-h);
	    for (w=border_width;w<r1;w++) 
        img_out(w,h) = img_in(w-border_width, border_width-1-h); 
	    for (w=r1;w<r2;w++) 
        img_out(w,h) = img_in(2*nx-w+border_width-1, border_width-1-h);
	}

  //middle rows
 	r1 = ny+border_width; r2 = nx+border_width; r3 = nx+2*border_width;
	for (h=border_width;h<r1;h++) 
  {
	    for (w=0;w<border_width;w++) 
        img_out(w,h) = img_in(border_width-1-w, h-border_width);        
	    for (w=border_width;w<r2;w++) 
        img_out(w,h) = img_in(w-border_width, h-border_width);        
	    for (w=r2;w<r3;w++) 
        img_out(w,h) = img_in(2*nx-w+border_width-1, h-border_width);
	}

  //bottom rows
	r1 = ny+border_width; r2 = ny+2*border_width;
	r3 = nx+border_width; r4 = nx+2*border_width;
	for (h=r1;h<r2;h++) 
  {
    for (w=0;w<border_width;w++)
      img_out(w,h) = img_in(border_width-1-w, 2*ny-h+border_width-1);     
		for (w=border_width;w<r3;w++)
      img_out(w,h) = img_in(w-border_width, 2*ny-h+border_width-1);    
		for (w=r3;w<r4;w++)
      img_out(w,h) = img_in(2*nx-w+border_width-1, 2*ny-h+border_width-1);    	
  }

}



//pad the image by mirroring the border pixels
template<class T, class destT >
void image_padding_mirror(const vil_image_view<T> &img_in, 
                        vil_image_view< destT > &img_out, int border_width)
{
  int nx,ny, nx_new, ny_new;
  int h, w, r1, r2, r3, r4;

  nx = img_in.ni();
  ny = img_in.nj();
  nx_new = nx + border_width*2;
  ny_new = ny + border_width*2;

  //adjust output image size only when it's smaller than padded size
  if (img_out.ni()<nx_new || img_out.nj()<ny_new)
    img_out.set_size(nx_new, ny_new);
  else
    img_out.fill(0);

  //padding the boundaries
  //top rows
  r1 = nx+border_width; r2 = nx+2*border_width;
	for (h=0;h<border_width;h++) 
  {
	    for (w=0;w<border_width;w++) 
        img_out(w,h) = img_in(border_width-1-w, border_width-1-h);
	    for (w=border_width;w<r1;w++) 
        img_out(w,h) = img_in(w-border_width, border_width-1-h); 
	    for (w=r1;w<r2;w++) 
        img_out(w,h) = img_in(2*nx-w+border_width-1, border_width-1-h);
	}

  //middle rows
 	r1 = ny+border_width; r2 = nx+border_width; r3 = nx+2*border_width;
	for (h=border_width;h<r1;h++) 
  {
	    for (w=0;w<border_width;w++) 
        img_out(w,h) = img_in(border_width-1-w, h-border_width);        
	    for (w=border_width;w<r2;w++) 
        img_out(w,h) = img_in(w-border_width, h-border_width);        
	    for (w=r2;w<r3;w++) 
        img_out(w,h) = img_in(2*nx-w+border_width-1, h-border_width);
	}

  //bottom rows
	r1 = ny+border_width; r2 = ny+2*border_width;
	r3 = nx+border_width; r4 = nx+2*border_width;
	for (h=r1;h<r2;h++) 
  {
    for (w=0;w<border_width;w++)
      img_out(w,h) = img_in(border_width-1-w, 2*ny-h+border_width-1);     
		for (w=border_width;w<r3;w++)
      img_out(w,h) = img_in(w-border_width, 2*ny-h+border_width-1);    
		for (w=r3;w<r4;w++)
      img_out(w,h) = img_in(2*nx-w+border_width-1, 2*ny-h+border_width-1);    	
  }

}



//pad the image by copying the border pixels
template<class T>
void image_padding_copy(const vil_image_view<T> &img_in, 
                        vil_image_view<T> &img_out, int border_width)
{
  int nx,ny, nx_new, ny_new;
  int h, w, r1, r2, r3, r4;

  nx = img_in.ni();
  ny = img_in.nj();
  nx_new = nx + border_width*2;
  ny_new = ny + border_width*2;

  //adjust output image size only when it's smaller than padded size
  if (img_out.ni()<nx_new || img_out.nj()<ny_new)
    img_out.set_size(nx_new, ny_new);
  else
    img_out.fill(0);

  //padding the boundaries by copying
  //top rows
  r1 = nx+border_width; r2 = nx+2*border_width;
	for (h=0;h<border_width;h++) 
  {
	    for (w=0;w<border_width;w++) 
        img_out(w,h) = img_in(0,0);
	    for (w=border_width;w<r1;w++) 
        img_out(w,h) = img_in(w-border_width, 0); 
	    for (w=r1;w<r2;w++) 
        img_out(w,h) = img_in(nx-1, 0);
	}

  //middle rows
 	r1 = ny+border_width; r2 = nx+border_width; r3 = nx+2*border_width;
	for (h=border_width;h<r1;h++) 
  {
	    for (w=0;w<border_width;w++) 
        img_out(w,h) = img_in(0, h-border_width);        
	    for (w=border_width;w<r2;w++) 
        img_out(w,h) = img_in(w-border_width, h-border_width);        
	    for (w=r2;w<r3;w++) 
        img_out(w,h) = img_in(nx-1, h-border_width);
	}

  //bottom rows
	r1 = ny+border_width; r2 = ny+2*border_width;
	r3 = nx+border_width; r4 = nx+2*border_width;
	for (h=r1;h<r2;h++) 
  {
    for (w=0;w<border_width;w++)
      img_out(w,h) = img_in(0, ny-1);     
		for (w=border_width;w<r3;w++)
      img_out(w,h) = img_in(w-border_width, ny-1);    
		for (w=r3;w<r4;w++)
      img_out(w,h) = img_in(nx-1, ny-1);    	
  }

}
*/


//pad an image matrix by zeros
template<class T>
void imagePaddingZero(const vnl_matrix<T> &img_in, 
                        vnl_matrix<T> &img_out, int border_width)
{
  int nx,ny, nx_new, ny_new;
  nx = img_in.columns();
  ny = img_in.rows();
  nx_new = nx + border_width*2;
  ny_new = ny + border_width*2;

  //adjust output image size only when it's smaller than padded size
  if (img_out.columns()<nx_new || img_out.rows()<ny_new)
    img_out.set_size(ny_new, nx_new);
  img_out.fill(0);


  img_out.set_size(ny_new, nx_new);
  img_out.fill(0);

  unsigned x,x_new, y, y_new;
  for (y=0,y_new=border_width; y<ny; y++, y_new++)
    for (x=0,x_new=border_width; x<nx; x++, x_new++)
      img_out(y_new, x_new) = img_in(y,x);  
}


//pad an image matrix by mirroring the border pixels
template<class T>
void imagePaddingMirror(const vnl_matrix<T> &img_in, 
                        vnl_matrix<T> &img_out, int border_width)
{
  int nx,ny, nx_new, ny_new;
  int h, w, r1, r2, r3, r4;

  nx = img_in.columns();
  ny = img_in.rows();
  nx_new = nx + border_width*2;
  ny_new = ny + border_width*2;

  //adjust output image size only when it's smaller than padded size
  if (img_out.columns()<(unsigned)nx_new || img_out.rows()<(unsigned)ny_new)
    img_out.set_size(ny_new, nx_new);
  else
    img_out.fill(0);

  //padding the boundaries
  //top rows
  r1 = nx+border_width; r2 = nx+2*border_width;
	for (h=0;h<border_width;h++) 
  {
	    for (w=0;w<border_width;w++) 
        img_out(h,w) = img_in(border_width-1-h, border_width-1-w);
	    for (w=border_width;w<r1;w++) 
        img_out(h,w) = img_in(border_width-1-h, w-border_width); 
	    for (w=r1;w<r2;w++) 
        img_out(h,w) = img_in(border_width-1-h, 2*nx-w+border_width-1);
	}

  //middle rows
 	r1 = ny+border_width; r2 = nx+border_width; r3 = nx+2*border_width;
	for (h=border_width;h<r1;h++) 
  {
	    for (w=0;w<border_width;w++) 
        img_out(h,w) = img_in(h-border_width, border_width-1-w);        
	    for (w=border_width;w<r2;w++) 
        img_out(h,w) = img_in(h-border_width, w-border_width);        
	    for (w=r2;w<r3;w++) 
        img_out(h,w) = img_in(h-border_width, 2*nx-w+border_width-1);
	}

  //bottom rows
	r1 = ny+border_width; r2 = ny+2*border_width;
	r3 = nx+border_width; r4 = nx+2*border_width;
	for (h=r1;h<r2;h++) 
  {
    for (w=0;w<border_width;w++)
      img_out(h,w) = img_in(2*ny-h+border_width-1, border_width-1-w);     
		for (w=border_width;w<r3;w++)
      img_out(h,w) = img_in(2*ny-h+border_width-1, w-border_width);    
		for (w=r3;w<r4;w++)
      img_out(h,w) = img_in(2*ny-h+border_width-1, 2*nx-w+border_width-1);    	
  }

}



//pad an image matrix by mirroring the border pixels
template<class T, class destT >
void imagePaddingMirror(const vnl_matrix<T> &img_in, 
                        vnl_matrix< destT > &img_out, int border_width)
{
  int nx,ny, nx_new, ny_new;
  int h, w, r1, r2, r3, r4;

  nx = img_in.columns();
  ny = img_in.rows();
  nx_new = nx + border_width*2;
  ny_new = ny + border_width*2;

  //adjust output image size only when it's smaller than padded size
  if (img_out.columns()<(unsigned)nx_new || img_out.rows()<(unsigned)ny_new)
    img_out.set_size(ny_new, nx_new);
  else
    img_out.fill(0);

  //padding the boundaries
  //top rows
  r1 = nx+border_width; r2 = nx+2*border_width;
	for (h=0;h<border_width;h++) 
  {
	    for (w=0;w<border_width;w++) 
        img_out(h,w) = img_in(border_width-1-h, border_width-1-w);
	    for (w=border_width;w<r1;w++) 
        img_out(h,w) = img_in(border_width-1-h, w-border_width); 
	    for (w=r1;w<r2;w++) 
        img_out(h,w) = img_in(border_width-1-h, 2*nx-w+border_width-1);
	}

  //middle rows
 	r1 = ny+border_width; r2 = nx+border_width; r3 = nx+2*border_width;
	for (h=border_width;h<r1;h++) 
  {
	    for (w=0;w<border_width;w++) 
        img_out(h,w) = img_in(h-border_width, border_width-1-w);        
	    for (w=border_width;w<r2;w++) 
        img_out(h,w) = img_in(h-border_width, w-border_width);        
	    for (w=r2;w<r3;w++) 
        img_out(h,w) = img_in(h-border_width, 2*nx-w+border_width-1);
	}

  //bottom rows
	r1 = ny+border_width; r2 = ny+2*border_width;
	r3 = nx+border_width; r4 = nx+2*border_width;
	for (h=r1;h<r2;h++) 
  {
    for (w=0;w<border_width;w++)
      img_out(h,w) = img_in(2*ny-h+border_width-1, border_width-1-w);     
		for (w=border_width;w<r3;w++)
      img_out(h,w) = img_in(2*ny-h+border_width-1, w-border_width);    
		for (w=r3;w<r4;w++)
      img_out(h,w) = img_in(2*ny-h+border_width-1, 2*nx-w+border_width-1);    	
  }

}


//pad an image matrix by copying the border pixels
template<class T>
void imagePaddingCopy(const vnl_matrix<T> &img_in, 
                        vnl_matrix<T> &img_out, int border_width)
{
  int nx,ny, nx_new, ny_new;
  int h, w, r1, r2, r3, r4;

  nx = img_in.columns();
  ny = img_in.rows();
  nx_new = nx + border_width*2;
  ny_new = ny + border_width*2;

  //adjust output image size only when it's smaller than padded size
  if (img_out.columns()<nx_new || img_out.rows()<ny_new)
    img_out.set_size(ny_new, nx_new);
  else
    img_out.fill(0);

  //padding the boundaries by copying
  //top rows
  r1 = nx+border_width; r2 = nx+2*border_width;
	for (h=0;h<border_width;h++) 
  {
	    for (w=0;w<border_width;w++) 
        img_out(h,w) = img_in(0,0);
	    for (w=border_width;w<r1;w++) 
        img_out(h,w) = img_in(0, w-border_width); 
	    for (w=r1;w<r2;w++) 
        img_out(h,w) = img_in(0, nx-1);
	}

  //middle rows
 	r1 = ny+border_width; r2 = nx+border_width; r3 = nx+2*border_width;
	for (h=border_width;h<r1;h++) 
  {
	    for (w=0;w<border_width;w++) 
        img_out(h,w) = img_in(h-border_width,0);        
	    for (w=border_width;w<r2;w++) 
        img_out(h,w) = img_in(h-border_width, w-border_width);        
	    for (w=r2;w<r3;w++) 
        img_out(h,w) = img_in(h-border_width, nx-1);
	}

  //bottom rows
	r1 = ny+border_width; r2 = ny+2*border_width;
	r3 = nx+border_width; r4 = nx+2*border_width;
	for (h=r1;h<r2;h++) 
  {
    for (w=0;w<border_width;w++)
      img_out(h,w) = img_in(ny-1, 0);     
		for (w=border_width;w<r3;w++)
      img_out(h,w) = img_in(ny-1, w-border_width);    
		for (w=r3;w<r4;w++)
      img_out(h,w) = img_in(ny-1, nx-1);    	
  }

}



	
//  }}} //end of namespaces

#endif

