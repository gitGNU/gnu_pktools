/**********************************************************************
pkfilter_bin.cc: program to filter raster images: median, min/max, morphological, filtering
Copyright (C) 2008-2016 Pieter Kempeneers

This file is part of pktools

pktools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pktools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with pktools.  If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/
#include <iostream>
#include <string>
#include "base/Optionpk.h"
#include "imageclasses/ImgRaster.h"
#include "AppFactory.h"

/******************************************************************************/
/*! \page pkfilter pkfilter
  program to filter raster images
  ## SYNOPSIS

  <code>
  Usage: pkfilter -i input -o output [-f filter | -perc value | -srf file [-srf file]* -win wavelength [-win wavelength]* | -wout wavelength -fwhm value [-wout wavelength -fwhm value]* -win wavelength [-win wavelength]*]
  </code>

  <code>
  Options: [-dx value [-dy value] | -dz value] [-nodata value]

  Advanced options: check table
  </code>

  \section pkfilter_description Description

  This utility implements spatial and spectral filtering for raster data. In the spatial domain (X, Y), the filter typically involves a rectangular convolution kernel (moving window). To avoid image shifting, the size of the window should be odd (3, 5, 7, ...). You can set the window sizes in X and Y directions separately with the options -dx and -dy. A circular kernel (disc) is applied if option -circ is set. An overview of the supported filters (option -f|--filter) is given below. You can create customized filters by defining your own filter taps (multiplicative elements of the filter kernel) via an ascii file (option -tap). In the spectral/temporal domain (Z) you can filter multi-band raster inputs. The kernel filter size can be set with the option -dz (use odd values only).

  \anchor pkfilter_functions

  \subsection pkfilter_functions_1 Filters in spatial (dx, dy) and spectral/temporal (dz) domain

  \subsubsection pkfilter_functions_1_1 Implemented as moving window: choose dx, dy or dz > 1 and odd (3, 5, 7, etc.)

  The number of output bands equals number of input bands

  |filter|description|
  |-------------|-------------|
  |dilate|morphological dilation|
  |erode|morphological erosion|
  |close|morpholigical closing (dilate+erode)|
  |open|morpholigical opening (erode+dilate)|
  |smoothnodata values|smooth nodata values (set nodata option!)|

  Example: "Smooth" (interpolate) nodata in spectral/temporal domain (-dz 1), using a linear interpolation
  \code
  pkfilter -i input.tif -o smoothed.tif -dz 1 -f smoothnodata -interp linear
  \endcode

  Example: Filter input.tif in spatial domain with morphological dilation filter with kernel size 3x3.

  \code
  pkfilter -i input.tif -o dilated.tif -dx 3 -dy 3 -f dilate
  \endcode

  \subsubsection pkfilter_functions_1_2 Implemented as either moving window or statistical function in spectral/temporal domain (choose dz=1). 

  In case of moving window, the number of output bands equals number of input bands. In case dz=1, the single output band is calculated as the result of the statistical function applied to all bands.

  |filter | description|
  |-------------|-------------|
  nvalid | report number of valid (not nodata) values in window
  median | perform a median filter in spatial (dx, dy) or spectral/temporal (dz) domain
  var | calculate variance in window
  min | calculate minimum in window
  max | calculate maximum in window
  sum | calculate sum in window
  mean | calculate mean in window
  stdev | calculate standard deviation in window
  savgolay | Savitzky-Golay filter (check examples page!)
  percentile | calculate percentile value in window
  proportion | calculate proportion in windoww

  Example: Median filter in spatial domain

  \code
  pkfilter -i input.tif -o median.tif -dx 3 -dy 3 -f median
  \endcode

  Example: Calculate statistical variance in spectral/temporal domain (single output band)

  \code
  pkfilter -i input.tif -o var.tif -dz 1 -f var
  \endcode

  \subsection pkfilter_functions_2 Wavelet filters 

  \subsubsection pkfilter_functions_2_1 Wavelet filter in in spatial or spectral/temporal (set dz = 1) domain. 

  The number of output bands equals number of input bands

  |filter | description|
  |-------------|-------------|
  dwt | discrete wavelet transform
  dwti | discrete inverse wavelet transform
  dwt_cut | discrete wavelet + inverse transform, using threshold option to cut percentile of coefficients 

  Example: Calculate discrete wavelet in spatial domain

  \code
  pkfilter -i lena.tif -o lena_dwt.tif -f dwt
  \endcode

  Example: Calculate discrete wavelet in spectral/temporal domain

  \code
  pkfilter -i timeseries.tif -o dwt.tif -f dwt -dz 1
  \endcode

  \subsubsection pkfilter_functions_2_2 Wavelet filter implemented in spectral/temporal domain only. 

  The number of output bands equals number of input bands

  |filter | description|
  |-------------|-------------|
  dwt_cut_from | discrete wavelet + inverse transform, setting all high frequence coefficients to zero (scale >= threshold)


  Example: Calculate low frequency time series based on discrete wavelet + inverse transform in spectral/temporal domain, retaining only coefficients until scale 3

  \code
  pkfilter -i timeseries.tif -o lowfrequency.tif -f dwt_cut_from -dz 1 -t 4
  \endcode

  \subsection pkfilter_functions_3 Filters in spatial domain only (dx, dy > 1 and odd). 

  The number of output bands equals number of input bands. 

  |filter | description|
  |-------------|-------------|
  mrf | Markov random field
  ismin | pixel is minimum?
  ismax | pixel is maximum?
  shift | perform a pixel shift in spatial window
  scramble | scramble pixels in a spatial window
  mode (majority voting) | perform a majority voring (set class option)
  sobelx | horizontal edge detection
  sobely | vertical edge detection 
  sobelxy | diagonal edge detection (NE-SW)
  sobelyx | diagonal edge detection (NW-SE)
  countid | count digital numbers in window
  order | rank pixels in order
  density | calculated the density
  homog | central pixel must be identical to all other pixels within window
  heterog | central pixel must be different than all other pixels within window
  sauvola | Sauvola's thresholding method

  Example: Sobel edge detection in horizontal direction

  \code
  pkfilter -i lena.tif -o sobelx.tif -f sobelx -dx 5 -dy 5
  \endcode

  \section pkfilter_options Options
  - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
  - short option `-h` shows basic options only, long option `--help` shows all options
  |short|long|type|default|description|
  |-----|----|----|-------|-----------|
  | i      | input                | std::string |       |input image file | 
  | o      | output               | std::string |       |Output image file | 
  | f      | filter               | std::string |       |filter function (nvalid, median, var, min, max, sum, mean, dilate, erode, close, open, homog (central pixel must be identical to all other pixels within window), heterog (central pixel must be different than all other pixels within window), sobelx (horizontal edge detection), sobely (vertical edge detection), sobelxy (diagonal edge detection NE-SW),sobelyx (diagonal edge detection NW-SE), density, countid, mode (majority voting), only for classes), smoothnodata (smooth nodata values only) values, ismin, ismax, order (rank pixels in order), stdev, mrf, dwt, dwti, dwt_cut, dwt_cut_from, scramble, shift, savgolay, percentile, proportion) | 
  | srf    | srf                  | std::string |       |list of ASCII files containing spectral response functions (two columns: wavelength response) | 
  | fwhm   | fwhm                 | double |       |list of full width half to apply spectral filtering (-fwhm band1 -fwhm band2 ...) | 
  | dx     | dx                   | double | 3     |filter kernel size in x, use odd values only | 
  | dy     | dy                   | double | 3     |filter kernel size in y, use odd values only | 
  | dz     | dz                   | int  |       |filter kernel size in z (spectral/temporal dimension), must be odd (example: 3).. Set dz>0 if 1-D filter must be used in band domain | 
  | nodata | nodata               | double |       |nodata value(s) (used for smoothnodata filter) | 
  | r      | resampling-method    | std::string | near  |Resampling method for shifting operation (near: nearest neighbour, bilinear: bi-linear interpolation). | 
  | co     | co                   | std::string |       |Creation option for output file. Multiple options can be specified. | 
  | wt     | wavelet              | std::string | daubechies |wavelet type: daubechies,daubechies_centered, haar, haar_centered, bspline, bspline_centered | 
  | wf     | family               | int  | 4     |wavelet family (vanishing moment, see also http://www.gnu.org/software/gsl/manual/html_node/DWT-Initialization.html) | 
  | nl     | nl                   | int  | 2     |Number of leftward (past) data points used in Savitzky-Golay filter) | 
  | nr     | nr                   | int  | 2     |Number of rightward (future) data points used in Savitzky-Golay filter) | 
  | ld     | ld                   | int  | 0     |order of the derivative desired in Savitzky-Golay filter (e.g., ld=0 for smoothed function) | 
  | m      | m                    | int  | 2     |order of the smoothing polynomial in Savitzky-Golay filter, also equal to the highest conserved moment; usual values are m = 2 or m = 4) | 
  | class  | class                | short |       |class value(s) to use for density, erosion, dilation, openening and closing, thresholding | 
  | t      | threshold            | double | 0     |threshold value(s) to use for threshold filter (one for each class), or threshold to cut for dwt_cut (use 0 to keep all) or dwt_cut_from, or sigma for shift | 
  | tap    | tap                  | std::string |       |text file containing taps used for spatial filtering (from ul to lr). Use dimX and dimY to specify tap dimensions in x and y. Leave empty for not using taps | 
  | tapz   | tapz                 | double |       |taps used for spectral filtering | 
  | pad    | pad                  | std::string | symmetric |Padding method for filtering (how to handle edge effects). Choose between: symmetric, replicate, circular, zero (pad with 0). | 
  | win    | wavelengthIn         | double |       |list of wavelengths in input spectrum (-win band1 -win band2 ...) | 
  | wout   | wavelengthOut        | double |       |list of wavelengths in output spectrum (-wout band1 -wout band2 ...) | 
  | d      | down                 | short | 1     |down sampling factor. Use value 1 for no downsampling). Use value n>1 for downsampling (aggregation) | 
  | beta   | beta                 | std::string |       |ASCII file with beta for each class transition in Markov Random Field | 
  | interp | interp               | std::string | akima |type of interpolation for spectral filtering (see http://www.gnu.org/software/gsl/manual/html_node/Interpolation-Types.html) | 
  | ot     | otype                | std::string |       |Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image | 
  | of     | oformat              | std::string | GTiff |Output image format (see also gdal_translate).| 
  | ct     | ct                   | std::string |       |color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid). Use none to ommit color table | 
  | circ   | circular             | bool | false |circular disc kernel for dilation and erosion | 
  | mem    | mem                  | unsigned long int | 0 |Buffer size (in MB) to read image data blocks in memory | 

  Usage: pkfilter -i input -o output [-f filter | -perc value | -srf file [-srf file]* -win wavelength [-win wavelength]* | -wout wavelength -fwhm value [-wout wavelength -fwhm value]* -win wavelength [-win wavelength]*]

  Examples
  ========
  Some examples how to use pkfilter can be found \ref examples_pkfilter "here"
**/

using namespace std;
/*------------------
  Main procedure
  ----------------*/
int main(int argc,char **argv) {
  Optionpk<std::string> input_opt("i","input","input image file");
  Optionpk<std::string> output_opt("o", "output", "Output image file");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate).","GTiff");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<unsigned long int>  memory_opt("mem", "mem", "Buffer size (in MB) to read image data blocks in memory",0,1);

  option_opt.setHide(1);
  oformat_opt.setHide(1);
  memory_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    memory_opt.retrieveOption(argc,argv);

    app::AppFactory app(argc,argv);

    if(doProcess&&input_opt.empty()){
      std::ostringstream errorStream;
      errorStream << "Error: no input file provided (use option -i). Use --help for help information" << std::endl;
      throw(errorStream.str());
    }
    if(doProcess&&output_opt.empty()){
      std::ostringstream errorStream;
      errorStream << "Error: no output file provided (use option -o). Use --help for help information" << std::endl;
      throw(errorStream.str());
    }

    ImgRaster input;
    string imageType;//=input.getImageType();
    if(input_opt.size()){
      input.open(input_opt[0],memory_opt[0]);
      if(oformat_opt.size())//default
        imageType=oformat_opt[0];
      else
        imageType=input.getImageType();
    }
    ImgRaster imgWriter;
    imgWriter.setFile(output_opt[0],imageType,memory_opt[0],option_opt);
    input.filter(imgWriter,app);

    input.close();
    imgWriter.close();
  }
  catch(string helpString){
    cerr << helpString << endl;
    cout << "Usage: pkfilter -i input -o output [-f filter | -perc value | -srf file [-srf file]* -win wavelength [-win wavelength]* | -wout wavelength -fwhm value [-wout wavelength -fwhm value]* -win wavelength [-win wavelength]*]" << endl;
    return(1);
  }
  return(0);
}
