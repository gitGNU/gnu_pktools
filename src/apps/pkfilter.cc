/**********************************************************************
pkfilter.cc: program to filter raster images: median, min/max, morphological, filtering
Copyright (C) 2008-2014 Pieter Kempeneers

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
#include <assert.h>
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <sys/types.h>
#include <stdio.h>
#include "base/Optionpk.h"
#include "base/Vector2d.h"
#include "algorithms/Filter2d.h"
#include "algorithms/Filter.h"
#include "fileclasses/FileReaderAscii.h"
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"
#include "algorithms/StatFactory.h"

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
  // Optionpk<std::string> tmpdir_opt("tmp", "tmp", "Temporary directory","/tmp",2);
  Optionpk<bool> disc_opt("circ", "circular", "circular disc kernel for dilation and erosion", false);
  // Optionpk<double> angle_opt("a", "angle", "angle used for directional filtering in dilation (North=0, East=90, South=180, West=270).");
  Optionpk<std::string> method_opt("f", "filter", "filter function (nvalid, median, var, min, max, sum, mean, dilate, erode, close, open, homog (central pixel must be identical to all other pixels within window), heterog (central pixel must be different than all other pixels within window), sauvola, sobelx (horizontal edge detection), sobely (vertical edge detection), sobelxy (diagonal edge detection NE-SW),sobelyx (diagonal edge detection NW-SE), density, countid, mode (majority voting), only for classes), smoothnodata (smooth nodata values only) values, ismin, ismax, order (rank pixels in order), stdev, mrf, dwt, dwti, dwt_cut, dwt_cut_from, scramble, shift, savgolay, percentile, proportion)");
  Optionpk<std::string> resample_opt("r", "resampling-method", "Resampling method for shifting operation (near: nearest neighbour, bilinear: bi-linear interpolation).", "near");
  Optionpk<double> dimX_opt("dx", "dx", "filter kernel size in x, use odd values only", 3);
  Optionpk<double> dimY_opt("dy", "dy", "filter kernel size in y, use odd values only", 3);
  Optionpk<int> dimZ_opt("dz", "dz", "filter kernel size in z (spectral/temporal dimension), must be odd (example: 3).. Set dz>0 if 1-D filter must be used in band domain");
  Optionpk<std::string> wavelet_type_opt("wt", "wavelet", "wavelet type: daubechies,daubechies_centered, haar, haar_centered, bspline, bspline_centered", "daubechies");
  Optionpk<int> family_opt("wf", "family", "wavelet family (vanishing moment, see also http://www.gnu.org/software/gsl/manual/html_node/DWT-Initialization.html)", 4);
  Optionpk<int> savgolay_nl_opt("nl", "nl", "Number of leftward (past) data points used in Savitzky-Golay filter)", 2);
  Optionpk<int> savgolay_nr_opt("nr", "nr", "Number of rightward (future) data points used in Savitzky-Golay filter)", 2);
  Optionpk<int> savgolay_ld_opt("ld", "ld", "order of the derivative desired in Savitzky-Golay filter (e.g., ld=0 for smoothed function)", 0);
  Optionpk<int> savgolay_m_opt("m", "m", "order of the smoothing polynomial in Savitzky-Golay filter, also equal to the highest conserved moment; usual values are m = 2 or m = 4)", 2);
  Optionpk<short> class_opt("class", "class", "class value(s) to use for density, erosion, dilation, openening and closing, thresholding");
  Optionpk<double> threshold_opt("t", "threshold", "threshold value(s) to use for threshold filter (one for each class), or threshold to cut for dwt_cut (use 0 to keep all) or dwt_cut_from, or sigma for shift", 0);
  Optionpk<double> nodata_opt("nodata", "nodata", "nodata value(s) (used for smoothnodata filter)");
  Optionpk<std::string> tap_opt("tap", "tap", "text file containing taps used for spatial filtering (from ul to lr). Use dimX and dimY to specify tap dimensions in x and y. Leave empty for not using taps");
  Optionpk<double> tapz_opt("tapz", "tapz", "taps used for spectral filtering");
  Optionpk<string> padding_opt("pad","pad", "Padding method for filtering (how to handle edge effects). Choose between: symmetric, replicate, circular, zero (pad with 0).", "symmetric");
  Optionpk<double> fwhm_opt("fwhm", "fwhm", "list of full width half to apply spectral filtering (-fwhm band1 -fwhm band2 ...)");
  Optionpk<std::string> srf_opt("srf", "srf", "list of ASCII files containing spectral response functions (two columns: wavelength response)");
  Optionpk<double> wavelengthIn_opt("win", "wavelengthIn", "list of wavelengths in input spectrum (-win band1 -win band2 ...)");
  Optionpk<double> wavelengthOut_opt("wout", "wavelengthOut", "list of wavelengths in output spectrum (-wout band1 -wout band2 ...)");
  Optionpk<std::string> interpolationType_opt("interp", "interp", "type of interpolation for spectral filtering (see http://www.gnu.org/software/gsl/manual/html_node/Interpolation-Types.html)","akima");
  Optionpk<std::string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image","");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate).","GTiff");
  Optionpk<string>  colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid). Use none to ommit color table");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<short> down_opt("d", "down", "down sampling factor. Use value 1 for no downsampling). Use value n>1 for downsampling (aggregation)", 1);
  Optionpk<string> beta_opt("beta", "beta", "ASCII file with beta for each class transition in Markov Random Field");
  // Optionpk<double> eps_opt("eps","eps", "error marging for linear feature",0);
  // Optionpk<bool> l1_opt("l1","l1", "obtain longest object length for linear feature",false);
  // Optionpk<bool> l2_opt("l2","l2", "obtain shortest object length for linear feature",false,2);
  // Optionpk<bool> a1_opt("a1","a1", "obtain angle found for longest object length for linear feature",false);
  // Optionpk<bool> a2_opt("a2","a2", "obtain angle found for shortest object length for linear feature",false);
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0,2);

  resample_opt.setHide(1);
  option_opt.setHide(1);
  wavelet_type_opt.setHide(1);
  family_opt.setHide(1);
  savgolay_nl_opt.setHide(1);
  savgolay_nr_opt.setHide(1);
  savgolay_ld_opt.setHide(1);
  savgolay_m_opt.setHide(1);
  class_opt.setHide(1);
  threshold_opt.setHide(1);
  tap_opt.setHide(1);
  tapz_opt.setHide(1);
  padding_opt.setHide(1);
  wavelengthIn_opt.setHide(1);
  wavelengthOut_opt.setHide(1);
  down_opt.setHide(1);
  beta_opt.setHide(1);
  // eps_opt.setHide(1);
  // l1_opt.setHide(1);
  // l2_opt.setHide(1);
  // a1_opt.setHide(1);
  // a2_opt.setHide(1);
  interpolationType_opt.setHide(1);
  otype_opt.setHide(1);
  oformat_opt.setHide(1);
  colorTable_opt.setHide(1);
  disc_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    // tmpdir_opt.retrieveOption(argc,argv);
    // angle_opt.retrieveOption(argc,argv);
    method_opt.retrieveOption(argc,argv);
    srf_opt.retrieveOption(argc,argv);
    fwhm_opt.retrieveOption(argc,argv);
    dimX_opt.retrieveOption(argc,argv);
    dimY_opt.retrieveOption(argc,argv);
    dimZ_opt.retrieveOption(argc,argv);
    nodata_opt.retrieveOption(argc,argv);
    resample_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    wavelet_type_opt.retrieveOption(argc,argv);
    family_opt.retrieveOption(argc,argv);
    savgolay_nl_opt.retrieveOption(argc,argv);
    savgolay_nr_opt.retrieveOption(argc,argv);
    savgolay_ld_opt.retrieveOption(argc,argv);
    savgolay_m_opt.retrieveOption(argc,argv);
    class_opt.retrieveOption(argc,argv);
    threshold_opt.retrieveOption(argc,argv);
    tap_opt.retrieveOption(argc,argv);
    tapz_opt.retrieveOption(argc,argv);
    padding_opt.retrieveOption(argc,argv);
    wavelengthIn_opt.retrieveOption(argc,argv);
    wavelengthOut_opt.retrieveOption(argc,argv);
    down_opt.retrieveOption(argc,argv);
    beta_opt.retrieveOption(argc,argv);
    // eps_opt.retrieveOption(argc,argv);
    // l1_opt.retrieveOption(argc,argv);
    // l2_opt.retrieveOption(argc,argv);
    // a1_opt.retrieveOption(argc,argv);
    // a2_opt.retrieveOption(argc,argv);
    interpolationType_opt.retrieveOption(argc,argv);
    otype_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
    disc_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    cout << endl;
    cout << "Usage: pkfilter -i input -o output [-f filter | -perc value | -srf file [-srf file]* -win wavelength [-win wavelength]* | -wout wavelength -fwhm value [-wout wavelength -fwhm value]* -win wavelength [-win wavelength]*]" << endl;
    cout << endl;
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }

  //not implemented yet, must debug first...
  vector<double> angle_opt;

  ImgReaderGdal input;
  ImgWriterGdal output;
  if(input_opt.empty()){
    cerr << "Error: no input file selected, use option -i" << endl;
    exit(1);
  }
  if(output_opt.empty()){
    cerr << "Error: no output file selected, use option -o" << endl;
    exit(1);
  }
  input.open(input_opt[0]);
  GDALDataType theType=GDT_Unknown;
  if(verbose_opt[0])
    cout << "possible output data types: ";
  for(int iType = 0; iType < GDT_TypeCount; ++iType){
    if(verbose_opt[0])
      cout << " " << GDALGetDataTypeName((GDALDataType)iType);
    if( GDALGetDataTypeName((GDALDataType)iType) != NULL
        && EQUAL(GDALGetDataTypeName((GDALDataType)iType),
                 otype_opt[0].c_str()))
      theType=(GDALDataType) iType;
  }
  if(theType==GDT_Unknown)
    theType=input.getDataType();

  if(verbose_opt[0])
    std::cout << std::endl << "Output pixel type:  " << GDALGetDataTypeName(theType) << endl;

  string imageType;//=input.getImageType();
  if(oformat_opt.size())
    imageType=oformat_opt[0];

  if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
    string theInterleave="INTERLEAVE=";
    theInterleave+=input.getInterleave();
    option_opt.push_back(theInterleave);
  }

  try{
    string errorString;
    int nband=input.nrOfBand();

    if(fwhm_opt.size())
      nband=fwhm_opt.size();
    else if(srf_opt.size())
      nband=srf_opt.size();
    else if(tap_opt.size()||tapz_opt.size())
      nband=input.nrOfBand();
    else{
      if(method_opt.empty()){
        errorString="Error: no filter selected, use option -f";
        throw(errorString);
      }
      switch(filter2d::Filter2d::getFilterType(method_opt[0])){
      case(filter2d::dilate):
      case(filter2d::erode):
      case(filter2d::close):
      case(filter2d::open):
      case(filter2d::smooth):
	//implemented in spectral/temporal domain (dimZ>1) and spatial domain
	if(dimZ_opt.size())
	  assert(dimZ_opt[0]>1);
      nband=input.nrOfBand();
      break;
      case(filter2d::dwt):
      case(filter2d::dwti):
      case(filter2d::dwt_cut):
      case(filter2d::smoothnodata):
	//implemented in spectral/temporal/spatial domain and nband always input.nrOfBand()
	nband=input.nrOfBand();
	break;
      case(filter2d::savgolay):
	nband=input.nrOfBand();
	if(dimZ_opt.empty())
	  dimZ_opt.push_back(1);
      case(filter2d::dwt_cut_from):
	//only implemented in spectral/temporal domain
	if(dimZ_opt.size()){
	  nband=input.nrOfBand();
	  assert(threshold_opt.size());
	}
	else{
          errorString="filter not implemented in spatial domain";
          throw(errorString);
	}
	break;
      case(filter2d::mrf)://deliberate fall through
	assert(class_opt.size()>1);
	if(verbose_opt[0])
	  std::cout << "opening output image " << output_opt[0] << std::endl;
	nband=class_opt.size();
      case(filter2d::ismin):
      case(filter2d::ismax):
      case(filter2d::shift):
      case(filter2d::scramble):
      case(filter2d::mode):
      case(filter2d::sobelx):
      case(filter2d::sobely):
      case(filter2d::sobelxy):
      case(filter2d::countid):
      case(filter2d::order):
      case(filter2d::density):
      case(filter2d::homog):
      case(filter2d::heterog):
      case(filter2d::sauvola):
	//only implemented in spatial domain
	if(dimZ_opt.size()){
          errorString="filter not implemented in spectral/temporal domain";
          throw(errorString);
	}
	break;
      // case(filter2d::percentile):
      // 	//implemented in spectral/temporal/spatial domain and nband 1 if dimZ>0
      // 	if(dimZ_opt.size()){
      // 	  dimZ_opt[0]=1;
      // 	  nband=1;
      // 	}
      // 	else
      // 	  nband=input.nrOfBand();
      // 	break;
      case(filter2d::sum):
      case(filter2d::mean):
      case(filter2d::min):
      case(filter2d::max):
      case(filter2d::var):
      case(filter2d::stdev):
      case(filter2d::nvalid):
      case(filter2d::median):
      case(filter2d::percentile):
      case(filter2d::proportion):
	//implemented in spectral/temporal/spatial domain and nband 1 if dimZ==1
	if(dimZ_opt.size()==1)
	  if(dimZ_opt[0]==1)
	    nband=1;
	else
	  nband=input.nrOfBand();
	break;
      default:
	errorString="filter not implemented";
        throw(errorString);
	break;
      }
    }
    std::cout << "opening output image " << output_opt[0] << " with " << nband << " bands" << std::endl;
    output.open(output_opt[0],(input.nrOfCol()+down_opt[0]-1)/down_opt[0],(input.nrOfRow()+down_opt[0]-1)/down_opt[0],nband,theType,imageType,option_opt);
  }
  catch(string errorstring){
    cerr << errorstring << endl;
    exit(1);
  }
  output.setProjection(input.getProjection());
  double gt[6];
  input.getGeoTransform(gt);
  gt[1]*=down_opt[0];//dx
  gt[5]*=down_opt[0];//dy
  output.setGeoTransform(gt);
  
  if(colorTable_opt.size()){
    if(colorTable_opt[0]!="none"){
      if(verbose_opt[0])
	cout << "set colortable " << colorTable_opt[0] << endl;
      assert(output.getDataType()==GDT_Byte);
      output.setColorTable(colorTable_opt[0]);
    }
  }
  else if(input.getColorTable()!=NULL)
    output.setColorTable(input.getColorTable());
  
  if(nodata_opt.size()){
      for(int iband=0;iband<output.nrOfBand();++iband)
	output.GDALSetNoDataValue(nodata_opt[0],iband);
  }

  filter2d::Filter2d filter2d;
  filter::Filter filter1d;
  if(verbose_opt[0])
    cout << "Set padding to " << padding_opt[0] << endl;
  filter1d.setPadding(padding_opt[0]);
  if(class_opt.size()){
    if(verbose_opt[0])
      std::cout<< "class values: ";
    for(int iclass=0;iclass<class_opt.size();++iclass){
      if(!dimZ_opt.size())
        filter2d.pushClass(class_opt[iclass]);
      else
        filter1d.pushClass(class_opt[iclass]);
      if(verbose_opt[0])
        std::cout<< class_opt[iclass] << " ";
    }
    if(verbose_opt[0])
      std::cout<< std::endl;
  }

  if(nodata_opt.size()){
    if(verbose_opt[0])
      std::cout<< "mask values: ";
    for(int imask=0;imask<nodata_opt.size();++imask){
      if(verbose_opt[0])
        std::cout<< nodata_opt[imask] << " ";
      filter1d.pushNoDataValue(nodata_opt[imask]);
      filter2d.pushNoDataValue(nodata_opt[imask]);
    }
    if(verbose_opt[0])
      std::cout<< std::endl;
  }
  filter1d.setThresholds(threshold_opt);
  filter2d.setThresholds(threshold_opt);

  if(tap_opt.size()){
    ifstream tapfile(tap_opt[0].c_str());
    assert(tapfile);
    Vector2d<double> taps(dimY_opt[0],dimX_opt[0]);

    for(int j=0;j<dimY_opt[0];++j){
      for(int i=0;i<dimX_opt[0];++i){
        tapfile >> taps[j][i];
      }
    }
    if(verbose_opt[0]){
      std::cout << "taps: ";
      for(int j=0;j<dimY_opt[0];++j){
        for(int i=0;i<dimX_opt[0];++i){
          std::cout<< taps[j][i] << " ";
        }
        std::cout<< std::endl;
      }
    }
    filter2d.setTaps(taps);
    try{
      filter2d.filter(input,output);
    }
    catch(string errorstring){
      cerr << errorstring << endl;
    }
    tapfile.close();
  }
  else if(tapz_opt.size()){
    if(verbose_opt[0]){
      std::cout << "taps: ";
      for(int itap=0;itap<tapz_opt.size();++itap)
	std::cout<< tapz_opt[itap] << " ";
      std::cout<< std::endl;
    }
    filter1d.setTaps(tapz_opt);
    filter1d.filter(input,output);
  }
  else if(fwhm_opt.size()){
    if(verbose_opt[0])
      std::cout << "spectral filtering to " << fwhm_opt.size() << " bands with provided fwhm " << std::endl;
    assert(wavelengthOut_opt.size()==fwhm_opt.size());
    assert(wavelengthIn_opt.size());

    Vector2d<double> lineInput(input.nrOfBand(),input.nrOfCol());
    Vector2d<double> lineOutput(wavelengthOut_opt.size(),input.nrOfCol());
    const char* pszMessage;
    void* pProgressArg=NULL;
    GDALProgressFunc pfnProgress=GDALTermProgress;
    double progress=0;
    pfnProgress(progress,pszMessage,pProgressArg);
    for(int y=0;y<input.nrOfRow();++y){
      if((y+1+down_opt[0]/2)%down_opt[0])
        continue;
      for(int iband=0;iband<input.nrOfBand();++iband)
        input.readData(lineInput[iband],y,iband);
      filter1d.applyFwhm<double>(wavelengthIn_opt,lineInput,wavelengthOut_opt,fwhm_opt, interpolationType_opt[0], lineOutput, down_opt[0], verbose_opt[0]);
      for(int iband=0;iband<output.nrOfBand();++iband){
        try{
          output.writeData(lineOutput[iband],y/down_opt[0],iband);
        }
        catch(string errorstring){
          cerr << errorstring << "in band " << iband << ", line " << y << endl;
        }
      }
      progress=(1.0+y)/output.nrOfRow();
      pfnProgress(progress,pszMessage,pProgressArg);
    }
  }
  else if(srf_opt.size()){
    if(verbose_opt[0])
      std::cout << "spectral filtering to " << srf_opt.size() << " bands with provided SRF " << std::endl;
    assert(wavelengthIn_opt.size());
    vector< Vector2d<double> > srf(srf_opt.size());//[0] srf_nr, [1]: wavelength, [2]: response
    ifstream srfFile;
    for(int isrf=0;isrf<srf_opt.size();++isrf){
      srf[isrf].resize(2);
      srfFile.open(srf_opt[isrf].c_str());
      double v;
      //add 0 to make sure srf is 0 at boundaries after interpolation step
      srf[isrf][0].push_back(0);
      srf[isrf][1].push_back(0);
      srf[isrf][0].push_back(1);
      srf[isrf][1].push_back(0);
      while(srfFile >> v){
        srf[isrf][0].push_back(v);
        srfFile >> v;
        srf[isrf][1].push_back(v);
      }
      srfFile.close();
      //add 0 to make sure srf[isrf] is 0 at boundaries after interpolation step
      srf[isrf][0].push_back(srf[isrf][0].back()+1);
      srf[isrf][1].push_back(0);    
      srf[isrf][0].push_back(srf[isrf][0].back()+1);
      srf[isrf][1].push_back(0);
      if(verbose_opt[0])
        cout << "srf file details: " << srf[isrf][0].size() << " wavelengths defined" << endl;    
    }
    assert(output.nrOfBand()==srf.size());
    double centreWavelength=0;
    Vector2d<double> lineInput(input.nrOfBand(),input.nrOfCol());
    const char* pszMessage;
    void* pProgressArg=NULL;
    GDALProgressFunc pfnProgress=GDALTermProgress;
    double progress=0;
    pfnProgress(progress,pszMessage,pProgressArg);
    for(int y=0;y<input.nrOfRow();++y){
      if((y+1+down_opt[0]/2)%down_opt[0])
        continue;
      for(int iband=0;iband<input.nrOfBand();++iband)
        input.readData(lineInput[iband],y,iband);
      for(int isrf=0;isrf<srf.size();++isrf){
        vector<double> lineOutput(output.nrOfCol());
        double delta=1.0;
        bool normalize=true;
	centreWavelength=filter1d.applySrf<double>(wavelengthIn_opt,lineInput,srf[isrf], interpolationType_opt[0], lineOutput, delta, normalize);
        if(verbose_opt[0])
          std::cout << "centre wavelength srf " << isrf << ": " << centreWavelength << std::endl;
        try{
          output.writeData(lineOutput,GDT_Float64,y/down_opt[0],isrf);
        }
        catch(string errorstring){
          cerr << errorstring << "in srf " << srf_opt[isrf] << ", line " << y << endl;
        }

      }
      progress=(1.0+y)/output.nrOfRow();
      pfnProgress(progress,pszMessage,pProgressArg);
    }

  }
  else{
    switch(filter2d::Filter2d::getFilterType(method_opt[0])){
    case(filter2d::dilate):
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for morphological operator" << std::endl;
	exit(1);
      }
      try{
	if(dimZ_opt.size()){
	  if(verbose_opt[0])
          std::cout<< "1-D filtering: dilate" << std::endl;
	  filter1d.morphology(input,output,"dilate",dimZ_opt[0],verbose_opt[0]);
	}
	else
	  filter2d.morphology(input,output,"dilate",dimX_opt[0],dimY_opt[0],angle_opt,disc_opt[0]);
      }
      catch(string errorstring){
	cerr << errorstring << endl;
      }
      break;
    case(filter2d::erode):
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for morphological operator" << std::endl;
	exit(1);
      }
      try{
	if(dimZ_opt.size()>0){
	  if(verbose_opt[0])
	    std::cout<< "1-D filtering: dilate" << std::endl;
	  filter1d.morphology(input,output,"erode",dimZ_opt[0]);
	}
	else{
	  filter2d.morphology(input,output,"erode",dimX_opt[0],dimY_opt[0],angle_opt,disc_opt[0]);
	}
      }
      catch(string errorstring){
	cerr << errorstring << endl;
      }
    break;
    case(filter2d::close):{//closing
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for morphological operator" << std::endl;
	exit(1);
      }

      ImgWriterGdal tmpout;
      tmpout.open("/vsimem/dilation.tif",input.nrOfCol(),input.nrOfRow(),input.nrOfBand(),input.getDataType(),input.getImageType());
      try{
        if(dimZ_opt.size()){
          filter1d.morphology(input,tmpout,"dilate",dimZ_opt[0]);
        }
        else{
	  filter2d.morphology(input,tmpout,"dilate",dimX_opt[0],dimY_opt[0],angle_opt,disc_opt[0]);
        }
      }
      catch(std::string errorString){
	std::cout<< errorString;
	exit(1);
      }
      tmpout.close();
      ImgReaderGdal tmpin;
      tmpin.open("/vsimem/dilation.tif");
      try{
	if(dimZ_opt.size()){
	  filter1d.morphology(tmpin,output,"erode",dimZ_opt[0]);
	}
	else{
	  filter2d.morphology(tmpin,output,"erode",dimX_opt[0],dimY_opt[0],angle_opt,disc_opt[0]);
	}
      }
      catch(string errorstring){
	cerr << errorstring << endl;
      }
      tmpin.close();
      break;
    }
    case(filter2d::open):{//opening
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for morphological operator" << std::endl;
	exit(1);
      }
      ImgWriterGdal tmpout;
      tmpout.open("/vsimem/erosion.tif",input.nrOfCol(),input.nrOfRow(),input.nrOfBand(),input.getDataType(),input.getImageType());
      try{
	if(dimZ_opt.size()){
	  filter1d.morphology(input,tmpout,"erode",dimZ_opt[0]);
	}
	else{
	  filter2d.morphology(input,tmpout,"erode",dimX_opt[0],dimY_opt[0],angle_opt,disc_opt[0]);
	}
      }
      catch(std::string errorString){
	std::cout<< errorString;
	exit(1);
      }
      tmpout.close();
      ImgReaderGdal tmpin;
      try{
	tmpin.open("/vsimem/erosion.tif");
	if(dimZ_opt.size()){
	  filter1d.morphology(tmpin,output,"dilate",dimZ_opt[0]);
	}
	else{
	  filter2d.morphology(tmpin,output,"dilate",dimX_opt[0],dimY_opt[0],angle_opt,disc_opt[0]);
	}
	tmpin.close();
	tmpout.close();
      }
      catch(string errorstring){
	cerr << errorstring << endl;
      }
      break;
    }
    case(filter2d::homog):{//spatially homogeneous
      try{
	filter2d.doit(input,output,"homog",dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      }
      catch(string errorstring){
	cerr << errorstring << endl;
      }
	break;
    }
    case(filter2d::heterog):{//spatially heterogeneous
      try{
	filter2d.doit(input,output,"heterog",dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      }
      catch(string errorstring){
	cerr << errorstring << endl;
      }
      break;
    }
    case(filter2d::sauvola):{//Implements Sauvola's thresholding method (http://fiji.sc/Auto_Local_Threshold)

  //test
      Vector2d<unsigned short> inBuffer;
      for(int iband=0;iband<input.nrOfBand();++iband){
        input.readDataBlock(inBuffer,0,input.nrOfCol()-1,0,input.nrOfRow()-1,iband);
      }
      try{
	filter2d.doit(input,output,"sauvola",dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      }
      catch(string errorstring){
	cerr << errorstring << endl;
      }
      break;
    }
    case(filter2d::shift):{//shift
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for shift operator" << std::endl;
	exit(1);
      }
      assert(input.nrOfBand());
      assert(input.nrOfCol());
      assert(input.nrOfRow());
      try{
        filter2d.shift(input,output,dimX_opt[0],dimY_opt[0],threshold_opt[0],filter2d::Filter2d::getResampleType(resample_opt[0]));
      }
      catch(string errorstring){
        cerr << errorstring << endl;
      }
      break;
    }
    // case(filter2d::linearfeature):{
    //   if(down_opt[0]!=1){
    // 	std::cerr << "Error: down option not supported for linear feature" << std::endl;
    // 	exit(1);
    //   }
    //   assert(input.nrOfBand());
    //   assert(input.nrOfCol());
    //   assert(input.nrOfRow());
    //   float theAngle=361;
    //   if(angle_opt.size())
    // 	theAngle=angle_opt[0];
    //   if(verbose_opt[0])
    // 	std::cout << "using angle " << theAngle << std::endl;
    //   try{
    // 	//using an angle step of 5 degrees and no maximum distance
    //     filter2d.linearFeature(input,output,theAngle,5,0,eps_opt[0],l1_opt[0],a1_opt[0],l2_opt[0],a2_opt[0],0,verbose_opt[0]);
    //   }
    //   catch(string errorstring){
    //     cerr << errorstring << endl;
    //   }
    //   break;
    // }
    case(filter2d::mrf):{//Markov Random Field
      if(verbose_opt[0])
	std::cout << "Markov Random Field filtering" << std::endl;
      try{
	if(beta_opt.size()){
	  //in file: classFrom classTo
	  //in variable: beta[classTo][classFrom]
	  FileReaderAscii betaReader(beta_opt[0]);
	  Vector2d<double> beta(class_opt.size(),class_opt.size());
	  vector<int> cols(class_opt.size());
	  for(int iclass=0;iclass<class_opt.size();++iclass)
	    cols[iclass]=iclass;
	  betaReader.readData(beta,cols);
	  if(verbose_opt[0]){
	    std::cout << "using values for beta:" << std::endl;
	    for(int iclass1=0;iclass1<class_opt.size();++iclass1)
	      std::cout << "      " << iclass1 << " (" << class_opt[iclass1] << ")";
	    std::cout << std::endl;
	    for(int iclass1=0;iclass1<class_opt.size();++iclass1){
	      std::cout << iclass1 << " (" << class_opt[iclass1] << ")";
	      for(int iclass2=0;iclass2<class_opt.size();++iclass2)
		std::cout << " " << beta[iclass2][iclass1] << " (" << class_opt[iclass2] << ")";
	      std::cout << std::endl;
	    }
	  }
	  filter2d.mrf(input, output, dimX_opt[0], dimY_opt[0], beta, true, down_opt[0], verbose_opt[0]);
	}
	else
	  filter2d.mrf(input, output, dimX_opt[0], dimY_opt[0], 1, true, down_opt[0], verbose_opt[0]);
      }
      catch(string errorstring){
	cerr << errorstring << endl;
      }
      break;
    }
    case(filter2d::sobelx):{//Sobel edge detection in X
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for sobel edge detection" << std::endl;
	exit(1);
      }
      Vector2d<double> theTaps(3,3);
      theTaps[0][0]=-1.0;
      theTaps[0][1]=0.0;
      theTaps[0][2]=1.0;
      theTaps[1][0]=-2.0;
      theTaps[1][1]=0.0;
      theTaps[1][2]=2.0;
      theTaps[2][0]=-1.0;
      theTaps[2][1]=0.0;
      theTaps[2][2]=1.0;
      filter2d.setTaps(theTaps);
      try{
	filter2d.filter(input,output,true,true);//absolute and normalize
      }
      catch(string errorstring){
	cerr << errorstring << endl;
      }
      break;
    }
    case(filter2d::sobely):{//Sobel edge detection in Y
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for sobel edge detection" << std::endl;
	exit(1);
      }
      Vector2d<double> theTaps(3,3);
      theTaps[0][0]=1.0;
      theTaps[0][1]=2.0;
      theTaps[0][2]=1.0;
      theTaps[1][0]=0.0;
      theTaps[1][1]=0.0;
      theTaps[1][2]=0.0;
      theTaps[2][0]=-1.0;
      theTaps[2][1]=-2.0;
      theTaps[2][2]=-1.0;
      filter2d.setTaps(theTaps);
      try{
	filter2d.filter(input,output,true,true);//absolute and normalize
      }
      catch(string errorstring){
	cerr << errorstring << endl;
      }
      break;
    }
    case(filter2d::sobelxy):{//Sobel edge detection in XY
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for sobel edge detection" << std::endl;
	exit(1);
      }
      Vector2d<double> theTaps(3,3);
      theTaps[0][0]=0.0;
      theTaps[0][1]=1.0;
      theTaps[0][2]=2.0;
      theTaps[1][0]=-1.0;
      theTaps[1][1]=0.0;
      theTaps[1][2]=1.0;
      theTaps[2][0]=-2.0;
      theTaps[2][1]=-1.0;
      theTaps[2][2]=0.0;
      filter2d.setTaps(theTaps);
      try{
	filter2d.filter(input,output,true,true);//absolute and normalize
      }
      catch(string errorstring){
	cerr << errorstring << endl;
      }
      break;
    }
    case(filter2d::sobelyx):{//Sobel edge detection in XY
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for sobel edge detection" << std::endl;
	exit(1);
      }
      Vector2d<double> theTaps(3,3);
      theTaps[0][0]=2.0;
      theTaps[0][1]=1.0;
      theTaps[0][2]=0.0;
      theTaps[1][0]=1.0;
      theTaps[1][1]=0.0;
      theTaps[1][2]=-1.0;
      theTaps[2][0]=0.0;
      theTaps[2][1]=-1.0;
      theTaps[2][2]=-2.0;
      filter2d.setTaps(theTaps);
      try{
	filter2d.filter(input,output,true,true);//absolute and normalize
      }
      catch(string errorstring){
	cerr << errorstring << endl;
      }
      break;
    }
    case(filter2d::smooth):{//Smoothing filter
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for this filter" << std::endl;
	exit(1);
      }
      try{
	if(dimZ_opt.size()){
	  if(verbose_opt[0])
	    std::cout<< "1-D filtering: smooth" << std::endl;
	  filter1d.smooth(input,output,dimZ_opt[0]);
	}
	else{
	  filter2d.smooth(input,output,dimX_opt[0],dimY_opt[0]);
	}
      }
      catch(string errorstring){
	cerr << errorstring << endl;
      }
      break;
    }
    case(filter2d::smoothnodata):{//Smoothing filter
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for this filter" << std::endl;
	exit(1);
      }
      try{
	if(dimZ_opt.size()){
	  if(verbose_opt[0])
	    std::cout<< "1-D filtering: smooth" << std::endl;
	  filter1d.smoothNoData(input,interpolationType_opt[0],output);
	}
	else{
	  if(verbose_opt[0])
	    std::cout<< "2-D filtering: smooth" << std::endl;
	  filter2d.smoothNoData(input,output,dimX_opt[0],dimY_opt[0]);
	}
      }
      catch(string errorstring){
	cerr << errorstring << endl;
      }
      break;
    }
    case(filter2d::dwt):
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for this filter" << std::endl;
	exit(1);
      }
      try{
	if(dimZ_opt.size()){
	  if(verbose_opt[0])
	    std::cout<< "DWT in spectral domain" << std::endl;
	  filter1d.dwtForward(input, output, wavelet_type_opt[0], family_opt[0]);
	}
	else
	  filter2d.dwtForward(input, output, wavelet_type_opt[0], family_opt[0]);
      }
      catch(string errorstring){
	cerr << errorstring << endl;
      }
      break;
    case(filter2d::dwti):
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for this filter" << std::endl;
	exit(1);
      }
      try{
	if(dimZ_opt.size()){
	  if(verbose_opt[0])
	    std::cout<< "inverse DWT in spectral domain" << std::endl;
	  filter1d.dwtInverse(input, output, wavelet_type_opt[0], family_opt[0]);
	}
	else
	  filter2d.dwtInverse(input, output, wavelet_type_opt[0], family_opt[0]);
      }
      catch(string errorstring){
	cerr << errorstring << endl;
      }
      break;
    case(filter2d::dwt_cut):
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for this filter" << std::endl;
	exit(1);
      }
      if(dimZ_opt.size()){
        if(verbose_opt[0])
          std::cout<< "DWT approximation in spectral domain" << std::endl;
	filter1d.dwtCut(input, output, wavelet_type_opt[0], family_opt[0], threshold_opt[0]);
      }
      else
	filter2d.dwtCut(input, output, wavelet_type_opt[0], family_opt[0], threshold_opt[0]);
      break;
    case(filter2d::dwt_cut_from):
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for this filter" << std::endl;
	exit(1);
      }
      try{
	if(dimZ_opt.size()){
	  if(verbose_opt[0])
	    std::cout<< "DWT approximation in spectral domain" << std::endl;
	  filter1d.dwtCutFrom(input, output, wavelet_type_opt[0], family_opt[0], static_cast<int>(threshold_opt[0]));
	}
	else{
	  string errorString="Error: this filter is not supported in 2D";
	  throw(errorString);
	}
      }
      catch(string errorstring){
	cerr << errorstring << endl;
      }
      break;
    case(filter2d::savgolay):{
      assert(savgolay_nl_opt.size());
      assert(savgolay_nr_opt.size());
      assert(savgolay_ld_opt.size());
      assert(savgolay_m_opt.size());
      if(verbose_opt[0])
      	std::cout << "Calculating Savitzky-Golay coefficients: " << endl;
      filter1d.getSavGolayCoefficients(tapz_opt, input.nrOfBand(), savgolay_nl_opt[0], savgolay_nr_opt[0], savgolay_ld_opt[0], savgolay_m_opt[0]);
      if(verbose_opt[0]){
      	std::cout << "taps (size is " << tapz_opt.size() << "): ";
      	for(int itap=0;itap<tapz_opt.size();++itap)
      	  std::cout<< tapz_opt[itap] << " ";
      	std::cout<< std::endl;
      }
      filter1d.setTaps(tapz_opt);
      filter1d.filter(input,output);
      break;
    }
    case(filter2d::percentile)://deliberate fall through
    case(filter2d::threshold)://deliberate fall through
      assert(threshold_opt.size());
      if(dimZ_opt.size())
	filter1d.setThresholds(threshold_opt);
      else
	filter2d.setThresholds(threshold_opt);
    case(filter2d::density)://deliberate fall through
      filter2d.setClasses(class_opt);
      if(verbose_opt[0])
	std::cout << "classes set" << std::endl;
    default:
      try{
	if(dimZ_opt.size()){
	  if(dimZ_opt[0]==1)
	    filter1d.stat(input,output,method_opt[0]);
	  else{
	    assert(down_opt[0]==1);//not implemented yet...
	    filter1d.filter(input,output,method_opt[0],dimZ_opt[0]);
	  }
	}
	else
	  filter2d.doit(input,output,method_opt[0],dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      }
      catch(string errorstring){
	cerr << errorstring << endl;
      }
      break;
    }
  }
  input.close();
  output.close();
  return 0;
}
