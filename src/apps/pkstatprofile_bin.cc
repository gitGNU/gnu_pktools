/**********************************************************************
pkstatprofile_bin.cc: program to calculate statistics in temporal or spectral profile
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
/*! \page pkstatprofile pkstatprofile
  program to calculate statistics in temporal or spectral profile
  ## SYNOPSIS

  <code>
  Usage: pkstatprofile -i input -o output [-f function]*
  </code>

  <code>
  Options: [-nodata value]

  Advanced options: check table
  </code>

  \section pkstatprofile_description Description

  This utility calculates statistics for a temporal (time series) or spectral profile

  \anchor pkstatprofile_functions

  |function | description|
  |-------------|-------------|
  mean | calculate mean in window
  median | perform a median filter in spatial (dx, dy) or spectral/temporal (dz) domain
  var | calculate variance in window
  stdev | calculate standard deviation in window
  min | calculate minimum in window
  max | calculate maximum in window
  sum | calculate sum in window
  mode | calculate mode of all values
  ismin | 1 if value is minimum, else 0
  ismax | 1 if value is maximum, else 0
  per | calculate percentile in time series (provide percentage value as argument)
  prop  |calculate proportion
  nvalid | report number of valid observations

  Example: Calculate min and max NDVI in time series

  \code
  pkstatprofile -i modis_ndvi_2010.tif -o modis_stats_2010.tif -f min -f max
  \endcode

  \section pkfilter_options Options
  - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
  - short option `-h` shows basic options only, long option `--help` shows all options
  |short|long|type|default|description|
  |-----|----|----|-------|-----------|
  | i      | input                | std::string |       |input image file | 
  | o      | output               | std::string |       |Output image file | 
  | f      | function             | std::string |       |statistics function (see table) | 
  | perc   | percentile           | double |  |percentile value(s) for percentile function | 
  |class | class | std::string | | class value(s) to use for mode, proportion |
  |nodata | nodata | double | | nodata value(s) |
  | mem    | mem                  | unsigned long int | 0 |Buffer size (in MB) to read image data blocks in memory | 

  Usage: pkstatprofile -i input -o output [-f function]*

**/

using namespace std;
/*------------------
  Main procedure
  ----------------*/
int main(int argc,char **argv) {
  Optionpk<std::string> input_opt("i","input","input image file");
  Optionpk<std::string> output_opt("o", "output", "Output image file");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate)","GTiff");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  // Optionpk<short> down_opt("d", "down", "down sampling factor. Use value 1 for no downsampling). Use value n>1 for downsampling (aggregation)", 1);
  Optionpk<unsigned long int>  memory_opt("mem", "mem", "Buffer size (in MB) to read image data blocks in memory",0,1);

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
    input.open(input_opt[0],memory_opt[0]);

    std::shared_ptr<ImgRaster> imgWriter = std::make_shared<ImgRaster>();

    string imageType;//=input.getImageType();
    if(oformat_opt.size())//default
      imageType=oformat_opt[0];
    else
      imageType=input.getImageType();

    imgWriter->setFile(output_opt[0],imageType,memory_opt[0],option_opt);
    input.statProfile(imgWriter,app);

    input.close();
    imgWriter->close();
  }
  catch(string helpString){
    cerr << helpString << endl;
    cout << "Usage: pkstatprofile -i input [-f function]* -o output" << endl;
    return(1);
  }
  return(0);
}
