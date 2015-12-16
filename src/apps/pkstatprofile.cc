/**********************************************************************
pkstatprofile.cc: program to calculate statistics in temporal or spectral profile
Copyright (C) 2008-2015 Pieter Kempeneers

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
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"
#include "algorithms/StatFactory.h"

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

Usage: pkstatprofile -i input -o ouptut [-f function]*

**/

using namespace std;
/*------------------
  Main procedure
  ----------------*/
int main(int argc,char **argv) {
  Optionpk<std::string> input_opt("i","input","input image file");
  Optionpk<std::string> output_opt("o", "output", "Output image file");
  Optionpk<std::string> function_opt("f", "function", "Statistics function (mean, median, var, stdev, min, max, sum, mode (provide classes), ismin, ismax, proportion (provide classes), percentile, nvalid");
  Optionpk<double> percentile_opt("perc","perc","Percentile value(s) used for rule percentile");
  Optionpk<short> class_opt("class", "class", "class value(s) to use for mode, proportion");
  Optionpk<double> nodata_opt("nodata", "nodata", "nodata value(s)");
  Optionpk<std::string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image","");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate)","GTiff");
  Optionpk<string>  colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid). Use none to ommit color table");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  // Optionpk<short> down_opt("d", "down", "down sampling factor. Use value 1 for no downsampling). Use value n>1 for downsampling (aggregation)", 1);
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0,2);

  percentile_opt.setHide(1);
  class_opt.setHide(1);
  otype_opt.setHide(1);
  oformat_opt.setHide(1);
  colorTable_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    function_opt.retrieveOption(argc,argv);
    percentile_opt.retrieveOption(argc,argv);
    nodata_opt.retrieveOption(argc,argv);
    otype_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    cout << endl;
    cout << "Usage: pkstatprofile -i input -o ouptut [-function]*" << endl;
    cout << endl;
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }

  if(function_opt.empty()){
    cerr << "Error: no function selected, use option -f" << endl;
    exit(1);
  }

  filter::Filter filter1d;
  filter1d.setNoDataValues(nodata_opt);
  for(int iclass=0;iclass<class_opt.size();++iclass)
    filter1d.pushClass(class_opt[iclass]);

  filter1d.setThresholds(percentile_opt);

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
  try{
    input.open(input_opt[0]);
  }
  catch(string errorstring){
    cout << errorstring << endl;
    exit(1);
  }

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
  if(verbose_opt[0])
    cout << "Calculating statistic metrics: " << function_opt.size() << endl;
  try{
    output.open(output_opt[0],input.nrOfCol(),input.nrOfRow(),function_opt.size(),theType,imageType,option_opt);
  }
  catch(string errorstring){
    cout << errorstring << endl;
    exit(4);
  }

  output.setProjection(input.getProjection());
  double gt[6];
  input.getGeoTransform(gt);
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

  try{
    filter1d.stats(input,output,function_opt);
  }
  catch(string errorstring){
    cerr << errorstring << endl;
  }
  input.close();
  output.close();
  return 0;
}
