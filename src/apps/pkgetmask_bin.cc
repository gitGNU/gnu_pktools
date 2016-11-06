/**********************************************************************
pkgetmask_bin.cc: program to create mask image based on values in input raster image
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
#include <vector>
#include <iostream>
#include <string>
#include "imageclasses/ImgRasterGdal.h"
#include "base/Optionpk.h"
#include "AppFactory.h"

/******************************************************************************/
/*! \page pkgetmask pkgetmask
 program to create mask image based on values in input raster image
## SYNOPSIS

<code>
  Usage: pkgetmask -i input -o output
</code>

<code>
  
  Options: [-min value]* [-max value]* [-data value]* [-nodata value]*
 
  Advanced options: [-b band]* [--operator AND|OR] [-ot type] [-of format] [-co option]* [-ct table] 

</code>

\section pkgetmask_description Description

The utility pkgetmask creates a mask raster dataset from an input raster dataset. Values smaller than the minimum value (-min) or larger than the maximum value (-max) will result in a -nodata value in the mask.\section pkgetmask_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |Input image file | 
 | o      | output               | std::string |       |Output mask file | 
 | min    | min                  | double |       |Values smaller than min threshold(s) are masked as invalid. Use one threshold for each band | 
 | max    | max                  | double |       |Values greater than max threshold(s) are masked as invalid. Use one threshold for each band | 
 | data   | data                 | unsigned short | 1     |value(s) for valid pixels: between min and max | 
 | nodata | nodata               | unsigned short | 0     |value(s) for invalid pixels: not between min and max | 
 | b      | band                 | short | 0     |band(s) used for mask | 
 | p      | operator             | std::string | OR    |Operator: [AND,OR]. | 
 | ot     | otype                | std::string | Byte  |Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image | 
 | of     | oformat              | std::string | GTiff |Output image format (see also gdal_translate).| 
 | co     | co                   | std::string |       |Creation option for output file. Multiple options can be specified. | 
 | ct     | ct                   | std::string |       |color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid) | 

Usage: pkgetmask -i input -o output


Examples
========
Some examples how to use pkgetmask can be found \ref examples_pkgetmask "here"
**/

using namespace std;
int main(int argc,char **argv) {
  Optionpk<string> input_opt("i", "input", "Input image file");
  Optionpk<string> output_opt("o", "output", "Output mask file");
  Optionpk<string> otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image", "Byte");
  Optionpk<string> oformat_opt("of", "oformat", "Output image format (see also gdal_translate).","GTiff");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");

  oformat_opt.setHide(1);
  option_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(1);
  }
  app::AppFactory app(argc,argv);
  if(doProcess&&input_opt.empty()){
    std::cerr << "Error: no input file provided (use option -i). Use --help for help information" << std::endl;
    exit(1);
  }
  if(doProcess&&output_opt.empty()){
    std::cerr << "Error: no output file provided (use option -o). Use --help for help information" << std::endl;
    exit(1);
  }
  try{
    ImgRasterGdal imgRaster;
    if(input_opt.size())
      imgRaster.open(input_opt[0]);
    string imageType;
    if(oformat_opt.size())//default
      imageType=oformat_opt[0];
    else
      imageType=imgRaster.getImageType();
    ImgRasterGdal imgWriter;
    if(output_opt.size())
      imgWriter.setFile(output_opt[0],imageType,option_opt);
    imgRaster.getMask(imgWriter,app);
    imgWriter.close();
  }
  catch(string helpString){//help was invoked
    std::cout << helpString << std::endl;
    cout << endl;
    cout << "Usage: pkgetmask -i input -o output" << endl;
    return(1);
  }
  return(0);
}
