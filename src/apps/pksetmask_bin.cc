/**********************************************************************
pksetmask_bin.cc: program to apply mask image (set invalid values) to raster image
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
/*! \page pksetmask pksetmask
 program to apply mask image (set invalid values) to raster image
## SYNOPSIS

<code>
  Usage: pksetmask -i input -m mask [-m mask]* -o output
</code>

<code>
  
  Options: [-mskband value]* [-msknodata value -nodata value]*
 
  Advanced options: [--operator '<'|'='|'<'] [-ot type] [-of format] [-co option]* [-ct table] 

</code>

\section pksetmask_description Description

The utility pksetmask sets a mask provided with option -m to an input raster dataset. The default operator is '='. Values in the input raster data where the mask has a nodata value (set with the option -msknodata) will then be set to nodata (set with -nodata). Other operators are less than (--operator '<') and larger than (--operator '>').\section pksetmask_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |Input image | 
 | m      | mask                 | std::string |       |Mask image(s) | 
 | msknodata | msknodata            | int  | 1     |Mask value(s) where image has nodata. Use one value for each mask, or multiple values for a single mask. | 
 | mskband | mskband              | short | 0     |Mask band to read (0 indexed). Provide band for each mask. | 
 | o      | output               | std::string |       |Output mask file | 
 | nodata | nodata               | int  | 0     |nodata value to put in image if not valid | 
 | p      | operator             | char | =     |Operator: < = > !. Use operator for each msknodata option | 
 | ot     | otype                | std::string |       |Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image | 
 | of     | oformat              | std::string | GTiff |Output image format (see also gdal_translate)| 
 | co     | co                   | std::string |       |Creation option for output file. Multiple options can be specified. | 
 | ct     | ct                   | std::string |       |color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid) | 

Usage: pksetmask -i input -m mask [-m mask]* [-msknodata value -nodata value]* -o output

Examples
========
Some examples how to use pksetmask can be found \ref examples_pksetmask "here"
FAQ
========
Frequently asked questions on pksetmask can be found \ref faq_pksetmask "here"
**/

using namespace std;

int main(int argc, char *argv[])
{
  //command line options
  Optionpk<string> input_opt("i", "input", "Input image");
  Optionpk<string> output_opt("o", "output", "Output mask file");
  Optionpk<string> oformat_opt("of", "oformat", "Output image format (see also gdal_translate)","GTiff");
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
    imgRaster.setMask(imgWriter,app);
    imgWriter.close();
  }
  catch(string helpString){//help was invoked
    std::cout << helpString << std::endl;
    cout << "Usage: pksetmask -i input -m mask [-m mask]* [-msknodata value -nodata value]* -o output" << endl;
    return(1);
  }
  return(0);
}
