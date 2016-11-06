/**********************************************************************
pkreclass_bin.cc: program to replace pixel values in raster image
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
#include <assert.h>
#include <map>
#include "base/Optionpk.h"
#include "imageclasses/ImgReaderOgr.h"
#include "imageclasses/ImgWriterOgr.h"
#include "imageclasses/ImgRasterGdal.h"
#include "base/Optionpk.h"
#include "AppFactory.h"
/******************************************************************************/
/*! \page pkreclass pkreclass
  program to replace pixel values in raster image
  ## SYNOPSIS

  <code>
  Usage: pkreclass -i input [-c from -r to]* -o output
  </code>

  \section pkreclass_options Options
  - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
  - short option `-h` shows basic options only, long option `--help` shows all options
  |short|long|type|default|description|
  |-----|----|----|-------|-----------|
  | i      | input                | std::string |       |Input image |
  | m      | mask                 | std::string |       |Mask image(s) |
  | msknodata | msknodata            | unsigned short | 1     |Mask value(s) where image has nodata. Use one value for each mask, or multiple values for a single mask. |
  | nodata | nodata               | int  | 0     |nodata value to put in image if not valid (0) |
  | code   | code                 | std::string |       |Recode text file (2 colums: from to) |
  | c      | class                | std::string |       |list of classes to reclass (in combination with reclass option) |
  | r      | reclass              | std::string |       |list of recoded classes (in combination with class option) |
  | ct     | ct                   | std::string |       |color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid) |
  | o      | output               | std::string |       |Output mask file |
  | ot     | otype                | std::string |       |Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image |
  | of     | oformat              | std::string | GTiff |Output image format (see also gdal_translate)|
  | b      | band                 | unsigned short | 0     |band index(es) to replace (other bands are copied to output) |
  | n      | fname                | std::string | label |field name of the shape file to be replaced |
  | co     | co                   | std::string |       |Creation option for output file. Multiple options can be specified. |
  | d      | description          | std::string |       |Set image description |
  | v      | verbose              | short | 0     |verbose |

  Usage: pkreclass -i input [-c from -r to]* -o output


  Examples
  ========
  Some examples how to use pkreclass can be found \ref examples_pkreclass "here"
**/

using namespace std;

int main(int argc, char *argv[])
{
  Optionpk<string> input_opt("i", "input", "Input image");
  Optionpk<string> output_opt("o", "output", "Output mask file");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate).","GTiff");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    app::AppFactory app(argc,argv);

    if(doProcess&&input_opt.empty()){
      std::cerr << "Error: no input file provided (use option -i). Use --help for help information" << std::endl;
      exit(1);
    }
    if(doProcess&&output_opt.empty()){
      std::cerr << "Error: no output file provided (use option -o). Use --help for help information" << std::endl;
      exit(1);
    }

    if(input_opt.empty()){
      std::cerr << "No input file provided (use option -i). Use --help for help information" << std::endl;
      exit(0);
    }
    if(output_opt.empty()){
      std::cerr << "No output file provided (use option -o). Use --help for help information" << std::endl;
      exit(0);
    }
    ImgRasterGdal imgRaster;
    //test
    std::cout << "opening: " << input_opt[0] << std::endl;
    if(input_opt.size())
      imgRaster.open(input_opt[0]);
    ImgRasterGdal imgWriter;
    string imageType;
    if(oformat_opt.size())//default
      imageType=oformat_opt[0];
    else
      imageType=imgRaster.getImageType();

    if(output_opt.size())
      imgWriter.setFile(output_opt[0],imageType,option_opt);
    imgRaster.reclass(imgWriter,app);
    imgRaster.close();
    imgWriter.close();
  }
  catch(string helpString){//help was invoked
    std::cout << helpString << std::endl;
    cout << "Usage: pkreclass -i input [-c from -r to]* -o output" << endl;
    return(1);
  }
  return(0);
}
