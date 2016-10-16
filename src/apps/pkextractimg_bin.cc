/**********************************************************************
pkextractimg_bin.cc: extract pixel values from raster image using a raster sample
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
#include "imageclasses/ImgRasterGdal.h"
#include "imageclasses/ImgWriterOgr.h"
#include "base/Optionpk.h"
#include "AppFactory.h"

/******************************************************************************/
/*! \page pkextractimg pkextractimg
 extract pixel values from raster image using a raster sample
## SYNOPSIS

<code>
  Usage: pkextractimg -i input -s sample -o output
</code>

<code>

  Options: [-c class]* [-t threshold]* [-f format] [-ft fieldType] [-lt labelType] [-b band]*

  Advanced options:
  [-sband band -eband band]* [-bndnodata band [-srcnodata value]*] [-bn attribute] [-cn attribute] [-down value]
</code>

\section pkextractimg_description Description

The utility pkextractimg extracts pixel values from an input raster dataset, based on the locations you provide via a sample file. The sample should be a raster dataset with categorical (integer) values. The typical use case is a land cover map that overlaps the input raster dataset. The utility then extracts pixels from the input raster for the respective land cover classes. To select a random subset of the sample raster dataset you can set the threshold option -t with a percentage value. You can provide a threshold value for each class (e.g. -t 80 -t 60). Use value 100 to select all pixels for selected class(es). As output, a new copy of the vector file is created with an extra attribute for the extracted pixel value. For each raster band in the input image, a separate attribute is created. For instance, if the raster dataset contains three bands, three attributes are created (b0, b1 and b2). 

\section pkextractimg_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |Raster input dataset containing band information | 
 | s      | sample               | std::string |       |Raster dataset with categorical values to sample the input raster dataset. Output will contain features with input band information included | 
 | o      | output               | std::string |       |Output sample dataset | 
 | c      | class                | int  |       |Class(es) to extract from input sample image. Leave empty to extract all valid data pixels from sample dataset | 
 | t      | threshold            | float | 100   |Probability threshold for selecting samples (randomly). Provide probability in percentage (>0) or absolute (<0). You can provide a threshold value for each class (e.g. -t 80 -t 60). Use value 100 to select all pixels for selected class(es) | 
 | f      | f                    | std::string | SQLite |Output sample dataset format | 
 | ft     | ftype                | std::string | Real  |Field type (only Real or Integer) | 
 | lt     | ltype                | std::string | Integer |Label type: In16 or String | 
 | b      | band                 | int  |       |Band index(es) to extract (0 based). Leave empty to use all bands | 
 | sband  | startband            | unsigned short |      |Start band sequence number | 
 | eband  | endband              | unsigned short |      |End band sequence number   | 
 | bndnodata | bndnodata            | int  | 0     |Band in input image to check if pixel is valid (used for srcnodata) | 
 | srcnodata | srcnodata            | double |       |Invalid value(s) for input image | 
 | bn     | bname                | std::string | b     |For single band input data, this extra attribute name will correspond to the raster values. For multi-band input data, multiple attributes with this prefix will be added (e.g. b0, b1, b2, etc.) | 
 | cn     | cname                | std::string | label |Name of the class label in the output vector dataset | 
 | down   | down                 | short | 1     |Down sampling factor | 

Usage: pkextractimg -i input [-s sample] -o output


Examples
========
Some examples how to use pkextractimg can be found \ref examples_pkextractimg "here"
**/

using namespace std;
using namespace app;

int main(int argc, char *argv[])
{
  Optionpk<string> image_opt("i", "input", "Raster input dataset containing band information");

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=image_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
  }

  app::AppFactory app(argc,argv);

  if(doProcess&&image_opt.empty()){
    std::cerr << "Error: no input file provided (use option -i). Use --help for help information" << std::endl;
    exit(1);
  }
  try{
    ImgRasterGdal imgRaster;
    if(image_opt.size()){
      imgRaster.open(image_opt[0]);
    }
    imgRaster.extractImg(app);
    imgRaster.close();
  }
  catch(string helpString){
    cerr << helpString << endl;
    cout << "Usage: pkextractImg -i input -s sample -o output" << endl;
    return(1);
  }
  return(0);
}
