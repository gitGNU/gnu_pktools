/**********************************************************************
pkcrop_bin.cc: perform raster data operations on image such as crop, extract and stack bands
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
#include <algorithm>
#include <vector>
#include <iostream>
#include <string>
#include "imageclasses/ImgRaster.h"
#include "imageclasses/ImgCollection.h"
#include "imageclasses/ImgReaderOgr.h"
#include "base/Vector2d.h"
#include "base/Optionpk.h"
#include "AppFactory.h"
/******************************************************************************/
/*! \page pkcrop pkcrop
 perform raster data operations on image such as crop, extract and stack bands
## SYNOPSIS

<code>
  Usage: pkcrop -i input -o output
</code>

<code>

  Options: [-of out_format] [-ot {Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}] [-b band]* [-ulx ULX -uly ULY -lrx LRX -lry LRY] [-dx xres] [-dy yres] [-r resampling_method] [-a_srs epsg:number] [-nodata value] 

  Advanced options:
  	   [-e vector [-cut] [-eo option]*] [-sband band -eband band]* [-co NAME=VALUE]* [-x center_x -y center_y] [-nx size_x -ny size_y] [-ns nsample -nl nlines] [-as min -as max] [-scale value]* [-off offset]* [-ct colortable] [-d description] [-align]
</code>

\section pkcrop_description Description

The utility pkcrop can subset and stack raster images. In the spatial domain it can crop a bounding box from a larger image. The output bounding box is selected by setting the new corner coordinates using the options -ulx -uly -lrx -lry. Alternatively you can set the new image center (-x -y) and size. This can be done either in projected coordinates (using the options -nx -ny) or in image coordinates (using the options -ns -nl). You can also use a vector file to set the new bounding box (option -e). In the spectral domain, pkcrop allows you to select individual bands from one or more input image(s). Bands are stored in the same order as provided on the command line, using the option -b. Band numbers start with index 0 (indicating the first band). The default is to select all input bands. If more input images are provided, the bands are stacked into a multi-band image. If the bounding boxes or spatial resolution are not identical for all input images, you should explicitly set them via the options. The pkcrop utility is not suitable to mosaic or composite images. Consider the utility pkcomposite instead.

\section pkcrop_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |Input image file(s). If input contains multiple images, a multi-band output is created | 
 | o      | output               | std::string |       |Output image file | 
 | a_srs  | a_srs                | std::string |       |Override the projection for the output file (leave blank to copy from input file, use epsg:3035 to use European projection and force to European grid | 
 | ulx    | ulx                  | double | 0     |Upper left x value bounding box | 
 | uly    | uly                  | double | 0     |Upper left y value bounding box | 
 | lrx    | lrx                  | double | 0     |Lower right x value bounding box | 
 | lry    | lry                  | double | 0     |Lower right y value bounding box | 
 | b      | band                 | unsigned int |       |band index to crop (leave empty to retain all bands) | 
 | sband  | startband            | unsigned int |      |Start band sequence number | 
 | eband  | endband              | unsigned int |      |End band sequence number   | 
 | as     | autoscale            | double |       |scale output to min and max, e.g., --autoscale 0 --autoscale 255 | 
 | ot     | otype                | std::string |       |Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image | 
 | of     | oformat              | std::string | GTiff |Output image format (see also gdal_translate)| 
 | ct     | ct                   | std::string |       |color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid) | 
 | dx     | dx                   | double |       |Output resolution in x (in meter) (empty: keep original resolution) | 
 | dy     | dy                   | double |       |Output resolution in y (in meter) (empty: keep original resolution) | 
 | r      | resampling-method    | std::string | near  |Resampling method (near: nearest neighbor, bilinear: bi-linear interpolation). | 
 | e      | extent               | std::string |       |get boundary from extent from polygons in vector file | 
 | cut      | crop_to_cutline    | bool | false |Crop the extent of the target dataset to the extent of the cutline | 
 | eo       | eo                 | std::string |       |special extent options controlling rasterization: ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG, e.g., -eo ATTRIBUTE=fieldname |
 | m      | mask                 | std::string |       |Use the specified file as a validity mask (0 is nodata) | 
 | msknodata | msknodata            | float | 0     |Mask value not to consider for crop
 | mskband | mskband              | short | 0     |Mask band to read (0 indexed) | 
 | co     | co                   | std::string |       |Creation option for output file. Multiple options can be specified. | 
 | x      | x                    | double |       |x-coordinate of image center to crop (in meter) | 
 | y      | y                    | double |       |y-coordinate of image center to crop (in meter) | 
 | nx     | nx                   | double |       |image size in x to crop (in meter) | 
 | ny     | ny                   | double |       |image size in y to crop (in meter) | 
 | ns     | ns                   | int  |       |number of samples  to crop (in pixels) | 
 | nl     | nl                   | int  |       |number of lines to crop (in pixels) | 
 | scale  | scale                | double |       |output=scale*input+offset | 
 | off    | offset               | double |       |output=scale*input+offset | 
 | nodata | nodata               | float |       |Nodata value to put in image if out of bounds. | 
 | align  | align                | bool  |       |Align output bounding box to input image | 
 | mem    | mem                  | unsigned long int | 0 |Buffer size (in MB) to read image data blocks in memory | 
 | d      | description          | std::string |       |Set image description | 

Examples
========
Some examples how to use pkcrop can be found \ref examples_pkcrop "here"
**/

using namespace std;

int main(int argc, char *argv[])
{
  Optionpk<string>  input_opt("i", "input", "Input image file(s). If input contains multiple images, a multi-band output is created");
  Optionpk<string>  output_opt("o", "output", "Output image file");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate).","GTiff");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<string>  projection_opt("a_srs", "a_srs", "Override the spatial reference for the output file (leave blank to copy from input file, use epsg:3035 to use European projection and force to European grid");
  Optionpk<double> scale_opt("scale", "scale", "output=scale*input+offset");
  Optionpk<double> offset_opt("offset", "offset", "output=scale*input+offset");
  Optionpk<unsigned long int>  memory_opt("mem", "mem", "Buffer size (in MB) to read image data blocks in memory",0,1);

  option_opt.setHide(1);
  scale_opt.setHide(1);
  offset_opt.setHide(1);
  memory_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    projection_opt.retrieveOption(argc,argv);
    scale_opt.retrieveOption(argc,argv);
    offset_opt.retrieveOption(argc,argv);
    memory_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
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
    ImgCollection imgCollection(input_opt.size());
    std::shared_ptr<ImgRaster> imgWriter = std::make_shared<ImgRaster>();
    if(imgCollection.size()){
      for(int ifile=0;ifile<input_opt.size();++ifile){
        imgCollection[ifile]->open(input_opt[ifile],memory_opt[0]);
        for(int iband=0;iband<scale_opt.size();++iband)
          imgCollection[ifile]->setScale(scale_opt[iband],iband);
        for(int iband=0;iband<offset_opt.size();++iband)
          imgCollection[ifile]->setOffset(offset_opt[iband],iband);
      }

      string imageType;
      if(oformat_opt.size())//default
        imageType=oformat_opt[0];
      else
        imageType=imgCollection[0]->getImageType();
      if(output_opt.size())
        imgWriter->setFile(output_opt[0],imageType,memory_opt[0],option_opt);
    }
    imgCollection.crop(imgWriter,app);

    for(int ifile=0;ifile<imgCollection.size();++ifile)
      imgCollection[ifile]->close();
    imgWriter->close();
  }
  catch(string helpString){//help was invoked
    std::cout << helpString << std::endl;
    cout << "Usage: pkcrop -i input [-i input]* -o output" << endl;
    return(1);
  }
  return(0);
}
  
