/**********************************************************************
pkdumpimg.cc: program to dump image content to ascii or std out
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
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <assert.h>
#include "base/Optionpk.h"
#include "imageclasses/ImgReaderOgr.h"
#include "imageclasses/ImgRasterGdal.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/******************************************************************************/
/*! \page pkdumpimg pkdumpimg
 program to dump image content to ascii or std out
## SYNOPSIS

<code>
  Usage: pkdumpimg -i input.txt [-o output]
</code>

<code>


  Options: [-of matrix | line] [-b band] [-e vector | -ulx value -uly value -lrx value -lry value]

  Advanced options: [-dx value -dy value] [-r resampling] -srcnodata value -dstnodata value

</code>

\section pkdumpimg_description Description

The utility pkdumpimg dumps the content of a raster dataset to (standard) output (screen or filename). The default is to dump the output in matrix format. Use -of line to dump each pixel value on a separate line, preceded by its position (x and y value). You can specify a bounding box to dump with either the extent of an OGR vector dataset or via the options -ulx -uly -lrx and -lry.

\section pkdumpimg_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |input image file |
 | o      | output               | std::string |       |Output ascii file (Default is empty: use stdout |
 | of     | oformat              | std::string | matrix |Output format (matrix form or list (x,y,z) form. Default is matrix form |
 | b      | band                 | int  |       |band index to crop |
 | e      | extent               | std::string |       |get boundary from extent from polygons in vector file |
 | ulx    | ulx                  | double | 0     |Upper left x value bounding box (in geocoordinates if georef is true) |
 | uly    | uly                  | double | 0     |Upper left y value bounding box (in geocoordinates if georef is true) |
 | lrx    | lrx                  | double | 0     |Lower left x value bounding box (in geocoordinates if georef is true) |
 | lry    | lry                  | double | 0     |Lower left y value bounding box (in geocoordinates if georef is true) |
 | dx     | dx                   | double | 0     |Output resolution in x (in meter) (0.0: keep original resolution) |
 | dy     | dy                   | double | 0     |Output resolution in y (in meter) (0.0: keep original resolution) |
 | r      | resampling-method    | std::string | GRIORA_NearestNeighbour |resample: GRIORA_NearestNeighbour|GRIORA_Bilinear|GRIORA_Cubic|GRIORA_CubicSpline|GRIORA_Lanczos|GRIORA_Average|GRIORA_Average|GRIORA_Gauss (check http://www.gdal.org/gdal_8h.html#a640ada511cbddeefac67c548e009d5a. |
 | srcnodata | srcnodata            | double |       |set no data value(s) for input image |
 | dstnodata | dstnodata            | short | 0     |nodata value for output if out of bounds. |

Usage: pkdumpimg -i input.txt [-o output]


Examples
========
Some examples how to use pkdumpimg can be found \ref examples_pkdumpimg "here"
**/

using namespace std;

int main(int argc, char *argv[])
{
  Optionpk<std::string> input_opt("i","input","input image file","");
  Optionpk<string> extent_opt("e", "extent", "get boundary from extent from polygons in vector file");
  Optionpk<double> ulx_opt("ulx", "ulx", "Upper left x value bounding box (in geocoordinates if georef is true)");
  Optionpk<double> uly_opt("uly", "uly", "Upper left y value bounding box (in geocoordinates if georef is true)");
  Optionpk<double> lrx_opt("lrx", "lrx", "Lower left x value bounding box (in geocoordinates if georef is true)");
  Optionpk<double> lry_opt("lry", "lry", "Lower left y value bounding box (in geocoordinates if georef is true)");
  Optionpk<double> dx_opt("dx", "dx", "Output resolution in x (in meter) (0.0: keep original resolution)");
  Optionpk<double> dy_opt("dy", "dy", "Output resolution in y (in meter) (0.0: keep original resolution)");
  Optionpk<std::string> resample_opt("r", "r", "resample: GRIORA_NearestNeighbour|GRIORA_Bilinear|GRIORA_Cubic|GRIORA_CubicSpline|GRIORA_Lanczos|GRIORA_Average|GRIORA_Average|GRIORA_Gauss (check http://www.gdal.org/gdal_8h.html#a640ada511cbddeefac67c548e009d5a)","GRIORA_NearestNeighbour");


  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    extent_opt.retrieveOption(argc,argv);
    ulx_opt.retrieveOption(argc,argv);
    uly_opt.retrieveOption(argc,argv);
    lrx_opt.retrieveOption(argc,argv);
    lry_opt.retrieveOption(argc,argv);
    dx_opt.retrieveOption(argc,argv);
    dy_opt.retrieveOption(argc,argv);
    resample_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
  }

  app::AppFactory app(argc,argv);

  if(doProcess&&input_opt.empty()){
    std::cerr << "No input file provided (use option -i). Use --help for help information" << std::endl;
    exit(1);
  }

  try{
    //get bounding box from extentReader if defined
    double cropulx=0;
    double cropuly=0;
    double croplrx=0;
    double croplry=0;
    if(ulx_opt.size())
      cropulx=ulx_opt[0];
    if(uly_opt.size())
      cropuly=uly_opt[0];
    if(lrx_opt.size())
      croplrx=lrx_opt[0];
    if(lry_opt.size())
      croplry=lry_opt[0];
    if(extent_opt.size()){
      ImgReaderOgr extentReader;
      for(int iextent=0;iextent<extent_opt.size();++iextent){
        extentReader.open(extent_opt[iextent]);
        if(!(extentReader.getExtent(ulx_opt[0],uly_opt[0],lrx_opt[0],lry_opt[0]))){
          cerr << "Error: could not get extent from " << extent_opt[0] << endl;
          exit(1);
        }
        if(ulx_opt.size())
          if(ulx_opt[0]<cropulx)
            cropulx=ulx_opt[0];
        if(uly_opt.size())
          if(uly_opt[0]>cropuly)
            cropuly=uly_opt[0];
        if(lry_opt.size())
          if(lry_opt[0]<croplry)
            croplry=lry_opt[0];
        if(lrx_opt.size())
          if(lrx_opt[0]>croplrx)
            croplrx=lrx_opt[0];
      }
      extentReader.close();
    }
    if(ulx_opt.size()&&uly_opt.size()&&lrx_opt.size()&&lry_opt.size()){
      app.setOption("ulx",cropulx);
      app.setOption("uly",cropuly);
      app.setOption("lrx",croplrx);
      app.setOption("lry",croplry);
    }
    ImgRasterGdal imgReader(input_opt[0]);
    imgReader.dumpimg(app);
    imgReader.close();
  }
  catch(string helpString){//help was invoked
    std::cout << helpString << std::endl;
    cout << "Usage: pkdumpimg -i input.txt [-o output]" << endl;
    return(1);
  }
  return(0);
}
