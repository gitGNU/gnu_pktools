/**********************************************************************
pkdiff_bin.cc: cprogram to compare two raster image file
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
#include <memory>
#include "imageclasses/ImgRasterGdal.h"
#include "base/Optionpk.h"
#include "AppFactory.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

using namespace std;
/******************************************************************************/
/*! \page pkdiff pkdiff
 program to compare two raster image files
## SYNOPSIS

<code>
  Usage: pkdiff -i input -ref reference
</code>

<code>

  Options: [-nodata value]* [-m mask] [-msknodata value]*

  Advanced options:
       [-o output] [-bnd value [-ct colortable] [-co NAME=VALUE]* 

</code>

\section pkdiff_description Description

The utility pkdiff compares two raster datasets, performing a pixel by pixel comparison. With no further options, the utility reports if the rasters are identical or different. The root mean squared error (-rmse) or regression (-reg) can also be calculated. If required, an output raster dataset can be written with a qualitative information per pixel: 0 (input=reference), 1 (input>reference) or 2 (input<reference). 

A particular use of the utility is to assess the accuracy of an input raster land cover map, based on a reference raster dataset (use pkvalidate to use a vector dataset as reference). A confusion matrix is produced if the option -cm|--confusion is set.
\section pkdiff_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |Input raster dataset. | 
 | ref    | reference            | std::string |       |Reference (raster or vector) dataset | 
 | b      | band                 | short | 0     |Input (reference) raster band. Optionally, you can define different bands for input and reference bands respectively: -b 1 -b 0. | 
 | rmse   | rmse                 | bool | false |Report root mean squared error | 
 | reg    | reg                  | bool | false |Report linear regression (Input = c0+c1*Reference) | 
 | cm     | confusion            | bool | false |Create confusion matrix (to std out) | 
 | c      | class                | std::string |       |List of class names. | 
 | r      | reclass              | short |       |List of class values (use same order as in classname option). | 
 | nodata | nodata               | double |       |No data value(s) in input or reference dataset are ignored | 
 | m      | mask                 | std::string |       |Use the first band of the specified file as a validity mask. Nodata values can be set with the option msknodata. | 
 | msknodata | msknodata            | double | 0     |Mask value(s) where image is invalid. Use negative value for valid data (example: use -t -1: if only -1 is valid value) | 
 | o      | output               | std::string |       |Output dataset (optional) | 
 | of     | oformat              | std::string | GTiff |Output image format (see also gdal_translate).| 
 | cmf    | cmf                  | std::string | ascii |Format for confusion matrix (ascii or latex) | 
 | cmo    | cmo                  | std::string |       |Output file for confusion matrix | 
 | se95   | se95                 | bool | false |Report standard error for 95 confidence interval | 
 | ct     | ct                   | std::string |       |Color table in ASCII format having 5 columns: id R G B ALFA (0: transparent, 255: solid). | 
 | co     | co                   | std::string |       |Creation option for output file. Multiple options can be specified. | 
 |        | commission           | short | 2     |Value for commission errors: input label < reference label | 

Usage: pkdiff -i input -ref reference


Examples
========
Some examples how to use pkdiff can be found \ref examples_pkdiff "here"
**/

int main(int argc, char *argv[])
{
  vector<double> priors;

  //--------------------------- command line options ------------------------------------
  Optionpk<string> input_opt("i", "input", "input image");
  Optionpk<string> reference_opt("ref", "reference", "Reference (raster or vector) dataset");

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    reference_opt.retrieveOption(argc,argv);

    app::AppFactory app(argc,argv);

    if(doProcess&&input_opt.empty()){
      std::ostringstream errorStream;
      errorStream << "Error: no input file provided (use option -i). Use --help for help information" << std::endl;
      throw(errorStream.str());
    }
    if(doProcess&&reference_opt.empty()){
      std::ostringstream errorStream;
      errorStream << "Error: no reference file provided (use option -ref). Use --help for help information" << std::endl;
      throw(errorStream.str());
    }
    ImgRasterGdal imgRaster;
    if(input_opt.size())
      imgRaster.open(input_opt[0]);
    ImgRasterGdal imgReference;
    if(reference_opt.size())
      imgReference.open(reference_opt[0]);

    imgRaster.diff(imgReference,app);
    imgRaster.close();
    imgReference.close();
  }
  catch(string helpString){
    cerr << helpString << endl;
    cout << "Usage: pkdiff -r reference [-cm] [-o output]" << endl;
    return(1);
  }
  return(0);
}
