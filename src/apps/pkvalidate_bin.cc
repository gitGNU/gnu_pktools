/**********************************************************************
pkvalidate.cc: program to validate raster dataset based on reference vector dataset
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
#include "imageclasses/ImgRaster.h"
#include "imageclasses/ImgCollection.h"
#include "imageclasses/ImgReaderOgr.h"
#include "imageclasses/ImgWriterOgr.h"
#include "base/Optionpk.h"
#include "algorithms/ConfusionMatrix.h"
#include "apps/AppFactory.h"

/******************************************************************************/
/*! \page pkvalidate pkvalidate
  program to validate raster dataset based on reference vector dataset
## SYNOPSIS

<code>
  Usage: pkvalidate -i input -ref reference
</code>

<code>

  Options: [-ln layer] [-b band] [-cm] [-lr attribute] [-c name -r value]* [-nodata value]* [-m mask] [-msknodata value]*

  Advanced options:
       [-o output] [-f OGR format] [-lc attribute] [-bnd value [-hom] [-circ]]

</code>

\section pkvalidate_description Description

The utility pkvalidate validates a raster dataset based on a reference vector dataset. 
The reference vector dataset must consist of point features. Polygon features are automatically converted to the centroid points before analyzing. 

A typical use of the utility is to assess the accuracy of an input raster land cover map, based on a reference vector dataset. The reference dataset must contain an attribute (label) for each class. A confusion matrix is produced if the option -cm|--confusion is set. An output vector dataset can be written that contains the reference feature points with the extracted data value of the raster input dataset as a new attribute.

\section pkvalidate_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |Input raster dataset. | 
 | ref    | reference            | std::string |       |Reference (raster or vector) dataset | 
 | ln     | ln                   | std::string |       |Layer name(s) in sample. Leave empty to select all (for vector reference datasets only) | 
 | b      | band                 | short | 0     |Input (reference) raster band. Optionally, you can define different bands for input and reference bands respectively: -b 1 -b 0. | 
 | cm     | confusion            | bool | true |Create confusion matrix (to std out) | 
 | lr     | lref                 | std::string | label |Attribute name of the reference label | 
 | c      | class                | std::string |       |List of class names. | 
 | r      | reclass              | short |       |List of class values (use same order as in classname option). | 
 | nodata | nodata               | double |       |No data value(s) in input or reference dataset are ignored | 
 | m      | mask                 | std::string |       |Use the first band of the specified file as a validity mask. Nodata values can be set with the option msknodata. | 
 | msknodata | msknodata            | double | 0     |Mask value(s) where image is invalid. Use negative value for valid data (example: use -t -1: if only -1 is valid value) | 
 | o      | output               | std::string |       |Output dataset (optional) | 
 | f      | f                    | std::string | SQLite |OGR format for output vector | 
 | lc     | lclass               | std::string | class |Attribute name of the classified label | 
 | cmf    | cmf                  | std::string | ascii |Format for confusion matrix (ascii or latex) | 
 | cmo    | cmo                  | std::string |       |Output file for confusion matrix | 
 | se95   | se95                 | bool | false |Report standard error for 95 confidence interval | 
 | bnd    | boundary             | short | 1     |Boundary for selecting the sample (for vector reference datasets only) | 
 | hom    | homogeneous          | bool | false |Only take regions with homogeneous boundary into account (for reference datasets only) | 
 | circ   | circular             | bool | false |Use circular boundary (for vector reference datasets only) | 
 | mem    | mem                  | unsigned long int | 0 |Buffer size (in MB) to read image data blocks in memory | 

Usage: pkvalidate -i input -ref reference


Examples
========
Some examples how to use pkvalidate can be found \ref examples_pkvalidate "here"
**/

using namespace std;
using namespace app;

int main(int argc, char *argv[])
{
  Optionpk<string> input_opt("i", "input", "Input raster dataset.");
  Optionpk<unsigned long int>  memory_opt("mem", "mem", "Buffer size (in MB) to read image data blocks in memory",0,1);

 memory_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    memory_opt.retrieveOption(argc,argv);

    app::AppFactory app(argc,argv);

    if(doProcess&&input_opt.empty()){
      std::ostringstream errorStream;
      errorStream << "Error: no input file provided (use option -i). Use --help for help information" << std::endl;
      throw(errorStream.str());
    }

    ImgCollection imgCollection(input_opt.size());
    // std::shared_ptr<ImgRaster> imgWriter = std::make_shared<ImgRaster>();
    if(imgCollection.size()){
      for(int ifile=0;ifile<input_opt.size();++ifile)
        imgCollection[ifile]->open(input_opt[ifile],memory_opt[0]);

      // string imageType;
      // if(oformat_opt.size())//default
      //   imageType=oformat_opt[0];
      // else
      //   imageType=imgCollection[0]->getImageType();
      // imgWriter->setFile(output_opt[0],imageType,memory_opt[0],option_opt);
    }
    imgCollection.validate(app);

    for(int ifile=0;ifile<imgCollection.size();++ifile)
      imgCollection[ifile]->close();
    // imgWriter->close();
  }
  catch(string helpString){
    cerr << helpString << endl;
    cout << "Usage: pkvalidate -i input -ref reference" << endl;
    return(1);
  }
  return(0);
}
