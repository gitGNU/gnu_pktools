/**********************************************************************
pkextractogr_bin.cc: extract pixel values from raster image from a vector sample
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
/*! \page pkextractogr pkextractogr
 extract pixel values from raster image using a vector dataset sample
## SYNOPSIS

<code>
  Usage: pkextractogr -i input [-s sample | -rand number | -grid size] -o output
</code>

<code>

  Options: [-ln layer]* [-c class]* [-t threshold]* [-f format] [-ft fieldType] [-b band]* [-r rule]*

  Advanced options:
  [-sband band -eband band]* [-bndnodata band -srcnodata value]* [-tp threshold] [-buf value [-circ]]
</code>

\section pkextractogr_description Description

The utility pkextractogr extracts pixel values from an input raster dataset, based on the locations you provide via a sample file. Alternatively, a random sample or systematic grid of points can also be extracted. The sample can be a vector file with points or polygons. In the case of polygons, you can either extract the values for all raster pixels that are covered by the polygons, or extract a single value for each polygon such as the centroid, mean, median, etc. As output, a new copy of the vector file is created with an extra attribute for the extracted pixel value. For each raster band in the input image, a separate attribute is created. For instance, if the raster dataset contains three bands, three attributes are created (b0, b1 and b2). 

A typical usage of pkextractogr is to prepare a training sample for one of the classifiers implemented in pktools.

\anchor pkextractogr_rules 

Overview of the possible extraction rules:

\section pkextractogr_rules Extraction rules:

extraction rule | output features
--------------- | ---------------
point | extract a single pixel within the polygon or on each point feature
allpoints | Extract all pixel values covered by the polygon
centroid | Extract pixel value at the centroid of the polygon
mean | Extract average of all pixel values within the polygon
stdev | Extract standard deviation of all pixel values within the polygon
median | Extract median of all pixel values within the polygon
min | Extract minimum value of all pixels within the polygon
max | Extract maximum value of all pixels within the polygon
sum | Extract sum of the values of all pixels within the polygon
mode | Extract the mode of classes within the polygon (classes must be set with the option class)
proportion | Extract proportion of class(es) within the polygon (classes must be set with the option class)
count | Extract count of class(es) within the polygon (classes must be set with the option class).
percentile | Extract percentile as defined by option perc (e.g, 95th percentile of values covered by polygon)

\section pkextractogr_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |Raster input dataset containing band information | 
 | s      | sample               | std::string |       |OGR vector dataset with features to be extracted from input data. Output will contain features with input band information included | 
 | ln     | ln                   | std::string |       |Layer name(s) in sample (leave empty to select all) | 
 | rand   | random               | unsigned int |       |Create simple random sample of points. Provide number of points to generate | 
 | grid   | grid                 | double |       |Create systematic grid of points. Provide cell grid size (in projected units, e.g,. m) | 
 | o      | output               | std::string |       |Output sample dataset | 
 | c      | class                | int  |       |Class(es) to extract from input sample image. Leave empty to extract all valid data pixels from sample dataset. Make sure to set classes if rule is set to mode, proportion or count | 
 | t      | threshold            | float | 100   |Probability threshold for selecting samples (randomly). Provide probability in percentage (>0) or absolute (<0). Use a single threshold per vector sample layer | 
 | perc   | perc                 | double | 95    |Percentile value(s) used for rule percentile | 
 | f      | f                    | std::string | SQLite |Output sample dataset format | 
 | ft     | ftype                | std::string | Real  |Field type (only Real or Integer) | 
 | b      | band                 | int  |       |Band index(es) to extract (0 based). Leave empty to use all bands | 
 | sband  | startband            | unsigned short |      |Start band sequence number | 
 | eband  | endband              | unsigned short |      |End band sequence number   | 
 | r      | rule                 | std::string | centroid |Rule how to report image information per feature (only for vector sample). point (single point or at centroid if polygon), allpoints (within polygon), centroid, mean, stdev, median, proportion, count, min, max, mode, sum, percentile. | 
 | bndnodata | bndnodata            | int  | 0     |Band(s) in input image to check if pixel is valid (used for srcnodata) | 
 | srcnodata | srcnodata            | double |       |Invalid value(s) for input image | 
 | tp     | thresholdPolygon     | float |       |(absolute) threshold for selecting samples in each polygon | 
 | buf    | buffer               | short | 0     |Buffer for calculating statistics for point features (in number of pixels)  | 
 | circ   | circular             | bool | false |Use a circular disc kernel buffer (for vector point sample datasets only, use in combination with buffer option) | 

Usage: pkextractogr -i input [-s sample | -rand number | -grid size] -o output -r rule


Examples
========
Some examples how to use pkextractogr can be found \ref examples_pkextractogr "here"
**/

namespace rule{
  enum RULE_TYPE {point=0, mean=1, proportion=2, custom=3, min=4, max=5, mode=6, centroid=7, sum=8, median=9, stdev=10, percentile=11, count=12, allpoints=13};
}

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
    imgRaster.extractOgr(app);
    imgRaster.close();
  }
  catch(string helpString){
    cerr << helpString << endl;
    cout << "Usage: pkextractogr -i input [-s sample | -rand number | -grid size] -o output" << endl;
    return(1);
  }
  return(0);
}
