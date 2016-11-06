/**********************************************************************
pkann_bin.cc: classify raster image using Support Vector Machine
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
/******************************************************************************/
/*! \page pkann pkann
 classify raster image using Artificial Neural Network
## SYNOPSIS

<code>
  Usage: pkann -t training [-i input -o output] [-cv value]
</code>

<code>

  
  Options: [-tln layer]* [-c name -r value]* [-of GDALformat|-f OGRformat] [-co NAME=VALUE]* [-ct filename] [-label attribute] [-prior value]* [-nn number]* [-m filename [-msknodata value]*] [-nodata value]

  Advanced options:
       [-b band] [-sband band -eband band]* [-bal size]* [-min] [-bag value] [-bs value] [-comb rule] [-cb filename] [-prob filename] [-pim priorimage] [--offset value] [--scale value] [--connection 0|1] [-w weights]* [--learning rate] [--maxit number]
</code>

\section pkann_description Description

The utility pkann implements an artificial neural network (ANN) to solve a supervised classification problem. The implementation is based on the open source C++ library <a href="http://leenissen.dk/fann/wp/">fann</a>). Both raster and vector files are supported as input. The output will contain the classification result, either in raster or vector format, corresponding to the format of the input. A training sample must be provided as an OGR vector dataset that contains the class labels and the features for each training point. The point locations are not considered in the training step. You can use the same training sample for classifying different images, provided the number of bands of the images are identical. Use the utility \ref pkextract "pkextract" to create a suitable training sample, based on a sample of points or polygons. For raster output maps you can attach a color table using the option -ct.

\section pkann_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |input image | 
 | t      | training             | std::string |       |training vector file. A single vector file contains all training features (must be set as: B0, B1, B2,...) for all classes (class numbers identified by label option). Use multiple training files for bootstrap aggregation (alternative to the bag and bsize options, where a random subset is taken from a single training file) | 
 | tln    | tln                  | std::string |       |training layer name(s) | 
 | label  | label                | std::string | label |identifier for class label in training vector file. | 
 | bal    | balance              | unsigned int | 0     |balance the input data to this number of samples for each class | 
 | min    | min                  | int  | 0     |if number of training pixels is less then min, do not take this class into account (0: consider all classes) | 
 | b      | band                 | short |       |band index (starting from 0, either use band option or use start to end) | 
 | sband  | startband            | unsigned short |      |Start band sequence number | 
 | eband  | endband              | unsigned short |      |End band sequence number   | 
 |        | offset               | double | 0     |offset value for each spectral band input features: refl[band]=(DN[band]-offset[band])/scale[band] | 
 | scale  | scale                | double | 0     |scale value for each spectral band input features: refl=(DN[band]-offset[band])/scale[band] (use 0 if scale min and max in each band to -1.0 and 1.0) | 
 | a      | aggreg               | unsigned short | 1     |how to combine aggregated classifiers, see also rc option (1: sum rule, 2: max rule). | 
 | prior  | prior                | double | 0     |prior probabilities for each class (e.g., -p 0.3 -p 0.3 -p 0.2 ) | 
 | pim    | priorimg             | std::string |       |prior probability image (multi-band img with band for each class | 
 | cv     | cv                   | unsigned short | 0     |n-fold cross validation mode | 
 | cmf    | cmf                  | std::string | ascii |Format for confusion matrix (ascii or latex) | 
 | nn     | nneuron              | unsigned int | 5     |number of neurons in hidden layers in neural network (multiple hidden layers are set by defining multiple number of neurons: -n 15 -n 1, default is one hidden layer with 5 neurons) | 
 |        | connection           | float | 1     |connection reate (default: 1.0 for a fully connected network) | 
 | w      | weights              | float | 0     |weights for neural network. Apply to fully connected network only, starting from first input neuron to last output neuron, including the bias neurons (last neuron in each but last layer) | 
 | l      | learning             | float | 0.7   |learning rate (default: 0.7) | 
 |        | maxit                | unsigned int | 500   |number of maximum iterations (epoch) (default: 500) | 
 | comb   | comb                 | unsigned short | 0     |how to combine bootstrap aggregation classifiers (0: sum rule, 1: product rule, 2: max rule). Also used to aggregate classes with rc option. Default is sum rule (0) | 
 | bag    | bag                  | unsigned short | 1     |Number of bootstrap aggregations (default is no bagging: 1) | 
 | bs     | bsize                | int  | 100   |Percentage of features used from available training features for each bootstrap aggregation (one size for all classes, or a different size for each class respectively | 
 | cb     | classbag             | std::string |       |output for each individual bootstrap aggregation (default is blank) | 
 | m      | mask                 | std::string |       |Only classify within specified mask (vector or raster). For raster mask, set nodata values with the option msknodata. | 
 | msknodata | msknodata            | short | 0     |mask value(s) not to consider for classification. Values will be taken over in classification image. Default is 0 | 
 | nodata | nodata               | unsigned short | 0     |nodata value to put where image is masked as nodata | 
 | o      | output               | std::string |       |output classification image | 
 | ot     | otype                | std::string |       |Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image | 
 | of     | oformat              | std::string |       |Output image format (see also gdal_translate). Empty string: inherit from input image | 
 | ct     | ct                   | std::string |       |colour table in ASCII format having 5 columns: id R G B ALFA (0: transparent, 255: solid) | 
 | co     | co                   | std::string |       |Creation option for output file. Multiple options can be specified. | 
 |        | prob                 | std::string |       |probability image. Default is no probability image | 
 | f      | f                    | std::string | SQLite |Output ogr format for active training sample | 
 | na     | nactive              | unsigned int | 1     |number of active training points | 
 | c      | class                | std::string |       |list of class names. | 
 | r      | reclass              | short |       |list of class values (use same order as in class opt). | 

Usage: pkann -t training [-i input -o output] [-cv value]


Examples
========
Some examples how to use pkann can be found \ref examples_pkann "here"
**/

using namespace std;

int main(int argc, char *argv[])
{
  vector<double> priors;

  //--------------------------- command line options ------------------------------------
  Optionpk<string> input_opt("i", "input", "input image");
  Optionpk<string> output_opt("o", "output", "Output classification image");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate).","GTiff");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");

  oformat_opt.setHide(1);
  option_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    app::AppFactory app(argc,argv);

    ImgRasterGdal imgRaster;
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

    imgRaster.ann(imgWriter,app);
    imgRaster.close();
    imgWriter.close();
  }
  catch(string helpString){
    cerr << helpString << endl;
    cout << "Usage: pkann -t training [-i input [-o output]] [-cv value]" << endl;
    return(1);
  }
  return(0);
}
