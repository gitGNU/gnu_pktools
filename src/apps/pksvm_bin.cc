/**********************************************************************
pksvm_bin.cc: classify raster image using Support Vector Machine
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
#include "imageclasses/ImgRaster.h"
#include "base/Optionpk.h"
#include "AppFactory.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
/******************************************************************************/
/*! \page pksvm pksvm
  classify raster image using Support Vector Machine
  ## SYNOPSIS

  <code>
  Usage: pksvm -t training [-i input -o output] [-cv value]
  </code>

  <code>

  Options: [-tln layer]* [-c name -r value]* [-of GDALformat|-f OGRformat] [-co NAME=VALUE]* [-ct filename] [-label attribute] [-prior value]* [-g gamma] [-cc cost] [-m filename [-msknodata value]*] [-nodata value]

  Advanced options:
  [-b band] [-sband band -eband band]* [-bal size]* [-min] [-bag value] [-bs value] [-comb rule] [-cb filename] [-prob filename] [-pim priorimage] [--offset value] [--scale value] [-svmt type] [-kt type] [-kd value]  [-c0 value] [-nu value] [-eloss value] [-cache value] [-etol value] [-shrink]
  </code>

  \section pksvm_description Description

  The utility pksvm implements a support vector machine (SVM) to solve a supervised classification problem. The implementation is based on the open source C++ library libSVM (http://www.csie.ntu.edu.tw/~cjlin/libsvm).
  Both raster and vector files are supported as input. The output will contain the classification result, either in raster or vector format, corresponding to the format of the input. A training sample must be provided as an OGR vector dataset that contains the class labels and the features for each training point. The point locations are not considered in the training step. You can use the same training sample for classifying different images, provided the number of bands of the images are identical. Use the utility \ref pkextract "pkextract" to create a suitable training sample, based on a sample of points or polygons. For raster output maps you can attach a color table using the option -ct.

  \section pksvm_options Options
  - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
  - short option `-h` shows basic options only, long option `--help` shows all options
  |short|long|type|default|description|
  |-----|----|----|-------|-----------|
  | t      | training             | std::string |       |Training vector file. A single vector file contains all training features (must be set as: b0, b1, b2,...) for all classes (class numbers identified by label option). Use multiple training files for bootstrap aggregation (alternative to the bag and bsize options, where a random subset is taken from a single training file) |
  | i      | input                | std::string |       |input image |
  | o      | output               | std::string |       |Output classification image |
  | cv     | cv                   | unsigned short | 0     |N-fold cross validation mode |
  | cmf    | cmf                  | std::string | ascii |Format for confusion matrix (ascii or latex) |
  | tln    | tln                  | std::string |       |Training layer name(s) |
  | c      | class                | std::string |       |List of class names. |
  | r      | reclass              | short |       |List of class values (use same order as in class opt). |
  | of     | oformat              | std::string | GTiff |Output image format (see also gdal_translate).|
  | f      | f                    | std::string | SQLite |Output ogr format for active training sample |
  | co     | co                   | std::string |       |Creation option for output file. Multiple options can be specified. |
  | ct     | ct                   | std::string |       |Color table in ASCII format having 5 columns: id R G B ALFA (0: transparent, 255: solid) |
  | label  | label                | std::string | label |Attribute name for class label in training vector file. |
  | prior  | prior                | double | 0     |Prior probabilities for each class (e.g., -p 0.3 -p 0.3 -p 0.2 ). Used for input only (ignored for cross validation) |
  | g      | gamma                | float | 1     |Gamma in kernel function |
  | cc     | ccost                | float | 1000  |The parameter C of C_SVC, epsilon_SVR, and nu_SVR |
  | m      | mask                 | std::string |       |Only classify within specified mask (raster dataset). Set nodata values with the option msknodata. |
  | e      | extent               | std::string |       |get boundary from extent from polygons in vector file | 
  | eo       | eo                 | std::string |       |special extent options controlling rasterization: ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG, e.g., -eo ATTRIBUTE=fieldname |
  | msknodata | msknodata            | short | 0     |Mask value(s) not to consider for classification. Values will be taken over in classification image. |
  | nodata | nodata               | unsigned short | 0     |Nodata value to put where image is masked as nodata |
  | b      | band                 | unsigened int |       |Band index (starting from 0, either use band option or use start to end) |
  | sband  | startband            | unsigned int |      |Start band sequence number |
  | eband  | endband              | unsigned int |      |End band sequence number   |
  | bal    | balance              | unsigned int | 0     |Balance the input data to this number of samples for each class |
  | min    | min                  | int  | 0     |If number of training pixels is less then min, do not take this class into account (0: consider all classes) |
  | bag    | bag                  | unsigned short | 1     |Number of bootstrap aggregations |
  | bagsize | bagsize              | int  | 100   |Percentage of features used from available training features for each bootstrap aggregation (one size for all classes, or a different size for each class respectively |
  | comb   | comb                 | unsigned short | 0     |How to combine bootstrap aggregation classifiers (0: sum rule, 1: product rule, 2: max rule). Also used to aggregate classes with rc option. |
  | cb     | classbag             | std::string |       |Output for each individual bootstrap aggregation |
  | prob   | prob                 | std::string |       |Probability image. |
  | pim    | priorimg             | std::string |       |Prior probability image (multi-band img with band for each class |
  | offset | offset               | double | 0     |Offset value for each spectral band input features: refl[band]=(DN[band]-offset[band])/scale[band] |
  | scale  | scale                | double | 0     |Scale value for each spectral band input features: refl=(DN[band]-offset[band])/scale[band] (use 0 if scale min and max in each band to -1.0 and 1.0) |
  | svmt   | svmtype              | std::string | C_SVC |Type of SVM (C_SVC, nu_SVC,one_class, epsilon_SVR, nu_SVR) |
  | kt     | kerneltype           | std::string | radial |Type of kernel function (linear,polynomial,radial,sigmoid)  |
  | kd     | kd                   | unsigned short | 3     |Degree in kernel function |
  | c0     | coef0                | float | 0     |Coef0 in kernel function |
  | nu     | nu                   | float | 0.5   |The parameter nu of nu_SVC, one_class SVM, and nu_SVR |
  | eloss  | eloss                | float | 0.1   |The epsilon in loss function of epsilon_SVR |
  | cache  | cache                | int  | 100   |Cache memory size in MB |
  | etol   | etol                 | float | 0.001 |The tolerance of termination criterion |
  | shrink | shrink               | bool | false |Whether to use the shrinking heuristics |
  | pe     | probest              | bool | true  |Whether to train a SVC or SVR model for probability estimates |
  | entropy | entropy              | std::string |       |Entropy image (measure for uncertainty of classifier output |
  | active | active               | std::string |       |Ogr output for active training sample. |
  | na     | nactive              | unsigned int | 1     |Number of active training points |
  | random | random               | bool | true  |Randomize training data for balancing and bagging |
  | mem    | mem                  | unsigned long int | 0 |Buffer size (in MB) to read image data blocks in memory |

  Usage: pksvm -t training [-i input -o output] [-cv value]


  Examples
  ========
  Some examples how to use pksvm can be found \ref examples_pksvm "here"
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
  Optionpk<unsigned long int> memory_opt("mem", "mem", "Buffer size (in MB) to read image data blocks in memory",0,1);

  oformat_opt.setHide(1);
  option_opt.setHide(1);
  memory_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    memory_opt.retrieveOption(argc,argv);

    app::AppFactory app(argc,argv);

    if(doProcess&&input_opt.empty()){
      std::ostringstream errorStream;
      errorStream << "Error: no input file provided (use option -i). Use --help for help information" << std::endl;
      throw(errorStream.str());
    }
    if(doProcess&&output_opt.empty()){
      std::ostringstream errorStream;
      errorStream << "Error: no output file provided (use option -o). Use --help for help information" << std::endl;
      throw(errorStream.str());
    }
    ImgRaster imgRaster;
    imgRaster.open(input_opt[0],memory_opt[0]);
    std::shared_ptr<ImgRaster> imgWriter = std::make_shared<ImgRaster>();
    string imageType;
    if(oformat_opt.size())//default
      imageType=oformat_opt[0];
    else
      imageType=imgRaster.getImageType();
    imgWriter->setFile(output_opt[0],imageType,memory_opt[0],option_opt);

    imgRaster.svm(imgWriter,app);
    imgRaster.close();
    imgWriter->close();
  }
  catch(string helpString){
    cerr << helpString << endl;
    cout << "Usage: pksvm -t training [-i input -o output] [-cv value]" << endl;
    return(1);
  }
  return(0);
}
