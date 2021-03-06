/**********************************************************************
pkann.cc: classify raster image using Artificial Neural Network
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
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>
#include "imageclasses/ImgRasterGdal.h"
#include "imageclasses/ImgReaderOgr.h"
#include "imageclasses/ImgWriterOgr.h"
#include "base/Optionpk.h"
#include "base/PosValue.h"
#include "algorithms/ConfusionMatrix.h"
#include "floatfann.h"
#include "algorithms/myfann_cpp.h"

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
  Optionpk<string> training_opt("t", "training", "training vector file. A single vector file contains all training features (must be set as: B0, B1, B2,...) for all classes (class numbers identified by label option). Use multiple training files for bootstrap aggregation (alternative to the bag and bsize options, where a random subset is taken from a single training file)"); 
  Optionpk<string> tlayer_opt("tln", "tln", "training layer name(s)");
  Optionpk<string> label_opt("label", "label", "identifier for class label in training vector file.","label"); 
  Optionpk<unsigned int> balance_opt("bal", "balance", "balance the input data to this number of samples for each class", 0);
  Optionpk<bool> random_opt("random", "random", "in case of balance, randomize input data", true,2);
  Optionpk<int> minSize_opt("min", "min", "if number of training pixels is less then min, do not take this class into account (0: consider all classes)", 0);
  Optionpk<unsigned int> band_opt("b", "band", "band index (starting from 0, either use band option or use start to end)");
  Optionpk<unsigned int> bstart_opt("sband", "startband", "Start band sequence number"); 
  Optionpk<unsigned int> bend_opt("eband", "endband", "End band sequence number"); 
  Optionpk<double> offset_opt("offset", "offset", "offset value for each spectral band input features: refl[band]=(DN[band]-offset[band])/scale[band]", 0.0);
  Optionpk<double> scale_opt("scale", "scale", "scale value for each spectral band input features: refl=(DN[band]-offset[band])/scale[band] (use 0 if scale min and max in each band to -1.0 and 1.0)", 0.0);
  Optionpk<unsigned short> aggreg_opt("a", "aggreg", "how to combine aggregated classifiers, see also rc option (1: sum rule, 2: max rule).",1);
  Optionpk<double> priors_opt("prior", "prior", "prior probabilities for each class (e.g., -p 0.3 -p 0.3 -p 0.2 )", 0.0); 
  Optionpk<string> priorimg_opt("pim", "priorimg", "prior probability image (multi-band img with band for each class","",2); 
  Optionpk<unsigned short> cv_opt("cv", "cv", "n-fold cross validation mode",0);
  Optionpk<string> cmformat_opt("cmf","cmf","Format for confusion matrix (ascii or latex)","ascii");
  Optionpk<unsigned int> nneuron_opt("nn", "nneuron", "number of neurons in hidden layers in neural network (multiple hidden layers are set by defining multiple number of neurons: -n 15 -n 1, default is one hidden layer with 5 neurons)", 5); 
  Optionpk<float> connection_opt("\0", "connection", "connection reate (default: 1.0 for a fully connected network)", 1.0); 
  Optionpk<float> learning_opt("l", "learning", "learning rate (default: 0.7)", 0.7); 
  Optionpk<float> weights_opt("w", "weights", "weights for neural network. Apply to fully connected network only, starting from first input neuron to last output neuron, including the bias neurons (last neuron in each but last layer)", 0.0); 
  Optionpk<unsigned int> maxit_opt("\0", "maxit", "number of maximum iterations (epoch) (default: 500)", 500); 
  Optionpk<unsigned short> comb_opt("comb", "comb", "how to combine bootstrap aggregation classifiers (0: sum rule, 1: product rule, 2: max rule). Also used to aggregate classes with rc option. Default is sum rule (0)",0); 
  Optionpk<unsigned short> bag_opt("bag", "bag", "Number of bootstrap aggregations (default is no bagging: 1)", 1);
  Optionpk<int> bagSize_opt("bs", "bsize", "Percentage of features used from available training features for each bootstrap aggregation (one size for all classes, or a different size for each class respectively", 100);
  Optionpk<string> classBag_opt("cb", "classbag", "output for each individual bootstrap aggregation (default is blank)"); 
  Optionpk<string> mask_opt("m", "mask", "Only classify within specified mask (vector or raster). For raster mask, set nodata values with the option msknodata.");
  Optionpk<short> msknodata_opt("msknodata", "msknodata", "mask value(s) not to consider for classification. Values will be taken over in classification image. Default is 0", 0);
  Optionpk<unsigned short> nodata_opt("nodata", "nodata", "nodata value to put where image is masked as nodata", 0);
  Optionpk<string> output_opt("o", "output", "output classification image"); 
  Optionpk<string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate).","GTiff");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<string> colorTable_opt("ct", "ct", "colour table in ASCII format having 5 columns: id R G B ALFA (0: transparent, 255: solid)"); 
  Optionpk<string> prob_opt("\0", "prob", "probability image. Default is no probability image"); 
  Optionpk<string> entropy_opt("entropy", "entropy", "entropy image (measure for uncertainty of classifier output","",2); 
  Optionpk<string> active_opt("active", "active", "ogr output for active training sample.","",2); 
  Optionpk<string> ogrformat_opt("f", "f", "Output ogr format for active training sample","SQLite");
  Optionpk<unsigned int> nactive_opt("na", "nactive", "number of active training points",1);
  Optionpk<string> classname_opt("c", "class", "list of class names."); 
  Optionpk<short> classvalue_opt("r", "reclass", "list of class values (use same order as in class opt)."); 
  Optionpk<short> verbose_opt("v", "verbose", "set to: 0 (results only), 1 (confusion matrix), 2 (debug)",0,2);

  option_opt.setHide(1);
  oformat_opt.setHide(1);
  band_opt.setHide(1);
  bstart_opt.setHide(1);
  bend_opt.setHide(1);
  balance_opt.setHide(1);
  minSize_opt.setHide(1);
  bag_opt.setHide(1);
  bagSize_opt.setHide(1);
  comb_opt.setHide(1);
  classBag_opt.setHide(1);
  minSize_opt.setHide(1);
  prob_opt.setHide(1);
  priorimg_opt.setHide(1);
  minSize_opt.setHide(1);
  offset_opt.setHide(1);
  scale_opt.setHide(1);
  connection_opt.setHide(1);
  weights_opt.setHide(1);
  maxit_opt.setHide(1);
  learning_opt.setHide(1);
  verbose_opt.setHide(2);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    training_opt.retrieveOption(argc,argv);
    tlayer_opt.retrieveOption(argc,argv);
    label_opt.retrieveOption(argc,argv);
    balance_opt.retrieveOption(argc,argv);
    random_opt.retrieveOption(argc,argv);
    minSize_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    bstart_opt.retrieveOption(argc,argv);
    bend_opt.retrieveOption(argc,argv);
    offset_opt.retrieveOption(argc,argv);
    scale_opt.retrieveOption(argc,argv);
    aggreg_opt.retrieveOption(argc,argv);
    priors_opt.retrieveOption(argc,argv);
    priorimg_opt.retrieveOption(argc,argv);
    cv_opt.retrieveOption(argc,argv);
    cmformat_opt.retrieveOption(argc,argv);
    nneuron_opt.retrieveOption(argc,argv);
    connection_opt.retrieveOption(argc,argv);
    weights_opt.retrieveOption(argc,argv);
    learning_opt.retrieveOption(argc,argv);
    maxit_opt.retrieveOption(argc,argv);
    comb_opt.retrieveOption(argc,argv);
    bag_opt.retrieveOption(argc,argv);
    bagSize_opt.retrieveOption(argc,argv);
    classBag_opt.retrieveOption(argc,argv);
    mask_opt.retrieveOption(argc,argv);
    msknodata_opt.retrieveOption(argc,argv);
    nodata_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    otype_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    prob_opt.retrieveOption(argc,argv);
    entropy_opt.retrieveOption(argc,argv);
    active_opt.retrieveOption(argc,argv);
    ogrformat_opt.retrieveOption(argc,argv);
    nactive_opt.retrieveOption(argc,argv);
    classname_opt.retrieveOption(argc,argv);
    classvalue_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    cout << endl;
    cout << "Usage: pkann -t training [-i input -o output] [-cv value]" << endl;
    cout << endl;
    cout << "short option -h shows basic options only, use long option --help to show all options" << endl;
    exit(0);//help was invoked, stop processing
  }

  if(entropy_opt[0]=="")
    entropy_opt.clear();
  if(active_opt[0]=="")
    active_opt.clear();
  if(priorimg_opt[0]=="")
    priorimg_opt.clear();

  if(verbose_opt[0]>=1){
    if(input_opt.size())
      cout << "image filename: " << input_opt[0] << endl;
    if(mask_opt.size())
      cout << "mask filename: " << mask_opt[0] << endl;
    if(training_opt.size()){
      cout << "training vector file: " << endl;
      for(unsigned int ifile=0;ifile<training_opt.size();++ifile)
        cout << training_opt[ifile] << endl;
    }
    else
      cerr << "no training file set!" << endl;
    cout << "verbose: " << verbose_opt[0] << endl;
  }
  unsigned short nbag=(training_opt.size()>1)?training_opt.size():bag_opt[0];
  if(verbose_opt[0]>=1)
    cout << "number of bootstrap aggregations: " << nbag << endl;
  
  ImgReaderOgr extentReader;
  // OGRLayer  *readLayer;

  double ulx=0;
  double uly=0;
  double lrx=0;
  double lry=0;

  bool maskIsVector=false;
  if(mask_opt.size()){
    try{
      extentReader.open(mask_opt[0]);
      maskIsVector=true;
      // readLayer = extentReader.getDataSource()->GetLayer(0);
      if(!(extentReader.getExtent(ulx,uly,lrx,lry))){
	cerr << "Error: could not get extent from " << mask_opt[0] << endl;
	exit(1);
      }
    }
    catch(string errorString){
      maskIsVector=false;
    }
  }

  ImgWriterOgr activeWriter;
  if(active_opt.size()){
    ImgReaderOgr trainingReader(training_opt[0]);
    activeWriter.open(active_opt[0],ogrformat_opt[0]);
    activeWriter.createLayer(active_opt[0],trainingReader.getProjection(),wkbPoint,NULL);
    activeWriter.copyFields(trainingReader);
  }
  vector<PosValue> activePoints(nactive_opt[0]);
  for(unsigned int iactive=0;iactive<activePoints.size();++iactive){
    activePoints[iactive].value=1.0;
    activePoints[iactive].posx=0.0;
    activePoints[iactive].posy=0.0;
  }

  unsigned int totalSamples=0;
  unsigned int nactive=0;
  vector<FANN::neural_net> net(nbag);//the neural network

  unsigned int nclass=0;
  unsigned int nband=0;
  unsigned int startBand=2;//first two bands represent X and Y pos

  if(priors_opt.size()>1){//priors from argument list
    priors.resize(priors_opt.size());
    double normPrior=0;
    for(unsigned int iclass=0;iclass<priors_opt.size();++iclass){
      priors[iclass]=priors_opt[iclass];
      normPrior+=priors[iclass];
    }
    //normalize
    for(unsigned int iclass=0;iclass<priors_opt.size();++iclass)
      priors[iclass]/=normPrior;
  }

  //convert start and end band options to vector of band indexes
  try{
    if(bstart_opt.size()){
      if(bend_opt.size()!=bstart_opt.size()){
	string errorstring="Error: options for start and end band indexes must be provided as pairs, missing end band";
	throw(errorstring);
      }
      band_opt.clear();
      for(unsigned int ipair=0;ipair<bstart_opt.size();++ipair){
	if(bend_opt[ipair]<=bstart_opt[ipair]){
	  string errorstring="Error: index for end band must be smaller then start band";
	  throw(errorstring);
	}
	for(unsigned int iband=bstart_opt[ipair];iband<=bend_opt[ipair];++iband)
	  band_opt.push_back(iband);
      }
    }
  }
  catch(string error){
    cerr << error << std::endl;
    exit(1);
  }
  //sort bands
  if(band_opt.size())
    std::sort(band_opt.begin(),band_opt.end());

  map<string,short> classValueMap;
  vector<std::string> nameVector;
  if(classname_opt.size()){
    assert(classname_opt.size()==classvalue_opt.size());
    for(unsigned int iclass=0;iclass<classname_opt.size();++iclass)
      classValueMap[classname_opt[iclass]]=classvalue_opt[iclass];
  }
  //----------------------------------- Training -------------------------------
  confusionmatrix::ConfusionMatrix cm;
  vector< vector<double> > offset(nbag);
  vector< vector<double> > scale(nbag);
  map<string,Vector2d<float> > trainingMap;
  vector< Vector2d<float> > trainingPixels;//[class][sample][band]
  vector<string> fields;
  for(unsigned int ibag=0;ibag<nbag;++ibag){
    //organize training data
    if(ibag<training_opt.size()){//if bag contains new training pixels
      trainingMap.clear();
      trainingPixels.clear();
      if(verbose_opt[0]>=1)
        cout << "reading imageVector file " << training_opt[0] << endl;
      try{
	ImgReaderOgr trainingReaderBag(training_opt[ibag]);
        if(band_opt.size())
          totalSamples=trainingReaderBag.readDataImageOgr(trainingMap,fields,band_opt,label_opt[0],tlayer_opt);
        else
          totalSamples=trainingReaderBag.readDataImageOgr(trainingMap,fields,0,0,label_opt[0],tlayer_opt);
        if(trainingMap.size()<2){
          string errorstring="Error: could not read at least two classes from training file, did you provide class labels in training sample (see option label)?";
          throw(errorstring);
        }
	trainingReaderBag.close();
      }
      catch(string error){
        cerr << error << std::endl;
        exit(1);
      }
      catch(...){
        cerr << "error caught" << std::endl;
        exit(1);
      }
      //delete class 0 ?
      // if(verbose_opt[0]>=1)
      //   std::cout << "erasing class 0 from training set (" << trainingMap[0].size() << " from " << totalSamples << ") samples" << std::endl;
      // totalSamples-=trainingMap[0].size();
      // trainingMap.erase(0);
      //convert map to vector
      if(verbose_opt[0]>1)
        std::cout << "training pixels: " << std::endl;
      map<string,Vector2d<float> >::iterator mapit=trainingMap.begin();
      while(mapit!=trainingMap.end()){
        //delete small classes
        if((mapit->second).size()<minSize_opt[0]){
          trainingMap.erase(mapit);
          continue;
        }
        trainingPixels.push_back(mapit->second);
        if(verbose_opt[0]>1)
          std::cout << mapit->first << ": " << (mapit->second).size() << " samples" << std::endl;
        ++mapit;
        }
      if(!ibag){
        nclass=trainingPixels.size();
	if(classname_opt.size())
	  assert(nclass==classname_opt.size());
        nband=(training_opt.size())?trainingPixels[0][0].size()-2:trainingPixels[0][0].size();//X and Y
      }
      else{
        assert(nclass==trainingPixels.size());
        assert(nband==(training_opt.size())?trainingPixels[0][0].size()-2:trainingPixels[0][0].size());
      }

      //do not remove outliers here: could easily be obtained through ogr2ogr -where 'B2<110' output.shp input.shp
      //balance training data
      if(balance_opt[0]>0){
        while(balance_opt.size()<nclass)
          balance_opt.push_back(balance_opt.back());
        if(random_opt[0])
          srand(time(NULL));
        totalSamples=0;
        for(unsigned int iclass=0;iclass<nclass;++iclass){
          if(trainingPixels[iclass].size()>balance_opt[iclass]){
            while(trainingPixels[iclass].size()>balance_opt[iclass]){
              unsigned int index=rand()%trainingPixels[iclass].size();
              trainingPixels[iclass].erase(trainingPixels[iclass].begin()+index);
            }
          }
          else{
            unsigned int oldsize=trainingPixels[iclass].size();
            for(unsigned int isample=trainingPixels[iclass].size();isample<balance_opt[iclass];++isample){
              unsigned int index = rand()%oldsize;
              trainingPixels[iclass].push_back(trainingPixels[iclass][index]);
            }
          }
          totalSamples+=trainingPixels[iclass].size();
        }
      }
    
      //set scale and offset
      offset[ibag].resize(nband);
      scale[ibag].resize(nband);
      if(offset_opt.size()>1)
        assert(offset_opt.size()==nband);
      if(scale_opt.size()>1)
        assert(scale_opt.size()==nband);
      for(unsigned int iband=0;iband<nband;++iband){
        if(verbose_opt[0]>=1)
          cout << "scaling for band" << iband << endl;
        offset[ibag][iband]=(offset_opt.size()==1)?offset_opt[0]:offset_opt[iband];
        scale[ibag][iband]=(scale_opt.size()==1)?scale_opt[0]:scale_opt[iband];
        //search for min and maximum
        if(scale[ibag][iband]<=0){
          float theMin=trainingPixels[0][0][iband+startBand];
          float theMax=trainingPixels[0][0][iband+startBand];
          for(unsigned int iclass=0;iclass<nclass;++iclass){
            for(unsigned int isample=0;isample<trainingPixels[iclass].size();++isample){
              if(theMin>trainingPixels[iclass][isample][iband+startBand])
                theMin=trainingPixels[iclass][isample][iband+startBand];
              if(theMax<trainingPixels[iclass][isample][iband+startBand])
                theMax=trainingPixels[iclass][isample][iband+startBand];
            }
          }
          offset[ibag][iband]=theMin+(theMax-theMin)/2.0;
          scale[ibag][iband]=(theMax-theMin)/2.0;
          if(verbose_opt[0]>=1){
            std::cout << "Extreme image values for band " << iband << ": [" << theMin << "," << theMax << "]" << std::endl;
            std::cout << "Using offset, scale: " << offset[ibag][iband] << ", " << scale[ibag][iband] << std::endl;
            std::cout << "scaled values for band " << iband << ": [" << (theMin-offset[ibag][iband])/scale[ibag][iband] << "," << (theMax-offset[ibag][iband])/scale[ibag][iband] << "]" << std::endl;
          }
        }
      }
    }
    else{//use same offset and scale 
      offset[ibag].resize(nband);
      scale[ibag].resize(nband);
      for(unsigned int iband=0;iband<nband;++iband){
        offset[ibag][iband]=offset[0][iband];
        scale[ibag][iband]=scale[0][iband];
      }
    }
      
    if(!ibag){
      if(priors_opt.size()==1){//default: equal priors for each class
        priors.resize(nclass);
        for(unsigned int iclass=0;iclass<nclass;++iclass)
          priors[iclass]=1.0/nclass;
      }
      assert(priors_opt.size()==1||priors_opt.size()==nclass);
    
      //set bagsize for each class if not done already via command line
      while(bagSize_opt.size()<nclass)
        bagSize_opt.push_back(bagSize_opt.back());

      if(verbose_opt[0]>=1){
        std::cout << "number of bands: " << nband << std::endl;
        std::cout << "number of classes: " << nclass << std::endl;
        std::cout << "priors:";
	if(priorimg_opt.empty()){
          for(unsigned int iclass=0;iclass<nclass;++iclass)
            std::cout << " " << priors[iclass];
          std::cout << std::endl;
        }
      }
      map<string,Vector2d<float> >::iterator mapit=trainingMap.begin();
      bool doSort=true;
      while(mapit!=trainingMap.end()){
	nameVector.push_back(mapit->first);
	if(classValueMap.size()){
	  //check if name in training is covered by classname_opt (values can not be 0)
	  if(classValueMap[mapit->first]>0){
	    if(cm.getClassIndex(type2string<short>(classValueMap[mapit->first]))<0)
	      cm.pushBackClassName(type2string<short>(classValueMap[mapit->first]),doSort);
	  }
	  else{
	    std::cerr << "Error: names in classname option are not complete, please check names in training vector and make sure classvalue is > 0" << std::endl;
	    exit(1);
	  }
	}
	else
	  cm.pushBackClassName(mapit->first,doSort);
	++mapit;
      }
      if(classname_opt.empty()){
        //std::cerr << "Warning: no class name and value pair provided for all " << nclass << " classes, using string2type<int> instead!" << std::endl;
        for(unsigned int iclass=0;iclass<nclass;++iclass){
          if(verbose_opt[0])
            std::cout << iclass << " " << cm.getClass(iclass) << " -> " << string2type<short>(cm.getClass(iclass)) << std::endl;
          classValueMap[cm.getClass(iclass)]=string2type<short>(cm.getClass(iclass));
        }
      }
      if(priors_opt.size()==nameVector.size()){
	std::cerr << "Warning: please check if priors are provided in correct order!!!" << std::endl;
	for(unsigned int iclass=0;iclass<nameVector.size();++iclass)
	  std::cerr << nameVector[iclass] << " " << priors_opt[iclass] << std::endl;
      }
    }//if(!ibag)

    //Calculate features of training set
    vector< Vector2d<float> > trainingFeatures(nclass);
    for(unsigned int iclass=0;iclass<nclass;++iclass){
      unsigned int nctraining=0;
      if(verbose_opt[0]>=1)
        cout << "calculating features for class " << iclass << endl;
      if(random_opt[0])
        srand(time(NULL));
      nctraining=(bagSize_opt[iclass]<100)? trainingPixels[iclass].size()/100.0*bagSize_opt[iclass] : trainingPixels[iclass].size();//bagSize_opt[iclass] given in % of training size
      if(nctraining<=0)
        nctraining=1;
      assert(nctraining<=trainingPixels[iclass].size());
      if(bagSize_opt[iclass]<100)
        random_shuffle(trainingPixels[iclass].begin(),trainingPixels[iclass].end());
      
      trainingFeatures[iclass].resize(nctraining);
      for(unsigned int isample=0;isample<nctraining;++isample){
        //scale pixel values according to scale and offset!!!
        for(unsigned int iband=0;iband<nband;++iband){
          float value=trainingPixels[iclass][isample][iband+startBand];
          trainingFeatures[iclass][isample].push_back((value-offset[ibag][iband])/scale[ibag][iband]);
        }
      }
      assert(trainingFeatures[iclass].size()==nctraining);
    }
    
    unsigned int nFeatures=trainingFeatures[0][0].size();
    unsigned int ntraining=0;
    for(unsigned int iclass=0;iclass<nclass;++iclass)
      ntraining+=trainingFeatures[iclass].size();

    const unsigned int num_layers = nneuron_opt.size()+2;
    const float desired_error = 0.0003;
    const unsigned int iterations_between_reports = (verbose_opt[0])? maxit_opt[0]+1:0;
    if(verbose_opt[0]>=1){
      cout << "number of features: " << nFeatures << endl;
      cout << "creating artificial neural network with " << nneuron_opt.size() << " hidden layer, having " << endl;
      for(unsigned int ilayer=0;ilayer<nneuron_opt.size();++ilayer)
        cout << nneuron_opt[ilayer] << " ";
      cout << "neurons" << endl;
      cout << "connection_opt[0]: " << connection_opt[0] << std::endl;
      cout << "num_layers: " << num_layers << std::endl;
      cout << "nFeatures: " << nFeatures << std::endl;
      cout << "nneuron_opt[0]: " << nneuron_opt[0] << std::endl;
      cout << "number of classes (nclass): " << nclass << std::endl;
    }
    switch(num_layers){
    case(3):{
      // net[ibag].create_sparse(connection_opt[0],num_layers, nFeatures, nneuron_opt[0], nclass);//replace all create_sparse with create_sparse_array due to bug in FANN!
      unsigned int layers[3];
      layers[0]=nFeatures;
      layers[1]=nneuron_opt[0];
      layers[2]=nclass;
      net[ibag].create_sparse_array(connection_opt[0],num_layers,layers);
      break;
    }
    case(4):{
      unsigned int layers[4];
      layers[0]=nFeatures;
      layers[1]=nneuron_opt[0];
      layers[2]=nneuron_opt[1];
      layers[3]=nclass;
      // layers.push_back(nFeatures);
      // for(int ihidden=0;ihidden<nneuron_opt.size();++ihidden)
      // 	layers.push_back(nneuron_opt[ihidden]);
      // layers.push_back(nclass);
      net[ibag].create_sparse_array(connection_opt[0],num_layers,layers);
      break;
    }
    default:
      cerr << "Only 1 or 2 hidden layers are supported!" << endl;
      exit(1);
      break;
    }
    if(verbose_opt[0]>=1)
      cout << "network created" << endl;
  
    net[ibag].set_learning_rate(learning_opt[0]);

    //   net.set_activation_steepness_hidden(1.0);
    //   net.set_activation_steepness_output(1.0);
    
    net[ibag].set_activation_function_hidden(FANN::SIGMOID_SYMMETRIC_STEPWISE);
    net[ibag].set_activation_function_output(FANN::SIGMOID_SYMMETRIC_STEPWISE);

    // Set additional properties such as the training algorithm
    //   net.set_training_algorithm(FANN::TRAIN_QUICKPROP);

    // Output network type and parameters
    if(verbose_opt[0]>=1){
      cout << endl << "Network Type                         :  ";
      switch (net[ibag].get_network_type())
        {
        case FANN::LAYER:
          cout << "LAYER" << endl;
          break;
        case FANN::SHORTCUT:
          cout << "SHORTCUT" << endl;
          break;
        default:
          cout << "UNKNOWN" << endl;
          break;
        }
      net[ibag].print_parameters();
    }
      
    if(cv_opt[0]>1){
      if(verbose_opt[0])
        std::cout << "cross validation" << std::endl;
      vector<unsigned short> referenceVector;
      vector<unsigned short> outputVector;
      float rmse=net[ibag].cross_validation(trainingFeatures,
                                            ntraining,
                                            cv_opt[0],
                                            maxit_opt[0],
                                            desired_error,
                                            referenceVector,
                                            outputVector,
                                            verbose_opt[0]);
      map<string,Vector2d<float> >::iterator mapit=trainingMap.begin();
      for(unsigned int isample=0;isample<referenceVector.size();++isample){
	string refClassName=nameVector[referenceVector[isample]];
	string className=nameVector[outputVector[isample]];
	if(classValueMap.size())
	  cm.incrementResult(type2string<short>(classValueMap[refClassName]),type2string<short>(classValueMap[className]),1.0/nbag);
	else
	  cm.incrementResult(cm.getClass(referenceVector[isample]),cm.getClass(outputVector[isample]),1.0/nbag);
      }        
    }
  
    if(verbose_opt[0]>=1)
      cout << endl << "Set training data" << endl;

    if(verbose_opt[0]>=1)
      cout << endl << "Training network" << endl;
    
    if(verbose_opt[0]>=1){
      cout << "Max Epochs " << setw(8) << maxit_opt[0] << ". "
           << "Desired Error: " << left << desired_error << right << endl;
    }
    if(weights_opt.size()==net[ibag].get_total_connections()){//no new training needed (same training sample)
      vector<fann_connection> convector;
      net[ibag].get_connection_array(convector);
      for(unsigned int i_connection=0;i_connection<net[ibag].get_total_connections();++i_connection)
        convector[i_connection].weight=weights_opt[i_connection];
      net[ibag].set_weight_array(convector);
    }
    else{
      bool initWeights=true;
      net[ibag].train_on_data(trainingFeatures,ntraining,initWeights, maxit_opt[0],
                              iterations_between_reports, desired_error);
    }


    if(verbose_opt[0]>=2){
      net[ibag].print_connections();
      vector<fann_connection> convector;
      net[ibag].get_connection_array(convector);
      for(unsigned int i_connection=0;i_connection<net[ibag].get_total_connections();++i_connection)
        cout << "connection " << i_connection << ": " << convector[i_connection].weight << endl;

    }
  }//for ibag
  if(cv_opt[0]>1){
    assert(cm.nReference());
    cm.setFormat(cmformat_opt[0]);
    cm.reportSE95(false);
    std::cout << cm << std::endl;
    cout << "class #samples userAcc prodAcc" << endl;
    double se95_ua=0;
    double se95_pa=0;
    double se95_oa=0;
    double dua=0;
    double dpa=0;
    double doa=0;
    for(unsigned int iclass=0;iclass<cm.nClasses();++iclass){
      dua=cm.ua_pct(cm.getClass(iclass),&se95_ua);
      dpa=cm.pa_pct(cm.getClass(iclass),&se95_pa);
      cout << cm.getClass(iclass) << " " << cm.nReference(cm.getClass(iclass)) << " " << dua << " (" << se95_ua << ")" << " " << dpa << " (" << se95_pa << ")" << endl;
    }
    std::cout << "Kappa: " << cm.kappa() << std::endl;
    doa=cm.oa_pct(&se95_oa);
    std::cout << "Overall Accuracy: " << doa << " (" << se95_oa << ")"  << std::endl;
  }
  //--------------------------------- end of training -----------------------------------
  if(input_opt.empty())
    exit(0);

    const char* pszMessage;
    void* pProgressArg=NULL;
    GDALProgressFunc pfnProgress=GDALTermProgress;
    float progress=0;
  //-------------------------------- open image file ------------------------------------
  bool inputIsRaster=false;
  ImgReaderOgr imgReaderOgr;
  try{
    imgReaderOgr.open(input_opt[0]);
    imgReaderOgr.close();
  }
  catch(string errorString){
    inputIsRaster=true;
  }
  if(inputIsRaster){
  // if(input_opt[0].find(".shp")==string::npos){
    ImgRasterGdal testImage;
    try{
      if(verbose_opt[0]>=1)
        cout << "opening image " << input_opt[0] << endl; 
      testImage.open(input_opt[0]);
    }
    catch(string error){
      cerr << error << endl;
      exit(2);
    }
    ImgRasterGdal priorReader;
    if(priorimg_opt.size()){
      try{
	if(verbose_opt[0]>=1)
          std::cout << "opening prior image " << priorimg_opt[0] << std::endl;
        priorReader.open(priorimg_opt[0]);
        assert(priorReader.nrOfCol()==testImage.nrOfCol());
        assert(priorReader.nrOfRow()==testImage.nrOfRow());
      }
      catch(string error){
        cerr << error << std::endl;
        exit(2);
      }
      catch(...){
        cerr << "error caught" << std::endl;
        exit(1);
      }
    }

    unsigned int nrow=testImage.nrOfRow();
    unsigned int ncol=testImage.nrOfCol();
    if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
      string theInterleave="INTERLEAVE=";
      theInterleave+=testImage.getInterleave();
      option_opt.push_back(theInterleave);
    }
    vector<char> classOut(ncol);//classified line for writing to image file

    //   assert(nband==testImage.nrOfBand());
    ImgRasterGdal classImageBag;
    ImgRasterGdal classImageOut;
    ImgRasterGdal probImage;
    ImgRasterGdal entropyImage;

    string imageType;//testImage.getImageType();
    if(oformat_opt.size())//default
      imageType=oformat_opt[0];
    try{

      if(verbose_opt[0]>=1)
        cout << "opening class image for writing output " << output_opt[0] << endl;
      if(classBag_opt.size()){
        classImageBag.open(classBag_opt[0],ncol,nrow,nbag,GDT_Byte,imageType,option_opt);
	classImageBag.GDALSetNoDataValue(nodata_opt[0]);
        classImageBag.copyGeoTransform(testImage);
        classImageBag.setProjection(testImage.getProjection());
      }
      classImageOut.open(output_opt[0],ncol,nrow,1,GDT_Byte,imageType,option_opt);
      classImageOut.GDALSetNoDataValue(nodata_opt[0]);
      classImageOut.copyGeoTransform(testImage);
      classImageOut.setProjection(testImage.getProjection());
      if(colorTable_opt.size())
        classImageOut.setColorTable(colorTable_opt[0],0);
      if(prob_opt.size()){
        probImage.open(prob_opt[0],ncol,nrow,nclass,GDT_Byte,imageType,option_opt);
	probImage.GDALSetNoDataValue(nodata_opt[0]);
        probImage.copyGeoTransform(testImage);
        probImage.setProjection(testImage.getProjection());
      }
      if(entropy_opt.size()){
        entropyImage.open(entropy_opt[0],ncol,nrow,1,GDT_Byte,imageType,option_opt);
	entropyImage.GDALSetNoDataValue(nodata_opt[0]);
        entropyImage.copyGeoTransform(testImage);
        entropyImage.setProjection(testImage.getProjection());
      }
    }
    catch(string error){
      cerr << error << endl;
    }
  
    ImgRasterGdal maskWriter;

    if(maskIsVector){
      try{
        //todo: produce unique name or perform in memory solving issue on flush memory buffer (check gdal development list on how to retrieve gdal mem buffer)
	maskWriter.open("/vsimem/mask.tif",ncol,nrow,1,GDT_Float32,imageType,option_opt);
	maskWriter.GDALSetNoDataValue(nodata_opt[0]);
        maskWriter.copyGeoTransform(testImage);
        maskWriter.setProjection(testImage.getProjection());
	vector<double> burnValues(1,1);//burn value is 1 (single band)
	maskWriter.rasterizeOgr(extentReader,burnValues);
	extentReader.close();
	maskWriter.close();
      }
      catch(string error){
        cerr << error << std::endl;
        exit(2);
      }
      catch(...){
        cerr << "error caught" << std::endl;
        exit(1);
      }
      mask_opt.clear();
      mask_opt.push_back("/vsimem/mask.tif");
    }
    ImgRasterGdal maskReader;
    if(mask_opt.size()){
      try{
        if(verbose_opt[0]>=1)
          std::cout << "opening mask image file " << mask_opt[0] << std::endl;
        maskReader.open(mask_opt[0],imageType);
      }
      catch(string error){
        cerr << error << std::endl;
        exit(2);
      }
      catch(...){
        cerr << "error caught" << std::endl;
        exit(1);
      }
    }

    for(unsigned int iline=0;iline<nrow;++iline){
      vector<float> buffer(ncol);
      vector<short> lineMask;
      if(mask_opt.size())
        lineMask.resize(maskReader.nrOfCol());
      Vector2d<float> linePrior;
      if(priorimg_opt.size())
	 linePrior.resize(nclass,ncol);//prior prob for each class
      Vector2d<float> hpixel(ncol);
      Vector2d<float> fpixel(ncol);
      Vector2d<float> probOut(nclass,ncol);//posterior prob for each (internal) class
      vector<float> entropy(ncol);
      Vector2d<char> classBag;//classified line for writing to image file
      if(classBag_opt.size())
        classBag.resize(nbag,ncol);
      //read all bands of all pixels in this line in hline
      try{
        if(band_opt.size()){
          for(unsigned int iband=0;iband<band_opt.size();++iband){
            if(verbose_opt[0]==2)
              std::cout << "reading band " << band_opt[iband] << std::endl;
            assert(band_opt[iband]>=0);
            assert(band_opt[iband]<testImage.nrOfBand());
            testImage.readData(buffer,iline,band_opt[iband]);
            for(unsigned int icol=0;icol<ncol;++icol)
              hpixel[icol].push_back(buffer[icol]);
          }
        }
        else{
          for(unsigned int iband=0;iband<nband;++iband){
            if(verbose_opt[0]==2)
              std::cout << "reading band " << iband << std::endl;
            assert(iband>=0);
            assert(iband<testImage.nrOfBand());
            testImage.readData(buffer,iline,iband);
            for(unsigned int icol=0;icol<ncol;++icol)
              hpixel[icol].push_back(buffer[icol]);
          }
        }
      }
      catch(string theError){
        cerr << "Error reading " << input_opt[0] << ": " << theError << std::endl;
        exit(3);
      }
      catch(...){
        cerr << "error caught" << std::endl;
        exit(3);
      }
      assert(nband==hpixel[0].size());
      if(verbose_opt[0]==2)
        cout << "used bands: " << nband << endl;
      //read prior
      if(priorimg_opt.size()){
        try{
	  for(short iclass=0;iclass<nclass;++iclass){
	    if(verbose_opt.size()>1)
	      std::cout << "Reading " << priorimg_opt[0] << " band " << iclass << " line " << iline << std::endl;
	    priorReader.readData(linePrior[iclass],iline,iclass);
	  }
        }
        catch(string theError){
	  std::cerr << "Error reading " << priorimg_opt[0] << ": " << theError << std::endl;
          exit(3);
        }
        catch(...){
          cerr << "error caught" << std::endl;
          exit(3);
        }
      }
      double oldRowMask=-1;//keep track of row mask to optimize number of line readings
      //process per pixel
      for(unsigned int icol=0;icol<ncol;++icol){
        assert(hpixel[icol].size()==nband);
	bool doClassify=true;
        bool masked=false;
	double geox=0;
	double geoy=0;
        if(maskIsVector){
	  doClassify=false;
	  testImage.image2geo(icol,iline,geox,geoy);
	  //check enveloppe first
	  if(uly>=geoy&&lry<=geoy&&ulx<=geox&&lrx>=geox){
	    doClassify=true;
	  }
	}
        if(mask_opt.size()){
	  //read mask
	  double colMask=0;
	  double rowMask=0;

	  testImage.image2geo(icol,iline,geox,geoy);
	  maskReader.geo2image(geox,geoy,colMask,rowMask);
	  colMask=static_cast<unsigned int>(colMask);
	  rowMask=static_cast<unsigned int>(rowMask);
	  if(rowMask>=0&&rowMask<maskReader.nrOfRow()&&colMask>=0&&colMask<maskReader.nrOfCol()){
	    if(static_cast<unsigned int>(rowMask)!=static_cast<unsigned int>(oldRowMask)){
	      assert(rowMask>=0&&rowMask<maskReader.nrOfRow());
	      try{
		// maskReader.readData(lineMask[imask],GDT_Int32,static_cast<unsigned int>(rowMask));
		maskReader.readData(lineMask,static_cast<unsigned int>(rowMask));
	      }
	      catch(string errorstring){
		cerr << errorstring << endl;
		exit(1);
	      }
	      catch(...){
		cerr << "error caught" << std::endl;
		exit(3);
	      }
	      oldRowMask=rowMask;
	    }
	    short theMask=0;
	    for(short ivalue=0;ivalue<msknodata_opt.size();++ivalue){
	      // if(msknodata_opt[ivalue]>=0){//values set in msknodata_opt are invalid
		if(lineMask[colMask]==msknodata_opt[ivalue]){
		  theMask=lineMask[colMask];
		  masked=true;
		  break;
		}
	      // }
	      // else{//only values set in msknodata_opt are valid
	      // 	if(lineMask[colMask]!=-msknodata_opt[ivalue]){
	      // 	  theMask=lineMask[colMask];
	      // 	  masked=true;
	      // 	}
	      // 	else{
	      // 	  masked=false;
	      // 	  break;
	      // 	}
	      // }
	    }
	    if(masked){
	      if(classBag_opt.size())
		for(unsigned int ibag=0;ibag<nbag;++ibag)
		  classBag[ibag][icol]=theMask;
	      classOut[icol]=theMask;
	      continue;
	    }
	  }
	  bool valid=false;
	  for(unsigned int iband=0;iband<hpixel[icol].size();++iband){
	    if(hpixel[icol][iband]){
	      valid=true;
	      break;
	    }
	  }
	  if(!valid)
	    doClassify=false;

	}
        for(short iclass=0;iclass<nclass;++iclass)
          probOut[iclass][icol]=0;
	if(!doClassify){
	  if(classBag_opt.size())
	    for(unsigned int ibag=0;ibag<nbag;++ibag)
	      classBag[ibag][icol]=nodata_opt[0];
	  classOut[icol]=nodata_opt[0];
	  continue;//next column
	}
	if(verbose_opt[0]>1)
	  std::cout << "begin classification " << std::endl;
        //----------------------------------- classification -------------------
        for(unsigned int ibag=0;ibag<nbag;++ibag){
          //calculate image features
          fpixel[icol].clear();
          for(unsigned int iband=0;iband<nband;++iband)
            fpixel[icol].push_back((hpixel[icol][iband]-offset[ibag][iband])/scale[ibag][iband]);
          vector<float> result(nclass);
          result=net[ibag].run(fpixel[icol]);
          vector<float> prValues(nclass);
          float maxP=0;

          //calculate posterior prob of bag 
          if(classBag_opt.size()){
            //search for max prob within bag
            maxP=0;
            classBag[ibag][icol]=0;
          }
	  double normPrior=0;
	  if(priorimg_opt.size()){
	    for(short iclass=0;iclass<nclass;++iclass)
	      normPrior+=linePrior[iclass][icol];
	  }
          for(unsigned int iclass=0;iclass<nclass;++iclass){
	    result[iclass]=(result[iclass]+1.0)/2.0;//bring back to scale [0,1]
	    if(priorimg_opt.size())
	      priors[iclass]=linePrior[iclass][icol]/normPrior;//todo: check if correct for all cases... (automatic classValueMap and manual input for names and values)
            switch(comb_opt[0]){
            default:
            case(0)://sum rule
              probOut[iclass][icol]+=result[iclass]*priors[iclass];//add probabilities for each bag
              break;
            case(1)://product rule
              probOut[iclass][icol]*=pow(static_cast<float>(priors[iclass]),static_cast<float>(1.0-nbag)/nbag)*result[iclass];//multiply probabilities for each bag
              break;
            case(2)://max rule
              if(priors[iclass]*result[iclass]>probOut[iclass][icol])
                probOut[iclass][icol]=priors[iclass]*result[iclass];
              break;
            }
            if(classBag_opt.size()){
              //search for max prob within bag
              // if(prValues[iclass]>maxP){
              //   maxP=prValues[iclass];
              //   classBag[ibag][icol]=vcode[iclass];
              // }
              if(result[iclass]>maxP){
                maxP=result[iclass];
                classBag[ibag][icol]=iclass;
              }
            }
          }
        }//ibag

        //search for max class prob
        float maxBag1=0;//max probability
        float maxBag2=0;//second max probability
        float normBag=0;
        for(short iclass=0;iclass<nclass;++iclass){
          if(probOut[iclass][icol]>maxBag1){
            maxBag1=probOut[iclass][icol];
            classOut[icol]=classValueMap[nameVector[iclass]];
          }
	  else if(probOut[iclass][icol]>maxBag2)
            maxBag2=probOut[iclass][icol];
          normBag+=probOut[iclass][icol];
        }
        //normalize probOut and convert to percentage
        entropy[icol]=0;
        for(short iclass=0;iclass<nclass;++iclass){
          float prv=probOut[iclass][icol];
          prv/=normBag;
          entropy[icol]-=prv*log(prv)/log(2.0);
          prv*=100.0;
            
          probOut[iclass][icol]=static_cast<short>(prv+0.5);
          // assert(classValueMap[nameVector[iclass]]<probOut.size());
          // assert(classValueMap[nameVector[iclass]]>=0);
          // probOut[classValueMap[nameVector[iclass]]][icol]=static_cast<short>(prv+0.5);
        }
        entropy[icol]/=log(static_cast<double>(nclass))/log(2.0);
        entropy[icol]=static_cast<short>(100*entropy[icol]+0.5);
	if(active_opt.size()){
	  if(entropy[icol]>activePoints.back().value){
	    activePoints.back().value=entropy[icol];//replace largest value (last)
	    activePoints.back().posx=icol;
	    activePoints.back().posy=iline;
	    std::sort(activePoints.begin(),activePoints.end(),Decrease_PosValue());//sort in descending order (largest first, smallest last)
	    if(verbose_opt[0])
	      std::cout << activePoints.back().posx << " " << activePoints.back().posy << " " << activePoints.back().value << std::endl;
	  }
        }
      }//icol
      //----------------------------------- write output ------------------------------------------
      if(classBag_opt.size())
        for(unsigned int ibag=0;ibag<nbag;++ibag)
          classImageBag.writeData(classBag[ibag],iline,ibag);
      if(prob_opt.size()){
        for(unsigned int iclass=0;iclass<nclass;++iclass)
          probImage.writeData(probOut[iclass],iline,iclass);
      }
      if(entropy_opt.size()){
        entropyImage.writeData(entropy,iline);
      }
      classImageOut.writeData(classOut,iline);
      if(!verbose_opt[0]){
        progress=static_cast<float>(iline+1.0)/classImageOut.nrOfRow();
        pfnProgress(progress,pszMessage,pProgressArg);
      }
    }
    //write active learning points
    if(active_opt.size()){
      for(unsigned int iactive=0;iactive<activePoints.size();++iactive){
	std::map<string,double> pointMap;
	for(unsigned int iband=0;iband<testImage.nrOfBand();++iband){
	  double value;
	  testImage.readData(value,static_cast<unsigned int>(activePoints[iactive].posx),static_cast<unsigned int>(activePoints[iactive].posy),iband);
	  ostringstream fs;
	  fs << "B" << iband;
	  pointMap[fs.str()]=value;
	}
	pointMap[label_opt[0]]=0;
	double x, y;
	testImage.image2geo(activePoints[iactive].posx,activePoints[iactive].posy,x,y);
        std::string fieldname="id";//number of the point
	activeWriter.addPoint(x,y,pointMap,fieldname,++nactive);
      }
    }

    testImage.close();
    if(mask_opt.size())
      maskReader.close();
    if(priorimg_opt.size())
      priorReader.close();
    if(prob_opt.size())
      probImage.close();
    if(entropy_opt.size())
      entropyImage.close();
    if(classBag_opt.size())
      classImageBag.close();
    classImageOut.close();
  }
  else{//classify vector file
    cm.clearResults();
    //notice that fields have already been set by readDataImageOgr (taking into account appropriate bands)
    for(unsigned int ivalidation=0;ivalidation<input_opt.size();++ivalidation){
      if(output_opt.size())
        assert(output_opt.size()==input_opt.size());
      if(verbose_opt[0])
        cout << "opening img reader " << input_opt[ivalidation] << endl;
      imgReaderOgr.open(input_opt[ivalidation]);
      ImgWriterOgr imgWriterOgr;

      if(output_opt.size()){
	if(verbose_opt[0])
	  std::cout << "opening img writer and copying fields from img reader" << output_opt[ivalidation] << std::endl;
	imgWriterOgr.open(output_opt[ivalidation],imgReaderOgr);
      }
      if(verbose_opt[0])
	cout << "number of layers in input ogr file: " << imgReaderOgr.getLayerCount() << endl;
      for(unsigned int ilayer=0;ilayer<imgReaderOgr.getLayerCount();++ilayer){
	if(verbose_opt[0])
	  cout << "processing input layer " << ilayer << endl;
	if(output_opt.size()){
	  if(verbose_opt[0])
	    std::cout << "creating field class" << std::endl;
	  if(classValueMap.size())
	    imgWriterOgr.createField("class",OFTInteger,ilayer);
	  else
	    imgWriterOgr.createField("class",OFTString,ilayer);
	}
	unsigned int nFeatures=imgReaderOgr.getFeatureCount(ilayer);
	unsigned int ifeature=0;
	progress=0;
	pfnProgress(progress,pszMessage,pProgressArg);
	OGRFeature *poFeature;
	while( (poFeature = imgReaderOgr.getLayer(ilayer)->GetNextFeature()) != NULL ){
	  if(verbose_opt[0]>1)
	    cout << "feature " << ifeature << endl;
	  if( poFeature == NULL ){
	    cout << "Warning: could not read feature " << ifeature << " in layer " << imgReaderOgr.getLayerName(ilayer) << endl;
	    continue;
	  }
	  OGRFeature *poDstFeature = NULL;
	  if(output_opt.size()){
	    poDstFeature=imgWriterOgr.createFeature(ilayer);
	    if( poDstFeature->SetFrom( poFeature, TRUE ) != OGRERR_NONE ){
	      CPLError( CE_Failure, CPLE_AppDefined,
			"Unable to translate feature %d from layer %s.\n",
			poFeature->GetFID(), imgWriterOgr.getLayerName(ilayer).c_str() );
	      OGRFeature::DestroyFeature( poFeature );
	      OGRFeature::DestroyFeature( poDstFeature );
	    }
	  }
	  vector<float> validationPixel;
	  vector<float> validationFeature;
        
	  imgReaderOgr.readData(validationPixel,OFTReal,fields,poFeature,ilayer);
	  assert(validationPixel.size()==nband);
	  vector<float> probOut(nclass);//posterior prob for each class
	  for(unsigned int iclass=0;iclass<nclass;++iclass)
	    probOut[iclass]=0;
	  for(unsigned int ibag=0;ibag<nbag;++ibag){
	    for(unsigned int iband=0;iband<nband;++iband){
	      validationFeature.push_back((validationPixel[iband]-offset[ibag][iband])/scale[ibag][iband]);
	      if(verbose_opt[0]==2)
		std:: cout << " " << validationFeature.back();
	    }
	    if(verbose_opt[0]==2)
	      std::cout << std:: endl;
	    vector<float> result(nclass);
	    result=net[ibag].run(validationFeature);

	    if(verbose_opt[0]>1){
	      for(unsigned int iclass=0;iclass<result.size();++iclass)
		std::cout << result[iclass] << " ";
	      std::cout << std::endl;
	    }
	    //calculate posterior prob of bag 
	    for(unsigned int iclass=0;iclass<nclass;++iclass){
	      result[iclass]=(result[iclass]+1.0)/2.0;//bring back to scale [0,1]
	      switch(comb_opt[0]){
	      default:
	      case(0)://sum rule
		probOut[iclass]+=result[iclass]*priors[iclass];//add probabilities for each bag
              break;
	      case(1)://product rule
		probOut[iclass]*=pow(static_cast<float>(priors[iclass]),static_cast<float>(1.0-nbag)/nbag)*result[iclass];//multiply probabilities for each bag
		break;
	      case(2)://max rule
		if(priors[iclass]*result[iclass]>probOut[iclass])
		  probOut[iclass]=priors[iclass]*result[iclass];
		break;
	      }
	    }
	  }//for ibag
	  //search for max class prob
	  float maxBag=0;
	  string classOut="Unclassified";
	  for(unsigned int iclass=0;iclass<nclass;++iclass){
	    if(verbose_opt[0]>1)
	      std::cout << probOut[iclass] << " ";
	    if(probOut[iclass]>maxBag){
	      maxBag=probOut[iclass];
	      classOut=nameVector[iclass];
	    }
	  }
	  //look for class name
	  if(verbose_opt[0]>1){
	    if(classValueMap.size())
	      std::cout << "->" << classValueMap[classOut] << std::endl;
	    else	    
	      std::cout << "->" << classOut << std::endl;
	  }
	  if(output_opt.size()){
	    if(classValueMap.size())
	      poDstFeature->SetField("class",classValueMap[classOut]);
	    else	    
	      poDstFeature->SetField("class",classOut.c_str());
	    poDstFeature->SetFID( poFeature->GetFID() );
	  }
	  int labelIndex=poFeature->GetFieldIndex(label_opt[0].c_str());
	  if(labelIndex>=0){
	    string classRef=poFeature->GetFieldAsString(labelIndex);
	    if(classRef!="0"){
	      if(classValueMap.size())
		cm.incrementResult(type2string<short>(classValueMap[classRef]),type2string<short>(classValueMap[classOut]),1);
	      else
		cm.incrementResult(classRef,classOut,1);
	    }
	  }
	  CPLErrorReset();
	  if(output_opt.size()){
	    if(imgWriterOgr.createFeature(poDstFeature,ilayer) != OGRERR_NONE){
	      CPLError( CE_Failure, CPLE_AppDefined,
			"Unable to translate feature %d from layer %s.\n",
			poFeature->GetFID(), imgWriterOgr.getLayerName(ilayer).c_str() );
	      OGRFeature::DestroyFeature( poDstFeature );
	      OGRFeature::DestroyFeature( poDstFeature );
	    }
	  }
	  ++ifeature;
	  if(!verbose_opt[0]){
	    progress=static_cast<float>(ifeature+1.0)/nFeatures;
	    pfnProgress(progress,pszMessage,pProgressArg);
	  }
	  OGRFeature::DestroyFeature( poFeature );
	  OGRFeature::DestroyFeature( poDstFeature );
	}//get next feature
      }//next layer
      imgReaderOgr.close();
      if(output_opt.size())
      imgWriterOgr.close();
    }
    if(cm.nReference()){
      std::cout << cm << std::endl;
      // cout << "class #samples userAcc prodAcc" << endl;
      // double se95_ua=0;
      // double se95_pa=0;
      // double se95_oa=0;
      // double dua=0;
      // double dpa=0;
      // double doa=0;
      // for(short iclass=0;iclass<cm.nClasses();++iclass){
      // 	dua=cm.ua_pct(cm.getClass(iclass),&se95_ua);
      // 	dpa=cm.pa_pct(cm.getClass(iclass),&se95_pa);
      // 	cout << cm.getClass(iclass) << " " << cm.nReference(cm.getClass(iclass)) << " " << dua << " (" << se95_ua << ")" << " " << dpa << " (" << se95_pa << ")" << endl;
      // }
      // std::cout << "Kappa: " << cm.kappa() << std::endl;
      // doa=cm.oa_pct(&se95_oa);
      // std::cout << "Overall Accuracy: " << doa << " (" << se95_oa << ")"  << std::endl;
    }
  }
  try{
    if(active_opt.size())
      activeWriter.close();
  }
  catch(string errorString){
    std::cerr << "Error: errorString" << std::endl;
  }
  return 0;
}
