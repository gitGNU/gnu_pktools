/**********************************************************************
pkoptsvm.cc: program to optimize parameters for support vector machine classifier pksvm
Copyright (C) 2008-2014 Pieter Kempeneers

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
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <math.h>
//#include <nlopt.hpp>
#include "base/Optionpk.h"
#include "algorithms/ConfusionMatrix.h"
#include "algorithms/FeatureSelector.h"
//#include "algorithms/OptFactory.h"
#include "algorithms/CostFactorySVM.h"
#include "algorithms/svm.h"
#include "imageclasses/ImgReaderOgr.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/******************************************************************************/
/*! \page pkoptsvm pkoptsvm
 program to optimize parameters for support vector machine classifier pksvm
## SYNOPSIS

<code>
  Usage: pkoptsvm -t training
</code>

<code>

  Options: [-cc startvalue -cc endvalue] [-g startvalue -g endvalue] [-stepcc stepsize] [-stepg stepsize]

  Advanced options:
</code>

\section pkoptsvm_description Description

The support vector machine depends on several parameters. Ideally, these parameters should be optimized for each classification problem. In case of a radial basis kernel function, two important parameters are \em{cost} and \em{gamma}. The utility pkoptsvm can optimize these two parameters, based on an accuracy assessment (the Kappa value). If an input test set (-i) is provided, it is used for the accuracy assessment. If not, the accuracy assessment is based on a cross validation (-cv) of the training sample.

The optimization routine uses a grid search. The initial and final values of the parameters can be set with -cc startvalue -cc endvalue and -g startvalue -g endvalue for cost and gamma respectively. The search uses a multiplicative step for iterating the parameters (set with the options -stepcc and -stepg). An often used approach is to define a relatively large multiplicative step first (e.g 10) to obtain an initial estimate for both parameters. The estimate can then be optimized by defining a smaller step (>1) with constrained start and end values for the parameters cost and gamma.

\section pkoptsvm_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | t      | training             | std::string |       |training vector file. A single vector file contains all training features (must be set as: b0, b1, b2,...) for all classes (class numbers identified by label option). | 
 | cc     | ccost                | float | 1     |min and max boundaries the parameter C of C-SVC, epsilon-SVR, and nu-SVR (optional: initial value) | 
 | g      | gamma                | float | 0     |min max boundaries for gamma in kernel function (optional: initial value) | 
 | stepcc | stepcc               | double | 2     |multiplicative step for ccost in GRID search | 
 | stepg  | stepg                | double | 2     |multiplicative step for gamma in GRID search | 
 | i      | input                | std::string |       |input test vector file | 
 | tln    | tln                  | std::string |       |training layer name(s) | 
 | label  | label                | std::string | label |identifier for class label in training vector file. | 
 | bal    | balance              | unsigned int | 0     |balance the input data to this number of samples for each class | 
 | random | random               | bool | true  |in case of balance, randomize input data | 
 | min    | min                  | int  | 0     |if number of training pixels is less then min, do not take this class into account | 
 | b      | band                 | unsigned short |      |band index (starting from 0, either use band option or use start to end) | 
 | sband  | startband            | unsigned short |      |Start band sequence number | 
 | eband  | endband              | unsigned short |      |End band sequence number   | 
 | offset | offset               | double | 0     |offset value for each spectral band input features: refl[band]=(DN[band]-offset[band])/scale[band] | 
 | scale  | scale                | double | 0     |scale value for each spectral band input features: refl=(DN[band]-offset[band])/scale[band] (use 0 if scale min and max in each band to -1.0 and 1.0) | 
 | svmt   | svmtype              | std::string | C_SVC |type of SVM (C_SVC, nu_SVC,one_class, epsilon_SVR, nu_SVR) | 
 | kt     | kerneltype           | std::string | radial |type of kernel function (linear,polynomial,radial,sigmoid)  | 
 | kd     | kd                   | unsigned short | 3     |degree in kernel function | 
 | c0     | coef0                | float | 0     |coef0 in kernel function | 
 | nu     | nu                   | float | 0.5   |the parameter nu of nu-SVC, one-class SVM, and nu-SVR | 
 | eloss  | eloss                | float | 0.1   |the epsilon in loss function of epsilon-SVR | 
 | cache  | cache                | int  | 100   |cache memory size in MB | 
 | etol   | etol                 | float | 0.001 |the tolerance of termination criterion | 
 | shrink | shrink               | bool | false |whether to use the shrinking heuristics | 
 | pe     | probest              | bool | true  |whether to train a SVC or SVR model for probability estimates | 
 | cv     | cv                   | unsigned short | 2     |n-fold cross validation mode | 
 | cf     | cf                   | bool | false |use Overall Accuracy instead of kappa | 
 | maxit  | maxit                | unsigned int | 500   |maximum number of iterations | 
 | tol    | tolerance            | double | 0.0001 |relative tolerance for stopping criterion | 
 | c      | class                | std::string |       |list of class names. | 
 | r      | reclass              | short |       |list of class values (use same order as in class opt). | 

Usage: pkoptsvm -t training


**/

using namespace std;

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
                                    //declare objective function
double objFunction(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data);

//global parameters used in objective function
map<string,short> classValueMap;
vector<std::string> nameVector;
vector<unsigned int> nctraining;
vector<unsigned int> nctest;
Optionpk<std::string> svm_type_opt("svmt", "svmtype", "type of SVM (C_SVC, nu_SVC,one_class, epsilon_SVR, nu_SVR)","C_SVC");
Optionpk<std::string> kernel_type_opt("kt", "kerneltype", "type of kernel function (linear,polynomial,radial,sigmoid) ","radial");
Optionpk<unsigned short> kernel_degree_opt("kd", "kd", "degree in kernel function",3);
Optionpk<float> coef0_opt("c0", "coef0", "coef0 in kernel function",0);
Optionpk<float> nu_opt("nu", "nu", "the parameter nu of nu-SVC, one-class SVM, and nu-SVR",0.5);
Optionpk<float> epsilon_loss_opt("eloss", "eloss", "the epsilon in loss function of epsilon-SVR",0.1);
Optionpk<int> cache_opt("cache", "cache", "cache memory size in MB",100);
Optionpk<float> epsilon_tol_opt("etol", "etol", "the tolerance of termination criterion",0.001);
Optionpk<bool> shrinking_opt("shrink", "shrink", "whether to use the shrinking heuristics",false);
Optionpk<bool> prob_est_opt("pe", "probest", "whether to train a SVC or SVR model for probability estimates",true,2);
Optionpk<bool> costfunction_opt("cf", "cf", "use Overall Accuracy instead of kappa",false);
// Optionpk<bool> weight_opt("wi", "wi", "set the parameter C of class i to weight*C, for C-SVC",true);
Optionpk<unsigned short> cv_opt("cv", "cv", "n-fold cross validation mode",2);
Optionpk<string> classname_opt("c", "class", "list of class names."); 
Optionpk<short> classvalue_opt("r", "reclass", "list of class values (use same order as in class opt)."); 
Optionpk<short> verbose_opt("v", "verbose", "use 1 to output intermediate results for plotting",0,2);

double objFunction(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data){

  assert(grad.empty());
  vector<Vector2d<float> > *tf=reinterpret_cast<vector<Vector2d<float> >*> (my_func_data);
  float ccost=x[0];
  float gamma=x[1];
  double error=1.0/epsilon_tol_opt[0];
  double kappa=1.0;
  double oa=1.0;

  CostFactorySVM costfactory(svm_type_opt[0], kernel_type_opt[0], kernel_degree_opt[0], gamma, coef0_opt[0], ccost, nu_opt[0],  epsilon_loss_opt[0], cache_opt[0], epsilon_tol_opt[0], shrinking_opt[0], prob_est_opt[0], cv_opt[0], verbose_opt[0]);

  assert(tf->size());
  // if(nctest>0)
  //   costfactory.setCv(0);

  costfactory.setCv(cv_opt[0]);

  if(classname_opt.size()){
    assert(classname_opt.size()==classvalue_opt.size());
    for(int iclass=0;iclass<classname_opt.size();++iclass)
      costfactory.setClassValueMap(classname_opt[iclass],classvalue_opt[iclass]);
  }
  //set names in confusion matrix using nameVector
  costfactory.setNameVector(nameVector);
  // vector<string> nameVector=costfactory.getNameVector();
  for(int iname=0;iname<nameVector.size();++iname){
    if(costfactory.getClassValueMap().empty()){
      costfactory.pushBackClassName(nameVector[iname]);
      // cm.pushBackClassName(nameVector[iname]);
    }
    else if(costfactory.getClassIndex(type2string<short>((costfactory.getClassValueMap())[nameVector[iname]]))<0)
      costfactory.pushBackClassName(type2string<short>((costfactory.getClassValueMap())[nameVector[iname]]));
  }

  costfactory.setNcTraining(nctraining);
  costfactory.setNcTest(nctest);

  kappa=costfactory.getCost(*tf);
  return(kappa);
}

int main(int argc, char *argv[])
{
  map<short,int> reclassMap;
  vector<int> vreclass;
  Optionpk<string> training_opt("t", "training", "training vector file. A single vector file contains all training features (must be set as: b0, b1, b2,...) for all classes (class numbers identified by label option)."); 
  Optionpk<float> ccost_opt("cc", "ccost", "min and max boundaries the parameter C of C-SVC, epsilon-SVR, and nu-SVR (optional: initial value)",1);
  Optionpk<float> gamma_opt("g", "gamma", "min max boundaries for gamma in kernel function (optional: initial value)",0);
  Optionpk<double> stepcc_opt("stepcc","stepcc","multiplicative step for ccost in GRID search",2);
  Optionpk<double> stepg_opt("stepg","stepg","multiplicative step for gamma in GRID search",2);
  Optionpk<string> input_opt("i", "input", "input test vector file"); 
  Optionpk<string> tlayer_opt("tln", "tln", "training layer name(s)");
  Optionpk<string> label_opt("label", "label", "identifier for class label in training vector file.","label"); 
  // Optionpk<unsigned short> reclass_opt("\0", "rc", "reclass code (e.g. --rc=12 --rc=23 to reclass first two classes to 12 and 23 resp.).", 0);
  Optionpk<unsigned int> balance_opt("bal", "balance", "balance the input data to this number of samples for each class", 0);
  Optionpk<bool> random_opt("random","random", "in case of balance, randomize input data", true);
  Optionpk<int> minSize_opt("min", "min", "if number of training pixels is less then min, do not take this class into account", 0);
  Optionpk<unsigned short> band_opt("b", "band", "band index (starting from 0, either use band option or use start to end)");
  Optionpk<unsigned short> bstart_opt("sband", "startband", "Start band sequence number"); 
  Optionpk<unsigned short> bend_opt("eband", "endband", "End band sequence number"); 
  Optionpk<double> offset_opt("offset", "offset", "offset value for each spectral band input features: refl[band]=(DN[band]-offset[band])/scale[band]", 0.0);
  Optionpk<double> scale_opt("scale", "scale", "scale value for each spectral band input features: refl=(DN[band]-offset[band])/scale[band] (use 0 if scale min and max in each band to -1.0 and 1.0)", 0.0);
  Optionpk<unsigned int> maxit_opt("maxit","maxit","maximum number of iterations",500);
//Optionpk<string> algorithm_opt("a", "algorithm", "GRID, or any optimization algorithm from http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms","GRID"); 
  Optionpk<double> tolerance_opt("tol","tolerance","relative tolerance for stopping criterion",0.0001);

  input_opt.setHide(1);
  tlayer_opt.setHide(1);
  label_opt.setHide(1);
  balance_opt.setHide(1);
  random_opt.setHide(1);
  minSize_opt.setHide(1);
  band_opt.setHide(1);
  bstart_opt.setHide(1);
  bend_opt.setHide(1);
  offset_opt.setHide(1);
  scale_opt.setHide(1);
  svm_type_opt.setHide(1);
  kernel_type_opt.setHide(1);
  kernel_degree_opt.setHide(1);
  coef0_opt.setHide(1);
  nu_opt.setHide(1);
  epsilon_loss_opt.setHide(1);
  cache_opt.setHide(1);
  epsilon_tol_opt.setHide(1);
  shrinking_opt.setHide(1);
  prob_est_opt.setHide(1);
  cv_opt.setHide(1);
  costfunction_opt.setHide(1);
  maxit_opt.setHide(1);
  tolerance_opt.setHide(1);
//  algorithm_opt.setHide(1);
  classname_opt.setHide(1);
  classvalue_opt.setHide(1);
  
  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=training_opt.retrieveOption(argc,argv);
    ccost_opt.retrieveOption(argc,argv);
    gamma_opt.retrieveOption(argc,argv);
    stepcc_opt.retrieveOption(argc,argv);
    stepg_opt.retrieveOption(argc,argv);
    input_opt.retrieveOption(argc,argv);
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
    svm_type_opt.retrieveOption(argc,argv);
    kernel_type_opt.retrieveOption(argc,argv);
    kernel_degree_opt.retrieveOption(argc,argv);
    coef0_opt.retrieveOption(argc,argv);
    nu_opt.retrieveOption(argc,argv);
    epsilon_loss_opt.retrieveOption(argc,argv);
    cache_opt.retrieveOption(argc,argv);
    epsilon_tol_opt.retrieveOption(argc,argv);
    shrinking_opt.retrieveOption(argc,argv);
    prob_est_opt.retrieveOption(argc,argv);
    cv_opt.retrieveOption(argc,argv);
    costfunction_opt.retrieveOption(argc,argv);
    maxit_opt.retrieveOption(argc,argv);
    tolerance_opt.retrieveOption(argc,argv);
//    algorithm_opt.retrieveOption(argc,argv);
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
    cout << "Usage: pkoptsvm -t training" << endl;
    cout << endl;
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }

  assert(training_opt.size());
  if(input_opt.size())
    cv_opt[0]=0;

  if(verbose_opt[0]>=1){
    if(input_opt.size())
      std::cout << "input filename: " << input_opt[0] << std::endl;
    std::cout << "training vector file: " << std::endl;
    for(int ifile=0;ifile<training_opt.size();++ifile)
      std::cout << training_opt[ifile] << std::endl;
    std::cout << "verbose: " << verbose_opt[0] << std::endl;
  }

  unsigned int totalSamples=0;
  unsigned int totalTestSamples=0;

  unsigned short nclass=0;
  int nband=0;
  int startBand=2;//first two bands represent X and Y pos

  vector<double> offset;
  vector<double> scale;
  vector< Vector2d<float> > trainingPixels;//[class][sample][band]
  vector< Vector2d<float> > testPixels;//[class][sample][band]

  // if(priors_opt.size()>1){//priors from argument list
  //   priors.resize(priors_opt.size());
  //   double normPrior=0;
  //   for(int iclass=0;iclass<priors_opt.size();++iclass){
  //     priors[iclass]=priors_opt[iclass];
  //     normPrior+=priors[iclass];
  //   }
  //   //normalize
  //   for(int iclass=0;iclass<priors_opt.size();++iclass)
  //     priors[iclass]/=normPrior;
  // }

  //convert start and end band options to vector of band indexes
  try{
    if(bstart_opt.size()){
      if(bend_opt.size()!=bstart_opt.size()){
	string errorstring="Error: options for start and end band indexes must be provided as pairs, missing end band";
	throw(errorstring);
      }
      band_opt.clear();
      for(int ipair=0;ipair<bstart_opt.size();++ipair){
	if(bend_opt[ipair]<=bstart_opt[ipair]){
	  string errorstring="Error: index for end band must be smaller then start band";
	  throw(errorstring);
	}
	for(int iband=bstart_opt[ipair];iband<=bend_opt[ipair];++iband)
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

  // map<string,short> classValueMap;//global variable for now (due to getCost)
  if(classname_opt.size()){
    assert(classname_opt.size()==classvalue_opt.size());
    for(int iclass=0;iclass<classname_opt.size();++iclass)
      classValueMap[classname_opt[iclass]]=classvalue_opt[iclass];
  }

  //----------------------------------- Training -------------------------------
  struct svm_problem prob;
  vector<string> fields;
  //organize training data
  trainingPixels.clear();
  testPixels.clear();
  map<string,Vector2d<float> > trainingMap;
  map<string,Vector2d<float> > testMap;
  if(verbose_opt[0]>=1)
    std::cout << "reading training file " << training_opt[0] << std::endl;
  try{
    ImgReaderOgr trainingReader(training_opt[0]);
    if(band_opt.size()){
      totalSamples=trainingReader.readDataImageOgr(trainingMap,fields,band_opt,label_opt[0],tlayer_opt,verbose_opt[0]);
      if(input_opt.size()){
	ImgReaderOgr inputReader(input_opt[0]);
	totalTestSamples=inputReader.readDataImageOgr(testMap,fields,band_opt,label_opt[0],tlayer_opt,verbose_opt[0]);
	inputReader.close();
      }
    }
    else{
      totalSamples=trainingReader.readDataImageOgr(trainingMap,fields,0,0,label_opt[0],tlayer_opt,verbose_opt[0]);
      if(input_opt.size()){
	ImgReaderOgr inputReader(input_opt[0]);
	totalTestSamples=inputReader.readDataImageOgr(testMap,fields,0,0,label_opt[0],tlayer_opt,verbose_opt[0]);
	inputReader.close();
      }
      trainingReader.close();
    }
    if(trainingMap.size()<2){
      // map<string,Vector2d<float> >::iterator mapit=trainingMap.begin();
      // while(mapit!=trainingMap.end())
      // 	cerr << mapit->first << " -> " << classValueMap[mapit->first] << std::endl;
      string errorstring="Error: could not read at least two classes from training input file";
      throw(errorstring);
    }
    if(input_opt.size()&&testMap.size()<2){
      string errorstring="Error: could not read at least two classes from test input file";
      throw(errorstring);
    }
  }
  catch(string error){
    cerr << error << std::endl;
    exit(1);
  }
  catch(...){
    cerr << "error caught" << std::endl;
    exit(1);
  }
  //todo delete class 0 ?
  // if(verbose_opt[0]>=1)
  //   std::cout << "erasing class 0 from training set (" << trainingMap[0].size() << " from " << totalSamples << ") samples" << std::endl;
  // totalSamples-=trainingMap[0].size();
  // trainingMap.erase(0);

  if(verbose_opt[0]>1)
    std::cout << "training pixels: " << std::endl;
  map<string,Vector2d<float> >::iterator mapit;
  mapit=trainingMap.begin();
  while(mapit!=trainingMap.end()){
    if(classValueMap.size()){
      //check if name in training is covered by classname_opt (values can not be 0)
      if(classValueMap[mapit->first]>0){
	if(verbose_opt[0])
	  std::cout << mapit->first << " -> " << classValueMap[mapit->first] << std::endl;
      }
      else{
	std::cerr << "Error: names in classname option are not complete, please check names in training vector and make sure classvalue is > 0" << std::endl;
	exit(1);
      }
    }    
    //delete small classes
    if((mapit->second).size()<minSize_opt[0]){
      trainingMap.erase(mapit);
      continue;
    }
    nameVector.push_back(mapit->first);
    trainingPixels.push_back(mapit->second);
    if(verbose_opt[0]>1)
      std::cout << mapit->first << ": " << (mapit->second).size() << " samples" << std::endl;
    // trainingPixels.push_back(mapit->second); ??
    // ++iclass;
    ++mapit;
  }
  nclass=trainingPixels.size();
  if(classname_opt.size())
    assert(nclass==classname_opt.size());
  nband=trainingPixels[0][0].size()-2;//X and Y//trainingPixels[0][0].size();

  mapit=testMap.begin();
  while(mapit!=testMap.end()){
    if(classValueMap.size()){
      //check if name in test is covered by classname_opt (values can not be 0)
      if(classValueMap[mapit->first]>0){
	;//ok, no need to print to std::cout 
      }
      else{
	std::cerr << "Error: names in classname option are not complete, please check names in test vector and make sure classvalue is > 0" << std::endl;
	exit(1);
      }
    }    
    //no need to delete small classes for test sample
    testPixels.push_back(mapit->second);
    if(verbose_opt[0]>1)
      std::cout << mapit->first << ": " << (mapit->second).size() << " samples" << std::endl;
    ++mapit;
  }
  if(input_opt.size()){
    assert(nclass==testPixels.size());
    assert(nband=testPixels[0][0].size()-2);//X and Y//testPixels[0][0].size();
    assert(!cv_opt[0]);
  }

  //do not remove outliers here: could easily be obtained through ogr2ogr -where 'B2<110' output.shp input.shp
  //balance training data
  if(balance_opt[0]>0){
    if(random_opt[0])
      srand(time(NULL));
    totalSamples=0;
    for(int iclass=0;iclass<nclass;++iclass){
      if(trainingPixels[iclass].size()>balance_opt[0]){
        while(trainingPixels[iclass].size()>balance_opt[0]){
          int index=rand()%trainingPixels[iclass].size();
          trainingPixels[iclass].erase(trainingPixels[iclass].begin()+index);
        }
      }
      else{
        int oldsize=trainingPixels[iclass].size();
        for(int isample=trainingPixels[iclass].size();isample<balance_opt[0];++isample){
          int index = rand()%oldsize;
          trainingPixels[iclass].push_back(trainingPixels[iclass][index]);
        }
      }
      totalSamples+=trainingPixels[iclass].size();
    }
    assert(totalSamples==nclass*balance_opt[0]);
  }

  //no need to balance test sample    
  //set scale and offset
  offset.resize(nband);
  scale.resize(nband);
  if(offset_opt.size()>1)
    assert(offset_opt.size()==nband);
  if(scale_opt.size()>1)
    assert(scale_opt.size()==nband);
  for(int iband=0;iband<nband;++iband){
    if(verbose_opt[0]>1)
      std::cout << "scaling for band" << iband << std::endl;
    offset[iband]=(offset_opt.size()==1)?offset_opt[0]:offset_opt[iband];
    scale[iband]=(scale_opt.size()==1)?scale_opt[0]:scale_opt[iband];
    //search for min and maximum
    if(scale[iband]<=0){
      float theMin=trainingPixels[0][0][iband+startBand];
      float theMax=trainingPixels[0][0][iband+startBand];
      for(int iclass=0;iclass<nclass;++iclass){
        for(int isample=0;isample<trainingPixels[iclass].size();++isample){
          if(theMin>trainingPixels[iclass][isample][iband+startBand])
            theMin=trainingPixels[iclass][isample][iband+startBand];
          if(theMax<trainingPixels[iclass][isample][iband+startBand])
            theMax=trainingPixels[iclass][isample][iband+startBand];
        }
      }
      offset[iband]=theMin+(theMax-theMin)/2.0;
      scale[iband]=(theMax-theMin)/2.0;
      if(verbose_opt[0]>1){
        std::cout << "Extreme image values for band " << iband << ": [" << theMin << "," << theMax << "]" << std::endl;
        std::cout << "Using offset, scale: " << offset[iband] << ", " << scale[iband] << std::endl;
        std::cout << "scaled values for band " << iband << ": [" << (theMin-offset[iband])/scale[iband] << "," << (theMax-offset[iband])/scale[iband] << "]" << std::endl;
      }
    }
  }

  // if(priors_opt.size()==1){//default: equal priors for each class
  //   priors.resize(nclass);
  //   for(int iclass=0;iclass<nclass;++iclass)
  //     priors[iclass]=1.0/nclass;
  // }
  // assert(priors_opt.size()==1||priors_opt.size()==nclass);
    
  if(verbose_opt[0]>=1){
    std::cout << "number of bands: " << nband << std::endl;
    std::cout << "number of classes: " << nclass << std::endl;
    // std::cout << "priors:";
    // for(int iclass=0;iclass<nclass;++iclass)
    //   std::cout << " " << priors[iclass];
    // std::cout << std::endl;
  }

  //Calculate features of training (and test) set
  nctraining.resize(nclass);
  nctest.resize(nclass);
  vector< Vector2d<float> > trainingFeatures(nclass);
  for(int iclass=0;iclass<nclass;++iclass){
    if(verbose_opt[0]>=1)
      std::cout << "calculating features for class " << iclass << std::endl;
    nctraining[iclass]=trainingPixels[iclass].size();
    if(verbose_opt[0]>=1)
      std::cout << "nctraining[" << iclass << "]: " << nctraining[iclass] << std::endl;
    if(testPixels.size()>iclass){
      nctest[iclass]=testPixels[iclass].size();
      if(verbose_opt[0]>=1){
	std::cout << "nctest[" << iclass << "]: " << nctest[iclass] << std::endl;
      }
    }
    else
      nctest[iclass]=0;
    // trainingFeatures[iclass].resize(nctraining[iclass]);
    trainingFeatures[iclass].resize(nctraining[iclass]+nctest[iclass]);
    for(int isample=0;isample<nctraining[iclass];++isample){
      //scale pixel values according to scale and offset!!!
      for(int iband=0;iband<nband;++iband){
        assert(trainingPixels[iclass].size()>isample);
        assert(trainingPixels[iclass][isample].size()>iband+startBand);
        assert(offset.size()>iband);
        assert(scale.size()>iband);
        float value=trainingPixels[iclass][isample][iband+startBand];
        trainingFeatures[iclass][isample].push_back((value-offset[iband])/scale[iband]);
      }
    }
    // assert(trainingFeatures[iclass].size()==nctraining[iclass]);
    for(int isample=0;isample<nctest[iclass];++isample){
      //scale pixel values according to scale and offset!!!
      for(int iband=0;iband<nband;++iband){
        assert(testPixels[iclass].size()>isample);
        assert(testPixels[iclass][isample].size()>iband+startBand);
        assert(offset.size()>iband);
        assert(scale.size()>iband);
        float value=testPixels[iclass][isample][iband+startBand];
        // testFeatures[iclass][isample].push_back((value-offset[iband])/scale[iband]);
        trainingFeatures[iclass][nctraining[iclass]+isample].push_back((value-offset[iband])/scale[iband]);
      }
    }
    assert(trainingFeatures[iclass].size()==nctraining[iclass]+nctest[iclass]);
  }

  assert(ccost_opt.size()>1);//must have boundaries at least (initial value is optional)
  if(ccost_opt.size()<3)//create initial value
    ccost_opt.push_back(sqrt(ccost_opt[0]*ccost_opt[1]));
  assert(gamma_opt.size()>1);//must have boundaries at least (initial value is optional)
  if(gamma_opt.size()<3)//create initial value
    gamma_opt.push_back(sqrt(gamma_opt[0]*gamma_opt[1]));//will be translated to 1.0/nFeatures
  assert(ccost_opt.size()==3);//min, init, max
  assert(gamma_opt.size()==3);//min, init, max
  assert(gamma_opt[0]<gamma_opt[1]);
  assert(gamma_opt[0]<gamma_opt[2]);
  assert(gamma_opt[2]<gamma_opt[1]);
  assert(ccost_opt[0]<ccost_opt[1]);
  assert(ccost_opt[0]<ccost_opt[2]);
  assert(ccost_opt[2]<ccost_opt[1]);

  std::vector<double> x(2);
//  if(algorithm_opt[0]=="GRID"){
  if (1){
    // double minError=1000;
    // double minCost=0;
    // double minGamma=0;
    double maxKappa=0;
    double maxCost=0;
    double maxGamma=0;
    const char* pszMessage;
    void* pProgressArg=NULL;
    GDALProgressFunc pfnProgress=GDALTermProgress;
    double progress=0;
    if(!verbose_opt[0])
      pfnProgress(progress,pszMessage,pProgressArg);
    double ncost=log(ccost_opt[1])/log(stepcc_opt[0])-log(ccost_opt[0])/log(stepcc_opt[0]);
    double ngamma=log(gamma_opt[1])/log(stepg_opt[0])-log(gamma_opt[0])/log(stepg_opt[0]);
    for(double ccost=ccost_opt[0];ccost<=ccost_opt[1];ccost*=stepcc_opt[0]){
      for(double gamma=gamma_opt[0];gamma<=gamma_opt[1];gamma*=stepg_opt[0]){
	x[0]=ccost;
	x[1]=gamma;
	std::vector<double> theGrad;
	double kappa=0;
	kappa=objFunction(x,theGrad,&trainingFeatures);
	if(kappa>maxKappa){
	  maxKappa=kappa;
	  maxCost=ccost;
	  maxGamma=gamma;
	}
	if(verbose_opt[0])
	  std::cout << ccost << " " << gamma << " " << kappa<< std::endl;
	progress+=1.0/ncost/ngamma;
	if(!verbose_opt[0])
	  pfnProgress(progress,pszMessage,pProgressArg);
      }
    }
    progress=1.0;
    if(!verbose_opt[0])
      pfnProgress(progress,pszMessage,pProgressArg);
    x[0]=maxCost;
    x[1]=maxGamma;
  }
  //else{
  //  nlopt::opt optimizer=OptFactory::getOptimizer(algorithm_opt[0],2);
  //  if(verbose_opt[0]>1)
  //    std::cout << "optimization algorithm: " << optimizer.get_algorithm_name() << "..." << std::endl;
  //  std::vector<double> lb(2);
  //  std::vector<double> init(2);
  //  std::vector<double> ub(2);

  //  lb[0]=ccost_opt[0];
  //  lb[1]=(gamma_opt[0]>0)? gamma_opt[0] : 1.0/trainingFeatures[0][0].size();
  //  init[0]=ccost_opt[2];
  //  init[1]=(gamma_opt[2]>0)? gamma_opt[1] : 1.0/trainingFeatures[0][0].size();
  //  ub[0]=ccost_opt[1];
  //  ub[1]=(gamma_opt[1]>0)? gamma_opt[1] : 1.0/trainingFeatures[0][0].size();
  //  // optimizer.set_min_objective(objFunction, &trainingFeatures);
  //  optimizer.set_max_objective(objFunction, &trainingFeatures);
  //  optimizer.set_lower_bounds(lb);
  //  optimizer.set_upper_bounds(ub);
  //  if(verbose_opt[0]>1)
  //    std::cout << "set stopping criteria" << std::endl;
  //  //set stopping criteria
  //  if(maxit_opt[0])
  //    optimizer.set_maxeval(maxit_opt[0]);
  //  else
  //    optimizer.set_xtol_rel(tolerance_opt[0]);
  //  double minf=0;
  //  x=init;
  //  try{
  //    optimizer.optimize(x, minf);
  //  }
  //  catch(string error){
  //    cerr << error << std::endl;
  //    exit(1);
  //  }
  //  catch (exception& e){
  //    cout << e.what() << endl;
  //  }
  //  catch(...){
  //    cerr << "error caught" << std::endl;
  //    exit(1);
  //  }

  //  double ccost=x[0];
  //  double gamma=x[1];
  //  if(verbose_opt[0])
  //    std::cout << "optimized with " << optimizer.get_algorithm_name() << "..." << std::endl;
  //}
  std::cout << " --ccost " << x[0];
  std::cout << " --gamma " << x[1];
  std::cout << std::endl;
}
