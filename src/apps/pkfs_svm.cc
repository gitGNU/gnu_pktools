/**********************************************************************
pkfs_svm.cc: feature selection for svm classifier
Copyright (C) 2008-2012 Pieter Kempeneers

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
#include <string>
#include <map>
#include <algorithm>
#include "base/Optionpk.h"
#include "algorithms/ConfusionMatrix.h"
#include "algorithms/FeatureSelector.h"
#include "algorithms/svm.h"
#include "imageclasses/ImgReaderOgr.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace svm{
  enum SVM_TYPE {C_SVC=0, nu_SVC=1,one_class=2, epsilon_SVR=3, nu_SVR=4};
  enum KERNEL_TYPE {linear=0,polynomial=1,radial=2,sigmoid=3};
}

enum SelectorValue  { NA=0, SFFS=1, SFS=2, SBS=3, BFS=4 };

using namespace std;

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

//global parameters used in cost function getCost
map<string,short> classValueMap;
vector<std::string> nameVector;
vector<unsigned int> nctraining;
vector<unsigned int> nctest;
Optionpk<std::string> svm_type_opt("svmt", "svmtype", "type of SVM (C_SVC, nu_SVC,one_class, epsilon_SVR, nu_SVR)","C_SVC");
Optionpk<std::string> kernel_type_opt("kt", "kerneltype", "type of kernel function (linear,polynomial,radial,sigmoid) ","radial");
Optionpk<unsigned short> kernel_degree_opt("kd", "kd", "degree in kernel function",3);
Optionpk<float> gamma_opt("g", "gamma", "gamma in kernel function",0);
Optionpk<float> coef0_opt("c0", "coef0", "coef0 in kernel function",0);
Optionpk<float> ccost_opt("cc", "ccost", "the parameter C of C-SVC, epsilon-SVR, and nu-SVR",1);
Optionpk<float> nu_opt("nu", "nu", "the parameter nu of nu-SVC, one-class SVM, and nu-SVR",0.5);
Optionpk<float> epsilon_loss_opt("eloss", "eloss", "the epsilon in loss function of epsilon-SVR",0.1);
Optionpk<int> cache_opt("cache", "cache", "cache memory size in MB",100);
Optionpk<float> epsilon_tol_opt("etol", "etol", "the tolerance of termination criterion",0.001);
Optionpk<bool> shrinking_opt("shrink", "shrink", "whether to use the shrinking heuristics",false);
Optionpk<bool> prob_est_opt("pe", "probest", "whether to train a SVC or SVR model for probability estimates",true,2);
// Optionpk<bool> weight_opt("wi", "wi", "set the parameter C of class i to weight*C, for C-SVC",true);
Optionpk<unsigned short> cv_opt("cv", "cv", "n-fold cross validation mode",2);
Optionpk<string> classname_opt("c", "class", "list of class names."); 
Optionpk<short> classvalue_opt("r", "reclass", "list of class values (use same order as in classname opt."); 
Optionpk<short> verbose_opt("v", "verbose", "set to: 0 (results only), 1 (confusion matrix), 2 (debug)",0);

double getCost(const vector<Vector2d<float> > &trainingFeatures)
{
  std::map<std::string, svm::SVM_TYPE> svmMap;

  svmMap["C_SVC"]=svm::C_SVC;
  svmMap["nu_SVC"]=svm::nu_SVC;
  svmMap["one_class"]=svm::one_class;
  svmMap["epsilon_SVR"]=svm::epsilon_SVR;
  svmMap["nu_SVR"]=svm::nu_SVR;

  std::map<std::string, svm::KERNEL_TYPE> kernelMap;

  kernelMap["linear"]=svm::linear;
  kernelMap["polynomial"]=svm::polynomial;
  kernelMap["radial"]=svm::radial;
  kernelMap["sigmoid;"]=svm::sigmoid;

  unsigned short nclass=trainingFeatures.size();
  unsigned int ntraining=0;
  unsigned int ntest=0;
  for(int iclass=0;iclass<nclass;++iclass){
    ntraining+=nctraining[iclass];
    ntest+=nctest[iclass];
  }
  if(ntest)
    assert(!cv_opt[0]);
  if(!cv_opt[0])
    assert(ntest);
  unsigned short nFeatures=trainingFeatures[0][0].size();

  struct svm_parameter param;
  param.svm_type = svmMap[svm_type_opt[0]];
  param.kernel_type = kernelMap[kernel_type_opt[0]];
  param.degree = kernel_degree_opt[0];
  param.gamma = (gamma_opt[0]>0)? gamma_opt[0] : 1.0/nFeatures;
  param.coef0 = coef0_opt[0];
  param.nu = nu_opt[0];
  param.cache_size = cache_opt[0];
  param.C = ccost_opt[0];
  param.eps = epsilon_tol_opt[0];
  param.p = epsilon_loss_opt[0];
  param.shrinking = (shrinking_opt[0])? 1 : 0;
  param.probability = (prob_est_opt[0])? 1 : 0;
  param.nr_weight = 0;//not used: I use priors and balancing
  param.weight_label = NULL;
  param.weight = NULL;
  param.verbose=(verbose_opt[0]>1)? true:false;
  struct svm_model* svm;
  struct svm_problem prob;
  struct svm_node* x_space;

  prob.l=ntraining;
  prob.y = Malloc(double,prob.l);
  prob.x = Malloc(struct svm_node *,prob.l);
  x_space = Malloc(struct svm_node,(nFeatures+1)*ntraining);
  unsigned long int spaceIndex=0;
  int lIndex=0;
  for(int iclass=0;iclass<nclass;++iclass){
    // for(int isample=0;isample<trainingFeatures[iclass].size();++isample){
    for(int isample=0;isample<nctraining[iclass];++isample){
      prob.x[lIndex]=&(x_space[spaceIndex]);
      for(int ifeature=0;ifeature<nFeatures;++ifeature){
        x_space[spaceIndex].index=ifeature+1;
        x_space[spaceIndex].value=trainingFeatures[iclass][isample][ifeature];
        ++spaceIndex;
      }
      x_space[spaceIndex++].index=-1;
      prob.y[lIndex]=iclass;
      ++lIndex;
    }
  }

  assert(lIndex==prob.l);
  if(verbose_opt[0]>2)
    std::cout << "checking parameters" << std::endl;
  svm_check_parameter(&prob,&param);
  if(verbose_opt[0]>2)
    std::cout << "parameters ok, training" << std::endl;
  svm=svm_train(&prob,&param);
  if(verbose_opt[0]>2)
    std::cout << "SVM is now trained" << std::endl;

  ConfusionMatrix cm;
  //set names in confusion matrix using nameVector
  for(int iname=0;iname<nameVector.size();++iname){
    if(classValueMap.empty())
      cm.pushBackClassName(nameVector[iname]);
    else if(cm.getClassIndex(type2string<short>(classValueMap[nameVector[iname]]))<0)
      cm.pushBackClassName(type2string<short>(classValueMap[nameVector[iname]]));
  }
  if(cv_opt[0]>1){
    double *target = Malloc(double,prob.l);
    svm_cross_validation(&prob,&param,cv_opt[0],target);
    assert(param.svm_type != EPSILON_SVR&&param.svm_type != NU_SVR);//only for regression
    for(int i=0;i<prob.l;i++){
      string refClassName=nameVector[prob.y[i]];
      string className=nameVector[target[i]];
      if(classValueMap.size())
	cm.incrementResult(type2string<short>(classValueMap[refClassName]),type2string<short>(classValueMap[className]),1.0);
      else
	cm.incrementResult(cm.getClass(prob.y[i]),cm.getClass(target[i]),1.0);
    }
    free(target);
  }
  else{
    struct svm_node *x_test;
    x_test = Malloc(struct svm_node,(nFeatures+1));
    for(int iclass=0;iclass<nclass;++iclass){
      for(int isample=0;isample<nctest[iclass];++isample){
	for(int ifeature=0;ifeature<nFeatures;++ifeature){
	  x_test[ifeature].index=ifeature+1;
	  x_test[ifeature].value=trainingFeatures[iclass][nctraining[iclass]+isample][ifeature];
	}
	x_test[nFeatures].index=-1;
	double predict_label=0;
	//todo: make distinction between svm_predict and svm_predict_probability?
	predict_label = svm_predict(svm,x_test);
	string refClassName=nameVector[iclass];
	string className=nameVector[static_cast<short>(predict_label)];
	if(classValueMap.size())
	  cm.incrementResult(type2string<short>(classValueMap[refClassName]),type2string<short>(classValueMap[className]),1.0);
	else
	  cm.incrementResult(refClassName,className,1.0);
      }
    }
    free(x_test);
  }
  if(verbose_opt[0]>1)
    std::cout << cm << std::endl;
  assert(cm.nReference());
  // if(verbose_opt[0])
  //   std::cout << cm << std::endl;
  // std::cout << "Kappa: " << cm.kappa() << std::endl;
  // double se95_oa=0;
  // double doa=0;
  // doa=cm.oa_pct(&se95_oa);
  // std::cout << "Overall Accuracy: " << doa << " (" << se95_oa << ")"  << std::endl;

  // *NOTE* Because svm_model contains pointers to svm_problem, you can
  // not free the memory used by svm_problem if you are still using the
  // svm_model produced by svm_train(). 
  // however, we will re-train the svm later on after the feature selection
  free(prob.y);
  free(prob.x);
  free(x_space);
  svm_free_and_destroy_model(&(svm));
  return(cm.kappa());
}

int main(int argc, char *argv[])
{
  // vector<double> priors;
  
  //--------------------------- command line options ------------------------------------
  Optionpk<string> input_opt("i", "input", "input test set (leave empty to perform a cross validation based on training only)"); 
  Optionpk<string> training_opt("t", "training", "training shape file. A single shape file contains all training features (must be set as: B0, B1, B2,...) for all classes (class numbers identified by label option)."); 
  Optionpk<string> label_opt("\0", "label", "identifier for class label in training shape file.","label"); 
  Optionpk<unsigned short> maxFeatures_opt("n", "nf", "number of features to select (0 to select optimal number, see also ecost option)", 0);
  Optionpk<unsigned int> balance_opt("\0", "balance", "balance the input data to this number of samples for each class", 0);
  Optionpk<int> minSize_opt("m", "min", "if number of training pixels is less then min, do not take this class into account", 0);
  Optionpk<double> start_opt("s", "start", "start band sequence number",0); 
  Optionpk<double> end_opt("e", "end", "end band sequence number (set to 0 to include all bands)", 0); 
  Optionpk<short> band_opt("b", "band", "band index (starting from 0, either use band option or use start to end)");
  Optionpk<double> offset_opt("\0", "offset", "offset value for each spectral band input features: refl[band]=(DN[band]-offset[band])/scale[band]", 0.0);
  Optionpk<double> scale_opt("\0", "scale", "scale value for each spectral band input features: refl=(DN[band]-offset[band])/scale[band] (use 0 if scale min and max in each band to -1.0 and 1.0)", 0.0);
  Optionpk<string> selector_opt("sm", "sm", "feature selection method (sffs=sequential floating forward search,sfs=sequential forward search, sbs, sequential backward search ,bfs=brute force search)","sffs"); 
  Optionpk<float> epsilon_cost_opt("ecost", "ecost", "epsilon for stopping criterion in cost function to determine optimal number of features",0.001);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    training_opt.retrieveOption(argc,argv);
    maxFeatures_opt.retrieveOption(argc,argv);
    label_opt.retrieveOption(argc,argv);
    balance_opt.retrieveOption(argc,argv);
    minSize_opt.retrieveOption(argc,argv);
    start_opt.retrieveOption(argc,argv);
    end_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    offset_opt.retrieveOption(argc,argv);
    scale_opt.retrieveOption(argc,argv);
    // priors_opt.retrieveOption(argc,argv);
    svm_type_opt.retrieveOption(argc,argv);
    kernel_type_opt.retrieveOption(argc,argv);
    kernel_degree_opt.retrieveOption(argc,argv);
    gamma_opt.retrieveOption(argc,argv);
    coef0_opt.retrieveOption(argc,argv);
    ccost_opt.retrieveOption(argc,argv);
    nu_opt.retrieveOption(argc,argv);
    epsilon_loss_opt.retrieveOption(argc,argv);
    cache_opt.retrieveOption(argc,argv);
    epsilon_tol_opt.retrieveOption(argc,argv);
    shrinking_opt.retrieveOption(argc,argv);
    prob_est_opt.retrieveOption(argc,argv);
    cv_opt.retrieveOption(argc,argv);
    selector_opt.retrieveOption(argc,argv);
    epsilon_cost_opt.retrieveOption(argc,argv);
    classname_opt.retrieveOption(argc,argv);
    classvalue_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }

  assert(training_opt.size());
  if(input_opt.size())
    cv_opt[0]=0;

  if(verbose_opt[0]>=1){
    if(input_opt.size())
      std::cout << "input filename: " << input_opt[0] << std::endl;
    std::cout << "training shape file: " << std::endl;
    for(int ifile=0;ifile<training_opt.size();++ifile)
      std::cout << training_opt[ifile] << std::endl;
    std::cout << "verbose: " << verbose_opt[0] << std::endl;
  }

  static std::map<std::string, SelectorValue> selMap;
  //initialize selMap
  selMap["sffs"]=SFFS;
  selMap["sfs"]=SFS;
  selMap["sbs"]=SBS;
  selMap["bfs"]=BFS;

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
      totalSamples=trainingReader.readDataImageShape(trainingMap,fields,band_opt,label_opt[0],verbose_opt[0]);
      if(input_opt.size()){
	ImgReaderOgr inputReader(input_opt[0]);
	totalTestSamples=inputReader.readDataImageShape(testMap,fields,band_opt,label_opt[0],verbose_opt[0]);
	inputReader.close();
      }
    }
    else{
      totalSamples=trainingReader.readDataImageShape(trainingMap,fields,start_opt[0],end_opt[0],label_opt[0],verbose_opt[0]);
      if(input_opt.size()){
	ImgReaderOgr inputReader(input_opt[0]);
	totalTestSamples=inputReader.readDataImageShape(testMap,fields,start_opt[0],end_opt[0],label_opt[0],verbose_opt[0]);
	inputReader.close();
      }
    }
    if(trainingMap.size()<2){
      string errorstring="Error: could not read at least two classes from training input file";
      throw(errorstring);
    }
    if(input_opt.size()&&testMap.size()<2){
      string errorstring="Error: could not read at least two classes from test input file";
      throw(errorstring);
    }
    trainingReader.close();
  }
  catch(string error){
    cerr << error << std::endl;
    exit(1);
  }
  catch(...){
    cerr << "error catched" << std::endl;
    exit(1);
  }
  //todo: delete class 0 ?
  // if(verbose_opt[0]>=1)
  //   std::cout << "erasing class 0 from training set (" << trainingMap[0].size() << " from " << totalSamples << ") samples" << std::endl;
  // totalSamples-=trainingMap[0].size();
  // trainingMap.erase(0);

  if(verbose_opt[0]>1)
    std::cout << "training pixels: " << std::endl;
  map<string,Vector2d<float> >::iterator mapit=trainingMap.begin();
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
    if(random)
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
    
  int nFeatures=trainingFeatures[0][0].size();
  int maxFeatures=(maxFeatures_opt[0])? maxFeatures_opt[0] : 1;
  double previousCost=-1;
  double cost=0;
  list<int> subset;//set of selected features (levels) for each class combination
  FeatureSelector selector;
  try{
    if(maxFeatures==nFeatures){
      subset.clear();
      for(int ifeature=0;ifeature<nFeatures;++ifeature)
        subset.push_back(ifeature);
      cost=getCost(trainingFeatures);
    }
    else{
      while(fabs(cost-previousCost)>epsilon_cost_opt[0]){
        previousCost=cost;
        switch(selMap[selector_opt[0]]){
        case(SFFS):
          subset.clear();//needed to clear in case of floating and brute force search
          cost=selector.floating(trainingFeatures,&getCost,subset,maxFeatures,verbose_opt[0]);
          break;
        case(SFS):
          cost=selector.forward(trainingFeatures,&getCost,subset,maxFeatures,verbose_opt[0]);
          break;
        case(SBS):
          cost=selector.backward(trainingFeatures,&getCost,subset,maxFeatures,verbose_opt[0]);
          break;
        case(BFS):
          subset.clear();//needed to clear in case of floating and brute force search
          cost=selector.bruteForce(trainingFeatures,&getCost,subset,maxFeatures,verbose_opt[0]);
          break;
        default:
          std::cout << "Error: selector not supported, please use sffs, sfs, sbs or bfs" << std::endl;
          exit(1);
          break;
        }
        if(verbose_opt[0]){
          std::cout << "cost: " << cost << std::endl;
          std::cout << "previousCost: " << previousCost << std::endl;
          std::cout << std::setprecision(12) << "cost-previousCost: " << cost - previousCost << " ( " << epsilon_cost_opt[0] << ")" << std::endl;
        }
        if(!maxFeatures_opt[0])
          ++maxFeatures;
        else
          break;
      }
    }
  }
  catch(...){
    std::cout << "catched feature selection" << std::endl;
    exit(1);
  }

  if(verbose_opt[0])
    cout <<"cost: " << cost << endl;
  for(list<int>::const_iterator lit=subset.begin();lit!=subset.end();++lit)
    std::cout << " -b " << *lit;
  std::cout << std::endl;
    // if((*(lit))!=subset.back())
    // else
    //   cout << endl;

  // *NOTE* Because svm_model contains pointers to svm_problem, you can
  // not free the memory used by svm_problem if you are still using the
  // svm_model produced by svm_train(). 

  // free(prob.y);
  // free(prob.x);
  // free(x_space);
  // svm_destroy_param(&param);
  return 0;
}
                            
