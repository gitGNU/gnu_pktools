/**********************************************************************
pkfsann.cc: feature selection for nn classifier
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
#include <stdlib.h>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include "base/Optionpk.h"
#include "imageclasses/ImgReaderOgr.h"
#include "algorithms/ConfusionMatrix.h"
#include "algorithms/FeatureSelector.h"
#include "floatfann.h"
#include "algorithms/myfann_cpp.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

enum SelectorValue  { NA=0, SFFS=1, SFS=2, SBS=3, BFS=4 };

using namespace std;

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))


//global parameters used in cost function getCost
ConfusionMatrix cm;
map<string,short> classValueMap;
vector<std::string> nameVector;
vector<unsigned int> nctraining;
vector<unsigned int> nctest;
Optionpk<unsigned int> nneuron_opt("\0", "nneuron", "number of neurons in hidden layers in neural network (multiple hidden layers are set by defining multiple number of neurons: -n 15 -n 1, default is one hidden layer with 5 neurons)", 5); 
Optionpk<float> connection_opt("\0", "connection", "connection reate (default: 1.0 for a fully connected network)", 1.0); 
Optionpk<float> weights_opt("w", "weights", "weights for neural network. Apply to fully connected network only, starting from first input neuron to last output neuron, including the bias neurons (last neuron in each but last layer)", 0.0); 
Optionpk<float> learning_opt("l", "learning", "learning rate (default: 0.7)", 0.7); 
Optionpk<unsigned int> maxit_opt("\0", "maxit", "number of maximum iterations (epoch) (default: 500)", 500); 
// Optionpk<bool> weight_opt("wi", "wi", "set the parameter C of class i to weight*C, for C-SVC",true);
Optionpk<unsigned short> cv_opt("cv", "cv", "n-fold cross validation mode",2);
Optionpk<string> classname_opt("c", "class", "list of class names."); 
Optionpk<short> classvalue_opt("r", "reclass", "list of class values (use same order as in classname opt."); 
Optionpk<short> verbose_opt("v", "verbose", "set to: 0 (results only), 1 (confusion matrix), 2 (debug)",0);

double getCost(const vector<Vector2d<float> > &trainingFeatures)
{
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

  FANN::neural_net net;//the neural network
  const unsigned int num_layers = nneuron_opt.size()+2;
  const float desired_error = 0.0003;
  const unsigned int iterations_between_reports = (verbose_opt[0])?maxit_opt[0]+1:0;
  if(verbose_opt[0]>1){
    cout << "creating artificial neural network with " << nneuron_opt.size() << " hidden layer, having " << endl;
    for(int ilayer=0;ilayer<nneuron_opt.size();++ilayer)
      cout << nneuron_opt[ilayer] << " ";
    cout << "neurons" << endl;
  }
  switch(num_layers){
  case(3):{
    unsigned int layers[3];
    layers[0]=nFeatures;
    layers[1]=nneuron_opt[0];
    layers[2]=nclass;
    net.create_sparse_array(connection_opt[0],num_layers,layers);
    // net.create_sparse(connection_opt[0],num_layers, nFeatures, nneuron_opt[0], nclass);
    break;
  }
  case(4):{
    unsigned int layers[4];
    layers[0]=nFeatures;
    layers[1]=nneuron_opt[0];
    layers[2]=nneuron_opt[1];
    layers[3]=nclass;
    net.create_sparse_array(connection_opt[0],num_layers,layers);
    // net.create_sparse(connection_opt[0],num_layers, nFeatures, nneuron_opt[0], nneuron_opt[1], nclass);
    break;
  }
  default:
    cerr << "Only 1 or 2 hidden layers are supported!" << endl;
    exit(1);
    break;
  }

  net.set_learning_rate(learning_opt[0]);
  //   net.set_activation_steepness_hidden(1.0);
  //   net.set_activation_steepness_output(1.0);
    
  net.set_activation_function_hidden(FANN::SIGMOID_SYMMETRIC_STEPWISE);
  net.set_activation_function_output(FANN::SIGMOID_SYMMETRIC_STEPWISE);

  // Set additional properties such as the training algorithm
  //   net.set_training_algorithm(FANN::TRAIN_QUICKPROP);

  vector<unsigned short> referenceVector;
  vector<unsigned short> outputVector;
  float rmse=0;
  // ConfusionMatrix cm;
  // //set names in confusion matrix using nameVector
  // for(int iname=0;iname<nameVector.size();++iname){
  //   if(classValueMap.empty())
  //     cm.pushBackClassName(nameVector[iname]);
  //   else if(cm.getClassIndex(type2string<short>(classValueMap[nameVector[iname]]))<0)
  //     cm.pushBackClassName(type2string<short>(classValueMap[nameVector[iname]]));
  // }
  vector<Vector2d<float> > tmpFeatures(nclass);
  for(int iclass=0;iclass<nclass;++iclass){
    tmpFeatures[iclass].resize(trainingFeatures[iclass].size(),nFeatures);
    for(unsigned int isample=0;isample<nctraining[iclass];++isample){
	for(int ifeature=0;ifeature<nFeatures;++ifeature){
          tmpFeatures[iclass][isample][ifeature]=trainingFeatures[iclass][isample][ifeature];
      }
    }
  }
  if(cv_opt[0]>0){
    rmse=net.cross_validation(tmpFeatures,
                              ntraining,
                              cv_opt[0],
                              maxit_opt[0],
                              desired_error,
                              referenceVector,
                              outputVector,
                              verbose_opt[0]);
    for(int isample=0;isample<referenceVector.size();++isample){
      string refClassName=nameVector[referenceVector[isample]];
      string className=nameVector[outputVector[isample]];
      if(classValueMap.size())
	cm.incrementResult(type2string<short>(classValueMap[refClassName]),type2string<short>(classValueMap[className]),1.0);
      else
	cm.incrementResult(cm.getClass(referenceVector[isample]),cm.getClass(outputVector[isample]),1.0);
    }
  }
  else{//not working yet. please repair...
    assert(cv_opt[0]>0);
    bool initWeights=true;
    net.train_on_data(tmpFeatures,ntraining,initWeights, maxit_opt[0],
                      iterations_between_reports, desired_error);
    vector<Vector2d<float> > testFeatures(nclass);
    vector<float> result(nclass);
    int maxClass=-1;
    for(int iclass=0;iclass<nclass;++iclass){
      testFeatures.resize(nctest[iclass],nFeatures);
      for(unsigned int isample=0;isample<nctraining[iclass];++isample){
	for(int ifeature=0;ifeature<nFeatures;++ifeature){
          testFeatures[iclass][isample][ifeature]=trainingFeatures[iclass][nctraining[iclass]+isample][ifeature];
        }
        result=net.run(testFeatures[iclass][isample]);
        string refClassName=nameVector[iclass];
        float maxP=-1;
        for(int ic=0;ic<nclass;++ic){
          float pv=(result[ic]+1.0)/2.0;//bring back to scale [0,1]
          if(pv>maxP){
            maxP=pv;
            maxClass=ic;
          }
        }
        string className=nameVector[maxClass];
        if(classValueMap.size())
          cm.incrementResult(type2string<short>(classValueMap[refClassName]),type2string<short>(classValueMap[className]),1.0);
        else
          cm.incrementResult(cm.getClass(referenceVector[isample]),cm.getClass(outputVector[isample]),1.0);
      }
    }
  }
  assert(cm.nReference());
  return(cm.kappa());
}

int main(int argc, char *argv[])
{
  // vector<double> priors;
  
  //--------------------------- command line options ------------------------------------
  Optionpk<string> input_opt("i", "input", "input test set (leave empty to perform a cross validation based on training only)"); 
  Optionpk<string> training_opt("t", "training", "training vector file. A single vector file contains all training features (must be set as: B0, B1, B2,...) for all classes (class numbers identified by label option). Use multiple training files for bootstrap aggregation (alternative to the bag and bsize options, where a random subset is taken from a single training file)"); 
  Optionpk<string> tlayer_opt("tln", "tln", "training layer name(s)");
  Optionpk<string> label_opt("\0", "label", "identifier for class label in training vector file.","label"); 
  Optionpk<unsigned short> maxFeatures_opt("n", "nf", "number of features to select (0 to select optimal number, see also ecost option)", 0);
  Optionpk<unsigned int> balance_opt("\0", "balance", "balance the input data to this number of samples for each class", 0);
  Optionpk<bool> random_opt("random","random", "in case of balance, randomize input data", true);
  Optionpk<int> minSize_opt("m", "min", "if number of training pixels is less then min, do not take this class into account", 0);
  Optionpk<double> start_opt("s", "start", "start band sequence number",0); 
  Optionpk<double> end_opt("e", "end", "end band sequence number (set to 0 to include all bands)", 0); 
  Optionpk<short> band_opt("b", "band", "band index (starting from 0, either use band option or use start to end)");
  Optionpk<double> offset_opt("\0", "offset", "offset value for each spectral band input features: refl[band]=(DN[band]-offset[band])/scale[band]", 0.0);
  Optionpk<double> scale_opt("\0", "scale", "scale value for each spectral band input features: refl=(DN[band]-offset[band])/scale[band] (use 0 if scale min and max in each band to -1.0 and 1.0)", 0.0);
  Optionpk<unsigned short> aggreg_opt("a", "aggreg", "how to combine aggregated classifiers, see also rc option (0: no aggregation, 1: sum rule, 2: max rule).",0);
  // Optionpk<double> priors_opt("p", "prior", "prior probabilities for each class (e.g., -p 0.3 -p 0.3 -p 0.2 )", 0.0);
  Optionpk<string> selector_opt("sm", "sm", "feature selection method (sffs=sequential floating forward search,sfs=sequential forward search, sbs, sequential backward search ,bfs=brute force search)","sffs"); 
  Optionpk<float> epsilon_cost_opt("ecost", "ecost", "epsilon for stopping criterion in cost function to determine optimal number of features",0.001);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    training_opt.retrieveOption(argc,argv);
    maxFeatures_opt.retrieveOption(argc,argv);
    tlayer_opt.retrieveOption(argc,argv);
    label_opt.retrieveOption(argc,argv);
    balance_opt.retrieveOption(argc,argv);
    random_opt.retrieveOption(argc,argv);
    minSize_opt.retrieveOption(argc,argv);
    start_opt.retrieveOption(argc,argv);
    end_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    offset_opt.retrieveOption(argc,argv);
    scale_opt.retrieveOption(argc,argv);
    aggreg_opt.retrieveOption(argc,argv);
    // priors_opt.retrieveOption(argc,argv);
    nneuron_opt.retrieveOption(argc,argv);
    connection_opt.retrieveOption(argc,argv);
    weights_opt.retrieveOption(argc,argv);
    learning_opt.retrieveOption(argc,argv);
    maxit_opt.retrieveOption(argc,argv);
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
  
  static std::map<std::string, SelectorValue> selMap;
  //initialize selMap
  selMap["sffs"]=SFFS;
  selMap["sfs"]=SFS;
  selMap["sbs"]=SBS;
  selMap["bfs"]=BFS;

  assert(training_opt.size());
  if(input_opt.size())
    cv_opt[0]=0;
  if(verbose_opt[0]>=1)
    std::cout << "training vector file: " << training_opt[0] << std::endl;

  unsigned int totalSamples=0;
  unsigned int totalTestSamples=0;

  unsigned short nclass=0;
  int nband=0;
  int startBand=2;//first two bands represent X and Y pos

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
  vector<double> offset;
  vector<double> scale;
  vector< Vector2d<float> > trainingPixels;//[class][sample][band]
  vector< Vector2d<float> > testPixels;//[class][sample][band]
  map<string,Vector2d<float> > trainingMap;
  map<string,Vector2d<float> > testMap;
  vector<string> fields;

  //organize training data
  trainingPixels.clear();
  if(verbose_opt[0]>=1)
    std::cout << "reading imageVector file " << training_opt[0] << std::endl;
  try{
    ImgReaderOgr trainingReader(training_opt[0]);
    if(band_opt.size()){
      totalSamples=trainingReader.readDataImageOgr(trainingMap,fields,band_opt,label_opt[0],tlayer_opt,verbose_opt[0]);
      if(input_opt.size()){
	ImgReaderOgr inputReader(input_opt[0]);
	totalTestSamples=trainingReader.readDataImageOgr(testMap,fields,band_opt,label_opt[0],tlayer_opt,verbose_opt[0]);
	inputReader.close();
      }
    }
    else{
      totalSamples=trainingReader.readDataImageOgr(trainingMap,fields,start_opt[0],end_opt[0],label_opt[0],tlayer_opt,verbose_opt[0]);
      if(input_opt.size()){
	ImgReaderOgr inputReader(input_opt[0]);
	totalTestSamples=trainingReader.readDataImageOgr(testMap,fields,start_opt[0],end_opt[0],label_opt[0],tlayer_opt,verbose_opt[0]);
	inputReader.close();
      }
    }
    if(trainingMap.size()<2){
      string errorstring="Error: could not read at least two classes from training file";
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
  catch(std::exception& e){
    std::cerr << "Error: ";
    std::cerr << e.what() << std::endl;
    std::cerr << CPLGetLastErrorMsg() << std::endl; 
    exit(1);
  }
  catch(...){
    cerr << "error catched" << std::endl;
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
    // if(classValueMap.size()){
    //   //check if name in training is covered by classname_opt (values can not be 0)
    //   if(classValueMap[mapit->first]>0){
    // 	if(verbose_opt[0])
    // 	  std::cout << mapit->first << " -> " << classValueMap[mapit->first] << std::endl;
    //   }
    //   else{
    // 	std::cerr << "Error: names in classname option are not complete, please check names in training vector and make sure classvalue is > 0" << std::endl;
    // 	exit(1);
    //   }
    // }    
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
  nclass=trainingPixels.size();
  if(classname_opt.size())
    assert(nclass==classname_opt.size());
  nband=trainingPixels[0][0].size()-2;//X and Y//trainingPixels[0][0].size();

  mapit=testMap.begin();
  while(mapit!=testMap.end()){
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

  //new
  mapit=trainingMap.begin();
  bool doSort=true;
  try{
    while(mapit!=trainingMap.end()){
      nameVector.push_back(mapit->first);
      if(classValueMap.size()){
	//check if name in training is covered by classname_opt (values can not be 0)
	if(classValueMap[mapit->first]>0){
	  if(cm.getClassIndex(type2string<short>(classValueMap[mapit->first]))<0){
	    cm.pushBackClassName(type2string<short>(classValueMap[mapit->first]),doSort);
	  }
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
  }
  catch(BadConversion conversionString){
    std::cerr << "Error: did you provide class pairs names (-c) and integer values (-r) for each class in training vector?" << std::endl;
    exit(1);
  }
  if(classname_opt.empty()){
    //std::cerr << "Warning: no class name and value pair provided for all " << nclass << " classes, using string2type<int> instead!" << std::endl;
    for(int iclass=0;iclass<nclass;++iclass){
      if(verbose_opt[0])
	std::cout << iclass << " " << cm.getClass(iclass) << " -> " << string2type<short>(cm.getClass(iclass)) << std::endl;
      classValueMap[cm.getClass(iclass)]=string2type<short>(cm.getClass(iclass));
    }
  }

  //Calculate features of trainig set
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
    if(maxFeatures>=nFeatures){
      subset.clear();
      for(int ifeature=0;ifeature<nFeatures;++ifeature)
        subset.push_back(ifeature);
      cost=getCost(trainingFeatures);
    }
    else{
      while(fabs(cost-previousCost)>=epsilon_cost_opt[0]){
        previousCost=cost;
        switch(selMap[selector_opt[0]]){
        case(SFFS):
          subset.clear();//needed to clear in case of floating and brute force search
          cost=selector.floating(trainingFeatures,&getCost,subset,maxFeatures,epsilon_cost_opt[0],verbose_opt[0]);
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
        if(verbose_opt[0]>1){
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
                            
