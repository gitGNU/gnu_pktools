/**********************************************************************
pkfs_nn.cc: feature selection for nn classifier
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
#include "floatfann.h"
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"
#include "imageclasses/ImgReaderOgr.h"
#include "imageclasses/ImgWriterOgr.h"
#include "base/Optionpk.h"
#include "algorithms/myfann_cpp.h"
#include "algorithms/ConfusionMatrix.h"
#include "algorithms/FeatureSelector.h"
#include "pkclassify_nn.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

                                    //static enum SelectorValue { NA, SFFS, SFS, SBS, BFS };
enum SelectorValue  { NA=0, SFFS=1, SFS=2, SBS=3, BFS=4 };

//global parameters used in cost function getCost
Optionpk<unsigned int> nneuron_opt("\0", "nneuron", "number of neurons in hidden layers in neural network (multiple hidden layers are set by defining multiple number of neurons: -n 15 -n 1, default is one hidden layer with 5 neurons)", 5); 
Optionpk<float> connection_opt("\0", "connection", "connection reate (default: 1.0 for a fully connected network)", 1.0); 
Optionpk<float> weights_opt("w", "weights", "weights for neural network. Apply to fully connected network only, starting from first input neuron to last output neuron, including the bias neurons (last neuron in each but last layer)", 0.0); 
Optionpk<float> learning_opt("l", "learning", "learning rate (default: 0.7)", 0.7); 
Optionpk<unsigned int> maxit_opt("\0", "maxit", "number of maximum iterations (epoch) (default: 500)", 500); 
// Optionpk<bool> weight_opt("wi", "wi", "set the parameter C of class i to weight*C, for C-SVC",true);
Optionpk<unsigned short> cv_opt("cv", "cv", "n-fold cross validation mode",2);
Optionpk<short> verbose_opt("v", "verbose", "set to: 0 (results only), 1 (confusion matrix), 2 (debug)",0);

double getCost(const vector<Vector2d<float> > &trainingFeatures)
{
  vector<Vector2d<float> > tmpFeatures=trainingFeatures;//deep copy is guaranteed through constructor?
  unsigned short nclass=tmpFeatures.size();
  unsigned int ntraining=0;
  for(int iclass=0;iclass<nclass;++iclass)
    ntraining+=tmpFeatures[iclass].size();
  unsigned short nFeatures=tmpFeatures[0][0].size();

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
  case(3):
    net.create_sparse(connection_opt[0],num_layers, nFeatures, nneuron_opt[0], nclass);
    break;
  case(4):
    net.create_sparse(connection_opt[0],num_layers, nFeatures, nneuron_opt[0], nneuron_opt[1], nclass);
    break;
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
  float rmse=net.cross_validation(tmpFeatures,
                                  ntraining,
                                  cv_opt[0],
                                  maxit_opt[0],
                                  0,
                                  desired_error,
                                  referenceVector,
                                  outputVector,
                                  verbose_opt[0]);
  ConfusionMatrix cm(nclass);
  for(int isample=0;isample<referenceVector.size();++isample)
    cm.incrementResult(cm.getClass(referenceVector[isample]),cm.getClass(outputVector[isample]),1);
  assert(cm.nReference());
  // std::cout << cm << std::endl;
  // std::cout << "Kappa: " << cm.kappa() << std::endl;
  // double se95_oa=0;
  // double doa=0;
  // doa=cm.oa_pct(&se95_oa);
  // std::cout << "Overall Accuracy: " << doa << " (" << se95_oa << ")"  << std::endl;
  // std::cout << "rmse cross-validation: " << rmse << std::endl;
  return(cm.kappa());
}

int main(int argc, char *argv[])
{
  map<short,int> reclassMap;
  vector<int> vreclass;
  // vector<double> priors;
  
  //--------------------------- command line options ------------------------------------
  Optionpk<string> training_opt("t", "training", "training shape file. A single shape file contains all training features (must be set as: B0, B1, B2,...) for all classes (class numbers identified by label option). Use multiple training files for bootstrap aggregation (alternative to the bag and bsize options, where a random subset is taken from a single training file)"); 
  Optionpk<string> label_opt("\0", "label", "identifier for class label in training shape file.","label"); 
  Optionpk<unsigned short> maxFeatures_opt("n", "nf", "number of features to select (0 to select optimal number, see also ecost option)", 0);
  Optionpk<unsigned short> reclass_opt("\0", "rc", "reclass code (e.g. --rc=12 --rc=23 to reclass first two classes to 12 and 23 resp.).", 0);
  Optionpk<unsigned int> balance_opt("\0", "balance", "balance the input data to this number of samples for each class", 0);
  Optionpk<int> minSize_opt("m", "min", "if number of training pixels is less then min, do not take this class into account", 0);
  Optionpk<double> start_opt("s", "start", "start band sequence number (set to 0)",0); 
  Optionpk<double> end_opt("e", "end", "end band sequence number (set to 0 for all bands)", 0); 
  Optionpk<short> band_opt("b", "band", "band index (starting from 0, either use band option or use start to end)");
  Optionpk<double> offset_opt("\0", "offset", "offset value for each spectral band input features: refl[band]=(DN[band]-offset[band])/scale[band]", 0.0);
  Optionpk<double> scale_opt("\0", "scale", "scale value for each spectral band input features: refl=(DN[band]-offset[band])/scale[band] (use 0 if scale min and max in each band to -1.0 and 1.0)", 0.0);
  Optionpk<unsigned short> aggreg_opt("a", "aggreg", "how to combine aggregated classifiers, see also rc option (0: no aggregation, 1: sum rule, 2: max rule).",0);
  // Optionpk<double> priors_opt("p", "prior", "prior probabilities for each class (e.g., -p 0.3 -p 0.3 -p 0.2 )", 0.0);
  Optionpk<string> selector_opt("sm", "sm", "feature selection method (sffs=sequential floating forward search,sfs=sequential forward search, sbs, sequential backward search ,bfs=brute force search)","sffs"); 
  Optionpk<float> epsilon_cost_opt("ecost", "ecost", "epsilon for stopping criterion in cost function to determine optimal number of features",0.001);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=training_opt.retrieveOption(argc,argv);
    maxFeatures_opt.retrieveOption(argc,argv);
    label_opt.retrieveOption(argc,argv);
    reclass_opt.retrieveOption(argc,argv);
    balance_opt.retrieveOption(argc,argv);
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

  assert(training_opt[0].size());
  if(verbose_opt[0]>=1)
    std::cout << "training shape file: " << training_opt[0] << std::endl;

  unsigned int totalSamples=0;
  int nreclass=0;
  vector<int> vcode;//unique class codes in recode string

  unsigned short nclass=0;
  int nband=0;
  int startBand=2;//first two bands represent X and Y pos

  vector<double> offset;
  vector<double> scale;
  vector< Vector2d<float> > trainingPixels;//[class][sample][band]

  if(reclass_opt.size()>1){
    vreclass.resize(reclass_opt.size());
    for(int iclass=0;iclass<reclass_opt.size();++iclass){
      reclassMap[iclass]=reclass_opt[iclass];
      vreclass[iclass]=reclass_opt[iclass];
    }
  }
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
  //----------------------------------- Training -------------------------------
  vector<string> fields;
  //organize training data
  trainingPixels.clear();
  map<int,Vector2d<float> > trainingMap;
  if(verbose_opt[0]>=1)
    std::cout << "reading imageShape file " << training_opt[0] << std::endl;
  try{
    if(band_opt.size())
      totalSamples=readDataImageShape(training_opt[0],trainingMap,fields,band_opt,label_opt[0],verbose_opt[0]);
    else
      totalSamples=readDataImageShape(training_opt[0],trainingMap,fields,start_opt[0],end_opt[0],label_opt[0],verbose_opt[0]);
    if(trainingMap.size()<2){
      string errorstring="Error: could not read at least two classes from training file";
      throw(errorstring);
    }
  }
  catch(string error){
    cerr << error << std::endl;
    exit(1);
  }
  catch(...){
    cerr << "error catched" << std::endl;
    exit(1);
  }
  //delete class 0
  if(verbose_opt[0]>=1)
    std::cout << "erasing class 0 from training set (" << trainingMap[0].size() << " from " << totalSamples << ") samples" << std::endl;
  totalSamples-=trainingMap[0].size();
  trainingMap.erase(0);
  //convert map to vector
  short iclass=0;
  if(reclass_opt.size()==1){//no reclass option, read classes from shape
    reclassMap.clear();
    vreclass.clear();
  }
  if(verbose_opt[0]>1)
    std::cout << "training pixels: " << std::endl;
  map<int,Vector2d<float> >::iterator mapit=trainingMap.begin();
  while(mapit!=trainingMap.end()){
    //       for(map<int,Vector2d<float> >::const_iterator mapit=trainingMap.begin();mapit!=trainingMap.end();++mapit){
    //delete small classes
    if((mapit->second).size()<minSize_opt[0]){
      trainingMap.erase(mapit);
      continue;
      //todo: beware of reclass option: delete this reclass if no samples are left in this classes!!
    }
    if(reclass_opt.size()==1){//no reclass option, read classes from shape
      reclassMap[iclass]=(mapit->first);
      vreclass.push_back(mapit->first);
    }
    trainingPixels.push_back(mapit->second);
    if(verbose_opt[0]>1)
      std::cout << mapit->first << ": " << (mapit->second).size() << " samples" << std::endl;
    ++iclass;
    ++mapit;
  }
  nclass=trainingPixels.size();
  nband=trainingPixels[0][0].size()-2;//X and Y//trainingPixels[0][0].size();
  assert(reclassMap.size()==nclass);

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

  //recode vreclass to ordered vector, starting from 0 to nreclass
  vcode.clear();
  if(verbose_opt[0]>=1){
    std::cout << "before recoding: " << std::endl;
    for(int iclass = 0; iclass < vreclass.size(); iclass++)
      std::cout << " " << vreclass[iclass];
    std::cout << std::endl; 
  }
  vector<int> vord=vreclass;//ordered vector, starting from 0 to nreclass
  map<short,int> mreclass;
  for(int ic=0;ic<vreclass.size();++ic){
    if(mreclass.find(vreclass[ic])==mreclass.end())
      mreclass[vreclass[ic]]=iclass++;
  }
  for(int ic=0;ic<vreclass.size();++ic)
    vord[ic]=mreclass[vreclass[ic]];
  //construct uniqe class codes
  while(!vreclass.empty()){
    vcode.push_back(*(vreclass.begin()));
    //delete all these entries from vreclass
    vector<int>::iterator vit;
    while((vit=find(vreclass.begin(),vreclass.end(),vcode.back()))!=vreclass.end())
      vreclass.erase(vit);
  }
  if(verbose_opt[0]>=1){
    std::cout << "recode values: " << std::endl;
    for(int icode=0;icode<vcode.size();++icode)
      std::cout << vcode[icode] << " ";
    std::cout << std::endl;
  }
  vreclass=vord;
  if(verbose_opt[0]>=1){
    std::cout << "after recoding: " << std::endl;
    for(int iclass = 0; iclass < vord.size(); iclass++)
      std::cout << " " << vord[iclass];
    std::cout << std::endl; 
  }
      
  vector<int> vuniqueclass=vreclass;
  //remove duplicate elements from vuniqueclass
  sort( vuniqueclass.begin(), vuniqueclass.end() );
  vuniqueclass.erase( unique( vuniqueclass.begin(), vuniqueclass.end() ), vuniqueclass.end() );
  nreclass=vuniqueclass.size();
  if(verbose_opt[0]>=1){
    std::cout << "unique classes: " << std::endl;
    for(int iclass = 0; iclass < vuniqueclass.size(); iclass++)
      std::cout << " " << vuniqueclass[iclass];
    std::cout << std::endl; 
    std::cout << "number of reclasses: " << nreclass << std::endl;
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

  //Calculate features of trainig set
  vector< Vector2d<float> > trainingFeatures(nclass);
  for(int iclass=0;iclass<nclass;++iclass){
    int nctraining=0;
    if(verbose_opt[0]>=1)
      std::cout << "calculating features for class " << iclass << std::endl;
    if(random)
      srand(time(NULL));
    nctraining=trainingPixels[iclass].size();//bagSize_opt[0] given in % of training size
    if(verbose_opt[0]>=1)
      std::cout << "nctraining class " << iclass << ": " << nctraining << std::endl;
    int index=0;
      
    trainingFeatures[iclass].resize(nctraining);
    for(int isample=0;isample<nctraining;++isample){
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
    assert(trainingFeatures[iclass].size()==nctraining);
  }
    
  unsigned int ntraining=0;
  for(int iclass=0;iclass<nclass;++iclass){
    if(verbose_opt[0]>1)
      std::cout << "training sample size for class " << vcode[iclass] << ": " << trainingFeatures[iclass].size() << std::endl;
    ntraining+=trainingFeatures[iclass].size();
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
      while(cost-previousCost>epsilon_cost_opt[0]){
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
                            
