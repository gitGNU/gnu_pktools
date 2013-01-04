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
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"
#include "imageclasses/ImgReaderOgr.h"
#include "imageclasses/ImgWriterOgr.h"
#include "base/Optionpk.h"
#include "algorithms/ConfusionMatrix.h"
#include "algorithms/FeatureSelector.h"
#include "algorithms/svm.h"
#include "pkclassify_nn.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

                                    //static enum SelectorValue { NA, SFFS, SFS, SBS, BFS };
enum SelectorValue  { NA=0, SFFS=1, SFS=2, SBS=3, BFS=4 };

//global parameters used in cost function getCost
Optionpk<unsigned short> svm_type_opt("svmt", "svmtype", "type of SVM (0: C-SVC, 1: nu-SVC, 2: one-class SVM, 3: epsilon-SVR,	4: nu-SVR)",0);
Optionpk<unsigned short> kernel_type_opt("kt", "kerneltype", "type of kernel function (0: linear: u'*v, 1: polynomial: (gamma*u'*v + coef0)^degree, 2: radial basis function: exp(-gamma*|u-v|^2), 3: sigmoid: tanh(gamma*u'*v + coef0), 4: precomputed kernel (kernel values in training_set_file)",2);
Optionpk<unsigned short> kernel_degree_opt("kd", "kd", "degree in kernel function",3);
Optionpk<float> gamma_opt("g", "gamma", "gamma in kernel function",0);
Optionpk<float> coef0_opt("c0", "coef0", "coef0 in kernel function",0);
Optionpk<float> ccost_opt("cc", "ccost", "the parameter C of C-SVC, epsilon-SVR, and nu-SVR",1);
Optionpk<float> nu_opt("nu", "nu", "the parameter nu of nu-SVC, one-class SVM, and nu-SVR",0.5);
Optionpk<float> epsilon_loss_opt("eloss", "eloss", "the epsilon in loss function of epsilon-SVR",0.1);
Optionpk<int> cache_opt("cache", "cache", "cache memory size in MB",100);
Optionpk<float> epsilon_tol_opt("etol", "etol", "the tolerance of termination criterion",0.001);
Optionpk<bool> shrinking_opt("shrink", "shrink", "whether to use the shrinking heuristics",false);
Optionpk<bool> prob_est_opt("pe", "probest", "whether to train a SVC or SVR model for probability estimates",false);
// Optionpk<bool> weight_opt("wi", "wi", "set the parameter C of class i to weight*C, for C-SVC",true);
Optionpk<unsigned short> cv_opt("cv", "cv", "n-fold cross validation mode",2);
Optionpk<short> verbose_opt("v", "verbose", "set to: 0 (results only), 1 (confusion matrix), 2 (debug)",0);

double getCost(const vector<Vector2d<float> > &trainingFeatures)
{
  unsigned short nclass=trainingFeatures.size();
  unsigned int ntraining=0;
  for(int iclass=0;iclass<nclass;++iclass)
    ntraining+=trainingFeatures[iclass].size();
  unsigned short nFeatures=trainingFeatures[0][0].size();

  struct svm_parameter param;
  param.svm_type = svm_type_opt[0];
  param.kernel_type = kernel_type_opt[0];
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
    for(int isample=0;isample<trainingFeatures[iclass].size();++isample){
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
  if(verbose_opt[0]>1)
    std::cout << "checking parameters" << std::endl;
  svm_check_parameter(&prob,&param);
  if(verbose_opt[0]>1)
    std::cout << "parameters ok, training" << std::endl;
  svm=svm_train(&prob,&param);
  if(verbose_opt[0]>1)
    std::cout << "SVM is now trained" << std::endl;

  ConfusionMatrix cm(nclass);
  double *target = Malloc(double,prob.l);
  svm_cross_validation(&prob,&param,cv_opt[0],target);
  assert(param.svm_type != EPSILON_SVR&&param.svm_type != NU_SVR);//only for regression
  int total_correct=0;
  for(int i=0;i<prob.l;i++)
    cm.incrementResult(cm.getClass(prob.y[i]),cm.getClass(target[i]),1);
  assert(cm.nReference());
  // std::cout << "Kappa: " << cm.kappa() << std::endl;
  // double se95_oa=0;
  // double doa=0;
  // doa=cm.oa_pct(&se95_oa);
  // std::cout << "Overall Accuracy: " << doa << " (" << se95_oa << ")"  << std::endl;

  // *NOTE* Because svm_model contains pointers to svm_problem, you can
  // not free the memory used by svm_problem if you are still using the
  // svm_model produced by svm_train(). 
  // however, we will re-train the svm later on after the feature selection
  free(target);
  free(prob.y);
  free(prob.x);
  free(x_space);
  svm_free_and_destroy_model(&(svm));
  return(cm.kappa());
}

int main(int argc, char *argv[])
{
  map<short,int> reclassMap;
  vector<int> vreclass;
  // vector<double> priors;
  
  //--------------------------- command line options ------------------------------------
  std::string versionString="version ";
  versionString+=static_cast<string>(VERSION);
  versionString+=", Copyright (C) 2008-2012 Pieter Kempeneers.\n\
   This program comes with ABSOLUTELY NO WARRANTY; for details type use option -h.\n\
   This is free software, and you are welcome to redistribute it\n\
   under certain conditions; use option --license for details.";
  Optionpk<bool> version_opt("\0","version",versionString,false);
  Optionpk<bool> license_opt("lic","license","show license information",false);
  Optionpk<bool> help_opt("h","help","shows this help info",false);
  Optionpk<bool> todo_opt("\0","todo","",false);
  Optionpk<string> training_opt("t", "training", "training shape file. A single shape file contains all training features (must be set as: B0, B1, B2,...) for all classes (class numbers identified by label option). Use multiple training files for bootstrap aggregation (alternative to the bag and bsize options, where a random subset is taken from a single training file)"); 
  Optionpk<string> label_opt("\0", "label", "identifier for class label in training shape file.","label"); 
  Optionpk<unsigned short> maxFeatures_opt("n", "nf", "number of features to select (0 to select all)", 0);
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

  version_opt.retrieveOption(argc,argv);
  license_opt.retrieveOption(argc,argv);
  help_opt.retrieveOption(argc,argv);
  todo_opt.retrieveOption(argc,argv);
  training_opt.retrieveOption(argc,argv);
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
  verbose_opt.retrieveOption(argc,argv);

  if(version_opt[0]||todo_opt[0]){
    std::cout << version_opt.getHelp() << std::endl;
    std::cout << "todo: " << todo_opt.getHelp() << std::endl;
    exit(0);
  }
  if(license_opt[0]){
    std::cout << Optionpk<bool>::getGPLv3License() << std::endl;
    exit(0);
  }
  if(help_opt[0]){
    std::cout << "usage: pkfs_svm -t training [OPTIONS]" << std::endl;
    exit(0);
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
  struct svm_problem prob;
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
  Histogram hist;
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
  int maxFeatures=maxFeatures_opt[0];
  double previousCost=(maxFeatures_opt[0])? 1 : 0;
  double currentCost=1;
  list<int> subset;//set of selected features (levels) for each class combination
  double cost=0;
  FeatureSelector selector;
  try{
    if(maxFeatures==nFeatures){
      subset.clear();
      for(int ifeature=0;ifeature<nFeatures;++ifeature)
        subset.push_back(ifeature);
      cost=getCost(trainingFeatures);
    }
    else if(!maxFeatures_opt[0]){
      while(currentCost-previousCost>epsilon_cost_opt[0]){
        ++maxFeatures;
        switch(selMap[selector_opt[0]]){
        case(SFFS):
          cost=selector.floating(trainingFeatures,&getCost,subset,maxFeatures,verbose_opt[0]);
          break;
        case(SFS):
          cost=selector.forward(trainingFeatures,&getCost,subset,maxFeatures,verbose_opt[0]);
          break;
        case(SBS):
          cost=selector.backward(trainingFeatures,&getCost,subset,maxFeatures,verbose_opt[0]);
          break;
        case(BFS):
          cost=selector.bruteForce(trainingFeatures,&getCost,subset,maxFeatures,verbose_opt[0]);
          break;
        default:
          std::cout << "Error: selector not supported, please use sffs, sfs, sbs or bfs" << std::endl;
          exit(1);
          break;
        }
        previousCost=currentCost;
        currentCost=cost;
      }
    }
    else{
      switch(selMap[selector_opt[0]]){
      case(SFFS):
        cost=selector.floating(trainingFeatures,&getCost,subset,maxFeatures,verbose_opt[0]);
        break;
      case(SFS):
        cost=selector.forward(trainingFeatures,&getCost,subset,maxFeatures,verbose_opt[0]);
        break;
      case(SBS):
        cost=selector.backward(trainingFeatures,&getCost,subset,maxFeatures,verbose_opt[0]);
        break;
      case(BFS):
        cost=selector.bruteForce(trainingFeatures,&getCost,subset,maxFeatures,verbose_opt[0]);
        break;
      default:
        std::cout << "Error: selector not supported, please use sffs, sfs, sbs or bfs" << std::endl;
        exit(1);
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
                            
