/**********************************************************************
pksvm.cc: classify raster image using Support Vector Machine
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
#include <map>
#include <algorithm>
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"
#include "imageclasses/ImgReaderOgr.h"
#include "imageclasses/ImgWriterOgr.h"
#include "base/Optionpk.h"
#include "base/PosValue.h"
#include "algorithms/ConfusionMatrix.h"
#include "algorithms/svm.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace svm{
  enum SVM_TYPE {C_SVC=0, nu_SVC=1,one_class=2, epsilon_SVR=3, nu_SVR=4};
  enum KERNEL_TYPE {linear=0,polynomial=1,radial=2,sigmoid=3};
}

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

using namespace std;

int main(int argc, char *argv[])
{
  vector<double> priors;
  
  //--------------------------- command line options ------------------------------------
  Optionpk<string> input_opt("i", "input", "input image"); 
  Optionpk<string> training_opt("t", "training", "training vector file. A single vector file contains all training features (must be set as: b0, b1, b2,...) for all classes (class numbers identified by label option). Use multiple training files for bootstrap aggregation (alternative to the bag and bsize options, where a random subset is taken from a single training file)");
  Optionpk<string> tlayer_opt("tln", "tln", "training layer name(s)");
  Optionpk<string> label_opt("label", "label", "attribute name for class label in training vector file.","label"); 
  Optionpk<unsigned int> balance_opt("bal", "balance", "balance the input data to this number of samples for each class", 0);
  Optionpk<bool> random_opt("random", "random", "randomize training data for balancing and bagging", true, 2);
  Optionpk<int> minSize_opt("min", "min", "if number of training pixels is less then min, do not take this class into account (0: consider all classes)", 0);
  Optionpk<double> start_opt("s", "start", "start band sequence number",0); 
  Optionpk<double> end_opt("e", "end", "end band sequence number (set to 0 to include all bands)", 0); 
  Optionpk<short> band_opt("b", "band", "band index (starting from 0, either use band option or use start to end)");
  Optionpk<double> offset_opt("\0", "offset", "offset value for each spectral band input features: refl[band]=(DN[band]-offset[band])/scale[band]", 0.0);
  Optionpk<double> scale_opt("\0", "scale", "scale value for each spectral band input features: refl=(DN[band]-offset[band])/scale[band] (use 0 if scale min and max in each band to -1.0 and 1.0)", 0.0);
  Optionpk<double> priors_opt("p", "prior", "prior probabilities for each class (e.g., -p 0.3 -p 0.3 -p 0.2 ). Used for input only (ignored for cross validation)", 0.0); 
  Optionpk<string> priorimg_opt("pim", "priorimg", "prior probability image (multi-band img with band for each class","",2); 
  Optionpk<unsigned short> cv_opt("cv", "cv", "n-fold cross validation mode",0);
  Optionpk<std::string> svm_type_opt("svmt", "svmtype", "type of SVM (C_SVC, nu_SVC,one_class, epsilon_SVR, nu_SVR)","C_SVC");
  Optionpk<std::string> kernel_type_opt("kt", "kerneltype", "type of kernel function (linear,polynomial,radial,sigmoid) ","radial");
  Optionpk<unsigned short> kernel_degree_opt("kd", "kd", "degree in kernel function",3);
  Optionpk<float> gamma_opt("g", "gamma", "gamma in kernel function",1.0);
  Optionpk<float> coef0_opt("c0", "coef0", "coef0 in kernel function",0);
  Optionpk<float> ccost_opt("cc", "ccost", "the parameter C of C_SVC, epsilon_SVR, and nu_SVR",1000);
  Optionpk<float> nu_opt("nu", "nu", "the parameter nu of nu_SVC, one_class SVM, and nu_SVR",0.5);
  Optionpk<float> epsilon_loss_opt("eloss", "eloss", "the epsilon in loss function of epsilon_SVR",0.1);
  Optionpk<int> cache_opt("cache", "cache", "cache memory size in MB",100);
  Optionpk<float> epsilon_tol_opt("etol", "etol", "the tolerance of termination criterion",0.001);
  Optionpk<bool> shrinking_opt("shrink", "shrink", "whether to use the shrinking heuristics",false);
  Optionpk<bool> prob_est_opt("pe", "probest", "whether to train a SVC or SVR model for probability estimates",true,2);
  // Optionpk<bool> weight_opt("wi", "wi", "set the parameter C of class i to weight*C, for C_SVC",true);
  Optionpk<unsigned short> comb_opt("comb", "comb", "how to combine bootstrap aggregation classifiers (0: sum rule, 1: product rule, 2: max rule). Also used to aggregate classes with rc option.",0); 
  Optionpk<unsigned short> bag_opt("bag", "bag", "Number of bootstrap aggregations", 1);
  Optionpk<int> bagSize_opt("bs", "bsize", "Percentage of features used from available training features for each bootstrap aggregation (one size for all classes, or a different size for each class respectively", 100);
  Optionpk<string> classBag_opt("cb", "classbag", "output for each individual bootstrap aggregation");
  Optionpk<string> mask_opt("m", "mask", "mask image (support for single mask only, see also msknodata option)"); 
  Optionpk<short> msknodata_opt("msknodata", "msknodata", "mask value(s) not to consider for classification (use negative values if only these values should be taken into account). Values will be taken over in classification image.", 0);
  Optionpk<unsigned short> nodata_opt("nodata", "nodata", "nodata value to put where image is masked as nodata", 0);
  Optionpk<string> output_opt("o", "output", "output classification image"); 
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<string> colorTable_opt("ct", "ct", "colour table in ascii format having 5 columns: id R G B ALFA (0: transparent, 255: solid)"); 
  Optionpk<string> prob_opt("prob", "prob", "probability image."); 
  Optionpk<string> entropy_opt("entropy", "entropy", "entropy image (measure for uncertainty of classifier output","",2); 
  Optionpk<string> active_opt("active", "active", "ogr output for active training sample.","",2); 
  Optionpk<string> ogrformat_opt("f", "f", "Output ogr format for active training sample","SQLite");
  Optionpk<unsigned int> nactive_opt("na", "nactive", "number of active training points",1);
  Optionpk<string> classname_opt("c", "class", "list of class names."); 
  Optionpk<short> classvalue_opt("r", "reclass", "list of class values (use same order as in class opt)."); 
  Optionpk<short> verbose_opt("v", "verbose", "set to: 0 (results only), 1 (confusion matrix), 2 (debug)",0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=training_opt.retrieveOption(argc,argv);
    tlayer_opt.retrieveOption(argc,argv);
    input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    cv_opt.retrieveOption(argc,argv);
    classname_opt.retrieveOption(argc,argv);
    classvalue_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    ogrformat_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
    label_opt.retrieveOption(argc,argv);
    priors_opt.retrieveOption(argc,argv);
    gamma_opt.retrieveOption(argc,argv);
    ccost_opt.retrieveOption(argc,argv);
    mask_opt.retrieveOption(argc,argv);
    msknodata_opt.retrieveOption(argc,argv);
    nodata_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    start_opt.retrieveOption(argc,argv);
    end_opt.retrieveOption(argc,argv);
    balance_opt.retrieveOption(argc,argv);
    minSize_opt.retrieveOption(argc,argv);
    bag_opt.retrieveOption(argc,argv);
    bagSize_opt.retrieveOption(argc,argv);
    comb_opt.retrieveOption(argc,argv);
    classBag_opt.retrieveOption(argc,argv);
    prob_opt.retrieveOption(argc,argv);
    priorimg_opt.retrieveOption(argc,argv);
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
    entropy_opt.retrieveOption(argc,argv);
    active_opt.retrieveOption(argc,argv);
    nactive_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
    random_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }

  if(entropy_opt[0]=="")
    entropy_opt.clear();
  if(active_opt[0]=="")
    active_opt.clear();
  if(priorimg_opt[0]=="")
    priorimg_opt.clear();


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

  assert(training_opt.size());

  if(verbose_opt[0]>=1){
    if(input_opt.size())
      std::cout << "input filename: " << input_opt[0] << std::endl;
    if(mask_opt.size())
      std::cout << "mask filename: " << mask_opt[0] << std::endl;
    std::cout << "training vector file: " << std::endl;
    for(int ifile=0;ifile<training_opt.size();++ifile)
      std::cout << training_opt[ifile] << std::endl;
    std::cout << "verbose: " << verbose_opt[0] << std::endl;
  }
  unsigned short nbag=(training_opt.size()>1)?training_opt.size():bag_opt[0];
  if(verbose_opt[0]>=1)
    std::cout << "number of bootstrap aggregations: " << nbag << std::endl;

  ImgWriterOgr activeWriter;
  if(active_opt.size()){
    prob_est_opt[0]=true;
    ImgReaderOgr trainingReader(training_opt[0]);
    activeWriter.open(active_opt[0],ogrformat_opt[0]);
    activeWriter.createLayer(active_opt[0],trainingReader.getProjection(),wkbPoint,NULL);
    activeWriter.copyFields(trainingReader);
  }
  vector<PosValue> activePoints(nactive_opt[0]);
  for(int iactive=0;iactive<activePoints.size();++iactive){
    activePoints[iactive].value=1.0;
    activePoints[iactive].posx=0.0;
    activePoints[iactive].posy=0.0;
  }
  
  unsigned int totalSamples=0;
  unsigned int nactive=0;
  vector<struct svm_model*> svm(nbag);
  vector<struct svm_parameter> param(nbag);

  short nclass=0;
  int nband=0;
  int startBand=2;//first two bands represent X and Y pos

  //normalize priors from command line
  if(priors_opt.size()>1){//priors from argument list
    priors.resize(priors_opt.size());
    double normPrior=0;
    for(short iclass=0;iclass<priors_opt.size();++iclass){
      priors[iclass]=priors_opt[iclass];
      normPrior+=priors[iclass];
    }
    //normalize
    for(short iclass=0;iclass<priors_opt.size();++iclass)
      priors[iclass]/=normPrior;
  }

  //sort bands
  if(band_opt.size())
    std::sort(band_opt.begin(),band_opt.end());

  map<string,short> classValueMap;
  vector<std::string> nameVector;
  if(classname_opt.size()){
    assert(classname_opt.size()==classvalue_opt.size());
    for(int iclass=0;iclass<classname_opt.size();++iclass)
      classValueMap[classname_opt[iclass]]=classvalue_opt[iclass];
  }

  //----------------------------------- Training -------------------------------
  ConfusionMatrix cm;
  vector< vector<double> > offset(nbag);
  vector< vector<double> > scale(nbag);
  map<string,Vector2d<float> > trainingMap;
  vector< Vector2d<float> > trainingPixels;//[class][sample][band]
  vector<string> fields;

  vector<struct svm_problem> prob(nbag);
  vector<struct svm_node *> x_space(nbag);

  for(int ibag=0;ibag<nbag;++ibag){
    //organize training data
    if(ibag<training_opt.size()){//if bag contains new training pixels
      trainingMap.clear();
      trainingPixels.clear();
      if(verbose_opt[0]>=1)
        std::cout << "reading imageVector file " << training_opt[0] << std::endl;
      try{
	ImgReaderOgr trainingReaderBag(training_opt[ibag]);
        if(band_opt.size())
          totalSamples=trainingReaderBag.readDataImageOgr(trainingMap,fields,band_opt,label_opt[0],tlayer_opt,verbose_opt[0]);
        else
          totalSamples=trainingReaderBag.readDataImageOgr(trainingMap,fields,start_opt[0],end_opt[0],label_opt[0],tlayer_opt,verbose_opt[0]);
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
      //todo: delete class 0 ?
      // if(verbose_opt[0]>=1)
      //   std::cout << "erasing class 0 from training set (" << trainingMap[0].size() << " from " << totalSamples << ") samples" << std::endl;
      // totalSamples-=trainingMap[0].size();
      // trainingMap.erase(0);

      //convert map to vector
      // short iclass=0;
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
        nband=trainingPixels[0][0].size()-2;//X and Y//trainingPixels[0][0].size();
      }
      else{
        assert(nclass==trainingPixels.size());
        assert(nband==trainingPixels[0][0].size()-2);
      }

      //do not remove outliers here: could easily be obtained through ogr2ogr -where 'B2<110' output.shp input.shp
      //balance training data
      if(balance_opt[0]>0){
        while(balance_opt.size()<nclass)
          balance_opt.push_back(balance_opt.back());
        if(random_opt[0])
          srand(time(NULL));
        totalSamples=0;
        for(short iclass=0;iclass<nclass;++iclass){
          if(trainingPixels[iclass].size()>balance_opt[iclass]){
            while(trainingPixels[iclass].size()>balance_opt[iclass]){
              int index=rand()%trainingPixels[iclass].size();
              trainingPixels[iclass].erase(trainingPixels[iclass].begin()+index);
            }
          }
          else{
            int oldsize=trainingPixels[iclass].size();
            for(int isample=trainingPixels[iclass].size();isample<balance_opt[iclass];++isample){
              int index = rand()%oldsize;
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
      for(int iband=0;iband<nband;++iband){
        if(verbose_opt[0]>=1)
          std::cout << "scaling for band" << iband << std::endl;
        offset[ibag][iband]=(offset_opt.size()==1)?offset_opt[0]:offset_opt[iband];
        scale[ibag][iband]=(scale_opt.size()==1)?scale_opt[0]:scale_opt[iband];
        //search for min and maximum
        if(scale[ibag][iband]<=0){
          float theMin=trainingPixels[0][0][iband+startBand];
          float theMax=trainingPixels[0][0][iband+startBand];
          for(short iclass=0;iclass<nclass;++iclass){
            for(int isample=0;isample<trainingPixels[iclass].size();++isample){
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
      for(int iband=0;iband<nband;++iband){
        offset[ibag][iband]=offset[0][iband];
        scale[ibag][iband]=scale[0][iband];
      }
    }
      
    if(!ibag){
      if(priors_opt.size()==1){//default: equal priors for each class
        priors.resize(nclass);
        for(short iclass=0;iclass<nclass;++iclass)
          priors[iclass]=1.0/nclass;
      }
      assert(priors_opt.size()==1||priors_opt.size()==nclass);

      //set bagsize for each class if not done already via command line
      while(bagSize_opt.size()<nclass)
        bagSize_opt.push_back(bagSize_opt.back());

      if(verbose_opt[0]>=1){
        std::cout << "number of bands: " << nband << std::endl;
        std::cout << "number of classes: " << nclass << std::endl;
	if(priorimg_opt.empty()){
	  std::cout << "priors:";
	  for(short iclass=0;iclass<nclass;++iclass)
	    std::cout << " " << priors[iclass];
	  std::cout << std::endl;
	}
      }
      map<string,Vector2d<float> >::iterator mapit=trainingMap.begin();
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

      // if(priors_opt.size()==nameVector.size()){
      // 	std::cerr << "Warning: please check if priors are provided in correct order!!!" << std::endl;
      // 	for(int iclass=0;iclass<nameVector.size();++iclass)
      // 	  std::cerr << nameVector[iclass] << " " << priors_opt[iclass] << std::endl;
      // }
    }//if(!ibag)

    //Calculate features of training set
    vector< Vector2d<float> > trainingFeatures(nclass);
    for(short iclass=0;iclass<nclass;++iclass){
      int nctraining=0;
      if(verbose_opt[0]>=1)
        std::cout << "calculating features for class " << iclass << std::endl;
      if(random_opt[0])
        srand(time(NULL));
      nctraining=(bagSize_opt[iclass]<100)? trainingPixels[iclass].size()/100.0*bagSize_opt[iclass] : trainingPixels[iclass].size();//bagSize_opt[iclass] given in % of training size
      if(nctraining<=0)
        nctraining=1;
      assert(nctraining<=trainingPixels[iclass].size());
      int index=0;
      if(bagSize_opt[iclass]<100)
        random_shuffle(trainingPixels[iclass].begin(),trainingPixels[iclass].end());
      
      trainingFeatures[iclass].resize(nctraining);
      for(int isample=0;isample<nctraining;++isample){
        //scale pixel values according to scale and offset!!!
        for(int iband=0;iband<nband;++iband){
          float value=trainingPixels[iclass][isample][iband+startBand];
          trainingFeatures[iclass][isample].push_back((value-offset[ibag][iband])/scale[ibag][iband]);
        }
      }
      assert(trainingFeatures[iclass].size()==nctraining);
    }
    
    unsigned int nFeatures=trainingFeatures[0][0].size();
    unsigned int ntraining=0;
    for(short iclass=0;iclass<nclass;++iclass)
      ntraining+=trainingFeatures[iclass].size();

    prob[ibag].l=ntraining;
    prob[ibag].y = Malloc(double,prob[ibag].l);
    prob[ibag].x = Malloc(struct svm_node *,prob[ibag].l);
    x_space[ibag] = Malloc(struct svm_node,(nFeatures+1)*ntraining);
    unsigned long int spaceIndex=0;
    int lIndex=0;
    for(short iclass=0;iclass<nclass;++iclass){
      for(int isample=0;isample<trainingFeatures[iclass].size();++isample){
        prob[ibag].x[lIndex]=&(x_space[ibag][spaceIndex]);
        for(int ifeature=0;ifeature<nFeatures;++ifeature){
          x_space[ibag][spaceIndex].index=ifeature+1;
          x_space[ibag][spaceIndex].value=trainingFeatures[iclass][isample][ifeature];
          ++spaceIndex;
        }
        x_space[ibag][spaceIndex++].index=-1;
        prob[ibag].y[lIndex]=iclass;
        ++lIndex;
      }
    }
    assert(lIndex==prob[ibag].l);

    //set SVM parameters through command line options
    param[ibag].svm_type = svmMap[svm_type_opt[0]];
    param[ibag].kernel_type = kernelMap[kernel_type_opt[0]];
    param[ibag].degree = kernel_degree_opt[0];
    param[ibag].gamma = (gamma_opt[0]>0)? gamma_opt[0] : 1.0/nFeatures;
    param[ibag].coef0 = coef0_opt[0];
    param[ibag].nu = nu_opt[0];
    param[ibag].cache_size = cache_opt[0];
    param[ibag].C = ccost_opt[0];
    param[ibag].eps = epsilon_tol_opt[0];
    param[ibag].p = epsilon_loss_opt[0];
    param[ibag].shrinking = (shrinking_opt[0])? 1 : 0;
    param[ibag].probability = (prob_est_opt[0])? 1 : 0;
    param[ibag].nr_weight = 0;//not used: I use priors and balancing
    param[ibag].weight_label = NULL;
    param[ibag].weight = NULL;
    param[ibag].verbose=(verbose_opt[0]>1)? true:false;

    if(verbose_opt[0]>1)
      std::cout << "checking parameters" << std::endl;
    svm_check_parameter(&prob[ibag],&param[ibag]);
    if(verbose_opt[0])
      std::cout << "parameters ok, training" << std::endl;
    svm[ibag]=svm_train(&prob[ibag],&param[ibag]);
    if(cv_opt[0]>1){
      double *target = Malloc(double,prob[ibag].l);
      svm_cross_validation(&prob[ibag],&param[ibag],cv_opt[0],target);
      assert(param[ibag].svm_type != EPSILON_SVR&&param[ibag].svm_type != NU_SVR);//only for regression

      for(int i=0;i<prob[ibag].l;i++){
	string refClassName=nameVector[prob[ibag].y[i]];
	string className=nameVector[target[i]];
	if(classValueMap.size())
	  cm.incrementResult(type2string<short>(classValueMap[refClassName]),type2string<short>(classValueMap[className]),1.0/nbag);
	else
	  cm.incrementResult(cm.getClass(prob[ibag].y[i]),cm.getClass(target[i]),1.0/nbag);
      }
      free(target);
    }
    if(verbose_opt[0]>1)
      std::cout << "SVM is now trained" << std::endl;
    // *NOTE* Because svm_model contains pointers to svm_problem, you can
    // not free the memory used by svm_problem if you are still using the
    // svm_model produced by svm_train(). 
  }//for ibag
  if(cv_opt[0]>1){
    assert(cm.nReference());
    std::cout << cm << std::endl;
    cout << "class #samples userAcc prodAcc" << endl;
    double se95_ua=0;
    double se95_pa=0;
    double se95_oa=0;
    double dua=0;
    double dpa=0;
    double doa=0;
    for(short iclass=0;iclass<cm.nClasses();++iclass){
      dua=cm.ua(cm.getClass(iclass),&se95_ua);
      dpa=cm.pa(cm.getClass(iclass),&se95_pa);
      cout << cm.getClass(iclass) << " " << cm.nReference(cm.getClass(iclass)) << " " << dua << " (" << se95_ua << ")" << " " << dpa << " (" << se95_pa << ")" << endl;
    }
    std::cout << "Kappa: " << cm.kappa() << std::endl;
    doa=cm.oa(&se95_oa);
    std::cout << "Overall Accuracy: " << 100*doa << " (" << 100*se95_oa << ")"  << std::endl;
  }

  //--------------------------------- end of training -----------------------------------
  if(input_opt.empty())
    exit(0);

  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  float progress=0;
  if(!verbose_opt[0])
    pfnProgress(progress,pszMessage,pProgressArg);
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
    ImgReaderGdal testImage;
    try{
      if(verbose_opt[0]>=1)
        std::cout << "opening image " << input_opt[0] << std::endl; 
      testImage.open(input_opt[0]);
    }
    catch(string error){
      cerr << error << std::endl;
      exit(2);
    }
    ImgReaderGdal maskReader;
    if(mask_opt.size()){
      try{
        if(verbose_opt[0]>=1)
          std::cout << "opening mask image file " << mask_opt[0] << std::endl;
        maskReader.open(mask_opt[0]);
        assert(maskReader.nrOfCol()==testImage.nrOfCol());
        assert(maskReader.nrOfRow()==testImage.nrOfRow());
      }
      catch(string error){
        cerr << error << std::endl;
        exit(2);
      }
      catch(...){
        cerr << "error catched" << std::endl;
        exit(1);
      }
    }
    ImgReaderGdal priorReader;
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
        cerr << "error catched" << std::endl;
        exit(1);
      }
    }

    int nrow=testImage.nrOfRow();
    int ncol=testImage.nrOfCol();
    if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
      string theInterleave="INTERLEAVE=";
      theInterleave+=testImage.getInterleave();
      option_opt.push_back(theInterleave);
    }
    vector<char> classOut(ncol);//classified line for writing to image file

    //   assert(nband==testImage.nrOfBand());
    ImgWriterGdal classImageBag;
    ImgWriterGdal classImageOut;
    ImgWriterGdal probImage;
    ImgWriterGdal entropyImage;

    string imageType=testImage.getImageType();
    if(oformat_opt.size())//default
      imageType=oformat_opt[0];
    try{
      assert(output_opt.size());
      if(verbose_opt[0]>=1)
        std::cout << "opening class image for writing output " << output_opt[0] << std::endl;
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
      cerr << error << std::endl;
    }
  
    for(int iline=0;iline<nrow;++iline){
      vector<float> buffer(ncol);
      vector<short> lineMask;
      if(mask_opt.size())
        lineMask.resize(maskReader.nrOfCol());
      Vector2d<float> linePrior;
      if(priorimg_opt.size())
	 linePrior.resize(nclass,ncol);//prior prob for each class
      Vector2d<float> hpixel(ncol);
      Vector2d<float> probOut(nclass,ncol);//posterior prob for each (internal) class
      vector<float> entropy(ncol);
      Vector2d<char> classBag;//classified line for writing to image file
      if(classBag_opt.size())
        classBag.resize(nbag,ncol);
      try{
        if(band_opt.size()){
          for(int iband=0;iband<band_opt.size();++iband){
            if(verbose_opt[0]==2)
              std::cout << "reading band " << band_opt[iband] << std::endl;
            assert(band_opt[iband]>=0);
            assert(band_opt[iband]<testImage.nrOfBand());
            testImage.readData(buffer,GDT_Float32,iline,band_opt[iband]);
            for(int icol=0;icol<ncol;++icol)
              hpixel[icol].push_back(buffer[icol]);
          }
        }
        else{
          for(int iband=start_opt[0];iband<start_opt[0]+nband;++iband){
            if(verbose_opt[0]==2)
              std::cout << "reading band " << iband << std::endl;
            assert(iband>=0);
            assert(iband<testImage.nrOfBand());
            testImage.readData(buffer,GDT_Float32,iline,iband);
            for(int icol=0;icol<ncol;++icol)
              hpixel[icol].push_back(buffer[icol]);
          }
        }
      }
      catch(string theError){
        cerr << "Error reading " << input_opt[0] << ": " << theError << std::endl;
        exit(3);
      }
      catch(...){
        cerr << "error catched" << std::endl;
        exit(3);
      }
      assert(nband==hpixel[0].size());
      if(verbose_opt[0]>1)
        std::cout << "used bands: " << nband << std::endl;
      //read mask
      if(!lineMask.empty()){
        try{
          maskReader.readData(lineMask,GDT_Int16,iline);
        }
        catch(string theError){
          cerr << "Error reading " << mask_opt[0] << ": " << theError << std::endl;
          exit(3);
        }
        catch(...){
          cerr << "error catched" << std::endl;
          exit(3);
        }
      }
      //read prior
      if(priorimg_opt.size()){
        try{
	  for(short iclass=0;iclass<nclass;++iclass){
	    if(verbose_opt.size()>1)
	      std::cout << "Reading " << priorimg_opt[0] << " band " << iclass << " line " << iline << std::endl;
	    priorReader.readData(linePrior[iclass],GDT_Float32,iline,iclass);
	  }
        }
        catch(string theError){
	  std::cerr << "Error reading " << priorimg_opt[0] << ": " << theError << std::endl;
          exit(3);
        }
        catch(...){
          cerr << "error catched" << std::endl;
          exit(3);
        }
      }
      //process per pixel
      for(int icol=0;icol<ncol;++icol){
        assert(hpixel[icol].size()==nband);
        bool masked=false;
        if(!lineMask.empty()){
          short theMask=0;
          for(short ivalue=0;ivalue<msknodata_opt.size();++ivalue){
            if(msknodata_opt[ivalue]>=0){//values set in msknodata_opt are invalid
              if(lineMask[icol]==msknodata_opt[ivalue]){
                theMask=lineMask[icol];
                masked=true;
                break;
              }
            }
            else{//only values set in msknodata_opt are valid
              if(lineMask[icol]!=-msknodata_opt[ivalue]){
                theMask=lineMask[icol];
                masked=true;
              }
              else{
                masked=false;
                break;
              }
            }
          }
          if(masked){
            if(classBag_opt.size())
              for(int ibag=0;ibag<nbag;++ibag)
                classBag[ibag][icol]=theMask;
            classOut[icol]=theMask;
            continue;
          }
        }
        bool valid=false;
        for(int iband=0;iband<hpixel[icol].size();++iband){
          if(hpixel[icol][iband]){
            valid=true;
            break;
          }
        }
        if(!valid){
          if(classBag_opt.size())
            for(int ibag=0;ibag<nbag;++ibag)
              classBag[ibag][icol]=nodata_opt[0];
          classOut[icol]=nodata_opt[0];
          continue;//next column
        }
        for(short iclass=0;iclass<nclass;++iclass)
          probOut[iclass][icol]=0;
	if(verbose_opt[0]>1)
	  std::cout << "begin classification " << std::endl;
        //----------------------------------- classification -------------------
        for(int ibag=0;ibag<nbag;++ibag){
          vector<double> result(nclass);
          struct svm_node *x;
          x = (struct svm_node *) malloc((nband+1)*sizeof(struct svm_node));
          for(int iband=0;iband<nband;++iband){
            x[iband].index=iband+1;
            x[iband].value=(hpixel[icol][iband]-offset[ibag][iband])/scale[ibag][iband];
          }
          x[nband].index=-1;//to end svm feature vector
          double predict_label=0;
          vector<float> prValues(nclass);
          float maxP=0;
          if(!prob_est_opt[0]){
            predict_label = svm_predict(svm[ibag],x);
            for(short iclass=0;iclass<nclass;++iclass){
              if(iclass==static_cast<short>(predict_label))
                result[iclass]=1;
              else
                result[iclass]=0;
            }
          }
          else{
            assert(svm_check_probability_model(svm[ibag]));
            predict_label = svm_predict_probability(svm[ibag],x,&(result[0]));
          }
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
          for(short iclass=0;iclass<nclass;++iclass){
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
              //   classBag[ibag][icol]=iclass;
              // }
              if(result[iclass]>maxP){
                maxP=result[iclass];
                classBag[ibag][icol]=iclass;
              }
            }
          }
          free(x);
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
        for(int ibag=0;ibag<nbag;++ibag)
          classImageBag.writeData(classBag[ibag],GDT_Byte,iline,ibag);
      if(prob_opt.size()){
        for(short iclass=0;iclass<nclass;++iclass)
          probImage.writeData(probOut[iclass],GDT_Float32,iline,iclass);
      }
      if(entropy_opt.size()){
        entropyImage.writeData(entropy,GDT_Float32,iline);
      }
      classImageOut.writeData(classOut,GDT_Byte,iline);
      if(!verbose_opt[0]){
        progress=static_cast<float>(iline+1.0)/classImageOut.nrOfRow();
        pfnProgress(progress,pszMessage,pProgressArg);
      }
    }
    //write active learning points
    if(active_opt.size()){
      for(int iactive=0;iactive<activePoints.size();++iactive){
	std::map<string,double> pointMap;
	for(int iband=0;iband<testImage.nrOfBand();++iband){
	  double value;
	  testImage.readData(value,GDT_Float64,static_cast<int>(activePoints[iactive].posx),static_cast<int>(activePoints[iactive].posy),iband);
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
    for(int ivalidation=0;ivalidation<input_opt.size();++ivalidation){
      if(output_opt.size())
	assert(output_opt.size()==input_opt.size());
      if(verbose_opt[0])
        std::cout << "opening img reader " << input_opt[ivalidation] << std::endl;
      imgReaderOgr.open(input_opt[ivalidation]);
      ImgWriterOgr imgWriterOgr;

      if(output_opt.size()){
	if(verbose_opt[0])
	  std::cout << "opening img writer and copying fields from img reader" << output_opt[ivalidation] << std::endl;
	imgWriterOgr.open(output_opt[ivalidation],imgReaderOgr);
      }
      if(verbose_opt[0])
	cout << "number of layers in input ogr file: " << imgReaderOgr.getLayerCount() << endl;
      for(int ilayer=0;ilayer<imgReaderOgr.getLayerCount();++ilayer){
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
	    std::cout << "feature " << ifeature << std::endl;
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
	  for(short iclass=0;iclass<nclass;++iclass)
	    probOut[iclass]=0;
	  for(int ibag=0;ibag<nbag;++ibag){
	    for(int iband=0;iband<nband;++iband){
	      validationFeature.push_back((validationPixel[iband]-offset[ibag][iband])/scale[ibag][iband]);
	      if(verbose_opt[0]==2)
		std::cout << " " << validationFeature.back();
	    }
	    if(verbose_opt[0]==2)
	      std::cout << std::endl;
	    vector<double> result(nclass);
	    struct svm_node *x;
	    x = (struct svm_node *) malloc((validationFeature.size()+1)*sizeof(struct svm_node));
	    for(int i=0;i<validationFeature.size();++i){
	      x[i].index=i+1;
	      x[i].value=validationFeature[i];
	    }

	    x[validationFeature.size()].index=-1;//to end svm feature vector
	    double predict_label=0;
	    if(!prob_est_opt[0]){
	      predict_label = svm_predict(svm[ibag],x);
	      for(short iclass=0;iclass<nclass;++iclass){
		if(iclass==static_cast<short>(predict_label))
		  result[iclass]=1;
		else
		  result[iclass]=0;
	      }
	    }
	    else{
	      assert(svm_check_probability_model(svm[ibag]));
	      predict_label = svm_predict_probability(svm[ibag],x,&(result[0]));
	    }
	    if(verbose_opt[0]>1){
	      std::cout << "predict_label: " << predict_label << std::endl;
	      for(int iclass=0;iclass<result.size();++iclass)
		std::cout << result[iclass] << " ";
	      std::cout << std::endl;
	    }

	    //calculate posterior prob of bag 
	    for(short iclass=0;iclass<nclass;++iclass){
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
	    free(x);
	  }//for ibag

	  //search for max class prob
	  float maxBag=0;
	  float normBag=0;
	  string classOut="Unclassified";
	  for(short iclass=0;iclass<nclass;++iclass){
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
      cout << "class #samples userAcc prodAcc" << endl;
      double se95_ua=0;
      double se95_pa=0;
      double se95_oa=0;
      double dua=0;
      double dpa=0;
      double doa=0;
      for(short iclass=0;iclass<cm.nClasses();++iclass){
	dua=cm.ua_pct(cm.getClass(iclass),&se95_ua);
	dpa=cm.pa_pct(cm.getClass(iclass),&se95_pa);
	cout << cm.getClass(iclass) << " " << cm.nReference(cm.getClass(iclass)) << " " << dua << " (" << se95_ua << ")" << " " << dpa << " (" << se95_pa << ")" << endl;
      }
      std::cout << "Kappa: " << cm.kappa() << std::endl;
      doa=cm.oa(&se95_oa);
      std::cout << "Overall Accuracy: " << 100*doa << " (" << 100*se95_oa << ")"  << std::endl;
    }
  }
  try{
    if(active_opt.size())
      activeWriter.close();
  }
  catch(string errorString){
    std::cerr << "Error: errorString" << std::endl;
  }

  for(int ibag=0;ibag<nbag;++ibag){
    // svm_destroy_param[ibag](&param[ibag]);
    svm_destroy_param(&param[ibag]);
    free(prob[ibag].y);
    free(prob[ibag].x);
    free(x_space[ibag]);
    svm_free_and_destroy_model(&(svm[ibag]));
  }
  return 0;
}
