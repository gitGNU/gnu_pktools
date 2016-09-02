/**********************************************************************
pksvmogr_lib.cc: classify vector dataset using Support Vector Machine
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
#include <memory>
#include <algorithm>
#include "ImgReaderOgr.h"
#include "ImgWriterOgr.h"
#include "base/Optionpk.h"
#include "base/PosValue.h"
#include "algorithms/ConfusionMatrix.h"
#include "algorithms/svm.h"
#include "apps/AppFactory.h"

namespace svm{
  enum SVM_TYPE {C_SVC=0, nu_SVC=1,one_class=2, epsilon_SVR=3, nu_SVR=4};
  enum KERNEL_TYPE {linear=0,polynomial=1,radial=2,sigmoid=3};
}

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

using namespace std;
using namespace app;

/**
 * @param imgWriterOgr output classified vector dataset
 * @param app application specific option arguments
 * @return CE_None if successful, CE_Failure if failed
 **/
CPLErr ImgReaderOgr::svm(shared_ptr<ImgWriterOgr> imgWriter, const AppFactory& app){
  vector<double> priors;

  //--------------------------- command line options ------------------------------------
  Optionpk<string> training_opt("t", "training", "Training vector file. A single vector file contains all training features (must be set as: b0, b1, b2,...) for all classes (class numbers identified by label option). Use multiple training files for bootstrap aggregation (alternative to the bag and bsize options, where a random subset is taken from a single training file)");
  Optionpk<string> tlayer_opt("tln", "tln", "Training layer name(s)");
  Optionpk<string> label_opt("label", "label", "Attribute name for class label in training vector file.","label");
  Optionpk<unsigned int> balance_opt("bal", "balance", "Balance the input data to this number of samples for each class", 0);
  Optionpk<bool> random_opt("random", "random", "Randomize training data for balancing and bagging", true, 2);
  Optionpk<unsigned int> minSize_opt("min", "min", "If number of training pixels is less then min, do not take this class into account (0: consider all classes)", 0);
  Optionpk<unsigned int> band_opt("b", "band", "Band index (starting from 0, either use band option or use start to end)");
  Optionpk<unsigned int> bstart_opt("sband", "startband", "Start band sequence number");
  Optionpk<unsigned int> bend_opt("eband", "endband", "End band sequence number");
  Optionpk<double> offset_opt("offset", "offset", "Offset value for each spectral band input features: refl[band]=(DN[band]-offset[band])/scale[band]", 0.0);
  Optionpk<double> scale_opt("scale", "scale", "Scale value for each spectral band input features: refl=(DN[band]-offset[band])/scale[band] (use 0 if scale min and max in each band to -1.0 and 1.0)", 0.0);
  Optionpk<double> priors_opt("prior", "prior", "Prior probabilities for each class (e.g., -p 0.3 -p 0.3 -p 0.2 ). Used for input only (ignored for cross validation)", 0.0);
  Optionpk<string> priorimg_opt("pim", "priorimg", "Prior probability image (multi-band img with band for each class","",2);
  Optionpk<unsigned short> cv_opt("cv", "cv", "N-fold cross validation mode",0);
  Optionpk<string> cmformat_opt("cmf","cmf","Format for confusion matrix (ascii or latex)","ascii");
  Optionpk<std::string> svm_type_opt("svmt", "svmtype", "Type of SVM (C_SVC, nu_SVC,one_class, epsilon_SVR, nu_SVR)","C_SVC");
  Optionpk<std::string> kernel_type_opt("kt", "kerneltype", "Type of kernel function (linear,polynomial,radial,sigmoid) ","radial");
  Optionpk<unsigned short> kernel_degree_opt("kd", "kd", "Degree in kernel function",3);
  Optionpk<float> gamma_opt("g", "gamma", "Gamma in kernel function",1.0);
  Optionpk<float> coef0_opt("c0", "coef0", "Coef0 in kernel function",0);
  Optionpk<float> ccost_opt("cc", "ccost", "The parameter C of C_SVC, epsilon_SVR, and nu_SVR",1000);
  Optionpk<float> nu_opt("nu", "nu", "The parameter nu of nu_SVC, one_class SVM, and nu_SVR",0.5);
  Optionpk<float> epsilon_loss_opt("eloss", "eloss", "The epsilon in loss function of epsilon_SVR",0.1);
  Optionpk<int> cache_opt("cache", "cache", "Cache memory size in MB",100);
  Optionpk<float> epsilon_tol_opt("etol", "etol", "The tolerance of termination criterion",0.001);
  Optionpk<bool> shrinking_opt("shrink", "shrink", "Whether to use the shrinking heuristics",false);
  Optionpk<bool> prob_est_opt("pe", "probest", "Whether to train a SVC or SVR model for probability estimates",true,2);
  // Optionpk<bool> weight_opt("wi", "wi", "Set the parameter C of class i to weight*C, for C_SVC",true);
  Optionpk<unsigned short> comb_opt("comb", "comb", "How to combine bootstrap aggregation classifiers (0: sum rule, 1: product rule, 2: max rule). Also used to aggregate classes with rc option.",0);
  Optionpk<unsigned short> bag_opt("bag", "bag", "Number of bootstrap aggregations", 1);
  Optionpk<int> bagSize_opt("bagsize", "bagsize", "Percentage of features used from available training features for each bootstrap aggregation (one size for all classes, or a different size for each class respectively", 100);
  Optionpk<string> classBag_opt("cb", "classbag", "Output for each individual bootstrap aggregation");
  Optionpk<string> prob_opt("prob", "prob", "Probability image.");
  Optionpk<string> entropy_opt("entropy", "entropy", "Entropy image (measure for uncertainty of classifier output","",2);
  Optionpk<string> active_opt("active", "active", "Ogr output for active training sample.","",2);
  Optionpk<string> ogrformat_opt("f", "f", "Output ogr format for active training sample","SQLite");
  Optionpk<unsigned int> nactive_opt("na", "nactive", "Number of active training points",1);
  Optionpk<string> classname_opt("c", "class", "List of class names.");
  Optionpk<short> classvalue_opt("r", "reclass", "List of class values (use same order as in class opt).");
  Optionpk<short> verbose_opt("v", "verbose", "Verbose level",0,2);

  // oformat_opt.setHide(1);
  // option_opt.setHide(1);
  band_opt.setHide(1);
  bstart_opt.setHide(1);
  bend_opt.setHide(1);
  balance_opt.setHide(1);
  minSize_opt.setHide(1);
  bag_opt.setHide(1);
  bagSize_opt.setHide(1);
  comb_opt.setHide(1);
  classBag_opt.setHide(1);
  prob_opt.setHide(1);
  priorimg_opt.setHide(1);
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
  entropy_opt.setHide(1);
  active_opt.setHide(1);
  nactive_opt.setHide(1);
  random_opt.setHide(1);
  verbose_opt.setHide(2);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=training_opt.retrieveOption(app.getArgc(),app.getArgv());
    cv_opt.retrieveOption(app.getArgc(),app.getArgv());
    cmformat_opt.retrieveOption(app.getArgc(),app.getArgv());
    tlayer_opt.retrieveOption(app.getArgc(),app.getArgv());
    classname_opt.retrieveOption(app.getArgc(),app.getArgv());
    classvalue_opt.retrieveOption(app.getArgc(),app.getArgv());
    ogrformat_opt.retrieveOption(app.getArgc(),app.getArgv());
    colorTable_opt.retrieveOption(app.getArgc(),app.getArgv());
    label_opt.retrieveOption(app.getArgc(),app.getArgv());
    priors_opt.retrieveOption(app.getArgc(),app.getArgv());
    gamma_opt.retrieveOption(app.getArgc(),app.getArgv());
    ccost_opt.retrieveOption(app.getArgc(),app.getArgv());
    nodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    // Advanced options
    band_opt.retrieveOption(app.getArgc(),app.getArgv());
    bstart_opt.retrieveOption(app.getArgc(),app.getArgv());
    bend_opt.retrieveOption(app.getArgc(),app.getArgv());
    balance_opt.retrieveOption(app.getArgc(),app.getArgv());
    minSize_opt.retrieveOption(app.getArgc(),app.getArgv());
    bag_opt.retrieveOption(app.getArgc(),app.getArgv());
    bagSize_opt.retrieveOption(app.getArgc(),app.getArgv());
    comb_opt.retrieveOption(app.getArgc(),app.getArgv());
    classBag_opt.retrieveOption(app.getArgc(),app.getArgv());
    prob_opt.retrieveOption(app.getArgc(),app.getArgv());
    priorimg_opt.retrieveOption(app.getArgc(),app.getArgv());
    offset_opt.retrieveOption(app.getArgc(),app.getArgv());
    scale_opt.retrieveOption(app.getArgc(),app.getArgv());
    svm_type_opt.retrieveOption(app.getArgc(),app.getArgv());
    kernel_type_opt.retrieveOption(app.getArgc(),app.getArgv());
    kernel_degree_opt.retrieveOption(app.getArgc(),app.getArgv());
    coef0_opt.retrieveOption(app.getArgc(),app.getArgv());
    nu_opt.retrieveOption(app.getArgc(),app.getArgv());
    epsilon_loss_opt.retrieveOption(app.getArgc(),app.getArgv());
    cache_opt.retrieveOption(app.getArgc(),app.getArgv());
    epsilon_tol_opt.retrieveOption(app.getArgc(),app.getArgv());
    shrinking_opt.retrieveOption(app.getArgc(),app.getArgv());
    prob_est_opt.retrieveOption(app.getArgc(),app.getArgv());
    entropy_opt.retrieveOption(app.getArgc(),app.getArgv());
    active_opt.retrieveOption(app.getArgc(),app.getArgv());
    nactive_opt.retrieveOption(app.getArgc(),app.getArgv());
    verbose_opt.retrieveOption(app.getArgc(),app.getArgv());
    random_opt.retrieveOption(app.getArgc(),app.getArgv());
    // memory_opt.retrieveOption(app.getArgc(),app.getArgv());

    if(!doProcess){
      cout << endl;
      std::ostringstream helpStream;
      helpStream << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
      throw(helpStream.str());//help was invoked, stop processing
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
      std::cout << "training vector file: " << std::endl;
      for(unsigned int ifile=0;ifile<training_opt.size();++ifile)
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
    for(unsigned int iactive=0;iactive<activePoints.size();++iactive){
      activePoints[iactive].value=1.0;
      activePoints[iactive].posx=0.0;
      activePoints[iactive].posy=0.0;
    }

    unsigned int totalSamples=0;
    unsigned int nactive=0;
    vector<struct svm_model*> svm(nbag);
    vector<struct svm_parameter> param(nbag);

    unsigned int nclass=0;
    unsigned int nband=0;
    unsigned int startBand=2;//first two bands represent X and Y pos

    //normalize priors from command line
    if(priors_opt.size()>1){//priors from argument list
      priors.resize(priors_opt.size());
      double normPrior=0;
      for(unsigned short iclass=0;iclass<priors_opt.size();++iclass){
        priors[iclass]=priors_opt[iclass];
        normPrior+=priors[iclass];
      }
      //normalize
      for(unsigned short iclass=0;iclass<priors_opt.size();++iclass)
        priors[iclass]/=normPrior;
    }

    //convert start and end band options to vector of band indexes
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
        for(unsigned int iband=bstart_opt[ipair];iband<=bend_opt[ipair];++iband)
          band_opt.push_back(iband);
      }
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

    vector<struct svm_problem> prob(nbag);
    vector<struct svm_node *> x_space(nbag);

    for(int ibag=0;ibag<nbag;++ibag){
      //organize training data
      if(ibag<training_opt.size()){//if bag contains new training pixels
        trainingMap.clear();
        trainingPixels.clear();
        if(verbose_opt[0]>=1)
          std::cout << "reading imageVector file " << training_opt[0] << std::endl;
        ImgReaderOgr trainingReaderBag(training_opt[ibag]);
        if(band_opt.size())
          //todo: when tlayer_opt is provided, readDataImageOgr does not read any layer
          totalSamples=trainingReaderBag.readDataImageOgr(trainingMap,fields,band_opt,label_opt[0],tlayer_opt,verbose_opt[0]);
        else
          totalSamples=trainingReaderBag.readDataImageOgr(trainingMap,fields,0,0,label_opt[0],tlayer_opt,verbose_opt[0]);
        if(trainingMap.size()<2){
          string errorstring="Error: could not read at least two classes from training file, did you provide class labels in training sample (see option label)?";
          throw(errorstring);
        }
        trainingReaderBag.close();
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
        if(classname_opt.empty()){
          //std::cerr << "Warning: no class name and value pair provided for all " << nclass << " classes, using string2type<int> instead!" << std::endl;
          for(int iclass=0;iclass<nclass;++iclass){
            if(verbose_opt[0])
              std::cout << iclass << " " << cm.getClass(iclass) << " -> " << string2type<short>(cm.getClass(iclass)) << std::endl;
            classValueMap[cm.getClass(iclass)]=string2type<short>(cm.getClass(iclass));
          }
        }

        // if(priors_opt.size()==nameVector.size()){
        //  std::cerr << "Warning: please check if priors are provided in correct order!!!" << std::endl;
        //  for(int iclass=0;iclass<nameVector.size();++iclass)
        //    std::cerr << nameVector[iclass] << " " << priors_opt[iclass] << std::endl;
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
        if(bagSize_opt[iclass]<100)
          random_shuffle(trainingPixels[iclass].begin(),trainingPixels[iclass].end());
        if(verbose_opt[0]>1)
          std::cout << "nctraining (class " << iclass << "): " << nctraining << std::endl;
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
      if(verbose_opt[0]>=1)
        std::cout << "number of features: " << nFeatures << std::endl;
      unsigned int ntraining=0;
      for(short iclass=0;iclass<nclass;++iclass)
        ntraining+=trainingFeatures[iclass].size();
      if(verbose_opt[0]>=1)
        std::cout << "training size over all classes: " << ntraining << std::endl;

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
      if(verbose_opt[0]>1)
        std::cout << "SVM is now trained" << std::endl;
      if(cv_opt[0]>1){
        if(verbose_opt[0]>1)
          std::cout << "Cross validating" << std::endl;
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
      // *NOTE* Because svm_model contains pointers to svm_problem, you can
      // not free the memory used by svm_problem if you are still using the
      // svm_model produced by svm_train().
    }//for ibag
    if(cv_opt[0]>1){
      assert(cm.nReference());
      cm.setFormat(cmformat_opt[0]);
      cm.reportSE95(false);
      std::cout << cm << std::endl;
      // cout << "class #samples userAcc prodAcc" << endl;
      // double se95_ua=0;
      // double se95_pa=0;
      // double se95_oa=0;
      // double dua=0;
      // double dpa=0;
      // double doa=0;
      // for(short iclass=0;iclass<cm.nClasses();++iclass){
      //   dua=cm.ua(cm.getClass(iclass),&se95_ua);
      //   dpa=cm.pa(cm.getClass(iclass),&se95_pa);
      //   cout << cm.getClass(iclass) << " " << cm.nReference(cm.getClass(iclass)) << " " << dua << " (" << se95_ua << ")" << " " << dpa << " (" << se95_pa << ")" << endl;
      // }
      // std::cout << "Kappa: " << cm.kappa() << std::endl;
      // doa=cm.oa(&se95_oa);
      // std::cout << "Overall Accuracy: " << 100*doa << " (" << 100*se95_oa << ")"  << std::endl;
    }

    //--------------------------------- end of training -----------------------------------

    const char* pszMessage;
    void* pProgressArg=NULL;
    GDALProgressFunc pfnProgress=GDALTermProgress;
    float progress=0;
    if(!verbose_opt[0])
      pfnProgress(progress,pszMessage,pProgressArg);

    cm.clearResults();
    //notice that fields have already been set by readDataImageOgr (taking into account appropriate bands)
    for(int ivalidation=0;ivalidation<vectorCollection.size();++ivalidation){
      if(output_opt.size())
        assert(output_opt.size()==input_opt.size());
      if(verbose_opt[0])
        std::cout << "opening img reader " << input_opt[ivalidation] << std::endl;
      // imgReaderOgr.open(input_opt[ivalidation]);
      // ImgWriterOgr imgWriterOgr;

      // if(output_opt.size()){
      //   if(verbose_opt[0])
      //     std::cout << "opening img writer and copying fields from img reader" << output_opt[ivalidation] << std::endl;
      //   imgWriterOgr.open(output_opt[ivalidation],imgReaderOgr);
      // }
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
      // imgReaderOgr.close();
      // if(output_opt.size())
      //   imgWriterOgr.close();

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
      if(active_opt.size())
        activeWriter.close();

      for(int ibag=0;ibag<nbag;++ibag){
        // svm_destroy_param[ibag](&param[ibag]);
        svm_destroy_param(&param[ibag]);
        free(prob[ibag].y);
        free(prob[ibag].x);
        free(x_space[ibag]);
        svm_free_and_destroy_model(&(svm[ibag]));
      }
      return(CE_None);
    }
  }
  catch(BadConversion conversionString){
    std::cerr << "Error: did you provide class pairs names (-c) and integer values (-r) for each class in training vector?" << std::endl;
    return(CE_Failure);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    return(CE_Failure);
  }
}
