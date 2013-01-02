/**********************************************************************
pkclassify_svm.cc: classify raster image using Support Vector Machine
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
#include <map>
#include <algorithm>
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"
#include "imageclasses/ImgReaderOgr.h"
#include "imageclasses/ImgWriterOgr.h"
#include "base/Optionpk.h"
#include "algorithms/ConfusionMatrix.h"
#include "algorithms/svm.h"
#include "pkclassify_nn.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

int main(int argc, char *argv[])
{
  map<short,int> reclassMap;
  vector<int> vreclass;
  vector<double> priors;
  
  //--------------------------- command line options ------------------------------------

  std::string versionString="version ";
  versionString+=VERSION;
  versionString+=", Copyright (C) 2008-2012 Pieter Kempeneers.\n\
   This program comes with ABSOLUTELY NO WARRANTY; for details type use option -h.\n\
   This is free software, and you are welcome to redistribute it\n\
   under certain conditions; use option --license for details.";
  Optionpk<bool> version_opt("\0","version",versionString,false);
  Optionpk<bool> license_opt("lic","license","show license information",false);
  Optionpk<bool> help_opt("h","help","shows this help info",false);
  Optionpk<bool> todo_opt("\0","todo","",false);
  Optionpk<string> input_opt("i", "input", "input image"); 
  Optionpk<string> training_opt("t", "training", "training shape file. A single shape file contains all training features (must be set as: B0, B1, B2,...) for all classes (class numbers identified by label option). Use multiple training files for bootstrap aggregation (alternative to the bag and bsize options, where a random subset is taken from a single training file)"); 
  Optionpk<string> label_opt("\0", "label", "identifier for class label in training shape file.","label"); 
  Optionpk<unsigned short> reclass_opt("\0", "rc", "reclass code (e.g. --rc=12 --rc=23 to reclass first two classes to 12 and 23 resp.).", 0);
  Optionpk<unsigned int> balance_opt("\0", "balance", "balance the input data to this number of samples for each class", 0);
  Optionpk<int> minSize_opt("m", "min", "if number of training pixels is less then min, do not take this class into account", 0);
  Optionpk<double> start_opt("s", "start", "start band sequence number (set to 0)",0); 
  Optionpk<double> end_opt("e", "end", "end band sequence number (set to 0 for all bands)", 0); 
  Optionpk<double> offset_opt("\0", "offset", "offset value for each spectral band input features: refl[band]=(DN[band]-offset[band])/scale[band]", 0.0);
  Optionpk<double> scale_opt("\0", "scale", "scale value for each spectral band input features: refl=(DN[band]-offset[band])/scale[band] (use 0 if scale min and max in each band to -1.0 and 1.0)", 0.0);
  Optionpk<unsigned short> aggreg_opt("a", "aggreg", "how to combine aggregated classifiers, see also rc option (0: no aggregation, 1: sum rule, 2: max rule).",0);
  Optionpk<double> priors_opt("p", "prior", "prior probabilities for each class (e.g., -p 0.3 -p 0.3 -p 0.2 )", 0.0); 


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
  Optionpk<unsigned int> cv_opt("cv", "cv", "n-fold cross validation mode",0);
  Optionpk<unsigned short> comb_opt("c", "comb", "how to combine bootstrap aggregation classifiers (0: sum rule, 1: product rule, 2: max rule). Also used to aggregate classes with rc option.",0); 
  Optionpk<unsigned short> bag_opt("\0", "bag", "Number of bootstrap aggregations", 1);
  Optionpk<int> bagSize_opt("\0", "bsize", "Percentage of features used from available training features for each bootstrap aggregation", 100);
  Optionpk<string> classBag_opt("\0", "class", "output for each individual bootstrap aggregation");
  Optionpk<string> mask_opt("\0", "mask", "mask image (see also mvalue option"); 
  Optionpk<short> maskValue_opt("\0", "mvalue", "mask value(s) not to consider for classification (use negative values if only these values should be taken into account). Values will be taken over in classification image.", 0);
  Optionpk<unsigned short> flag_opt("f", "flag", "flag to put where image is invalid.", 0);
  Optionpk<string> output_opt("o", "output", "output classification image"); 
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image");
  Optionpk<string> option_opt("co", "co", "options: NAME=VALUE [-co COMPRESS=LZW] [-co INTERLEAVE=BAND]");
  Optionpk<string> colorTable_opt("\0", "ct", "colour table in ascii format having 5 columns: id R G B ALFA (0: transparent, 255: solid)"); 
  Optionpk<string> prob_opt("\0", "prob", "probability image."); 
  Optionpk<short> verbose_opt("v", "verbose", "set to: 0 (results only), 1 (confusion matrix), 2 (debug)",0);

  version_opt.retrieveOption(argc,argv);
  license_opt.retrieveOption(argc,argv);
  help_opt.retrieveOption(argc,argv);
  todo_opt.retrieveOption(argc,argv);

  input_opt.retrieveOption(argc,argv);
  training_opt.retrieveOption(argc,argv);
  label_opt.retrieveOption(argc,argv);
  reclass_opt.retrieveOption(argc,argv);
  balance_opt.retrieveOption(argc,argv);
  minSize_opt.retrieveOption(argc,argv);
  start_opt.retrieveOption(argc,argv);
  end_opt.retrieveOption(argc,argv);
  offset_opt.retrieveOption(argc,argv);
  scale_opt.retrieveOption(argc,argv);
  aggreg_opt.retrieveOption(argc,argv);
  priors_opt.retrieveOption(argc,argv);
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
  comb_opt.retrieveOption(argc,argv);
  bag_opt.retrieveOption(argc,argv);
  bagSize_opt.retrieveOption(argc,argv);
  classBag_opt.retrieveOption(argc,argv);
  mask_opt.retrieveOption(argc,argv);
  maskValue_opt.retrieveOption(argc,argv);
  flag_opt.retrieveOption(argc,argv);
  output_opt.retrieveOption(argc,argv);
  oformat_opt.retrieveOption(argc,argv);
  colorTable_opt.retrieveOption(argc,argv);
  option_opt.retrieveOption(argc,argv);
  prob_opt.retrieveOption(argc,argv);
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
    std::cout << "usage: pkclassify_svm -i testimage -o outputimage -t training [OPTIONS]" << std::endl;
    exit(0);
  }

  if(verbose_opt[0]>=1){
    if(input_opt.size())
      std::cout << "image filename: " << input_opt[0] << std::endl;
    if(mask_opt.size())
      std::cout << "mask filename: " << mask_opt[0] << std::endl;
    if(training_opt[0].size()){
      std::cout << "training shape file: " << std::endl;
      for(int ifile=0;ifile<training_opt.size();++ifile)
        std::cout << training_opt[ifile] << std::endl;
    }
    else
      cerr << "no training file set!" << std::endl;
    std::cout << "verbose: " << verbose_opt[0] << std::endl;
  }
  unsigned short nbag=(training_opt.size()>1)?training_opt.size():bag_opt[0];
  if(verbose_opt[0]>=1)
    std::cout << "number of bootstrap aggregations: " << nbag << std::endl;
  
  unsigned int totalSamples=0;
  int nreclass=0;
  vector<int> vcode;//unique class codes in recode string
  vector<struct svm_model*> svm(nbag);
  vector<struct svm_parameter> param(nbag);

  unsigned int nclass=0;
  int nband=0;
  int startBand=2;//first two bands represent X and Y pos

  vector< vector<double> > offset(nbag);
  vector< vector<double> > scale(nbag);
  vector< Vector2d<float> > trainingPixels;//[class][sample][band]

  if(reclass_opt.size()>1){
    vreclass.resize(reclass_opt.size());
    for(int iclass=0;iclass<reclass_opt.size();++iclass){
      reclassMap[iclass]=reclass_opt[iclass];
      vreclass[iclass]=reclass_opt[iclass];
    }
  }
  if(priors_opt.size()>1){//priors from argument list
    priors.resize(priors_opt.size());
    double normPrior=0;
    for(int iclass=0;iclass<priors_opt.size();++iclass){
      priors[iclass]=priors_opt[iclass];
      normPrior+=priors[iclass];
    }
    //normalize
    for(int iclass=0;iclass<priors_opt.size();++iclass)
      priors[iclass]/=normPrior;
  }

  //----------------------------------- Training -------------------------------
  vector<struct svm_problem> prob(nbag);
  vector<struct svm_node *> x_space(nbag);
  // struct svm_node *x_space;
  vector<string> fields;
  for(int ibag=0;ibag<nbag;++ibag){
    //organize training data
    if(ibag<training_opt.size()){//if bag contains new training pixels
      trainingPixels.clear();
      map<int,Vector2d<float> > trainingMap;
      if(verbose_opt[0]>=1)
        std::cout << "reading imageShape file " << training_opt[0] << std::endl;
      try{
        totalSamples=readDataImageShape(training_opt[ibag],trainingMap,fields,start_opt[0],end_opt[0],label_opt[0],verbose_opt[0]);
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
      if(!ibag){
        nclass=trainingPixels.size();
        nband=trainingPixels[0][0].size()-2;//X and Y//trainingPixels[0][0].size();
      }
      else{
        assert(nclass==trainingPixels.size());
        assert(nband==trainingPixels[0][0].size()-2);
      }
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
      offset[ibag].resize(nband);
      scale[ibag].resize(nband);
      if(offset_opt.size()>1)
        assert(offset_opt.size()==nband);
      if(scale_opt.size()>1)
        assert(scale_opt.size()==nband);
      Histogram hist;
      for(int iband=0;iband<nband;++iband){
        if(verbose_opt[0]>=1)
          std::cout << "scaling for band" << iband << std::endl;
        offset[ibag][iband]=(offset_opt.size()==1)?offset_opt[0]:offset_opt[iband];
        scale[ibag][iband]=(scale_opt.size()==1)?scale_opt[0]:scale_opt[iband];
        //search for min and maximum
        if(scale[ibag][iband]<=0){
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
      //recode vreclass to ordered vector, starting from 0 to nreclass
      vcode.clear();
      if(verbose_opt[0]>=1){
        std::cout << "before recoding: " << std::endl;
        for(int iclass = 0; iclass < vreclass.size(); iclass++)
          std::cout << " " << vreclass[iclass];
        std::cout << std::endl; 
      }
      vector<int> vord=vreclass;//ordered vector, starting from 0 to nreclass
      int iclass=0;
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
    
      if(priors_opt.size()==1){//default: equal priors for each class
        priors.resize(nclass);
        for(int iclass=0;iclass<nclass;++iclass)
          priors[iclass]=1.0/nclass;
      }
      assert(priors_opt.size()==1||priors_opt.size()==nclass);
    
      if(verbose_opt[0]>=1){
        std::cout << "number of bands: " << nband << std::endl;
        std::cout << "number of classes: " << nclass << std::endl;
        std::cout << "priors:";
        for(int iclass=0;iclass<nclass;++iclass)
          std::cout << " " << priors[iclass];
        std::cout << std::endl;
      }
    }//if(!ibag)

    //Calculate features of trainig set
    vector< Vector2d<float> > trainingFeatures(nclass);
    for(int iclass=0;iclass<nclass;++iclass){
      int nctraining=0;
      if(verbose_opt[0]>=1)
        std::cout << "calculating features for class " << iclass << std::endl;
      if(random)
        srand(time(NULL));
      nctraining=(bagSize_opt[0]<100)? trainingPixels[iclass].size()/100.0*bagSize_opt[0] : trainingPixels[iclass].size();//bagSize_opt[0] given in % of training size
      if(nctraining<=0)
        nctraining=1;
      assert(nctraining<=trainingPixels[iclass].size());
      int index=0;
      if(bagSize_opt[0]<100)
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
    for(int iclass=0;iclass<nclass;++iclass){
      if(verbose_opt[0]>=1)
        std::cout << "training sample size for class " << vcode[iclass] << ": " << trainingFeatures[iclass].size() << std::endl;
      ntraining+=trainingFeatures[iclass].size();
    }
    // vector<struct svm_problem> prob(ibag);
    // vector<struct svm_node *> x_space(ibag);
    prob[ibag].l=ntraining;
    prob[ibag].y = Malloc(double,prob[ibag].l);
    prob[ibag].x = Malloc(struct svm_node *,prob[ibag].l);
    x_space[ibag] = Malloc(struct svm_node,(nFeatures+1)*ntraining);
    unsigned long int spaceIndex=0;
    int lIndex=0;
    for(int iclass=0;iclass<nclass;++iclass){
      for(int isample=0;isample<trainingFeatures[iclass].size();++isample){
        // //test
        // std::cout << iclass;
        prob[ibag].x[lIndex]=&(x_space[ibag][spaceIndex]);
        for(int ifeature=0;ifeature<nFeatures;++ifeature){
          x_space[ibag][spaceIndex].index=ifeature+1;
          x_space[ibag][spaceIndex].value=trainingFeatures[iclass][isample][ifeature];
          // //test
          // std::cout << " " << x_space[ibag][spaceIndex].index << ":" << x_space[ibag][spaceIndex].value;
          ++spaceIndex;
        }
        x_space[ibag][spaceIndex++].index=-1;
        prob[ibag].y[lIndex]=iclass;
        ++lIndex;
      }
    }
    assert(lIndex==prob[ibag].l);

    //set SVM parameters through command line options
    param[ibag].svm_type = svm_type_opt[0];
    param[ibag].kernel_type = kernel_type_opt[0];
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

    if(verbose_opt[0])
      std::cout << "checking parameters" << std::endl;
    svm_check_parameter(&prob[ibag],&param[ibag]);
    if(verbose_opt[0])
      std::cout << "parameters ok, training" << std::endl;
    svm[ibag]=svm_train(&prob[ibag],&param[ibag]);
    
    if(verbose_opt[0])
      std::cout << "SVM is now trained" << std::endl;
    if(cv_opt[0]>0){
      std::cout << "Confusion matrix" << std::endl;
      ConfusionMatrix cm(nclass);
      // for(int iclass=0;iclass<nclass;++iclass)
      //   cm.pushBackClassName(type2string(iclass));
      double *target = Malloc(double,prob[ibag].l);
      svm_cross_validation(&prob[ibag],&param[ibag],cv_opt[0],target);
      assert(param[ibag].svm_type != EPSILON_SVR&&param[ibag].svm_type != NU_SVR);//only for regression
      int total_correct=0;
      for(int i=0;i<prob[ibag].l;i++)
        cm.incrementResult(cm.getClass(prob[ibag].y[i]),cm.getClass(target[i]),1);
      assert(cm.nReference());
      std::cout << cm << std::endl;
      std::cout << "Kappa: " << cm.kappa() << std::endl;
      double se95_oa=0;
      double doa=0;
      doa=cm.oa_pct(&se95_oa);
      std::cout << "Overall Accuracy: " << doa << " (" << se95_oa << ")"  << std::endl;
      free(target);
    }

    // *NOTE* Because svm_model contains pointers to svm_problem, you can
    // not free the memory used by svm_problem if you are still using the
    // svm_model produced by svm_train(). 

    // free(prob.y);
    // free(prob.x);
    // free(x_space);
    // svm_destroy_param(&param);
  }//for ibag

  //--------------------------------- end of training -----------------------------------
  if(!output_opt.size())
    exit(0);


  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  float progress=0;
  if(!verbose_opt[0])
    pfnProgress(progress,pszMessage,pProgressArg);
  //-------------------------------- open image file ------------------------------------
  if(input_opt[0].find(".shp")==string::npos){
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
    string imageType=testImage.getImageType();
    if(oformat_opt.size())//default
      imageType=oformat_opt[0];
    try{
      assert(output_opt.size());
      if(verbose_opt[0]>=1)
        std::cout << "opening class image for writing output " << output_opt[0] << std::endl;
      if(classBag_opt.size()){
        classImageBag.open(output_opt[0],ncol,nrow,nbag,GDT_Byte,imageType,option_opt);
        classImageBag.copyGeoTransform(testImage);
        classImageBag.setProjection(testImage.getProjection());
      }
      classImageOut.open(output_opt[0],ncol,nrow,1,GDT_Byte,imageType,option_opt);
      classImageOut.copyGeoTransform(testImage);
      classImageOut.setProjection(testImage.getProjection());
      if(colorTable_opt.size())
        classImageOut.setColorTable(colorTable_opt[0],0);
      if(prob_opt.size()){
        probImage.open(prob_opt[0],ncol,nrow,nreclass,GDT_Byte,imageType,option_opt);
        probImage.copyGeoTransform(testImage);
        probImage.setProjection(testImage.getProjection());
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
      Vector2d<float> hpixel(ncol,nband);
      // Vector2d<float> fpixel(ncol);
      Vector2d<float> prOut(nreclass,ncol);//posterior prob for each reclass
      Vector2d<char> classBag;//classified line for writing to image file
      if(classBag_opt.size())
        classBag.resize(nbag,ncol);
      //read all bands of all pixels in this line in hline
      for(int iband=start_opt[0];iband<start_opt[0]+nband;++iband){
        if(verbose_opt[0]==2)
          std::cout << "reading band " << iband << std::endl;
        assert(iband>=0);
        assert(iband<testImage.nrOfBand());
        try{
          testImage.readData(buffer,GDT_Float32,iline,iband);
        }
        catch(string theError){
          cerr << "Error reading " << input_opt[0] << ": " << theError << std::endl;
          exit(3);
        }
        catch(...){
          cerr << "error catched" << std::endl;
          exit(3);
        }
        for(int icol=0;icol<ncol;++icol)
          hpixel[icol][iband-start_opt[0]]=buffer[icol];
      }

      assert(nband==hpixel[0].size());
      if(verbose_opt[0]==2)
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
    
      //process per pixel
      for(int icol=0;icol<ncol;++icol){

        bool masked=false;
        if(!lineMask.empty()){
          short theMask=0;
          for(short ivalue=0;ivalue<maskValue_opt.size();++ivalue){
            if(maskValue_opt[ivalue]>=0){//values set in maskValue_opt are invalid
              if(lineMask[icol]==maskValue_opt[ivalue]){
                theMask=(flag_opt.size()==maskValue_opt.size())? flag_opt[ivalue] : flag_opt[0];// lineMask[icol];
                masked=true;
                break;
              }
            }
            else{//only values set in maskValue_opt are valid
              if(lineMask[icol]!=-maskValue_opt[ivalue]){
                  theMask=(flag_opt.size()==maskValue_opt.size())? flag_opt[ivalue] : flag_opt[0];// lineMask[icol];
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
        for(int iband=0;iband<nband;++iband){
          if(hpixel[icol][iband]){
            valid=true;
            break;
          }
        }
        if(!valid){
          if(classBag_opt.size())
            for(int ibag=0;ibag<nbag;++ibag)
              classBag[ibag][icol]=flag_opt[0];
          classOut[icol]=flag_opt[0];
          continue;//next column
        }
        for(int iclass=0;iclass<nreclass;++iclass)
          prOut[iclass][icol]=0;
        //----------------------------------- classification -------------------
        for(int ibag=0;ibag<nbag;++ibag){
          //calculate image features
          // fpixel[icol].clear();
          // for(int iband=0;iband<nband;++iband)
          //   fpixel[icol].push_back((hpixel[icol][iband]-offset[ibag][iband])/scale[ibag][iband]);
          vector<double> result(nclass);
          // result=net[ibag].run(fpixel[icol]);
          struct svm_node *x;
          // x = (struct svm_node *) malloc((fpixel[icol].size()+1)*sizeof(struct svm_node));
          x = (struct svm_node *) malloc((nband+1)*sizeof(struct svm_node));
          // for(int i=0;i<fpixel[icol].size();++i){
          for(int iband=0;iband<nband;++iband){
            x[iband].index=iband+1;
            // x[i].value=fpixel[icol][i];
            x[iband].value=(hpixel[icol][iband]-offset[ibag][iband])/scale[ibag][iband];
          }
          // x[fpixel[icol].size()].index=-1;//to end svm feature vector
          x[nband].index=-1;//to end svm feature vector
          double predict_label=0;
          vector<float> pValues(nclass);
          vector<float> prValues(nreclass);
          vector<float> priorsReclass(nreclass);
          float maxP=0;
          if(!aggreg_opt[0]){
            predict_label = svm_predict(svm[ibag],x);
            for(int iclass=0;iclass<nclass;++iclass){
              if(iclass==static_cast<int>(predict_label))
                result[iclass]=1;
              else
                result[iclass]=0;
            }
          }
          else{
            assert(svm_check_probability_model(svm[ibag]));
            predict_label = svm_predict_probability(svm[ibag],x,&(result[0]));
          }
          for(int iclass=0;iclass<nclass;++iclass){
            float pv=result[iclass];
            assert(pv>=0);
            assert(pv<=1);
            pv*=priors[iclass];
            pValues[iclass]=pv;
          }
          float normReclass=0;
          for(int iclass=0;iclass<nreclass;++iclass){
            prValues[iclass]=0;
            priorsReclass[iclass]=0;
            float maxPaggreg=0;
            for(int ic=0;ic<nclass;++ic){
              if(vreclass[ic]==iclass){
                priorsReclass[iclass]+=priors[ic];
                switch(aggreg_opt[0]){
                default:
                case(1)://sum rule (sum posterior probabilities of aggregated individual classes)
                  prValues[iclass]+=pValues[ic];
                break;
                case(0):
                case(2)://max rule (look for maximum post probability of aggregated individual classes)
                  if(pValues[ic]>maxPaggreg){
                    maxPaggreg=pValues[ic];
                    prValues[iclass]=maxPaggreg;
                  }
                  break;
                }
              }
            }
          }
          for(int iclass=0;iclass<nreclass;++iclass)
            normReclass+=prValues[iclass];
        
          //calculate posterior prob of bag 
          if(classBag_opt.size()){
            //search for max prob within bag
            maxP=0;
            classBag[ibag][icol]=0;
          }
          for(int iclass=0;iclass<nreclass;++iclass){
            float prv=prValues[iclass];
            prv/=normReclass;
            //           prv*=100.0;
            prValues[iclass]=prv;
            switch(comb_opt[0]){
            default:
            case(0)://sum rule
              prOut[iclass][icol]+=prValues[iclass]+static_cast<float>(1.0-nbag)/nbag*priorsReclass[iclass];//add probabilities for each bag
              break;
            case(1)://product rule
              prOut[iclass][icol]*=pow(priorsReclass[iclass],static_cast<float>(1.0-nbag)/nbag)*prValues[iclass];//add probabilities for each bag
              break;
            case(2)://max rule
              if(prValues[iclass]>prOut[iclass][icol])
                prOut[iclass][icol]=prValues[iclass];
              break;
            }
            if(classBag_opt.size()){
              //search for max prob within bag
              if(prValues[iclass]>maxP){
                maxP=prValues[iclass];
                classBag[ibag][icol]=vcode[iclass];
              }
            }
          }
          free(x);
        }//ibag

        //search for max class prob
        float maxBag=0;
        float normBag=0;
        for(int iclass=0;iclass<nreclass;++iclass){
          if(prOut[iclass][icol]>maxBag){
            maxBag=prOut[iclass][icol];
            classOut[icol]=vcode[iclass];
          }
          normBag+=prOut[iclass][icol];
        }
        //normalize prOut and convert to percentage
        if(prob_opt.size()){
          for(int iclass=0;iclass<nreclass;++iclass){
            float prv=prOut[iclass][icol];
            prv/=normBag;
            prv*=100.0;
            prOut[iclass][icol]=static_cast<short>(prv+0.5);
          }
        }
      }//icol
      //----------------------------------- write output ------------------------------------------
      if(classBag_opt.size())
        for(int ibag=0;ibag<nbag;++ibag)
          classImageBag.writeData(classBag[ibag],GDT_Byte,iline,ibag);
      if(prob_opt.size()){
        for(int iclass=0;iclass<nreclass;++iclass)
          probImage.writeData(prOut[iclass],GDT_Float32,iline,iclass);
      }
      classImageOut.writeData(classOut,GDT_Byte,iline);
      if(!verbose_opt[0]){
        progress=static_cast<float>(iline+1.0)/classImageOut.nrOfRow();
        pfnProgress(progress,pszMessage,pProgressArg);
      }
    }
    testImage.close();
    if(prob_opt.size())
      probImage.close();
    if(classBag_opt.size())
      classImageBag.close();
    classImageOut.close();
  }
  else{//classify shape file
    for(int ivalidation=0;ivalidation<input_opt.size();++ivalidation){
      assert(output_opt.size()==input_opt.size());
      if(verbose_opt[0])
        std::cout << "opening img reader " << input_opt[ivalidation] << std::endl;
      ImgReaderOgr imgReaderOgr(input_opt[ivalidation]);
      if(verbose_opt[0])
        std::cout << "opening img writer and copying fields from img reader" << output_opt[ivalidation] << std::endl;
      ImgWriterOgr imgWriterOgr(output_opt[ivalidation],imgReaderOgr,false);
      if(verbose_opt[0])
        std::cout << "creating field class" << std::endl;
      imgWriterOgr.createField("class",OFTInteger);
      OGRFeature *poFeature;
      unsigned int ifeature=0;
      unsigned int nFeatures=imgReaderOgr.getFeatureCount();
      while( (poFeature = imgReaderOgr.getLayer()->GetNextFeature()) != NULL ){
        if(verbose_opt[0]>1)
          std::cout << "feature " << ifeature << std::endl;
        if( poFeature == NULL )
          break;
        OGRFeature *poDstFeature = NULL;
        poDstFeature=imgWriterOgr.createFeature();
        if( poDstFeature->SetFrom( poFeature, TRUE ) != OGRERR_NONE ){
          CPLError( CE_Failure, CPLE_AppDefined,
                    "Unable to translate feature %d from layer %s.\n",
                    poFeature->GetFID(), imgWriterOgr.getLayerName().c_str() );
          OGRFeature::DestroyFeature( poFeature );
          OGRFeature::DestroyFeature( poDstFeature );
        }
        vector<float> validationPixel;
        vector<float> validationFeature;
        
        imgReaderOgr.readData(validationPixel,OFTReal,fields,poFeature);
        OGRFeature::DestroyFeature( poFeature );
//         assert(validationPixel.size()>=start_opt[0]+nband);
        assert(validationPixel.size()==nband);
        vector<float> prOut(nreclass);//posterior prob for each reclass
        for(int iclass=0;iclass<nreclass;++iclass)
          prOut[iclass]=0;
        for(int ibag=0;ibag<nbag;++ibag){
//           for(int iband=start_opt[0];iband<start_opt[0]+nband;++iband){
          for(int iband=0;iband<nband;++iband){
//             validationFeature.push_back((validationPixel[iband]-offset[ibag][iband-start_opt[0]])/scale[ibag][iband-start_opt[0]]);
            validationFeature.push_back((validationPixel[iband]-offset[ibag][iband])/scale[ibag][iband]);
            if(verbose_opt[0]==2)
              std::cout << " " << validationFeature.back();
          }
          if(verbose_opt[0]==2)
            std::cout << std::endl;
          vector<double> result(nclass);

          // result=net[ibag].run(validationFeature);
          struct svm_node *x;
          x = (struct svm_node *) malloc((validationFeature.size()+1)*sizeof(struct svm_node));
          for(int i=0;i<validationFeature.size();++i){
            x[i].index=i+1;
            x[i].value=validationFeature[i];
          }
          x[validationFeature.size()].index=-1;//to end svm feature vector
          double predict_label=0;
          vector<float> pValues(nclass);
          vector<float> prValues(nreclass);
          vector<float> priorsReclass(nreclass);
          float maxP=0;
          if(!aggreg_opt[0]){
            predict_label = svm_predict(svm[ibag],x);
            for(int iclass=0;iclass<nclass;++iclass){
              if(predict_label==static_cast<int>(predict_label))
                result[iclass]=1;
              else
                result[iclass]=0;
            }
          }
          else{
            assert(svm_check_probability_model(svm[ibag]));
            predict_label = svm_predict_probability(svm[ibag],x,&(result[0]));
          }
          // int maxClass=0;
          for(int iclass=0;iclass<nclass;++iclass){
            float pv=(result[iclass]+1.0)/2.0;//bring back to scale [0,1]
            pv*=priors[iclass];
            pValues[iclass]=pv;
          }
          float normReclass=0;
          for(int iclass=0;iclass<nreclass;++iclass){
            prValues[iclass]=0;
            priorsReclass[iclass]=0;
            float maxPaggreg=0;
            for(int ic=0;ic<nclass;++ic){
              if(vreclass[ic]==iclass){
                priorsReclass[iclass]+=priors[ic];
                switch(aggreg_opt[0]){
                default:
                case(0)://sum rule (sum posterior probabilities of aggregated individual classes)
                  prValues[iclass]+=pValues[ic];
                  break;
                case(1)://max rule (look for maximum post probability of aggregated individual classes)
                  if(pValues[ic]>maxPaggreg){
                    maxPaggreg=pValues[ic];
                    prValues[iclass]=maxPaggreg;
                  }
                  break;
                }
              }
            }
          }
          for(int iclass=0;iclass<nreclass;++iclass)
            normReclass+=prValues[iclass];
          //calculate posterior prob of bag 
          for(int iclass=0;iclass<nreclass;++iclass){
            float prv=prValues[iclass];
            prv/=normReclass;
            //           prv*=100.0;
            prValues[iclass]=prv;
            switch(comb_opt[0]){
            default:
            case(0)://sum rule
              prOut[iclass]+=prValues[iclass]+static_cast<float>(1.0-nbag)/nbag*priorsReclass[iclass];//add probabilities for each bag
              break;
            case(1)://product rule
              prOut[iclass]*=pow(priorsReclass[iclass],static_cast<float>(1.0-nbag)/nbag)*prValues[iclass];//add probabilities for each bag
              break;
            case(2)://max rule
              if(prValues[iclass]>prOut[iclass])
                prOut[iclass]=prValues[iclass];
              break;
            }
          }
          free(x);
        }//for ibag
        //search for max class prob
        float maxBag=0;
        float normBag=0;
        char classOut=0;
        for(int iclass=0;iclass<nreclass;++iclass){
          if(prOut[iclass]>maxBag){
            maxBag=prOut[iclass];
            classOut=vcode[iclass];
          }
          normBag+=prOut[iclass];
        }
        //normalize prOut and convert to percentage
        for(int iclass=0;iclass<nreclass;++iclass){
          float prv=prOut[iclass];
          prv/=normBag;
          prv*=100.0;
          prOut[iclass]=static_cast<short>(prv+0.5);
        }
        poDstFeature->SetField("class",classOut);
        poDstFeature->SetFID( poFeature->GetFID() );
        CPLErrorReset();
        if(imgWriterOgr.createFeature( poDstFeature ) != OGRERR_NONE){
          CPLError( CE_Failure, CPLE_AppDefined,
                    "Unable to translate feature %d from layer %s.\n",
                    poFeature->GetFID(), imgWriterOgr.getLayerName().c_str() );
          OGRFeature::DestroyFeature( poDstFeature );
          OGRFeature::DestroyFeature( poDstFeature );
        }
        ++ifeature;
        if(!verbose_opt[0]){
          progress=static_cast<float>(ifeature+1.0)/nFeatures;
          pfnProgress(progress,pszMessage,pProgressArg);
        }
      }
      imgReaderOgr.close();
      imgWriterOgr.close();
    }
  }

  for(int ibag=0;ibag<nbag;++ibag){
    svm_destroy_param[ibag](&param[ibag]);
    free(prob[ibag].y);
    free(prob[ibag].x);
    free(x_space[ibag]);
    svm_free_and_destroy_model(&(svm[ibag]));
  }
  return 0;
}
