/**********************************************************************
pkclassify_nn.cc: classify raster image using Artificial Neural Network
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
#include "pkclassify_nn.h"
#include <vector>
#include <map>
#include <algorithm>
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"
#include "imageclasses/ImgReaderOgr.h"
#include "imageclasses/ImgWriterOgr.h"
#include "base/Optionpk.h"
#include "algorithms/ConfusionMatrix.h"
#include "floatfann.h"
#include "myfann_cpp.h"

void reclass(const vector<float>& result, const vector<int>& vreclass, const vector<double>& priors, unsigned short aggregation, vector<float>& theResultReclass);

void reclass(const vector<float>& result, const vector<int>& vreclass, const vector<double>& priors, unsigned short aggregation, vector<float>& theResultReclass){
  unsigned int nclass=result.size();
  assert(priors.size()==nclass);
  assert(theResultReclass.size()>1);//must have size nreclass!
  unsigned int nreclass=theResultReclass.size();
  vector<float> pValues(nclass);
  float normReclass=0;
  for(int iclass=0;iclass<nclass;++iclass){
    float pv=(result[iclass]+1.0)/2.0;//bring back to scale [0,1]
    assert(pv>=0);
    assert(pv<=1);
    pv*=priors[iclass];
    pValues[iclass]=pv;
  }
  for(int iclass=0;iclass<nreclass;++iclass){
    theResultReclass[iclass]=0;
    float maxPaggreg=0;
    for(int ic=0;ic<nclass;++ic){
      if(vreclass[ic]==iclass){
	switch(aggregation){
	default:
	case(1)://sum rule (sum posterior probabilities of aggregated individual classes)
	  theResultReclass[iclass]+=pValues[ic];
	break;
	case(0):
	case(2)://max rule (look for maximum post probability of aggregated individual classes)
	  if(pValues[ic]>maxPaggreg){
	    maxPaggreg=pValues[ic];
	    theResultReclass[iclass]=maxPaggreg;
	  }
	break;
	}
      }
    }
    normReclass+=theResultReclass[iclass];
  }
  for(int iclass=0;iclass<nreclass;++iclass){
    float prv=theResultReclass[iclass];
    prv/=normReclass;
    theResultReclass[iclass]=prv;
  }
}

int main(int argc, char *argv[])
{
  map<short,int> reclassMap;
  vector<int> vreclass;
  vector<double> priors;
  vector<double> priorsReclass;
  
  //--------------------------- command line options ------------------------------------
  Optionpk<string> input_opt("i", "input", "input image"); 
  Optionpk<string> training_opt("t", "training", "training shape file. A single shape file contains all training features (must be set as: B0, B1, B2,...) for all classes (class numbers identified by label option). Use multiple training files for bootstrap aggregation (alternative to the bag and bsize options, where a random subset is taken from a single training file)"); 
  Optionpk<string> label_opt("label", "label", "identifier for class label in training shape file.","label"); 
  Optionpk<unsigned short> reclass_opt("rc", "rc", "reclass code (e.g. --rc=12 --rc=23 to reclass first two classes to 12 and 23 resp.)");
  Optionpk<unsigned int> balance_opt("bal", "balance", "balance the input data to this number of samples for each class", 0);
  Optionpk<int> minSize_opt("min", "min", "if number of training pixels is less then min, do not take this class into account (0: consider all classes)", 0);
  Optionpk<double> start_opt("s", "start", "start band sequence number (set to 0)",0); 
  Optionpk<double> end_opt("e", "end", "end band sequence number (set to 0 for all bands)", 0); 
  Optionpk<short> band_opt("b", "band", "band index (starting from 0, either use band option or use start to end)");
  Optionpk<double> offset_opt("\0", "offset", "offset value for each spectral band input features: refl[band]=(DN[band]-offset[band])/scale[band]", 0.0);
  Optionpk<double> scale_opt("\0", "scale", "scale value for each spectral band input features: refl=(DN[band]-offset[band])/scale[band] (use 0 if scale min and max in each band to -1.0 and 1.0)", 0.0);
  Optionpk<unsigned short> aggreg_opt("a", "aggreg", "how to combine aggregated classifiers, see also rc option (1: sum rule, 2: max rule).",1);
  Optionpk<double> priors_opt("p", "prior", "prior probabilities for each class (e.g., -p 0.3 -p 0.3 -p 0.2 )", 0.0); 
  Optionpk<unsigned short> cv_opt("cv", "cv", "n-fold cross validation mode",0);
  Optionpk<unsigned int> nneuron_opt("\0", "nneuron", "number of neurons in hidden layers in neural network (multiple hidden layers are set by defining multiple number of neurons: -n 15 -n 1, default is one hidden layer with 5 neurons)", 5); 
  Optionpk<float> connection_opt("\0", "connection", "connection reate (default: 1.0 for a fully connected network)", 1.0); 
  Optionpk<float> weights_opt("w", "weights", "weights for neural network. Apply to fully connected network only, starting from first input neuron to last output neuron, including the bias neurons (last neuron in each but last layer)", 0.0); 
  Optionpk<float> learning_opt("l", "learning", "learning rate (default: 0.7)", 0.7); 
  Optionpk<unsigned int> maxit_opt("\0", "maxit", "number of maximum iterations (epoch) (default: 500)", 500); 
  Optionpk<unsigned short> comb_opt("comb", "comb", "how to combine bootstrap aggregation classifiers (0: sum rule, 1: product rule, 2: max rule). Also used to aggregate classes with rc option. Default is sum rule (0)",0); 
  Optionpk<unsigned short> bag_opt("bag", "bag", "Number of bootstrap aggregations (default is no bagging: 1)", 1);
  Optionpk<int> bagSize_opt("bs", "bsize", "Percentage of features used from available training features for each bootstrap aggregation (one size for all classes, or a different size for each class respectively", 100);
  Optionpk<string> classBag_opt("cb", "classbag", "output for each individual bootstrap aggregation (default is blank)"); 
  Optionpk<string> mask_opt("m", "mask", "mask image (see also mvalue option (default is no mask)"); 
  Optionpk<short> maskValue_opt("mv", "mvalue", "mask value(s) not to consider for classification (use negative values if only these values should be taken into account). Values will be taken over in classification image. Default is 0", 0);
  Optionpk<unsigned short> flag_opt("f", "flag", "flag to put where image is invalid. Default is 0", 0);
  Optionpk<string> output_opt("o", "output", "output classification image"); 
  Optionpk<string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image");
  Optionpk<string> option_opt("co", "co", "options: NAME=VALUE [-co COMPRESS=LZW] [-co INTERLEAVE=BAND]");
  Optionpk<string> colorTable_opt("ct", "ct", "colour table in ascii format having 5 columns: id R G B ALFA (0: transparent, 255: solid)"); 
  Optionpk<string> prob_opt("\0", "prob", "probability image. Default is no probability image"); 
  Optionpk<short> verbose_opt("v", "verbose", "set to: 0 (results only), 1 (confusion matrix), 2 (debug)",0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    training_opt.retrieveOption(argc,argv);
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
    priors_opt.retrieveOption(argc,argv);
    cv_opt.retrieveOption(argc,argv);
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
    maskValue_opt.retrieveOption(argc,argv);
    flag_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    otype_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    prob_opt.retrieveOption(argc,argv);
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

  if(verbose_opt[0]>=1){
    if(input_opt.size())
      cout << "image filename: " << input_opt[0] << endl;
    if(mask_opt.size())
      cout << "mask filename: " << mask_opt[0] << endl;
    if(training_opt.size()){
      cout << "training shape file: " << endl;
      for(int ifile=0;ifile<training_opt.size();++ifile)
        cout << training_opt[ifile] << endl;
    }
    else
      cerr << "no training file set!" << endl;
    cout << "verbose: " << verbose_opt[0] << endl;
  }
  unsigned short nbag=(training_opt.size()>1)?training_opt.size():bag_opt[0];
  if(verbose_opt[0]>=1)
    cout << "number of bootstrap aggregations: " << nbag << endl;
  
  unsigned int totalSamples=0;
  int nreclass=0;
  vector<int> vcode;//unique class codes in recode string
  vector<FANN::neural_net> net(nbag);//the neural network

  unsigned int nclass=0;
  int nband=0;
  int startBand=2;//first two bands represent X and Y pos

  if(reclass_opt.size()){
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

  //sort bands
  if(band_opt.size())
    std::sort(band_opt.begin(),band_opt.end());
  //----------------------------------- Training -------------------------------
  ConfusionMatrix cm;
  vector< vector<double> > offset(nbag);
  vector< vector<double> > scale(nbag);
  map<string,Vector2d<float> > trainingMap;
  vector< Vector2d<float> > trainingPixels;//[class][sample][band]
  vector<string> fields;
  for(int ibag=0;ibag<nbag;++ibag){
    //organize training data
    if(ibag<training_opt.size()){//if bag contains new training pixels
      trainingMap.clear();
      trainingPixels.clear();
      if(verbose_opt[0]>=1)
        cout << "reading imageShape file " << training_opt[0] << endl;
      try{
        if(band_opt.size())
          totalSamples=readDataImageShape(training_opt[ibag],trainingMap,fields,band_opt,label_opt[0],verbose_opt[0]);
        else
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
      //delete class 0 ?
      // if(verbose_opt[0]>=1)
      //   std::cout << "erasing class 0 from training set (" << trainingMap[0].size() << " from " << totalSamples << ") samples" << std::endl;
      // totalSamples-=trainingMap[0].size();
      // trainingMap.erase(0);
      //convert map to vector
      short iclass=0;
      if(reclass_opt.empty()){//no reclass option, read classes from shape
        reclassMap.clear();
        vreclass.clear();
      }
      if(verbose_opt[0]>1)
        std::cout << "training pixels: " << std::endl;
      map<string,Vector2d<float> >::iterator mapit=trainingMap.begin();
      while(mapit!=trainingMap.end()){
        //delete small classes
        if((mapit->second).size()<minSize_opt[0]){
          trainingMap.erase(mapit);
          continue;
          //todo: beware of reclass option: delete this reclass if no samples are left in this classes!!
        }
        if(reclass_opt.empty()){//no reclass option, read classes from shape
          vreclass.push_back(iclass);
        }
        trainingPixels.push_back(mapit->second);
        if(verbose_opt[0]>1)
          std::cout << mapit->first << ": " << (mapit->second).size() << " samples" << std::endl;
        ++iclass;
        ++mapit;
      }
      if(!ibag){
        nclass=trainingPixels.size();
        nband=(training_opt.size())?trainingPixels[0][0].size()-2:trainingPixels[0][0].size();//X and Y
      }
      else{
        assert(nclass==trainingPixels.size());
        assert(nband==(training_opt.size())?trainingPixels[0][0].size()-2:trainingPixels[0][0].size());
      }
      // assert(reclassMap.size()==nclass);
      assert(vreclass.size()==nclass);

      //do not remove outliers here: could easily be obtained through ogr2ogr -where 'B2<110' output.shp input.shp
      //balance training data
      if(balance_opt[0]>0){
        while(balance_opt.size()<nclass)
          balance_opt.push_back(balance_opt.back());
        if(random)
          srand(time(NULL));
        totalSamples=0;
        for(int iclass=0;iclass<nclass;++iclass){
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
      Histogram hist;
      for(int iband=0;iband<nband;++iband){
        if(verbose_opt[0]>=1)
          cout << "scaling for band" << iband << endl;
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
            cout << "Extreme image values for band " << iband << ": [" << theMin << "," << theMax << "]" << endl;
            cout << "Using offset, scale: " << offset[ibag][iband] << ", " << scale[ibag][iband] << endl;
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
        cout << "before recoding: " << endl;
        for(int iclass = 0; iclass < vreclass.size(); iclass++)
          cout << " " << vreclass[iclass];
        cout << endl; 
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
        cout << "recode values: " << endl;
        for(int icode=0;icode<vcode.size();++icode)
          cout << vcode[icode] << " ";
        cout << endl;
      }
      vreclass=vord;
      if(verbose_opt[0]>=1){
        cout << "after recoding: " << endl;
        for(int iclass = 0; iclass < vord.size(); iclass++)
          cout << " " << vord[iclass];
        cout << endl; 
      }
      
      vector<int> vuniqueclass=vreclass;
      //remove duplicate elements from vuniqueclass
      sort( vuniqueclass.begin(), vuniqueclass.end() );
      vuniqueclass.erase( unique( vuniqueclass.begin(), vuniqueclass.end() ), vuniqueclass.end() );
      nreclass=vuniqueclass.size();
      if(verbose_opt[0]>=1){
        cout << "unique classes: " << endl;
        for(int iclass = 0; iclass < vuniqueclass.size(); iclass++)
          cout << " " << vuniqueclass[iclass];
        cout << endl; 
        cout << "number of reclasses: " << nreclass << endl;
      }
    
      if(priors_opt.size()==1){//default: equal priors for each class
        priors.resize(nclass);
        for(int iclass=0;iclass<nclass;++iclass)
          priors[iclass]=1.0/nclass;
      }
      assert(priors_opt.size()==1||priors_opt.size()==nclass);
    
      //set priors
      priorsReclass.resize(nreclass);
      for(int iclass=0;iclass<nreclass;++iclass){
	priorsReclass[iclass]=0;
	for(int ic=0;ic<nclass;++ic){
	  if(vreclass[ic]==iclass)
	    priorsReclass[iclass]+=priors[ic];
	}
      }
      //set bagsize for each class if not done already via command line
      while(bagSize_opt.size()<nclass)
        bagSize_opt.push_back(bagSize_opt.back());

      if(verbose_opt[0]>=1){
        cout << "number of bands: " << nband << endl;
        cout << "number of classes: " << nclass << endl;
        cout << "priors:";
        for(int iclass=0;iclass<nclass;++iclass)
          cout << " " << priors[iclass];
        cout << endl;
      }
    }

    //Calculate features of trainig set
    vector< Vector2d<float> > trainingFeatures(nclass);
    for(int iclass=0;iclass<nclass;++iclass){
      int nctraining=0;
      if(verbose_opt[0]>=1)
        cout << "calculating features for class " << iclass << endl;
      if(random)
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
    for(int iclass=0;iclass<nclass;++iclass){
      //vcode has size nreclass???
      // if(verbose_opt[0]>=1)
      //   cout << "training sample size for class " << vcode[iclass] << ": " << trainingFeatures[iclass].size() << endl;
      ntraining+=trainingFeatures[iclass].size();
    }
    const unsigned int num_layers = nneuron_opt.size()+2;
    const float desired_error = 0.0003;
    const unsigned int iterations_between_reports = (verbose_opt[0])? maxit_opt[0]+1:0;
    if(verbose_opt[0]>=1){
      cout << "creating artificial neural network with " << nneuron_opt.size() << " hidden layer, having " << endl;
      for(int ilayer=0;ilayer<nneuron_opt.size();++ilayer)
        cout << nneuron_opt[ilayer] << " ";
      cout << "neurons" << endl;
    }
    switch(num_layers){
    case(3):
      net[ibag].create_sparse(connection_opt[0],num_layers, nFeatures, nneuron_opt[0], nclass);
      break;
    case(4):
      net[ibag].create_sparse(connection_opt[0],num_layers, nFeatures, nneuron_opt[0], nneuron_opt[1], nclass);
      break;
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
      
    if(cv_opt[0]){
      if(verbose_opt[0])
        std::cout << "cross validation" << std::endl;
      vector<unsigned short> referenceVector;
      vector<unsigned short> outputVector;
      float rmse=net[ibag].cross_validation(trainingFeatures,
                                            ntraining,
                                            cv_opt[0],
                                            maxit_opt[0],
                                            0,
                                            desired_error,
                                            referenceVector,
                                            outputVector,
                                            verbose_opt[0]);
      map<string,Vector2d<float> >::iterator mapit=trainingMap.begin();
      if(reclass_opt.empty()){
        while(mapit!=trainingMap.end()){
          cm.pushBackClassName(mapit->first);
          ++mapit;
        }
      }
      else{
        if(verbose_opt[0]>1)
          std::cout << "classes for confusion matrix: " << std::endl;
        for(int iclass=0;iclass<nreclass;++iclass){
          ostringstream os;
          os << vcode[iclass];
          if(verbose_opt[0]>1)
            std::cout << os.str() << " ";
          cm.pushBackClassName(os.str());
        }
        if(verbose_opt[0]>1)
          std::cout << std::endl;
      }
      assert(cm.size()==nreclass);

      for(int isample=0;isample<referenceVector.size();++isample)
        cm.incrementResult(cm.getClass(vreclass[referenceVector[isample]]),cm.getClass(vreclass[outputVector[isample]]),1.0/nbag);
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
      for(int i_connection=0;i_connection<net[ibag].get_total_connections();++i_connection)
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
      for(int i_connection=0;i_connection<net[ibag].get_total_connections();++i_connection)
        cout << "connection " << i_connection << ": " << convector[i_connection].weight << endl;

    }
  }//for ibag
  if(cv_opt[0]>0){
    assert(cm.nReference());
    std::cout << cm << std::endl;
    cout << "class #samples userAcc prodAcc" << endl;
    double se95_ua=0;
    double se95_pa=0;
    double se95_oa=0;
    double dua=0;
    double dpa=0;
    double doa=0;
    for(int iclass=0;iclass<cm.nClasses();++iclass){
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

  //-------------------------------- open image file ------------------------------------
  if(input_opt[0].find(".shp")==string::npos){
    ImgReaderGdal testImage;
    try{
      if(verbose_opt[0]>=1)
        cout << "opening image " << input_opt[0] << endl; 
      testImage.open(input_opt[0]);
    }
    catch(string error){
      cerr << error << endl;
      exit(2);
    }
    ImgReaderGdal maskReader;
    if(mask_opt.size()){
      try{
        if(verbose_opt[0]>=1)
          cout << "opening mask image file " << mask_opt[0] << endl;
        maskReader.open(mask_opt[0]);
        assert(maskReader.nrOfCol()==testImage.nrOfCol());
        assert(maskReader.nrOfRow()==testImage.nrOfRow());
      }
      catch(string error){
        cerr << error << endl;
        exit(2);
      }
      catch(...){
        cerr << "error catched" << endl;
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

      if(verbose_opt[0]>=1)
        cout << "opening class image for writing output " << output_opt[0] << endl;
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
      cerr << error << endl;
    }
  
    const char* pszMessage;
    void* pProgressArg=NULL;
    GDALProgressFunc pfnProgress=GDALTermProgress;
    float progress=0;
    if(!verbose_opt[0])
      pfnProgress(progress,pszMessage,pProgressArg);
    for(int iline=0;iline<nrow;++iline){
      vector<float> buffer(ncol);
      vector<short> lineMask;
      if(mask_opt.size())
        lineMask.resize(maskReader.nrOfCol());
      Vector2d<float> hpixel(ncol);
      Vector2d<float> fpixel(ncol);
      Vector2d<float> prOut(nreclass,ncol);//posterior prob for each reclass
      Vector2d<char> classBag;//classified line for writing to image file
      if(classBag_opt.size())
        classBag.resize(nbag,ncol);
      //read all bands of all pixels in this line in hline
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
      // for(int iband=start_opt[0];iband<start_opt[0]+nband;++iband){
      //   if(verbose_opt[0]==2)
      //     cout << "reading band " << iband << endl;
      //   assert(iband>=0);
      //   assert(iband<testImage.nrOfBand());
      //   try{
      //     testImage.readData(buffer,GDT_Float32,iline,iband);
      //   }
      //   catch(string theError){
      //     cerr << "Error reading " << input_opt[0] << ": " << theError << endl;
      //     exit(3);
      //   }
      //   catch(...){
      //     cerr << "error catched" << endl;
      //     exit(3);
      //   }
      //   for(int icol=0;icol<ncol;++icol)
      //     hpixel[icol][iband-start_opt[0]]=buffer[icol];
      // }

      assert(nband==hpixel[0].size());
      if(verbose_opt[0]==2)
        cout << "used bands: " << nband << endl;
      //read mask
      if(!lineMask.empty()){
        try{
          maskReader.readData(lineMask,GDT_Int16,iline);
        }
        catch(string theError){
          cerr << "Error reading " << mask_opt[0] << ": " << theError << endl;
          exit(3);
        }
        catch(...){
          cerr << "error catched" << endl;
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
          fpixel[icol].clear();
          for(int iband=0;iband<nband;++iband)
            fpixel[icol].push_back((hpixel[icol][iband]-offset[ibag][iband])/scale[ibag][iband]);
          vector<float> result(nclass);
          result=net[ibag].run(fpixel[icol]);
          int maxClass=0;
          float maxP=0;
          vector<float> pValues(nclass);
          vector<float> prValues(nreclass);

	  reclass(result,vreclass,priors,aggreg_opt[0],prValues);

          //calculate posterior prob of bag 
          if(classBag_opt.size()){
            //search for max prob within bag
            maxP=0;
            classBag[ibag][icol]=0;
          }
          for(int iclass=0;iclass<nreclass;++iclass){
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
        //       assert(classOut[icol]);
        //normalize prOut and convert to percentage
        for(int iclass=0;iclass<nreclass;++iclass){
          float prv=prOut[iclass][icol];
          prv/=normBag;
          prv*=100.0;
          prOut[iclass][icol]=static_cast<short>(prv+0.5);
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
    cm.clearResults();
    //notice that fields have already been set by readDataImageShape (taking into account appropriate bands)
    for(int ivalidation=0;ivalidation<input_opt.size();++ivalidation){
      assert(output_opt.size()==input_opt.size());
      if(verbose_opt[0])
        cout << "opening img reader " << input_opt[ivalidation] << endl;
      ImgReaderOgr imgReaderOgr(input_opt[ivalidation]);
      if(verbose_opt[0])
        cout << "opening img writer and copying fields from img reader" << output_opt[ivalidation] << endl;
      ImgWriterOgr imgWriterOgr(output_opt[ivalidation],imgReaderOgr,false);
      if(verbose_opt[0])
        cout << "creating field class" << endl;
      imgWriterOgr.createField("class",OFTInteger);
      OGRFeature *poFeature;
      unsigned int ifeature=0;
      while( (poFeature = imgReaderOgr.getLayer()->GetNextFeature()) != NULL ){
        if(verbose_opt[0]>1)
          cout << "feature " << ifeature << endl;
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
              cout << " " << validationFeature.back();
          }
          if(verbose_opt[0]==2)
            cout << endl;
          vector<float> result(nclass);
          result=net[ibag].run(validationFeature);
          int maxClass=0;
          float maxP=0;
          vector<float> pValues(nclass);
          vector<float> prValues(nreclass);

	  reclass(result,vreclass,priors,aggreg_opt[0],prValues);

          //calculate posterior prob of bag 
          for(int iclass=0;iclass<nreclass;++iclass){
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
        }//for ibag
        //search for max class prob
        float maxBag=0;
        float normBag=0;
        char classOut=0;
	short classOutIndex=0;
        for(int iclass=0;iclass<nreclass;++iclass){
          if(prOut[iclass]>maxBag){
            maxBag=prOut[iclass];
            classOut=vcode[iclass];
	    classOutIndex=iclass;
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
	int labelIndex=poDstFeature->GetFieldIndex(label_opt[0].c_str());
	if(labelIndex>=0){
	  string classRef=poDstFeature->GetFieldAsString(labelIndex);
	  // cm.incrementResult(classRef,cm.getClass(classOutIndex),1);
	}
        CPLErrorReset();
        if(imgWriterOgr.createFeature( poDstFeature ) != OGRERR_NONE){
          CPLError( CE_Failure, CPLE_AppDefined,
                    "Unable to translate feature %d from layer %s.\n",
                    poFeature->GetFID(), imgWriterOgr.getLayerName().c_str() );
          OGRFeature::DestroyFeature( poDstFeature );
          OGRFeature::DestroyFeature( poDstFeature );
        }
        ++ifeature;
      }
      imgReaderOgr.close();
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
      doa=cm.oa_pct(&se95_oa);
      std::cout << "Overall Accuracy: " << doa << " (" << se95_oa << ")"  << std::endl;
    }
  }
  return 0;
}
