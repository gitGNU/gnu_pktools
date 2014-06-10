/**********************************************************************
pkkalman.cc: program to kalman raster images: median, min/max, morphological, kalmaning
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
#include <sstream>
#include <vector>
#include "base/Optionpk.h"
#include "base/Vector2d.h"
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"
#include "algorithms/StatFactory.h"
#include "algorithms/ImgRegression.h"

using namespace std;
/*------------------
  Main procedure
  ----------------*/
int main(int argc,char **argv) {
  Optionpk<string> model_opt("mod","model","model input datasets, e.g., MODIS (use: -mod model1 -mod model2 etc.");
  Optionpk<string> observation_opt("obs","observation","observation input datasets, e.g., landsat (use: -obs obs1 -obs obs2 etc.");
  Optionpk<int> tobservation_opt("tobs","tobservation","time sequence of observation input (sequence must have exact same length as observation input)"); 
  Optionpk<string> output_opt("o", "output", "Output raster dataset");
  Optionpk<float> threshold_opt("t", "threshold", "threshold for selecting samples (randomly). Provide probability in percentage (>0) or absolute (<0).", 0);
  Optionpk<double> modnodata_opt("modnodata", "modnodata", "invalid value for model input", 0);
  Optionpk<double> obsnodata_opt("obsnodata", "obsnodata", "invalid value for observation input", 0);
  Optionpk<double> eps_opt("eps", "eps", "epsilon for non zero division", 0.00001);
  Optionpk<int> down_opt("down", "down", "Downsampling factor for reading model data to calculate regression", 9);
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image","GTiff",2);
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode when positive", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=model_opt.retrieveOption(argc,argv);
    observation_opt.retrieveOption(argc,argv);
    tobservation_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    threshold_opt.retrieveOption(argc,argv);
    modnodata_opt.retrieveOption(argc,argv);
    obsnodata_opt.retrieveOption(argc,argv);
    eps_opt.retrieveOption(argc,argv);
    down_opt.retrieveOption(argc,argv);
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

  try{
    ostringstream errorStream;
    if(model_opt.size()<2){
      errorStream << "Error: no model dataset selected, use option -mod" << endl;
      throw(errorStream.str());
    }
    if(observation_opt.size()<1){
      errorStream << "Error: no observation dataset selected, use option -obs" << endl;
      throw(errorStream.str());
    }
    if(model_opt.size()<observation_opt.size()){
      errorStream << "Error: sequence of models should be larger than observations" << endl;
      throw(errorStream.str());
    }
    if(tobservation_opt.size()!=observation_opt.size()){
      errorStream << "Error: time sequence for observation must match size of observation dataset" << endl;
      throw(errorStream.str());
    }
    if(output_opt.empty()){
      errorStream << "Error: suffix for output datasets is not provided" << endl;
      throw(errorStream.str());
    }
  }
  catch(string errorString){
    std::cout << errorString << std::endl;
    exit(1);
  }

  imgregression::ImgRegression imgreg;
  // vector<ImgReaderGdal> imgReaderModel(model_opt.size());
  // vector<ImgReaderGdal> imgReaderObs(observation_opt.size());
  ImgReaderGdal imgReaderModel1;
  ImgReaderGdal imgReaderModel2;
  ImgReaderGdal imgReaderObs;
  ImgWriterGdal imgWriterEst;

  imgReaderObs.open(observation_opt[0]);
  int ncol=imgReaderObs.nrOfCol();
  int nrow=imgReaderObs.nrOfRow();

  string imageType=imgReaderObs.getImageType();
  if(oformat_opt.size())//default
    imageType=oformat_opt[0];
  if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
    string theInterleave="INTERLEAVE=";
    theInterleave+=imgReaderObs.getInterleave();
    option_opt.push_back(theInterleave);
  }

  int obsindex=0;
  //forward step

  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;

  imgreg.setDown(down_opt[0]);
  imgreg.setThreshold(threshold_opt[0]);

  double c0obs=0;
  double c1obs=0;
  double errObs=0;
  for(int modindex=0;modindex<model_opt.size()-1;++modindex){
    string output;
    if(output_opt.size()==model_opt.size())
      output=output_opt[modindex+1];
    else{
      ostringstream outputstream;
      outputstream << output_opt[0] << "_" << modindex+1 << ".tif";
      output=outputstream.str();
    }
    
    //two band output band0=estimation, band1=uncertainty
    imgWriterEst.open(output,ncol,nrow,2,GDT_Float32,imageType,option_opt);

    //calculate regression between two subsequence model inputs
    imgReaderModel1.open(model_opt[modindex]);
    imgReaderModel2.open(model_opt[modindex+1]);
    //calculate regression
    //we could re-use the points from second image from last run, but
    //to keep it general, we must redo it (overlap might have changed)
    
    pfnProgress(progress,pszMessage,pProgressArg);
    double c0mod=0;
    double c1mod=0;

    double errMod=imgreg.getRMSE(imgReaderModel1,imgReaderModel2,c0mod,c1mod,verbose_opt[0]);

    bool update=(tobservation_opt[obsindex]==modindex);
    if(update){//update
      imgReaderObs.open(observation_opt[obsindex]);
      imgReaderObs.setNoData(obsnodata_opt);
      //calculate regression between model and observation
      errObs=imgreg.getRMSE(imgReaderModel1,imgReaderObs,c0obs,c1obs,verbose_opt[0]);
    }
    //prediction (also to fill cloudy pixels in update mode)
    string input;
    if(output_opt.size()==model_opt.size())
      //todo: initialize output_opt[0] with model[0]
      input=output_opt[modindex];
    else{
      ostringstream outputstring;
      outputstring << output_opt[0] << "_" << modindex << ".tif";
      input=outputstring.str();
    }
    ImgReaderGdal imgReaderEst(input);

    vector<double> obsBuffer;
    vector<double> estReadBuffer;
    vector<double> estWriteBuffer;
    vector<double> uncertWriteBuffer;
    for(int irow=0;irow<imgWriterEst.nrOfRow();++irow){
      assert(irow<imgReaderEst.nrOfRow());
      imgReaderEst.readData(estReadBuffer,GDT_Float64,irow);
      if(update){
	imgReaderObs.readData(estWriteBuffer,GDT_Float64,irow);
	//todo: write uncertainty image: 0 if observation, 
      }
      for(int icol=0;icol<imgWriterEst.nrOfCol();++icol){
	if(update){//check for nodata in observation
	  if(!imgReaderObs.isNoData(estWriteBuffer[icol])){
	    uncertWriteBuffer[icol]=0;
	    continue;//no need to estimate
	  }
	}
	//predict
	double estValue=estReadBuffer[icol];
	double certNorm=(errMod*errMod+errObs*errObs);
	double certMod=errObs*errObs/certNorm;
	double certObs=errMod*errMod/certNorm;
	estWriteBuffer[icol]=(c0mod+c1mod*estValue)*certMod;
	estWriteBuffer[icol]+=(c0obs+c1obs*estValue)*certObs;

	double totalUncertainty=0;
	if(errMod<eps_opt[0])
	  totalUncertainty=errObs;
	else if(errObs<eps_opt[0])
	  totalUncertainty=errMod;
	else{
	  totalUncertainty=1.0/errMod/errMod+1/errObs/errObs;
	  totalUncertainty=sqrt(1.0/totalUncertainty);
	}
	uncertWriteBuffer[icol]=totalUncertainty;
      }
      imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,0);
      imgWriterEst.writeData(uncertWriteBuffer,GDT_Float64,irow,1);
    }

    imgWriterEst.close();
    imgReaderEst.close();

    if(update){
      imgReaderObs.close();
      ++obsindex;
    }
    imgReaderModel1.close();
    imgReaderModel2.close();
  }
}

