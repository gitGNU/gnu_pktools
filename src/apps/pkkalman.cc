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
#include "algorithms/StatFactory.h"
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"

using namespace std;
/*------------------
  Main procedure
  ----------------*/
int main(int argc,char **argv) {
  Optionpk<string> model_opt("mod","model","model input datasets, e.g., MODIS (use: -mod model1 -mod model2 etc.");
  Optionpk<string> observation_opt("obs","observation","observation input datasets, e.g., landsat (use: -obs obs1 -obs obs2 etc.");
  Optionpk<int> tobservation_opt("tobs","tobservation","time sequence of observation input (sequence must have exact same length as observation input)"); 
  Optionpk<string> output_opt("o", "output", "Suffix for output image datasets");
  Optionpk<float> threshold_opt("t", "threshold", "threshold for selecting samples (randomly). Provide probability in percentage (>0) or absolute (<0).", 100);
  Optionpk<double> modnodata_opt("modnodata", "modnodata", "invalid value for model input", 0);
  Optionpk<int> bndmodnodata_opt("bmnodata", "bndmodnodata", "Bands in model input to check if pixel is valid (used for srcnodata, min and max options)", 0);
  Optionpk<double> obsnodata_opt("obsnodata", "obsnodata", "invalid value for observation input", 0);
  Optionpk<int> bndobsnodata_opt("bonodata", "bndobsnodata", "Bands in observation input to check if pixel is valid (used for srcnodata, min and max options)", 0);
  Optionpk<int> down_opt("down", "down", "Downsampling factor for reading model data to calculate regression", 90);
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode when positive", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=model_opt.retrieveOption(argc,argv);
    observation_opt.retrieveOption(argc,argv);
    tobservation_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    threshold_opt.retrieveOption(argc,argv);
    modnodata_opt.retrieveOption(argc,argv);
    bndmodnodata_opt.retrieveOption(argc,argv);
    obsnodata_opt.retrieveOption(argc,argv);
    bndobsnodata_opt.retrieveOption(argc,argv);
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

  while(modnodata_opt.size()<bndmodnodata_opt.size())
    modnodata_opt.push_back(modnodata_opt[0]);
  while(bndmodnodata_opt.size()<modnodata_opt.size())
    bndmodnodata_opt.push_back(bndmodnodata_opt[0]);
  while(obsnodata_opt.size()<bndobsnodata_opt.size())
    obsnodata_opt.push_back(obsnodata_opt[0]);
  while(bndobsnodata_opt.size()<obsnodata_opt.size())
    bndobsnodata_opt.push_back(bndobsnodata_opt[0]);

  statfactory::StatFactory stat;

  vector<ImgReaderGdal> imgReaderModel(model_opt.size());
  vector<ImgReaderGdal> imgReaderObs(observation_opt.size());
  vector<ImgWriterGdal> imgWriterPred(model_opt.size());

  int obsindex=0;
  //forward step

  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  srand(time(NULL));

  for(int modindex=0;modindex<model_opt.size()-1;++modindex){
    //calculate regression between two subsequence model inputs
    imgReaderModel[modindex].open(model_opt[modindex]);
    imgReaderModel[modindex+1].open(model_opt[modindex]);
    //calculate regression
    //we could re-use the points from second image from last run, but
    //to keep it general, we must redo it (overlap might have changed)
    pfnProgress(progress,pszMessage,pProgressArg);
    int icol1=0,irow1=0;
    vector<double> rowBuffer1(imgReaderModel[modindex].nrOfCol());
    vector<double> rowBuffer2(imgReaderModel[modindex+1].nrOfCol());
    vector<double> buffer1;
    vector<double> buffer2;

    for(irow1=0;irow1<imgReaderModel[modindex].nrOfRow();++irow1){
      if(irow1%down_opt[0])
	continue;
      icol1=0;
      double icol2=0,irow2=0;
      double geox=0,geoy=0;
      imgReaderModel[modindex].readData(rowBuffer1,GDT_Float64,irow1);
      imgReaderModel[modindex].image2geo(icol1,irow1,geox,geoy);
      imgReaderModel[modindex+1].geo2image(geox,geoy,icol2,irow2);
      icol2=static_cast<int>(icol2);
      irow2=static_cast<int>(irow2);
      imgReaderModel[modindex+1].readData(rowBuffer2,GDT_Float64,irow2);
      for(icol1=0;icol1<imgReaderModel[modindex].nrOfCol();++icol1){
	if(icol1%down_opt[0])
	  continue;
	if(threshold_opt[0]>0){//percentual value
	  double p=static_cast<double>(rand())/(RAND_MAX);
	  p*=100.0;
	  if(p>threshold_opt[0])
	    continue;//do not select for now, go to next column
	}
	else if(buffer1.size()>-threshold_opt[0])//absolute value
	  continue;//do not select any more pixels
	imgReaderModel[modindex].image2geo(icol1,irow1,geox,geoy);
	imgReaderModel[modindex+1].geo2image(geox,geoy,icol2,irow2);
	icol2=static_cast<int>(icol2);
	irow2=static_cast<int>(irow2);
	//check for nodata
	double valmod1=rowBuffer1[icol1];
	double valmod2=rowBuffer2[icol2];
	bool readValid=true;
	for(int vband=0;vband<bndmodnodata_opt.size();++vband){
	  if(modnodata_opt.size()>vband){
	    if(valmod1==modnodata_opt[vband] || valmod2==modnodata_opt[vband]){
	      readValid=false;
	      break;
	    }
	  }
	}
	buffer1.push_back(valmod1);
	buffer2.push_back(valmod2);
      }
    }
    double c0=0;
    double c1=0;
    double err=stat.linear_regression_err(buffer1,buffer2,c0,c1);
    if(verbose_opt[0])
      cout << "linear regression model-model based on " << buffer1.size() << " points: " << c0 << "+" << c1 << " * x " << " with rmse: " << err << endl;

    if(tobservation_opt[obsindex]==modindex){//update
      imgReaderObs[obsindex].open(observation_opt[obsindex]);
      //calculate regression 
      //set ImgWriterPred[modindex]=imgReaderObs[obindex]
      //calculate regression between model and observation

    }
    else{//prediction
    }
    ++obsindex;
  }
}
