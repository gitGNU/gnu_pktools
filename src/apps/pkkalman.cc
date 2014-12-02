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
#include <algorithm>
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
  Optionpk<string> direction_opt("dir","direction","direction to run model (forward|backward|smooth)","forward");
  Optionpk<string> model_opt("mod","model","model input datasets, e.g., MODIS (use: -mod model1 -mod model2 etc.");
  Optionpk<string> observation_opt("obs","observation","observation input datasets, e.g., landsat (use: -obs obs1 -obs obs2 etc.");
  Optionpk<int> tmodel_opt("tmod","tmodel","time sequence of model input. Sequence must have exact same length as model input. Leave empty to have default sequence 0,1,2,etc."); 
  Optionpk<int> tobservation_opt("tobs","tobservation","time sequence of observation input. Sequence must have exact same length as observation input)"); 
  Optionpk<string>  projection_opt("a_srs", "a_srs", "Override the projection for the output file (leave blank to copy from input file, use epsg:3035 to use European projection and force to European grid");
  Optionpk<string> outputfw_opt("ofw", "outputfw", "Output raster dataset for forward model");
  Optionpk<string> outputbw_opt("obw", "outputbw", "Output raster dataset for backward model");
  Optionpk<string> outputfb_opt("ofb", "outputfb", "Output raster dataset for smooth model");
  Optionpk<double> modnodata_opt("modnodata", "modnodata", "invalid value for model input", 0);
  Optionpk<double> obsnodata_opt("obsnodata", "obsnodata", "invalid value for observation input", 0);
  Optionpk<double> modoffset_opt("modoffset", "modoffset", "offset used to read model input dataset (value=offset+scale*readValue", 0);
  Optionpk<double> obsoffset_opt("obsoffset", "obsoffset", "offset used to read observation input dataset (value=offset+scale*readValue", 0);
  Optionpk<double> modscale_opt("modscale", "modscale", "scale used to read model input dataset (value=offset+scale*readValue", 1);
  Optionpk<double> obsscale_opt("obsscale", "obsscale", "scale used to read observation input dataset (value=offset+scale*readValue", 1);
  Optionpk<double> eps_opt("eps", "eps", "epsilon for non zero division", 0.00001);
  Optionpk<double> uncertModel_opt("um", "uncertmodel", "Multiply this value with std dev of first model image to obtain uncertainty of model",2);
  Optionpk<double> uncertObs_opt("uo", "uncertobs", "Uncertainty of valid observations",0);
  Optionpk<double> deltaObs_opt("dobs", "deltaobs", "Threshold for relative difference (in percentage) in observation and model pixel");
  Optionpk<double> uncertNodata_opt("unodata", "uncertnodata", "Uncertainty in case of no-data values in observation", 10000);
  // Optionpk<double> regTime_opt("rt", "regtime", "Relative Weight for regression in time series", 1.0);
  // Optionpk<double> regSensor_opt("rs", "regsensor", "Relative Weight for regression in sensor series", 1.0);
  Optionpk<int> down_opt("down", "down", "Downsampling factor for reading model data to calculate regression", 9);
  Optionpk<float> threshold_opt("th", "threshold", "threshold for selecting samples (randomly). Provide probability in percentage (>0) or absolute (<0).", 0);
  Optionpk<int> minreg_opt("minreg", "minreg", "Minimum number of pixels to take into account for regression", 5, 2);
  // Optionpk<bool> regObs_opt("regobs", "regobs", "Perform regression between modeled and observed value",false);
  // Optionpk<double> checkDiff_opt("diff", "diff", "Flag observation as invalid if difference with model is above uncertainty",false);
  Optionpk<unsigned short> window_opt("win", "window", "window size for calculating regression (use 0 for global)", 0);
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image","GTiff",2);
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode when positive", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=direction_opt.retrieveOption(argc,argv);
    model_opt.retrieveOption(argc,argv);
    observation_opt.retrieveOption(argc,argv);
    tmodel_opt.retrieveOption(argc,argv);
    tobservation_opt.retrieveOption(argc,argv);
    projection_opt.retrieveOption(argc,argv);
    outputfw_opt.retrieveOption(argc,argv);
    outputbw_opt.retrieveOption(argc,argv);
    outputfb_opt.retrieveOption(argc,argv);
    modnodata_opt.retrieveOption(argc,argv);
    obsnodata_opt.retrieveOption(argc,argv);
    modoffset_opt.retrieveOption(argc,argv);
    modscale_opt.retrieveOption(argc,argv);
    obsoffset_opt.retrieveOption(argc,argv);
    obsscale_opt.retrieveOption(argc,argv);
    eps_opt.retrieveOption(argc,argv);
    uncertModel_opt.retrieveOption(argc,argv);
    uncertObs_opt.retrieveOption(argc,argv);
    deltaObs_opt.retrieveOption(argc,argv);
    uncertNodata_opt.retrieveOption(argc,argv);
    // regTime_opt.retrieveOption(argc,argv);
    // regSensor_opt.retrieveOption(argc,argv);
    down_opt.retrieveOption(argc,argv);
    threshold_opt.retrieveOption(argc,argv);
    minreg_opt.retrieveOption(argc,argv);
    // regObs_opt.retrieveOption(argc,argv);
    // checkDiff_opt.retrieveOption(argc,argv);
    window_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
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
    if(direction_opt[0]=="smooth"){
      if(outputfw_opt.empty()){
	errorStream << "Error: output forward datasets is not provided, use option -ofw" << endl;
	throw(errorStream.str());
      }
      if(outputbw_opt.empty()){
	errorStream << "Error: output backward datasets is not provided, use option -obw" << endl;
	throw(errorStream.str());
      }
      if(outputfb_opt.empty()){
	errorStream << "Error: output smooth datasets is not provided, use option -ofb" << endl;
	throw(errorStream.str());
      }
    }
    else{
      if(direction_opt[0]=="forward"&&outputfw_opt.empty()){
	errorStream << "Error: output forward datasets is not provided, use option -ofw" << endl;
	throw(errorStream.str());
      }
      else if(direction_opt[0]=="backward"&&outputbw_opt.empty()){
	errorStream << "Error: output backward datasets is not provided, use option -obw" << endl;
	throw(errorStream.str());
      }

      if(model_opt.size()<observation_opt.size()){
	errorStream << "Error: sequence of models should be larger than observations" << endl;
	throw(errorStream.str());
      }
      if(tmodel_opt.size()!=model_opt.size()){
	if(tmodel_opt.empty())
	  cout << "Warning: time sequence is not provided, self generating time sequence from 0 to " << model_opt.size() << endl;
	else
	  cout << "Warning: time sequence provided (" << tmodel_opt.size() << ") does not match number of model raster datasets (" << model_opt.size() << ")" << endl;
	tmodel_opt.clear();
	for(int tindex=0;tindex<model_opt.size();++tindex)
	  tmodel_opt.push_back(tindex);
      }
      if(tobservation_opt.size()!=observation_opt.size()){
	errorStream << "Error: time sequence for observation must match size of observation dataset" << endl;
	throw(errorStream.str());
      }
    }
  }
  catch(string errorString){
    std::cout << errorString << std::endl;
    exit(1);
  }

  statfactory::StatFactory stat;
  stat.setNoDataValues(modnodata_opt);
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
  if(projection_opt.empty())
    projection_opt.push_back(imgReaderObs.getProjection());
  double geotransform[6];
  imgReaderObs.getGeoTransform(geotransform);

  string imageType=imgReaderObs.getImageType();
  if(oformat_opt.size())//default
    imageType=oformat_opt[0];
  if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
    string theInterleave="INTERLEAVE=";
    theInterleave+=imgReaderObs.getInterleave();
    option_opt.push_back(theInterleave);
  }

  imgReaderObs.close();

  int obsindex=0;

  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;

  imgreg.setDown(down_opt[0]);
  imgreg.setThreshold(threshold_opt[0]);

  double c0modGlobal=0;//time regression coefficient c0 calculated on entire image 
  double c1modGlobal=0;//time regression coefficient c1 calculated on entire image 
  double c0mod=0;//time regression coefficient c0 calculated on local window
  double c1mod=0;//time regression coefficient c1 calculated on local window

  double c0obs=0;
  double c1obs=0;
  double errObs=uncertNodata_opt[0];//start with high initial value in case we do not have first observation at time 0

  vector<int> relobsindex;
  // cout << tmodel_opt << endl;
  // cout << tobservation_opt << endl;

  for(int tindex=0;tindex<tobservation_opt.size();++tindex){
    vector<int>::iterator modit;
    modit=lower_bound(tmodel_opt.begin(),tmodel_opt.end(),tobservation_opt[tindex]);
    int relpos=modit-tmodel_opt.begin();
    // if(relpos<0)
    //   relpos=0;
    relobsindex.push_back(relpos);
    if(verbose_opt[0])
      cout << "tobservation_opt[tindex] " << tobservation_opt[tindex] << " " << relobsindex.back() << " " << observation_opt[tindex] << " " << model_opt[relpos] << endl;
    // if(verbose_opt[0])
    //   cout << "tobservation_opt[tindex] " << tobservation_opt[tindex] << " " << relobsindex.back() << endl;
  }

  if(find(direction_opt.begin(),direction_opt.end(),"forward")!=direction_opt.end()){
    ///////////////////////////// forward model /////////////////////////
    obsindex=0;
    if(verbose_opt[0])
      cout << "Running forward model" << endl;
    //initialization
    string output;
    if(outputfw_opt.size()==model_opt.size())
      output=outputfw_opt[0];
    else{
      ostringstream outputstream;
      outputstream << outputfw_opt[0] << "_" << tmodel_opt[0] << ".tif";
      output=outputstream.str();
    }
    if(verbose_opt[0])
      cout << "Opening image " << output << " for writing " << endl;
    imgWriterEst.open(output,ncol,nrow,2,GDT_Float32,imageType,option_opt);
    imgWriterEst.setProjectionProj4(projection_opt[0]);
    imgWriterEst.setGeoTransform(geotransform);
    imgWriterEst.GDALSetNoDataValue(obsnodata_opt[0]);

    if(verbose_opt[0]){
      cout << "processing time " << tmodel_opt[0] << endl;
      if(obsindex<relobsindex.size())
	cout << "next observation " << tmodel_opt[relobsindex[obsindex]] << endl;
      else
	cout << "There is no next observation" << endl;
    }

    try{
      imgReaderModel1.open(model_opt[0]);
      imgReaderModel1.setNoData(modnodata_opt);
      imgReaderModel1.setOffset(modoffset_opt[0]);
      imgReaderModel1.setScale(modscale_opt[0]);
    }
    catch(string errorString){
      cerr << errorString << endl;
    }
    catch(...){
      cerr << "Error opening file " << model_opt[0] << endl;
    }
      
    //calculate standard deviation of image to serve as model uncertainty
    GDALRasterBand* rasterBand;
    rasterBand=imgReaderModel1.getRasterBand(0);
    double minValue, maxValue, meanValue, stdDev;
    void* pProgressData;
    rasterBand->ComputeStatistics(0,&minValue,&maxValue,&meanValue,&stdDev,pfnProgress,pProgressData);
    double x=0;
    double y=0;
    double modRow=0;
    double modCol=0;
    if(relobsindex[0]>0){//initialize output_opt[0] as model[0]
      //write first model as output
      if(verbose_opt[0])
	cout << "write first model as output" << endl;
      for(int irow=0;irow<nrow;++irow){
	vector<double> estReadBuffer;
	vector<double> estWriteBuffer(ncol);
	vector<double> uncertWriteBuffer(ncol);
	imgWriterEst.image2geo(0,irow,x,y);
	imgReaderModel1.geo2image(x,y,modCol,modRow);
	assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	try{
	  imgReaderModel1.readData(estReadBuffer,GDT_Float64,modRow);
	  //simple nearest neighbor
	  //stat.nearUp(estReadBuffer,estWriteBuffer);

	  for(int icol=0;icol<ncol;++icol){
	    imgWriterEst.image2geo(icol,irow,x,y);
	    imgReaderModel1.geo2image(x,y,modCol,modRow);
	    estWriteBuffer[icol]=estReadBuffer[modCol];
	    uncertWriteBuffer[icol]=uncertModel_opt[0]*stdDev;
	  }
	  imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,0);
	  imgWriterEst.writeData(uncertWriteBuffer,GDT_Float64,irow,1);
	}
	catch(string errorString){
	  cerr << errorString << endl;
	}
	catch(...){
	  cerr << "Error writing file " << imgWriterEst.getFileName() << endl;
	}
      }
    }
    else{//we have an observation at time 0
      if(verbose_opt[0])
	cout << "we have an observation at time 0" << endl;
      imgReaderObs.open(observation_opt[0]);
      imgReaderObs.getGeoTransform(geotransform);
      imgReaderObs.setNoData(obsnodata_opt);
      imgReaderObs.setOffset(obsoffset_opt[0]);
      imgReaderObs.setScale(obsscale_opt[0]);
      for(int irow=0;irow<nrow;++irow){
	vector<double> estReadBuffer;
	imgWriterEst.image2geo(0,irow,x,y);
	imgReaderModel1.geo2image(x,y,modCol,modRow);
	assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	imgReaderModel1.readData(estReadBuffer,GDT_Float64,modRow);
	vector<double> estWriteBuffer(ncol);
	vector<double> uncertWriteBuffer(ncol);
	vector<double> uncertObsBuffer;
	imgReaderObs.readData(estWriteBuffer,GDT_Float64,irow,0);
	if(imgReaderObs.nrOfBand()>1)
	  imgReaderObs.readData(uncertObsBuffer,GDT_Float64,irow,1);
	for(int icol=0;icol<ncol;++icol){
	  if(imgReaderObs.isNoData(estWriteBuffer[icol])){
	    imgWriterEst.image2geo(icol,irow,x,y);
	    imgReaderModel1.geo2image(x,y,modCol,modRow);
	    assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	    if(imgReaderModel1.isNoData(estReadBuffer[modCol]))//if both obs and model are no-data, set obs to nodata
	      estWriteBuffer[icol]=obsnodata_opt[0];
	    else
	      estWriteBuffer[icol]=estReadBuffer[modCol];
	    uncertWriteBuffer[icol]=uncertModel_opt[0]*stdDev;
	  }
	  else{
	    double uncertObs=uncertObs_opt[0];
	    if(uncertObsBuffer.size()>icol)
	      uncertObs=uncertObsBuffer[icol];
	    if(uncertModel_opt[0]*stdDev+uncertObs>eps_opt[0]){
	      imgWriterEst.image2geo(icol,irow,x,y);
	      imgReaderModel1.geo2image(x,y,modCol,modRow);
	      assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	      // double noemer=uncertObs*uncertObs+stdDev*stdDev;//todo: multiply stdDev with uncertModel_opt[0]
	      // estWriteBuffer[icol]*=uncertModel_opt[0]*stdDev*stdDev/noemer;
	      // estWriteBuffer[icol]+=uncertModel_opt[0]*uncertObs*uncertObs/noemer;//todo:check! error?
	      double kalmanGain=1;
	      if(!imgReaderModel1.isNoData(estReadBuffer[modCol])){
		//if model is no-data, retain observation value
		kalmanGain=uncertModel_opt[0]*stdDev/(uncertModel_opt[0]*stdDev+uncertObs);
		estWriteBuffer[icol]=estReadBuffer[modCol]+kalmanGain*(estWriteBuffer[icol]-estReadBuffer[modCol]);
	      }
	      uncertWriteBuffer[icol]=uncertModel_opt[0]*stdDev*(1-kalmanGain);
	    }
	    else{
	      //no need to fill write buffer (already done in imgReaderObs.readData
	      uncertWriteBuffer[icol]=uncertObs;
	    }
	  }
	}
	imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,0);
	imgWriterEst.writeData(uncertWriteBuffer,GDT_Float64,irow,1);
      }
      imgReaderObs.close();
      ++obsindex;
    }
    imgReaderModel1.close();
    imgWriterEst.close();

    for(int modindex=1;modindex<model_opt.size();++modindex){
      if(verbose_opt[0]){
	cout << "processing time " << tmodel_opt[modindex] << endl;
	if(obsindex<relobsindex.size())
	  cout << "next observation " << tmodel_opt[relobsindex[obsindex]] << endl;
	else
	  cout << "There is no next observation" << endl;
      }
      string output;
      if(outputfw_opt.size()==model_opt.size())
	output=outputfw_opt[modindex];
      else{
	ostringstream outputstream;
	outputstream << outputfw_opt[0] << "_" << tmodel_opt[modindex] << ".tif";
	// outputstream << output_opt[0] << "_" << modindex+1 << ".tif";
	output=outputstream.str();
      }
    
      //two band output band0=estimation, band1=uncertainty
      imgWriterEst.open(output,ncol,nrow,2,GDT_Float32,imageType,option_opt);
      imgWriterEst.setProjectionProj4(projection_opt[0]);
      imgWriterEst.setGeoTransform(geotransform);
      imgWriterEst.GDALSetNoDataValue(obsnodata_opt[0]);

      //calculate regression between two subsequence model inputs
      imgReaderModel1.open(model_opt[modindex-1]);
      imgReaderModel1.setNoData(modnodata_opt);
      imgReaderModel1.setOffset(modoffset_opt[0]);
      imgReaderModel1.setScale(modscale_opt[0]);
      imgReaderModel2.open(model_opt[modindex]);
      imgReaderModel2.setNoData(modnodata_opt);
      imgReaderModel2.setOffset(modoffset_opt[0]);
      imgReaderModel2.setScale(modscale_opt[0]);
      //calculate regression
      //we could re-use the points from second image from last run, but
      //to keep it general, we must redo it (overlap might have changed)
    
      pfnProgress(progress,pszMessage,pProgressArg);

      if(verbose_opt[0])
	cout << "Calculating regression for " << imgReaderModel1.getFileName() << " " << imgReaderModel2.getFileName() << endl;
      
      double errMod=imgreg.getRMSE(imgReaderModel1,imgReaderModel2,c0modGlobal,c1modGlobal);
      // double errMod=imgreg.getRMSE(imgReaderModel1,imgReaderModel2,c0modGlobal,c1modGlobal,verbose_opt[0]);

      bool update=false;
      if(obsindex<relobsindex.size()){
	update=(relobsindex[obsindex]==modindex);
      }
      if(update){
	if(verbose_opt[0])
	  cout << "***update " << relobsindex[obsindex] << " = " << modindex << " " << observation_opt[obsindex] << " ***" << endl;

	imgReaderObs.open(observation_opt[obsindex]);
	imgReaderObs.getGeoTransform(geotransform);
	imgReaderObs.setNoData(obsnodata_opt);
	imgReaderObs.setOffset(obsoffset_opt[0]);
	imgReaderObs.setScale(obsscale_opt[0]);
	//calculate regression between model and observation
	if(verbose_opt[0])
	  cout << "Calculating regression for " << imgReaderModel2.getFileName() << " " << imgReaderObs.getFileName() << endl;
	errObs=imgreg.getRMSE(imgReaderModel2,imgReaderObs,c0obs,c1obs,verbose_opt[0]);
	if(verbose_opt[0])
	  cout << "c0obs, c1obs: " << c0obs << ", " << c1obs << endl;
      }
      //prediction (also to fill cloudy pixels in update mode)
      string input;
      if(outputfw_opt.size()==model_opt.size())
	input=outputfw_opt[modindex-1];
      else{
	ostringstream outputstream;
	outputstream << outputfw_opt[0] << "_" << tmodel_opt[modindex-1] << ".tif";
	input=outputstream.str();
      }
      ImgReaderGdal imgReaderEst(input);
      imgReaderEst.setNoData(obsnodata_opt);
      imgReaderEst.setOffset(obsoffset_opt[0]);
      imgReaderEst.setScale(obsscale_opt[0]);

      vector<double> obsLineBuffer;
      vector<double> obsWindowBuffer;//buffer for observation to calculate average corresponding to model pixel
      vector<double> model1LineBuffer;
      vector<double> model2LineBuffer;
      vector<double> model1buffer;//buffer for model 1 to calculate time regression based on window
      vector<double> model2buffer;//buffer for model 2 to calculate time regression based on window
      vector<double> uncertObsLineBuffer;
      vector<double> estReadBuffer;
      // vector<double> estWindowBuffer;//buffer for estimate to calculate average corresponding to model pixel
      vector<double> uncertReadBuffer;
      vector<double> estWriteBuffer(ncol);
      vector<double> uncertWriteBuffer(ncol);

      for(int irow=0;irow<imgWriterEst.nrOfRow();++irow){
	assert(irow<imgReaderEst.nrOfRow());
	//not needed here, because we read entire window for each pixel...
	imgReaderEst.readData(estReadBuffer,GDT_Float64,irow,0);
	imgReaderEst.readData(uncertReadBuffer,GDT_Float64,irow,1);
	//read model2 in case current estimate is nodata
	imgReaderEst.image2geo(0,irow,x,y);
	imgReaderModel2.geo2image(x,y,modCol,modRow);
	assert(modRow>=0&&modRow<imgReaderModel2.nrOfRow());
	imgReaderModel2.readData(model2LineBuffer,GDT_Float64,modRow,0);

	imgReaderModel1.geo2image(x,y,modCol,modRow);
	assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	imgReaderModel1.readData(model1LineBuffer,GDT_Float64,modRow,0);

	if(update){
	  imgReaderObs.readData(obsLineBuffer,GDT_Float64,irow,0);
	  if(imgReaderObs.nrOfBand()>1)
	    imgReaderObs.readData(uncertObsLineBuffer,GDT_Float64,irow,1);
	}
	for(int icol=0;icol<imgWriterEst.nrOfCol();++icol){
	  int minCol=(icol>down_opt[0]/2) ? icol-down_opt[0]/2 : 0;
	  int maxCol=(icol+down_opt[0]/2<imgReaderEst.nrOfCol()) ? icol+down_opt[0]/2 : imgReaderEst.nrOfCol()-1;
	  int minRow=(irow>down_opt[0]/2) ? irow-down_opt[0]/2 : 0;
	  int maxRow=(irow+down_opt[0]/2<imgReaderEst.nrOfRow()) ? irow+down_opt[0]/2 : imgReaderEst.nrOfRow()-1;
	  if(update)
	    imgReaderObs.readDataBlock(obsWindowBuffer,GDT_Float64,minCol,maxCol,minRow,maxRow,0);
	  // imgReaderEst.readDataBlock(estWindowBuffer,GDT_Float64,minCol,maxCol,minRow,maxRow,0);
	  double estValue=estReadBuffer[icol];
	  // double estValue=estWindowBuffer[estWindowBuffer.size()/2];
	  double estMeanValue=0;//stat.mean(estWindowBuffer);
	  double nvalid=0;
	  //time update
	  imgReaderEst.image2geo(icol,irow,x,y);
	  imgReaderModel2.geo2image(x,y,modCol,modRow);
	  assert(modRow>=0&&modRow<imgReaderModel2.nrOfRow());
	  double modValue=model2LineBuffer[modCol];
	  bool estNodata=imgReaderEst.isNoData(estValue);
	  //todo: debug checkDiff_opt 
	  // if(checkDiff_opt[0]){
	  //   vector<double>::iterator itwin=estWindowBuffer.begin();
	  //   while(itwin!=estWindowBuffer.end()){
	  //     if(!imgReaderEst.isNoData(*itwin)){
	  // 	estMeanValue+=*itwin;
	  // 	++nvalid;
	  //     }
	  //     ++itwin;
	  //   }
	  //   estMeanValue/=nvalid;
	  //   double difference=(estMeanValue-c0obs)/c1obs-model1LineBuffer[modCol];
	  //   difference*=difference;
	  //   estNodata=estNodata||(sqrt(difference)>uncertModel_opt[0]*stdDev);
	  // }
	  if(estNodata){
	    //we have not found any valid data yet, better here to take the current model value if valid
	    if(imgReaderModel2.isNoData(modValue)){//if both estimate and model are no-data, set obs to nodata
	      estWriteBuffer[icol]=obsnodata_opt[0];
	      uncertWriteBuffer[icol]=uncertNodata_opt[0];
	    }
	    else{
	      estWriteBuffer[icol]=modValue;
	      uncertWriteBuffer[icol]=uncertModel_opt[0]*stdDev;
	    }
	  }	  
	  else{
	    if(window_opt[0]>0){
	      try{
		// imgReaderModel2.geo2image(x,y,modCol,modRow);//did that already
		minCol=(modCol>window_opt[0]/2) ? modCol-window_opt[0]/2 : 0;
		maxCol=(modCol+window_opt[0]/2<imgReaderModel2.nrOfCol()) ? modCol+window_opt[0]/2 : imgReaderModel2.nrOfCol()-1;
		minRow=(modRow>window_opt[0]/2) ? modRow-window_opt[0]/2 : 0;
		maxRow=(modRow+window_opt[0]/2<imgReaderModel2.nrOfRow()) ? modRow+window_opt[0]/2 : imgReaderModel2.nrOfRow()-1;
		imgReaderModel2.readDataBlock(model2buffer,GDT_Float64,minCol,maxCol,minRow,maxRow,0);

		imgReaderModel1.geo2image(x,y,modCol,modRow);
		assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
		minCol=(modCol>window_opt[0]/2) ? modCol-window_opt[0]/2 : 0;
		maxCol=(modCol+window_opt[0]/2<imgReaderModel1.nrOfCol()) ? modCol+window_opt[0]/2 : imgReaderModel1.nrOfCol()-1;
		minRow=(modRow>window_opt[0]/2) ? modRow-window_opt[0]/2 : 0;
		maxRow=(modRow+window_opt[0]/2<imgReaderModel1.nrOfRow()) ? modRow+window_opt[0]/2 : imgReaderModel1.nrOfRow()-1;
		imgReaderModel1.readDataBlock(model1buffer,GDT_Float64,minCol,maxCol,minRow,maxRow,0);
		// imgReaderEst.image2geo(icol,irow,x,y);
	      }
	      catch(string errorString){
		cerr << "Error reading data block for " << minCol << "-" << maxCol << ", " << minRow << "-" << maxRow << endl;
	      }
	      //erase no-data from buffer
	      vector<double>::iterator it1=model1buffer.begin();
	      vector<double>::iterator it2=model2buffer.begin();
	      while(it1!=model1buffer.end()&&it2!=model2buffer.end()){
		//erase nodata
		bool modNodata=false;
		modNodata=modNodata||imgReaderModel1.isNoData(*it1);
		modNodata=modNodata||imgReaderModel2.isNoData(*it2);
		//todo: debug checkDiff_opt 
		// if(checkDiff_opt[0]){
		//   double difference=(estMeanValue-c0obs)/c1obs-*it1;
		//   difference*=difference;
		//   modNodata=modNodata||(sqrt(difference)>uncertModel_opt[0]*stdDev);
		// }
		if(modNodata){
		  model1buffer.erase(it1);
		  model2buffer.erase(it2);
		}
		else{
		  ++it1;
		  ++it2;
		}
	      }
	      if(model1buffer.size()>minreg_opt[0]&&model2buffer.size()>minreg_opt[0])
		errMod=stat.linear_regression_err(model1buffer,model2buffer,c0mod,c1mod);
	      else{//use global regression...
		c0mod=c0modGlobal;
		c1mod=c1modGlobal;
	      }
	    }
	    else{
	      c0mod=c0modGlobal;
	      c1mod=c1modGlobal;
	    }
	    double certNorm=(errMod*errMod+errObs*errObs);
	    double certMod=errObs*errObs/certNorm;
	    double certObs=errMod*errMod/certNorm;
	    double regTime=(c0mod+c1mod*estValue)*certMod;

	    // double regSensor=(c0obs+c1obs*estValue)*certObs;
	    double regSensor=(c0obs+c1obs*modValue)*certObs;
	    estWriteBuffer[icol]=regTime+regSensor;
	    double totalUncertainty=0;
	    if(errMod<eps_opt[0])
	      totalUncertainty=errObs;
	    else if(errObs<eps_opt[0])
	      totalUncertainty=errMod;
	    else{
	      totalUncertainty=1.0/errMod/errMod+1/errObs/errObs;
	      totalUncertainty=sqrt(1.0/totalUncertainty);
	    }
	    uncertWriteBuffer[icol]=totalUncertainty+uncertReadBuffer[icol];
	  }
	  //observation update
	  if(update&&!imgReaderObs.isNoData(obsLineBuffer[icol])){
	    bool doUpdate=true;
	    if(deltaObs_opt.size()){
	      statfactory::StatFactory statobs;
	      statobs.setNoDataValues(obsnodata_opt);
	      double obsMeanValue=statobs.mean(obsWindowBuffer);
	      double difference=(obsMeanValue-c0obs)/c1obs-modValue;
	      difference/=modValue;//make relative difference
	      if(deltaObs_opt[0]<0)
		doUpdate=((100*difference)<deltaObs_opt[0]);
	      else{
		difference*=difference;
		doUpdate=(100*sqrt(difference)<deltaObs_opt[0]);
	      }
	    }
	    if(doUpdate){
	      double kalmanGain=1;
	      double uncertObs=uncertObs_opt[0];
	      if(uncertObsLineBuffer.size()>icol)
		uncertObs=uncertObsLineBuffer[icol];
	      if((uncertWriteBuffer[icol]+uncertObs)>eps_opt[0])
		kalmanGain=uncertWriteBuffer[icol]/(uncertWriteBuffer[icol]+uncertObs);
	      assert(kalmanGain<=1);
	      estWriteBuffer[icol]+=kalmanGain*(obsLineBuffer[icol]-estWriteBuffer[icol]);
	      uncertWriteBuffer[icol]*=(1-kalmanGain);
	    }
	  }
	}
	imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,0);
	imgWriterEst.writeData(uncertWriteBuffer,GDT_Float64,irow,1);
	progress=static_cast<float>((irow+1.0)/imgWriterEst.nrOfRow());
	pfnProgress(progress,pszMessage,pProgressArg);
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
  if(find(direction_opt.begin(),direction_opt.end(),"backward")!=direction_opt.end()){
    ///////////////////////////// backward model /////////////////////////
    obsindex=relobsindex.size()-1;
    if(verbose_opt[0])
      cout << "Running backward model" << endl;
    //initialization
    string output;
    if(outputbw_opt.size()==model_opt.size())
      output=outputbw_opt.back();
    else{
      ostringstream outputstream;
      outputstream << outputbw_opt[0] << "_" << tmodel_opt.back() << ".tif";
      output=outputstream.str();
    }
    if(verbose_opt[0])
      cout << "Opening image " << output << " for writing " << endl;
    imgWriterEst.open(output,ncol,nrow,2,GDT_Float32,imageType,option_opt);
    imgWriterEst.setProjectionProj4(projection_opt[0]);
    imgWriterEst.setGeoTransform(geotransform);
    imgWriterEst.GDALSetNoDataValue(obsnodata_opt[0]);

    if(verbose_opt[0]){
      cout << "processing time " << tmodel_opt.back() << endl;
      if(obsindex<relobsindex.size())
	cout << "next observation " << tmodel_opt[relobsindex[obsindex]] << endl;
      else
	cout << "There is no next observation" << endl;
    }

    try{
      imgReaderModel1.open(model_opt.back());
      imgReaderModel1.setNoData(modnodata_opt);
      imgReaderModel1.setOffset(modoffset_opt[0]);
      imgReaderModel1.setScale(modscale_opt[0]);
    }
    catch(string errorString){
      cerr << errorString << endl;
    }
    catch(...){
      cerr << "Error opening file " << model_opt[0] << endl;
    }

    //calculate standard deviation of image to serve as model uncertainty
    GDALRasterBand* rasterBand;
    rasterBand=imgReaderModel1.getRasterBand(0);
    double minValue, maxValue, meanValue, stdDev;
    void* pProgressData;
    rasterBand->ComputeStatistics(0,&minValue,&maxValue,&meanValue,&stdDev,pfnProgress,pProgressData);
    double x=0;
    double y=0;
    double modRow=0;
    double modCol=0;
    if(relobsindex.back()<model_opt.size()-1){//initialize output_opt.back() as model[0]
      //write last model as output
      if(verbose_opt[0])
	cout << "write last model as output" << endl;
      for(int irow=0;irow<nrow;++irow){
	vector<double> estReadBuffer;
	vector<double> estWriteBuffer(ncol);
	vector<double> uncertWriteBuffer(ncol);
	imgWriterEst.image2geo(0,irow,x,y);
	imgReaderModel1.geo2image(x,y,modCol,modRow);
	assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	try{
	  imgReaderModel1.readData(estReadBuffer,GDT_Float64,modRow);
	  //simple nearest neighbor
	  //stat.nearUp(estReadBuffer,estWriteBuffer);

	  for(int icol=0;icol<imgWriterEst.nrOfCol();++icol){
	    imgWriterEst.image2geo(icol,irow,x,y);	    
	    imgReaderModel1.geo2image(x,y,modCol,modRow);
	    estWriteBuffer[icol]=estReadBuffer[modCol];
	    uncertWriteBuffer[icol]=uncertModel_opt[0]*stdDev;
	  }
	  imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,0);
	  imgWriterEst.writeData(uncertWriteBuffer,GDT_Float64,irow,1);
	}
	catch(string errorString){
	  cerr << errorString << endl;
	}
	catch(...){
	  cerr << "Error writing file " << imgWriterEst.getFileName() << endl;
	}
      }
    }
    else{//we have an observation at end time
      if(verbose_opt[0])
	cout << "we have an observation at end time" << endl;
      imgReaderObs.open(observation_opt.back());
      imgReaderObs.getGeoTransform(geotransform);
      imgReaderObs.setNoData(obsnodata_opt);
      imgReaderObs.setOffset(obsoffset_opt[0]);
      imgReaderObs.setScale(obsscale_opt[0]);
      for(int irow=0;irow<nrow;++irow){
	vector<double> estReadBuffer;
	imgWriterEst.image2geo(0,irow,x,y);
	imgReaderModel1.geo2image(x,y,modCol,modRow);
	assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	imgReaderModel1.readData(estReadBuffer,GDT_Float64,modRow);
	vector<double> estWriteBuffer(ncol);
	vector<double> uncertWriteBuffer(ncol);
	vector<double> uncertObsLineBuffer;
	imgReaderObs.readData(estWriteBuffer,GDT_Float64,irow,0);
	if(imgReaderObs.nrOfBand()>1)
	  imgReaderObs.readData(uncertObsLineBuffer,GDT_Float64,irow,1);
	for(int icol=0;icol<imgWriterEst.nrOfCol();++icol){
	  if(imgReaderObs.isNoData(estWriteBuffer[icol])){
	    imgWriterEst.image2geo(icol,irow,x,y);
	    imgReaderModel1.geo2image(x,y,modCol,modRow);
	    assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	    if(imgReaderModel1.isNoData(estReadBuffer[modCol]))//if both obs and model are no-data, set obs to nodata
	      estWriteBuffer[icol]=obsnodata_opt[0];
	    else
	      estWriteBuffer[icol]=estReadBuffer[modCol];
	    uncertWriteBuffer[icol]=uncertModel_opt[0]*stdDev;
	  }
	  else{
	    double uncertObs=uncertObs_opt[0];
	    if(uncertObsLineBuffer.size()>icol)
	      uncertObs=uncertObsLineBuffer[icol];
	    if(uncertModel_opt[0]*stdDev+uncertObs>eps_opt[0]){
	      imgWriterEst.image2geo(icol,irow,x,y);
	      imgReaderModel1.geo2image(x,y,modCol,modRow);
	      assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	      // double noemer=uncertObs*uncertObs+stdDev*stdDev;//todo: multiply stdDev with uncertModel_opt[0]
	      // estWriteBuffer[icol]*=uncertModel_opt[0]*stdDev*stdDev/noemer;
	      // estWriteBuffer[icol]+=uncertModel_opt[0]*uncertObs*uncertObs/noemer;//todo:check! error?
	      double kalmanGain=1;
	      if(!imgReaderModel1.isNoData(estReadBuffer[modCol])){
		//if model is no-data, retain observation value
		kalmanGain=uncertModel_opt[0]*stdDev/(uncertModel_opt[0]*stdDev+uncertObs);
		estWriteBuffer[icol]=estReadBuffer[modCol]+kalmanGain*(estWriteBuffer[icol]-estReadBuffer[modCol]);
	      }
	      uncertWriteBuffer[icol]=uncertModel_opt[0]*stdDev*(1-kalmanGain);
	    }
	    else{
	      //no need to fill write buffer (already done in imgReaderObs.readData
	      uncertWriteBuffer[icol]=uncertObs;
	    }
	  }
	}
	imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,0);
	imgWriterEst.writeData(uncertWriteBuffer,GDT_Float64,irow,1);
      }
      imgReaderObs.close();
      --obsindex;
    }
    imgReaderModel1.close();
    imgWriterEst.close();

    for(int modindex=model_opt.size()-2;modindex>=0;--modindex){
      if(verbose_opt[0]){
	cout << "processing time " << tmodel_opt[modindex] << endl;
	if(obsindex<relobsindex.size())
	  cout << "next observation " << tmodel_opt[relobsindex[obsindex]] << endl;
	else
	  cout << "There is no next observation" << endl;
      }
      string output;
      if(outputbw_opt.size()==model_opt.size())
	output=outputbw_opt[modindex];
      else{
	ostringstream outputstream;
	outputstream << outputbw_opt[0] << "_" << tmodel_opt[modindex] << ".tif";
	// outputstream << output_opt[0] << "_" << modindex+1 << ".tif";
	output=outputstream.str();
      }

      //two band output band0=estimation, band1=uncertainty
      imgWriterEst.open(output,ncol,nrow,2,GDT_Float32,imageType,option_opt);
      imgWriterEst.setProjectionProj4(projection_opt[0]);
      imgWriterEst.setGeoTransform(geotransform);
      imgWriterEst.GDALSetNoDataValue(obsnodata_opt[0]);

      //calculate regression between two subsequence model inputs
      imgReaderModel1.open(model_opt[modindex+1]);
      imgReaderModel1.setNoData(modnodata_opt);
      imgReaderModel1.setOffset(modoffset_opt[0]);
      imgReaderModel1.setScale(modscale_opt[0]);
      imgReaderModel2.open(model_opt[modindex]);
      imgReaderModel2.setNoData(modnodata_opt);
      imgReaderModel2.setOffset(modoffset_opt[0]);
      imgReaderModel2.setScale(modscale_opt[0]);
      //calculate regression
      //we could re-use the points from second image from last run, but
      //to keep it general, we must redo it (overlap might have changed)
    
      pfnProgress(progress,pszMessage,pProgressArg);

      if(verbose_opt[0])
	cout << "Calculating regression for " << imgReaderModel1.getFileName() << " " << imgReaderModel2.getFileName() << endl;

      double errMod=imgreg.getRMSE(imgReaderModel1,imgReaderModel2,c0modGlobal,c1modGlobal);
      // double errMod=imgreg.getRMSE(imgReaderModel1,imgReaderModel2,c0modGlobal,c1modGlobal,verbose_opt[0]);

      bool update=false;
      if(obsindex<relobsindex.size()){
	update=(relobsindex[obsindex]==modindex);
      }
      if(update){
	if(verbose_opt[0])
	  cout << "***update " << relobsindex[obsindex] << " = " << modindex << " " << observation_opt[obsindex] << " ***" << endl;

	imgReaderObs.open(observation_opt[obsindex]);
	imgReaderObs.getGeoTransform(geotransform);
	imgReaderObs.setNoData(obsnodata_opt);
	imgReaderObs.setOffset(obsoffset_opt[0]);
	imgReaderObs.setScale(obsscale_opt[0]);
	//calculate regression between model and observation
	if(verbose_opt[0])
	  cout << "Calculating regression for " << imgReaderModel2.getFileName() << " " << imgReaderObs.getFileName() << endl;
	errObs=imgreg.getRMSE(imgReaderModel2,imgReaderObs,c0obs,c1obs,verbose_opt[0]);
	if(verbose_opt[0])
	  cout << "c0obs, c1obs: " << c0obs << ", " << c1obs << endl;
      }
      //prediction (also to fill cloudy pixels in update mode)
      string input;
      if(outputbw_opt.size()==model_opt.size())
	input=outputbw_opt[modindex+1];
      else{
	ostringstream outputstream;
	outputstream << outputbw_opt[0] << "_" << tmodel_opt[modindex+1] << ".tif";
	input=outputstream.str();
      }
      ImgReaderGdal imgReaderEst(input);
      imgReaderEst.setNoData(obsnodata_opt);
      imgReaderEst.setOffset(obsoffset_opt[0]);
      imgReaderEst.setScale(obsscale_opt[0]);

      vector<double> obsLineBuffer;
      vector<double> obsWindowBuffer;//buffer for observation to calculate average corresponding to model pixel
      vector<double> model1LineBuffer;
      vector<double> model2LineBuffer;
      vector<double> model1buffer;//buffer for model 1 to calculate time regression based on window
      vector<double> model2buffer;//buffer for model 2 to calculate time regression based on window
      vector<double> uncertObsLineBuffer;
      vector<double> estReadBuffer;
      // vector<double> estWindowBuffer;//buffer for estimate to calculate average corresponding to model pixel
      vector<double> uncertReadBuffer;
      vector<double> estWriteBuffer(ncol);
      vector<double> uncertWriteBuffer(ncol);

      for(int irow=0;irow<imgWriterEst.nrOfRow();++irow){
	assert(irow<imgReaderEst.nrOfRow());
	//not needed here, because we read entire window for each pixel...
	imgReaderEst.readData(estReadBuffer,GDT_Float64,irow,0);
	imgReaderEst.readData(uncertReadBuffer,GDT_Float64,irow,1);
	//read model2 in case current estimate is nodata
	imgReaderEst.image2geo(0,irow,x,y);
	imgReaderModel2.geo2image(x,y,modCol,modRow);
	assert(modRow>=0&&modRow<imgReaderModel2.nrOfRow());
	imgReaderModel2.readData(model2LineBuffer,GDT_Float64,modRow,0);

	imgReaderModel1.geo2image(x,y,modCol,modRow);
	assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	imgReaderModel1.readData(model1LineBuffer,GDT_Float64,modRow,0);

	if(update){
	  imgReaderObs.readData(obsLineBuffer,GDT_Float64,irow,0);
	  if(imgReaderObs.nrOfBand()>1)
	    imgReaderObs.readData(uncertObsLineBuffer,GDT_Float64,irow,1);
	}
	for(int icol=0;icol<imgWriterEst.nrOfCol();++icol){
	  int minCol=(icol>down_opt[0]/2) ? icol-down_opt[0]/2 : 0;
	  int maxCol=(icol+down_opt[0]/2<imgReaderEst.nrOfCol()) ? icol+down_opt[0]/2 : imgReaderEst.nrOfCol()-1;
	  int minRow=(irow>down_opt[0]/2) ? irow-down_opt[0]/2 : 0;
	  int maxRow=(irow+down_opt[0]/2<imgReaderEst.nrOfRow()) ? irow+down_opt[0]/2 : imgReaderEst.nrOfRow()-1;
	  if(update)
	    imgReaderObs.readDataBlock(obsWindowBuffer,GDT_Float64,minCol,maxCol,minRow,maxRow,0);
	  // imgReaderEst.readDataBlock(estWindowBuffer,GDT_Float64,minCol,maxCol,minRow,maxRow,0);
	  double estValue=estReadBuffer[icol];
	  // double estValue=estWindowBuffer[estWindowBuffer.size()/2];
	  double estMeanValue=0;//stat.mean(estWindowBuffer);
	  double nvalid=0;
	  //time update
	  imgReaderEst.image2geo(icol,irow,x,y);
	  imgReaderModel2.geo2image(x,y,modCol,modRow);
	  assert(modRow>=0&&modRow<imgReaderModel2.nrOfRow());
	  double modValue=model2LineBuffer[modCol];
	  bool estNodata=imgReaderEst.isNoData(estValue);
	  //todo: debug checkDiff_opt 
	  // if(checkDiff_opt[0]){
	  //   vector<double>::iterator itwin=estWindowBuffer.begin();
	  //   while(itwin!=estWindowBuffer.end()){
	  //     if(!imgReaderEst.isNoData(*itwin)){
	  // 	estMeanValue+=*itwin;
	  // 	++nvalid;
	  //     }
	  //     ++itwin;
	  //   }
	  //   estMeanValue/=nvalid;
	  //   double difference=(estMeanValue-c0obs)/c1obs-model1LineBuffer[modCol];
	  //   difference*=difference;
	  //   estNodata=estNodata||(sqrt(difference)>uncertModel_opt[0]*stdDev);
	  // }
	  if(estNodata){
	    //we have not found any valid data yet, better here to take the current model value if valid
	    if(imgReaderModel2.isNoData(modValue)){//if both estimate and model are no-data, set obs to nodata
	      estWriteBuffer[icol]=obsnodata_opt[0];
	      uncertWriteBuffer[icol]=uncertNodata_opt[0];
	    }
	    else{
	      estWriteBuffer[icol]=modValue;
	      uncertWriteBuffer[icol]=uncertModel_opt[0]*stdDev;
	    }
	  }	  
	  else{
	    if(window_opt[0]>0){
	      try{
		// imgReaderModel2.geo2image(x,y,modCol,modRow);//did that already
		minCol=(modCol>window_opt[0]/2) ? modCol-window_opt[0]/2 : 0;
		maxCol=(modCol+window_opt[0]/2<imgReaderModel2.nrOfCol()) ? modCol+window_opt[0]/2 : imgReaderModel2.nrOfCol()-1;
		minRow=(modRow>window_opt[0]/2) ? modRow-window_opt[0]/2 : 0;
		maxRow=(modRow+window_opt[0]/2<imgReaderModel2.nrOfRow()) ? modRow+window_opt[0]/2 : imgReaderModel2.nrOfRow()-1;
		imgReaderModel2.readDataBlock(model2buffer,GDT_Float64,minCol,maxCol,minRow,maxRow,0);

		imgReaderModel1.geo2image(x,y,modCol,modRow);
		assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
		minCol=(modCol>window_opt[0]/2) ? modCol-window_opt[0]/2 : 0;
		maxCol=(modCol+window_opt[0]/2<imgReaderModel1.nrOfCol()) ? modCol+window_opt[0]/2 : imgReaderModel1.nrOfCol()-1;
		minRow=(modRow>window_opt[0]/2) ? modRow-window_opt[0]/2 : 0;
		maxRow=(modRow+window_opt[0]/2<imgReaderModel1.nrOfRow()) ? modRow+window_opt[0]/2 : imgReaderModel1.nrOfRow()-1;
		imgReaderModel1.readDataBlock(model1buffer,GDT_Float64,minCol,maxCol,minRow,maxRow,0);
		// imgReaderEst.image2geo(icol,irow,x,y);
	      }
	      catch(string errorString){
		cerr << "Error reading data block for " << minCol << "-" << maxCol << ", " << minRow << "-" << maxRow << endl;
	      }
	      //erase no-data from buffer
	      vector<double>::iterator it1=model1buffer.begin();
	      vector<double>::iterator it2=model2buffer.begin();
	      while(it1!=model1buffer.end()&&it2!=model2buffer.end()){
		//erase nodata
		bool modNodata=false;
		modNodata=modNodata||imgReaderModel1.isNoData(*it1);
		modNodata=modNodata||imgReaderModel2.isNoData(*it2);
		//todo: debug checkDiff_opt 
		// if(checkDiff_opt[0]){
		//   double difference=(estMeanValue-c0obs)/c1obs-*it1;
		//   difference*=difference;
		//   modNodata=modNodata||(sqrt(difference)>uncertModel_opt[0]*stdDev);
		// }
		if(modNodata){
		  model1buffer.erase(it1);
		  model2buffer.erase(it2);
		}
		else{
		  ++it1;
		  ++it2;
		}
	      }
	      if(model1buffer.size()>minreg_opt[0]&&model2buffer.size()>minreg_opt[5])
		errMod=stat.linear_regression_err(model1buffer,model2buffer,c0mod,c1mod);
	      else{//use global regression...
		c0mod=c0modGlobal;
		c1mod=c1modGlobal;
	      }
	    }
	    else{
	      c0mod=c0modGlobal;
	      c1mod=c1modGlobal;
	    }
	    double certNorm=(errMod*errMod+errObs*errObs);
	    double certMod=errObs*errObs/certNorm;
	    double certObs=errMod*errMod/certNorm;
	    double regTime=(c0mod+c1mod*estValue)*certMod;

	    // double regSensor=(c0obs+c1obs*estValue)*certObs;
	    double regSensor=(c0obs+c1obs*modValue)*certObs;
	    estWriteBuffer[icol]=regTime+regSensor;
	    double totalUncertainty=0;
	    if(errMod<eps_opt[0])
	      totalUncertainty=errObs;
	    else if(errObs<eps_opt[0])
	      totalUncertainty=errMod;
	    else{
	      totalUncertainty=1.0/errMod/errMod+1/errObs/errObs;
	      totalUncertainty=sqrt(1.0/totalUncertainty);
	    }
	    uncertWriteBuffer[icol]=totalUncertainty+uncertReadBuffer[icol];
	  }
	  //observation update
	  if(update&&!imgReaderObs.isNoData(obsLineBuffer[icol])){
	    bool doUpdate=true;
	    if(deltaObs_opt.size()){
	      statfactory::StatFactory statobs;
	      statobs.setNoDataValues(obsnodata_opt);
	      double obsMeanValue=statobs.mean(obsWindowBuffer);
	      double difference=(obsMeanValue-c0obs)/c1obs-modValue;
	      difference*=difference;
	      difference/=modValue;//make relative difference
	      difference*=difference;
	      if(deltaObs_opt[0]<0)
		doUpdate=((100*difference)<deltaObs_opt[0]);
	      else{
		difference*=difference;
		doUpdate=(100*sqrt(difference)<deltaObs_opt[0]);
	      }
	    }
	    if(doUpdate){
	      double kalmanGain=1;
	      double uncertObs=uncertObs_opt[0];
	      if(uncertObsLineBuffer.size()>icol)
		uncertObs=uncertObsLineBuffer[icol];
	      if((uncertWriteBuffer[icol]+uncertObs)>eps_opt[0])
		kalmanGain=uncertWriteBuffer[icol]/(uncertWriteBuffer[icol]+uncertObs);
	      assert(kalmanGain<=1);
	      estWriteBuffer[icol]+=kalmanGain*(obsLineBuffer[icol]-estWriteBuffer[icol]);
	      uncertWriteBuffer[icol]*=(1-kalmanGain);
	    }
	  }
	}
	imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,0);
	imgWriterEst.writeData(uncertWriteBuffer,GDT_Float64,irow,1);
	progress=static_cast<float>((irow+1.0)/imgWriterEst.nrOfRow());
	pfnProgress(progress,pszMessage,pProgressArg);
      }

      imgWriterEst.close();
      imgReaderEst.close();

      if(update){
	imgReaderObs.close();
	--obsindex;
      }
      imgReaderModel1.close();
      imgReaderModel2.close();
    }
  }
  if(find(direction_opt.begin(),direction_opt.end(),"smooth")!=direction_opt.end()){
    ///////////////////////////// smooth model /////////////////////////
    obsindex=0;
    if(verbose_opt[0])
      cout << "Running smooth model" << endl;
    for(int modindex=0;modindex<model_opt.size();++modindex){
      if(verbose_opt[0]){
	cout << "processing time " << tmodel_opt[modindex] << endl;
	if(obsindex<relobsindex.size())
	  cout << "next observation " << tmodel_opt[relobsindex[obsindex]] << endl;
	else
	  cout << "There is no next observation" << endl;
      }
      string output;
      if(outputfb_opt.size()==model_opt.size())
	output=outputfb_opt[modindex];
      else{
	ostringstream outputstream;
	outputstream << outputfb_opt[0] << "_" << tmodel_opt[modindex] << ".tif";
	output=outputstream.str();
      }
    
      //two band output band0=estimation, band1=uncertainty
      imgWriterEst.open(output,ncol,nrow,2,GDT_Float32,imageType,option_opt);
      imgWriterEst.setProjectionProj4(projection_opt[0]);
      imgWriterEst.setGeoTransform(geotransform);
      imgWriterEst.GDALSetNoDataValue(obsnodata_opt[0]);

      //open forward and backward estimates
      //we assume forward in model and backward in observation...

      string inputfw;
      string inputbw;
      if(outputfw_opt.size()==model_opt.size())
	inputfw=outputfw_opt[modindex];
      else{
	ostringstream outputstream;
	outputstream << outputfw_opt[0] << "_" << tmodel_opt[modindex] << ".tif";
	inputfw=outputstream.str();
      }
      if(outputbw_opt.size()==model_opt.size())
	inputbw=outputbw_opt[modindex];
      else{
	ostringstream outputstream;
	outputstream << outputbw_opt[0] << "_" << tmodel_opt[modindex] << ".tif";
	inputbw=outputstream.str();
      }
      ImgReaderGdal imgReaderForward(inputfw);
      ImgReaderGdal imgReaderBackward(inputbw);
      imgReaderForward.setNoData(obsnodata_opt);
      imgReaderForward.setOffset(obsoffset_opt[0]);
      imgReaderForward.setScale(obsscale_opt[0]);
      imgReaderBackward.setNoData(obsnodata_opt);
      imgReaderBackward.setOffset(obsoffset_opt[0]);
      imgReaderBackward.setScale(obsscale_opt[0]);
    
      vector<double> estForwardBuffer;
      vector<double> estBackwardBuffer;
      vector<double> uncertObsLineBuffer;
      vector<double> uncertForwardBuffer;
      vector<double> uncertBackwardBuffer;
      vector<double> uncertReadBuffer;
      vector<double> estWriteBuffer(ncol);
      vector<double> uncertWriteBuffer(ncol);

      bool update=false;
      if(obsindex<relobsindex.size()){
	update=(relobsindex[obsindex]==modindex);
      }

      if(update){
	if(verbose_opt[0])
	  cout << "***update " << relobsindex[obsindex] << " = " << modindex << " " << observation_opt[obsindex] << " ***" << endl;
	imgReaderObs.open(observation_opt[obsindex]);
	imgReaderObs.getGeoTransform(geotransform);
	imgReaderObs.setNoData(obsnodata_opt);
	imgReaderObs.setOffset(obsoffset_opt[0]);
	imgReaderObs.setScale(obsscale_opt[0]);
	//calculate regression between model and observation
      }

      pfnProgress(progress,pszMessage,pProgressArg);

      for(int irow=0;irow<imgWriterEst.nrOfRow();++irow){
	assert(irow<imgReaderForward.nrOfRow());
	assert(irow<imgReaderBackward.nrOfRow());
	imgReaderForward.readData(estForwardBuffer,GDT_Float64,irow,0);
	imgReaderBackward.readData(estBackwardBuffer,GDT_Float64,irow,0);
	imgReaderForward.readData(uncertForwardBuffer,GDT_Float64,irow,1);
	imgReaderBackward.readData(uncertBackwardBuffer,GDT_Float64,irow,1);

	if(update){
	  imgReaderObs.readData(estWriteBuffer,GDT_Float64,irow,0);
	  if(imgReaderObs.nrOfBand()>1)
	    imgReaderObs.readData(uncertObsLineBuffer,GDT_Float64,irow,1);
	}

	for(int icol=0;icol<imgWriterEst.nrOfCol();++icol){
	  double A=estForwardBuffer[icol];
	  double B=estBackwardBuffer[icol];
	  double C=uncertForwardBuffer[icol]*uncertForwardBuffer[icol];
	  double D=uncertBackwardBuffer[icol]*uncertBackwardBuffer[icol];
	  double uncertObs=uncertObs_opt[0];

	  // if(update){//check for nodata in observation
	  //   if(imgReaderObs.isNoData(estWriteBuffer[icol]))
	  //     uncertObs=uncertNodata_opt[0];
	  //   else if(uncertObsLineBuffer.size()>icol)
	  //     uncertObs=uncertObsLineBuffer[icol];
	  // }

	  double noemer=(C+D);
	  //todo: consistently check for division by zero...
	  if(imgReaderForward.isNoData(A)&&imgReaderBackward.isNoData(B)){
	    estWriteBuffer[icol]=obsnodata_opt[0];
	    uncertWriteBuffer[icol]=uncertNodata_opt[0];
	  }
	  else if(imgReaderForward.isNoData(A)){
	    estWriteBuffer[icol]=B;
	    uncertWriteBuffer[icol]=uncertBackwardBuffer[icol];
	  }
	  else if(imgReaderForward.isNoData(B)){
	    estWriteBuffer[icol]=A;
	    uncertWriteBuffer[icol]=uncertForwardBuffer[icol];
	  }
	  else{
	    if(noemer<eps_opt[0]){//simple average if both uncertainties are ~>0
	      estWriteBuffer[icol]=0.5*(A+B);
	      uncertWriteBuffer[icol]=uncertObs;
	    }
	    else{
	      estWriteBuffer[icol]=(A*D+B*C)/noemer;
	      double P=0;
	      if(C>eps_opt[0])
		P+=1.0/C;
	      if(D>eps_opt[0])
		P+=1.0/D;
	      if(uncertObs*uncertObs>eps_opt[0])
		P-=1.0/uncertObs/uncertObs;
	      if(P>eps_opt[0])
		P=sqrt(1.0/P);
	      else
		P=0;
	      uncertWriteBuffer[icol]=P;
	    }
	  }
	}
	imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,0);
	imgWriterEst.writeData(uncertWriteBuffer,GDT_Float64,irow,1);
	progress=static_cast<float>((irow+1.0)/imgWriterEst.nrOfRow());
	pfnProgress(progress,pszMessage,pProgressArg);
      }

      imgWriterEst.close();
      imgReaderForward.close();
      imgReaderBackward.close();
      if(update){
	imgReaderObs.close();
	++obsindex;
      }
    }
  }
}

