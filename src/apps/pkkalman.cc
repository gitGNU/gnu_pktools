/**********************************************************************
pkkalman.cc: produce kalman filtered raster time series
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

/******************************************************************************/
/*! \page pkkalman pkkalman
produce kalman filtered raster time series
## SYNOPSIS

<code>
  
</code>

\section pkkalman_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options

|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | dir    | direction            | std::string | forward |direction to run model (forward\|backward\|smooth) | 
 | mod    | model                | std::string |       |model input datasets, e.g., MODIS (use: -mod model1 -mod model2 etc. | 
 | obs    | observation          | std::string |       |observation input datasets, e.g., landsat (use: -obs obs1 -obs obs2 etc. | 
 | tmod   | tmodel               | int  |       |time sequence of model input. Sequence must have exact same length as model input. Leave empty to have default sequence 0,1,2,etc. | 
 | tobs   | tobservation         | int  |       |time sequence of observation input. Sequence must have exact same length as observation input) | 
 | a_srs  | a_srs                | std::string |       |Override the projection for the output file (leave blank to copy from input file, use epsg:3035 to use European projection and force to European grid | 
 | ofw    | outputfw             | std::string |       |Output raster dataset for forward model | 
 | obw    | outputbw             | std::string |       |Output raster dataset for backward model | 
 | ofb    | outputfb             | std::string |       |Output raster dataset for smooth model | 
 | modnodata | modnodata            | double | 0     |invalid value for model input | 
 | obsnodata | obsnodata         | double | 0     |invalid value for observation input | 
 | obsmin | obsmin               | double |      |Minimum value for observation data | 
 | obsmax | obsmax               | double |      |Maximum value for observation data | 
 | eps    | eps                  | double | 1e-05 |epsilon for non zero division | 
 | um     | uncertmodel          | double | 2     |Multiply this value with std dev of first model image to obtain uncertainty of model | 
 | uo     | uncertobs            | double | 0     |Uncertainty of valid observations | 
 | unodata | uncertnodata         | double | 10000 |Uncertainty in case of no-data values in observation | 
 | q      | q                    | double | 1     |Process noise: expresses instability (variance) of proportions of fine res pixels within a moderate resolution pixel | 
 | down   | down                 | int  |       |Downsampling factor for reading model data to calculate regression | 
 | co     | co                   | std::string |       |Creation option for output file. Multiple options can be specified. | 
 | v      | verbose              | short | 0     |verbose mode when positive | 

**/

using namespace std;
/*------------------
  Main procedure
  ----------------*/
int main(int argc,char **argv) {
  Optionpk<string> direction_opt("dir","direction","direction to run model (forward|backward|smooth)","forward");
  Optionpk<string> model_opt("mod","model","model input datasets, e.g., MODIS (use: -mod model1 -mod model2 etc.)");
  Optionpk<string> observation_opt("obs","observation","observation input datasets, e.g., landsat (use: -obs obs1 -obs obs2 etc.)");
  Optionpk<int> tmodel_opt("tmod","tmodel","time sequence of model input. Sequence must have exact same length as model input. Leave empty to have default sequence 0,1,2,etc."); 
  Optionpk<int> tobservation_opt("tobs","tobservation","time sequence of observation input. Sequence must have exact same length as observation input)"); 
  Optionpk<string>  projection_opt("a_srs", "a_srs", "Override the projection for the output file (leave blank to copy from input file, use epsg:3035 to use European projection and force to European grid");
  Optionpk<string> outputfw_opt("ofw", "outputfw", "Output raster dataset for forward model");
  Optionpk<string> outputbw_opt("obw", "outputbw", "Output raster dataset for backward model");
  Optionpk<string> outputfb_opt("ofb", "outputfb", "Output raster dataset for smooth model");
  Optionpk<string> gain_opt("gain", "gain", "Output raster dataset for gain");
  
  Optionpk<double> modnodata_opt("modnodata", "modnodata", "invalid value for model input", 0);
  Optionpk<double> obsnodata_opt("obsnodata", "obsnodata", "invalid value for observation input", 0);
  Optionpk<double> obsmin_opt("obsmin", "obsmin", "Minimum value for observation data");
  Optionpk<double> obsmax_opt("obsmax", "obsmax", "Maximum value for observation data");
  Optionpk<double> eps_opt("eps", "eps", "epsilon for non zero division", 0.00001);
  Optionpk<double> uncertModel_opt("um", "uncertmodel", "Multiplication factor for uncertainty of model",1,1);
  Optionpk<double> uncertObs_opt("uo", "uncertobs", "Multiplication factor for uncertainty of valid observations",1,1);
  Optionpk<double> processNoise_opt("q", "q", "Process noise: expresses instability (variance) of proportions of fine res pixels within a moderate resolution pixel",1);
  Optionpk<double> uncertNodata_opt("unodata", "uncertnodata", "Uncertainty in case of no-data values in observation", 10000);
  Optionpk<int> down_opt("down", "down", "Downsampling factor for reading model data to calculate regression");
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
    gain_opt.retrieveOption(argc,argv);
    modnodata_opt.retrieveOption(argc,argv);
    obsnodata_opt.retrieveOption(argc,argv);
    obsmin_opt.retrieveOption(argc,argv);
    obsmax_opt.retrieveOption(argc,argv);
    eps_opt.retrieveOption(argc,argv);
    uncertModel_opt.retrieveOption(argc,argv);
    uncertObs_opt.retrieveOption(argc,argv);
    processNoise_opt.retrieveOption(argc,argv);
    uncertNodata_opt.retrieveOption(argc,argv);
    down_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    std::cerr << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }

  if(down_opt.empty()){
    std::cerr << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
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
  ImgReaderGdal imgReaderModel1;
  ImgReaderGdal imgReaderModel2;
  ImgReaderGdal imgReaderObs;
  ImgWriterGdal imgWriterEst;
  //test
  ImgWriterGdal imgWriterGain;

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

  if(down_opt.empty()){
    imgReaderModel1.open(model_opt[0]);
    double resModel=imgReaderModel1.getDeltaX();
    double resObs=imgReaderObs.getDeltaX();
    int down=static_cast<int>(ceil(resModel/resObs));
    if(!(down%2))
      down+=1;
    down_opt.push_back(down);
    imgReaderModel1.close();
  }
  imgReaderObs.close();

  int obsindex=0;

  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;

  double errObs=uncertNodata_opt[0];//start with high initial value in case we do not have first observation at time 0

  vector<int> relobsindex;

  for(int tindex=0;tindex<tobservation_opt.size();++tindex){
    vector<int>::iterator modit;
    modit=upper_bound(tmodel_opt.begin(),tmodel_opt.end(),tobservation_opt[tindex]);
    int relpos=modit-tmodel_opt.begin()-1;
    assert(relpos>=0);//todo: for now, we assume model is available at time before first measurement
    relobsindex.push_back(relpos);
    if(verbose_opt[0])
      cout << "observation " << tindex << ": " << "relative position in model time series is " << relpos << ", date of observation is (tobservation_opt[tindex]): " << tobservation_opt[tindex] << ", relobsindex.back(): " << relobsindex.back() << ", filename observation: " << observation_opt[tindex] << ", filename of corresponding model: " << model_opt[relpos] << endl;
  }

  int ndigit=log(1.0*tmodel_opt.back())/log(10.0)+1;

  double geox=0;
  double geoy=0;

  if(find(direction_opt.begin(),direction_opt.end(),"forward")!=direction_opt.end()){
    ///////////////////////////// forward model /////////////////////////
    cout << "Running forward model" << endl;
    obsindex=0;
    //initialization
    string output;
    if(outputfw_opt.size()==model_opt.size()){
      output=outputfw_opt[0];
    }
    else{
      ostringstream outputstream;
      outputstream << outputfw_opt[0] << "_";
      outputstream << setfill('0') << setw(ndigit) << tmodel_opt[0];
      outputstream << ".tif";
      output=outputstream.str();
    }
    if(verbose_opt[0])
      cout << "Opening image " << output << " for writing " << endl;

    imgWriterEst.open(output,ncol,nrow,2,GDT_Float64,imageType,option_opt);
    imgWriterEst.setProjectionProj4(projection_opt[0]);
    imgWriterEst.setGeoTransform(geotransform);
    imgWriterEst.GDALSetNoDataValue(obsnodata_opt[0]);

    //test
    if(gain_opt.size()){
      imgWriterGain.open(gain_opt[0],ncol,nrow,model_opt.size(),GDT_Float64,imageType,option_opt);
      imgWriterGain.setProjectionProj4(projection_opt[0]);
      imgWriterGain.setGeoTransform(geotransform);
      imgWriterGain.GDALSetNoDataValue(obsnodata_opt[0]);
    }

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
    double modRow=0;
    double modCol=0;
    double lowerCol=0;
    double upperCol=0;
    RESAMPLE theResample=BILINEAR;

    if(relobsindex[0]>0){//initialize output_opt[0] as model[0]
      //write first model as output
      if(verbose_opt[0])
	cout << "write first model as output" << endl;
      for(int jrow=0;jrow<nrow;jrow+=down_opt[0]){
	vector<double> estReadBuffer;
	vector<double> estWriteBuffer(ncol);
	vector<double> uncertWriteBuffer(ncol);
	//test
	vector<double> gainWriteBuffer(ncol);
	try{
	  for(int irow=jrow;irow<jrow+down_opt[0];++irow){
	    imgWriterEst.image2geo(0,irow,geox,geoy);
	    imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
	    if(modRow<0||modRow>=imgReaderModel1.nrOfRow()){
	      cerr << "Error: geo coordinates (" << geox << "," << geoy << ") not covered in model image " << imgReaderModel1.getFileName() << endl;
	      assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	    }
	    imgReaderModel1.readData(estReadBuffer,GDT_Float64,modRow,0,theResample);
	    for(int jcol=0;jcol<ncol;jcol+=down_opt[0]){
	      for(int icol=icol;icol<icol+down_opt[0];++icol){
		imgWriterEst.image2geo(icol,irow,geox,geoy);
		imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
		lowerCol=modCol-0.5;
		lowerCol=static_cast<int>(lowerCol);
		upperCol=modCol+0.5;
		upperCol=static_cast<int>(upperCol);
		if(lowerCol<0)
		  lowerCol=0;
		if(upperCol>=imgReaderModel1.nrOfCol())
		  upperCol=imgReaderModel1.nrOfCol()-1;
		double modValue=(modCol-0.5-lowerCol)*estReadBuffer[upperCol]+(1-modCol+0.5+lowerCol)*estReadBuffer[lowerCol];
		// double modValue=estReadBuffer[modCol];
		if(imgReaderModel1.isNoData(modValue)){
		  estWriteBuffer[icol]=obsnodata_opt[0];
		  uncertWriteBuffer[icol]=uncertNodata_opt[0];
		  gainWriteBuffer[icol]=obsnodata_opt[0];
		}
		else{
		  estWriteBuffer[icol]=modValue;
		  if(obsmin_opt.size()){
		    if(estWriteBuffer[icol]<obsmin_opt[0])
		      estWriteBuffer[icol]=obsmin_opt[0];
		  }
		  if(obsmax_opt.size()){
		    if(estWriteBuffer[icol]>obsmax_opt[0])
		      estWriteBuffer[icol]=obsmax_opt[0];
		  }
		  uncertWriteBuffer[icol]=uncertModel_opt[0]*stdDev*stdDev;
		  gainWriteBuffer[icol]=0;
		}
	      }
	    }
	    imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,0);
	    imgWriterEst.writeData(uncertWriteBuffer,GDT_Float64,irow,1);
	    if(gain_opt.size())
	      imgWriterGain.writeData(gainWriteBuffer,GDT_Float64,irow,0);
	  }
	}
	catch(string errorString){
	  cerr << errorString << endl;
	}
	catch(...){
	  cerr << "Error writing file " << imgWriterEst.getFileName() << endl;
	}
      }
    }
    else{//we have a measurement
      if(verbose_opt[0])
	cout << "we have a measurement at initial time" << endl;
      imgReaderObs.open(observation_opt[0]);
      imgReaderObs.getGeoTransform(geotransform);
      imgReaderObs.setNoData(obsnodata_opt);

      vector< vector<double> > obsLineVector(down_opt[0]);
      vector<double> obsLineBuffer;
      vector<double> obsWindowBuffer;//buffer for observation to calculate average corresponding to model pixel
      vector<double> estReadBuffer;
      vector<double> estWriteBuffer(ncol);
      vector<double> uncertWriteBuffer(ncol);
      vector<double> uncertObsLineBuffer;
      //test
      vector<double> gainWriteBuffer(ncol);

      if(verbose_opt[0])
	cout << "initialize obsLineVector" << endl;
      assert(down_opt[0]%2);//window size must be odd 
      for(int iline=-down_opt[0]/2;iline<down_opt[0]/2+1;++iline){
	if(iline<0)//replicate line 0
	  imgReaderObs.readData(obsLineVector[iline+down_opt[0]/2],GDT_Float64,0,0);
	else
	  imgReaderObs.readData(obsLineVector[iline+down_opt[0]/2],GDT_Float64,iline,0);
      }
      for(int jrow=0;jrow<nrow;jrow+=down_opt[0]){
	for(int irow=jrow;irow<jrow+down_opt[0];++irow){
	  imgWriterEst.image2geo(0,irow,geox,geoy);
	  imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
	  assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	  imgReaderModel1.readData(estReadBuffer,GDT_Float64,modRow,0,theResample);
	  int maxRow=(irow+down_opt[0]/2<imgReaderObs.nrOfRow()) ? irow+down_opt[0]/2 : imgReaderObs.nrOfRow()-1;
	  obsLineVector.erase(obsLineVector.begin());
	  imgReaderObs.readData(obsLineBuffer,GDT_Float64,maxRow,0);
	  obsLineVector.push_back(obsLineBuffer);
	  // obsLineBuffer=obsLineVector[down_opt[0]/2];
	  if(imgReaderObs.nrOfBand()>1)
	    imgReaderObs.readData(uncertObsLineBuffer,GDT_Float64,irow,1);

	  for(int jcol=0;jcol<ncol;jcol+=down_opt[0]){
	    for(int icol=jcol;icol<jcol+down_opt[0];++icol){
	      imgWriterEst.image2geo(icol,irow,geox,geoy);
	      imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
	      assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	      lowerCol=modCol-0.5;
	      lowerCol=static_cast<int>(lowerCol);
	      upperCol=modCol+0.5;
	      upperCol=static_cast<int>(upperCol);
	      if(lowerCol<0)
		lowerCol=0;
	      if(upperCol>=imgReaderModel1.nrOfCol())
		upperCol=imgReaderModel1.nrOfCol()-1;
	      double modValue=(modCol-0.5-lowerCol)*estReadBuffer[upperCol]+(1-modCol+0.5+lowerCol)*estReadBuffer[lowerCol];
	      // double modValue=estReadBuffer[modCol];
	      double errMod=uncertModel_opt[0]*stdDev*stdDev;
	      if(imgReaderModel1.isNoData(modValue)){//model is nodata: retain observation 
		if(imgReaderObs.isNoData(obsLineBuffer[icol])){//both model and observation nodata
		  estWriteBuffer[icol]=obsnodata_opt[0];
		  uncertWriteBuffer[icol]=uncertNodata_opt[0];
		  gainWriteBuffer[icol]=obsnodata_opt[0];
		}
		else{
		  estWriteBuffer[icol]=obsLineBuffer[icol];
		  if(obsmin_opt.size()){
		    if(estWriteBuffer[icol]<obsmin_opt[0])
		      estWriteBuffer[icol]=obsmin_opt[0];
		  }
		  if(obsmax_opt.size()){
		    if(estWriteBuffer[icol]>obsmax_opt[0])
		      estWriteBuffer[icol]=obsmax_opt[0];
		  }
		  if(uncertObsLineBuffer.size()>icol)
		    uncertWriteBuffer[icol]=uncertObsLineBuffer[icol];
		  else
		    uncertWriteBuffer[icol]=uncertObs_opt[0];
		}
	      }
	      else{//model is valid: calculate estimate from model
		estWriteBuffer[icol]=modValue;
		uncertWriteBuffer[icol]=errMod;//in case observation is not valid
		gainWriteBuffer[icol]=0;
	      }
	      //measurement update
	      if(!imgReaderObs.isNoData(obsLineBuffer[icol])){
		// estWriteBuffer[icol]=estReadBuffer[icol]*modValue1/modValue2
		double kalmanGain=1;
		int minCol=(icol>down_opt[0]/2) ? icol-down_opt[0]/2 : 0;
		int maxCol=(icol+down_opt[0]/2<imgReaderObs.nrOfCol()) ? icol+down_opt[0]/2 : imgReaderObs.nrOfCol()-1;
		int minRow=(irow>down_opt[0]/2) ? irow-down_opt[0]/2 : 0;
		int maxRow=(irow+down_opt[0]/2<imgReaderObs.nrOfRow()) ? irow+down_opt[0]/2 : imgReaderObs.nrOfRow()-1;
		obsWindowBuffer.clear();
		for(int iline=0;iline<obsLineVector.size();++iline){
		  for(int isample=minCol;isample<=maxCol;++isample){
		    assert(isample<obsLineVector[iline].size());
		    obsWindowBuffer.push_back(obsLineVector[iline][isample]);
		  }
		}
		if(!imgReaderModel1.isNoData(modValue)){//model is valid
		  statfactory::StatFactory statobs;
		  statobs.setNoDataValues(obsnodata_opt);
		  double obsMeanValue=statobs.mean(obsWindowBuffer);
		  double difference=0;
		  difference=obsMeanValue-modValue;
		  errObs=uncertObs_opt[0]*difference*difference;//uncertainty of the observation (R in Kalman equations)
		  double errorCovariance=errMod;//assumed initial errorCovariance (P in Kalman equations)
		  if(errorCovariance+errObs>eps_opt[0])
		    kalmanGain=errorCovariance/(errorCovariance+errObs);
		  else 
		    kalmanGain=1;
		  estWriteBuffer[icol]+=kalmanGain*(obsLineBuffer[icol]-estWriteBuffer[icol]);
		  if(obsmin_opt.size()){
		    if(estWriteBuffer[icol]<obsmin_opt[0])
		      estWriteBuffer[icol]=obsmin_opt[0];
		  }
		  if(obsmax_opt.size()){
		    if(estWriteBuffer[icol]>obsmax_opt[0])
		      estWriteBuffer[icol]=obsmax_opt[0];
		  }
		}
		assert(kalmanGain<=1);
		gainWriteBuffer[icol]=kalmanGain;
	      }
	    }
	  }
	  if(gain_opt.size())
	    imgWriterGain.writeData(gainWriteBuffer,GDT_Float64,irow,0);
	  imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,0);
	  imgWriterEst.writeData(uncertWriteBuffer,GDT_Float64,irow,1);
	}
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
	outputstream << outputfw_opt[0] << "_";
	outputstream << setfill('0') << setw(ndigit) << tmodel_opt[modindex];
	outputstream << ".tif";
	// outputstream << outputfw_opt[0] << "_" << tmodel_opt[modindex] << ".tif";
	output=outputstream.str();
      }
    
      //two band output band0=estimation, band1=uncertainty
      imgWriterEst.open(output,ncol,nrow,2,GDT_Float64,imageType,option_opt);
      imgWriterEst.setProjectionProj4(projection_opt[0]);
      imgWriterEst.setGeoTransform(geotransform);
      imgWriterEst.GDALSetNoDataValue(obsnodata_opt[0]);

      //calculate regression between two subsequence model inputs
      imgReaderModel1.open(model_opt[modindex-1]);
      imgReaderModel1.setNoData(modnodata_opt);
      imgReaderModel2.open(model_opt[modindex]);
      imgReaderModel2.setNoData(modnodata_opt);
      //calculate regression
      //we could re-use the points from second image from last run, but
      //to keep it general, we must redo it (overlap might have changed)
    
      pfnProgress(progress,pszMessage,pProgressArg);

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
      }
      //prediction (also to fill cloudy pixels in measurement update mode)
      string input;
      if(outputfw_opt.size()==model_opt.size())
	input=outputfw_opt[modindex-1];
      else{
	ostringstream outputstream;
	outputstream << outputfw_opt[0] << "_";
	outputstream << setfill('0') << setw(ndigit) << tmodel_opt[modindex-1];
	outputstream << ".tif";
	// outputstream << outputfw_opt[0] << "_" << tmodel_opt[modindex-1] << ".tif";
	input=outputstream.str();
      }
      if(verbose_opt[0])
	cout << "opening " << input << endl;
      ImgReaderGdal imgReaderEst(input);
      imgReaderEst.setNoData(obsnodata_opt);

      vector< vector<double> > obsLineVector(down_opt[0]);
      vector<double> obsLineBuffer;
      vector<double> obsWindowBuffer;//buffer for observation to calculate average corresponding to model pixel
      vector<double> model1LineBuffer;
      vector<double> model2LineBuffer;
      vector<double> model1buffer;//buffer for model 1 to calculate time regression based on window
      vector<double> model2buffer;//buffer for model 2 to calculate time regression based on window
      vector<double> uncertObsLineBuffer;
      vector< vector<double> > estLineVector(down_opt[0]);
      vector<double> estLineBuffer;
      vector<double> estWindowBuffer;//buffer for estimate to calculate average corresponding to model pixel
      vector<double> uncertReadBuffer;
      vector<double> estWriteBuffer(ncol);
      vector<double> uncertWriteBuffer(ncol);
      vector<double> gainWriteBuffer(ncol);

      //initialize obsLineVector if update
      if(update){
	if(verbose_opt[0])
	  cout << "initialize obsLineVector" << endl;
	assert(down_opt[0]%2);//window size must be odd 
	for(int iline=-down_opt[0]/2;iline<down_opt[0]/2+1;++iline){
	  if(iline<0)//replicate line 0
	    imgReaderObs.readData(obsLineVector[iline+down_opt[0]/2],GDT_Float64,0,0);
	  else
	    imgReaderObs.readData(obsLineVector[iline+down_opt[0]/2],GDT_Float64,iline,0);
	}
      }
      //initialize estLineVector
      if(verbose_opt[0])
	cout << "initialize estLineVector" << endl;
      assert(down_opt[0]%2);//window size must be odd 
      for(int iline=-down_opt[0]/2;iline<down_opt[0]/2+1;++iline){
	if(iline<0)//replicate line 0
	  imgReaderEst.readData(estLineVector[iline+down_opt[0]/2],GDT_Float64,0,0);
	else
	  imgReaderEst.readData(estLineVector[iline+down_opt[0]/2],GDT_Float64,iline,0);
      }
      statfactory::StatFactory statobs;
      statobs.setNoDataValues(obsnodata_opt);

      for(int jrow=0;jrow<nrow;jrow+=down_opt[0]){
	//todo: read entire window for uncertReadBuffer...
	for(int irow=jrow;irow<jrow+down_opt[0];++irow){
	  imgReaderEst.readData(uncertReadBuffer,GDT_Float64,irow,1);
	  imgReaderEst.image2geo(0,irow,geox,geoy);
	  imgReaderModel2.geo2image(geox,geoy,modCol,modRow);
	  assert(modRow>=0&&modRow<imgReaderModel2.nrOfRow());
	  imgReaderModel2.readData(model2LineBuffer,GDT_Float64,modRow,0,theResample);

	  imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
	  assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	  imgReaderModel1.readData(model1LineBuffer,GDT_Float64,modRow,0,theResample);

	  int maxRow=(irow+down_opt[0]/2<imgReaderEst.nrOfRow()) ? irow+down_opt[0]/2 : imgReaderEst.nrOfRow()-1;
	  estLineVector.erase(estLineVector.begin());
	  imgReaderEst.readData(estLineBuffer,GDT_Float64,maxRow,0);
	  estLineVector.push_back(estLineBuffer);
	  estLineBuffer=estLineVector[down_opt[0]/2];

	  if(update){
	    int maxRow=(irow+down_opt[0]/2<imgReaderObs.nrOfRow()) ? irow+down_opt[0]/2 : imgReaderObs.nrOfRow()-1;
	    obsLineVector.erase(obsLineVector.begin());
	    imgReaderObs.readData(obsLineBuffer,GDT_Float64,maxRow,0);
	    obsLineVector.push_back(obsLineBuffer);
	    obsLineBuffer=obsLineVector[down_opt[0]/2];
	    // imgReaderObs.readData(obsLineBuffer,GDT_Float64,irow,0);
	    if(imgReaderObs.nrOfBand()>1)
	      imgReaderObs.readData(uncertObsLineBuffer,GDT_Float64,irow,1);
	  }
	  for(int jcol=0;jcol<ncol;jcol+=down_opt[0]){
	    for(int icol=jcol;icol<jcol+down_opt[0];++icol){
	      imgReaderEst.image2geo(icol,irow,geox,geoy);
	      int minCol=(icol>down_opt[0]/2) ? icol-down_opt[0]/2 : 0;
	      int maxCol=(icol+down_opt[0]/2<imgReaderEst.nrOfCol()) ? icol+down_opt[0]/2 : imgReaderEst.nrOfCol()-1;
	      int minRow=(irow>down_opt[0]/2) ? irow-down_opt[0]/2 : 0;
	      int maxRow=(irow+down_opt[0]/2<imgReaderEst.nrOfRow()) ? irow+down_opt[0]/2 : imgReaderEst.nrOfRow()-1;
	      estWindowBuffer.clear();
	      for(int iline=0;iline<estLineVector.size();++iline){
		for(int isample=minCol;isample<=maxCol;++isample){
		  assert(isample<estLineVector[iline].size());
		  estWindowBuffer.push_back(estLineVector[iline][isample]);
		}
	      }
	      if(update){
		obsWindowBuffer.clear();
		for(int iline=0;iline<obsLineVector.size();++iline){
		  for(int isample=minCol;isample<=maxCol;++isample){
		    assert(isample<obsLineVector[iline].size());
		    obsWindowBuffer.push_back(obsLineVector[iline][isample]);
		  }
		}
	      }
	      double estValue=estLineBuffer[icol];
	      imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
	      assert(modRow>=0&&modRow<imgReaderModel2.nrOfRow());
	      lowerCol=modCol-0.5;
	      lowerCol=static_cast<int>(lowerCol);
	      upperCol=modCol+0.5;
	      upperCol=static_cast<int>(upperCol);
	      if(lowerCol<0)
		lowerCol=0;
	      if(upperCol>=imgReaderModel1.nrOfCol())
		upperCol=imgReaderModel1.nrOfCol()-1;
	      double modValue1=(modCol-0.5-lowerCol)*model1LineBuffer[upperCol]+(1-modCol+0.5+lowerCol)*model1LineBuffer[lowerCol];
	      // double modValue1=model1LineBuffer[modCol];
	      imgReaderModel2.geo2image(geox,geoy,modCol,modRow);
	      assert(modRow>=0&&modRow<imgReaderModel2.nrOfRow());
	      lowerCol=modCol-0.5;
	      lowerCol=static_cast<int>(lowerCol);
	      upperCol=modCol+0.5;
	      upperCol=static_cast<int>(upperCol);
	      if(lowerCol<0)
		lowerCol=0;
	      if(upperCol>=imgReaderModel1.nrOfCol())
		upperCol=imgReaderModel1.nrOfCol()-1;
	      double modValue2=(modCol-0.5-lowerCol)*model2LineBuffer[upperCol]+(1-modCol+0.5+lowerCol)*model2LineBuffer[lowerCol];
	      // double modValue2=model2LineBuffer[modCol];
	      if(imgReaderEst.isNoData(estValue)){
		//we have not found any valid data yet, better here to take the current model value if valid
		if(imgReaderModel2.isNoData(modValue2)){//if both estimate and model are no-data, set obs to nodata
		  estWriteBuffer[icol]=obsnodata_opt[0];
		  uncertWriteBuffer[icol]=uncertNodata_opt[0];
		  gainWriteBuffer[icol]=0;
		}
		else{
		  estWriteBuffer[icol]=modValue2;
		  if(obsmin_opt.size()){
		    if(estWriteBuffer[icol]<obsmin_opt[0])
		      estWriteBuffer[icol]=obsmin_opt[0];
		  }
		  if(obsmax_opt.size()){
		    if(estWriteBuffer[icol]>obsmax_opt[0])
		      estWriteBuffer[icol]=obsmax_opt[0];
		  }
		  uncertWriteBuffer[icol]=uncertModel_opt[0]*stdDev*stdDev;
		  gainWriteBuffer[icol]=0;
		}
	      }
	      else{//previous estimate is valid
		double estMeanValue=statobs.mean(estWindowBuffer);
		double nvalid=0;
		//time update
		double processNoiseVariance=processNoise_opt[0];
		//todo: estimate process noise variance expressing instability of weights over time
		//estimate stability of weight distribution from model (low resolution) data in a window mod1 -> mod2 and assume distribution holds at fine spatial resolution. 

		if(imgReaderModel1.isNoData(modValue1)||imgReaderModel2.isNoData(modValue2)){
		  estWriteBuffer[icol]=estValue;
		  uncertWriteBuffer[icol]=uncertReadBuffer[icol]+processNoiseVariance;
		}
		else{//model is good
		  double modRatio=modValue2/modValue1;//transition matrix A in Kalman equations
		  estWriteBuffer[icol]=estValue*modRatio;
		  uncertWriteBuffer[icol]=uncertReadBuffer[icol]*modRatio*modRatio+processNoiseVariance;
		}
		if(obsmin_opt.size()){
		  if(estWriteBuffer[icol]<obsmin_opt[0])
		    estWriteBuffer[icol]=obsmin_opt[0];
		}
		if(obsmax_opt.size()){
		  if(estWriteBuffer[icol]>obsmax_opt[0])
		    estWriteBuffer[icol]=obsmax_opt[0];
		}
	      }
	      //measurement update
	      if(update&&!imgReaderObs.isNoData(obsLineBuffer[icol])){
		double kalmanGain=1;
		if(!imgReaderModel2.isNoData(modValue2)){//model is valid
		  statfactory::StatFactory statobs;
		  statobs.setNoDataValues(obsnodata_opt);
		  double obsMeanValue=statobs.mean(obsWindowBuffer);
		  double difference=0;
		  difference=obsMeanValue-modValue2;
		  errObs=uncertObs_opt[0]*difference*difference;//uncertainty of the observation (R in Kalman equations)

		  if(errObs<eps_opt[0])
		    errObs=eps_opt[0];
		  double errorCovariance=uncertWriteBuffer[icol];//P in Kalman equations

		  if(errorCovariance+errObs>eps_opt[0])
		    kalmanGain=errorCovariance/(errorCovariance+errObs);
		  else 
		    kalmanGain=1;
		  estWriteBuffer[icol]+=kalmanGain*(obsLineBuffer[icol]-estWriteBuffer[icol]);
		  if(obsmin_opt.size()){
		    if(estWriteBuffer[icol]<obsmin_opt[0])
		      estWriteBuffer[icol]=obsmin_opt[0];
		  }
		  if(obsmax_opt.size()){
		    if(estWriteBuffer[icol]>obsmax_opt[0])
		      estWriteBuffer[icol]=obsmax_opt[0];
		  }
		  uncertWriteBuffer[icol]*=(1-kalmanGain);
		}
		assert(kalmanGain<=1);
		gainWriteBuffer[icol]=kalmanGain;
	      }
	    }
	  }
	  //test
	  if(gain_opt.size()){
	    imgWriterGain.writeData(gainWriteBuffer,GDT_Float64,irow,modindex);
	  }
	  imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,0);
	  imgWriterEst.writeData(uncertWriteBuffer,GDT_Float64,irow,1);
	  progress=static_cast<float>((irow+1.0)/imgWriterEst.nrOfRow());
	  pfnProgress(progress,pszMessage,pProgressArg);
	}
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
    if(gain_opt.size())
      imgWriterGain.close();
  }
  if(find(direction_opt.begin(),direction_opt.end(),"backward")!=direction_opt.end()){
    ///////////////////////////// backward model /////////////////////////
    cout << "Running backward model" << endl;
    obsindex=relobsindex.size()-1;
    //initialization
    string output;
    if(outputbw_opt.size()==model_opt.size())
      output=outputbw_opt.back();
    else{
      ostringstream outputstream;
      outputstream << outputbw_opt[0] << "_";
      outputstream << setfill('0') << setw(ndigit) << tmodel_opt.back();
      outputstream << ".tif";
      // outputstream << outputbw_opt[0] << "_" << tmodel_opt.back() << ".tif";
      output=outputstream.str();
    }
    if(verbose_opt[0])
      cout << "Opening image " << output << " for writing " << endl;

    imgWriterEst.open(output,ncol,nrow,2,GDT_Float64,imageType,option_opt);
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
    double modRow=0;
    double modCol=0;
    double lowerCol=0;
    double upperCol=0;
    RESAMPLE theResample=BILINEAR;

    if(relobsindex.back()<model_opt.size()-1){//initialize output_opt.back() as model[0]
      //write last model as output
      if(verbose_opt[0])
	cout << "write last model as output" << endl;
      for(int jrow=0;jrow<nrow;jrow+=down_opt[0]){
	vector<double> estReadBuffer;
	vector<double> estWriteBuffer(ncol);
	vector<double> uncertWriteBuffer(ncol);
	try{
	  for(int irow=jrow;irow<jrow+down_opt[0];++irow){
	    imgWriterEst.image2geo(0,irow,geox,geoy);
	    imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
	    assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	    imgReaderModel1.readData(estReadBuffer,GDT_Float64,modRow,0,theResample);
	    for(int jcol=0;jcol<ncol;jcol+=down_opt[0]){
	      for(int icol=icol;icol<icol+down_opt[0];++icol){
		imgWriterEst.image2geo(icol,irow,geox,geoy);
		imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
		lowerCol=modCol-0.5;
		lowerCol=static_cast<int>(lowerCol);
		upperCol=modCol+0.5;
		upperCol=static_cast<int>(upperCol);
		if(lowerCol<0)
		  lowerCol=0;
		if(upperCol>=imgReaderModel1.nrOfCol())
		  upperCol=imgReaderModel1.nrOfCol()-1;
		double modValue=(modCol-0.5-lowerCol)*estReadBuffer[upperCol]+(1-modCol+0.5+lowerCol)*estReadBuffer[lowerCol];

		// double modValue=estReadBuffer[modCol];
		if(imgReaderModel1.isNoData(modValue)){
		  estWriteBuffer[icol]=obsnodata_opt[0];
		  uncertWriteBuffer[icol]=uncertNodata_opt[0];
		}
		else{
		  estWriteBuffer[icol]=modValue;
		  if(obsmin_opt.size()){
		    if(estWriteBuffer[icol]<obsmin_opt[0])
		      estWriteBuffer[icol]=obsmin_opt[0];
		  }
		  if(obsmax_opt.size()){
		    if(estWriteBuffer[icol]>obsmax_opt[0])
		      estWriteBuffer[icol]=obsmax_opt[0];
		  }
		  uncertWriteBuffer[icol]=uncertModel_opt[0]*stdDev*stdDev;
		}
	      }
	    }
	    imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,0);
	    imgWriterEst.writeData(uncertWriteBuffer,GDT_Float64,irow,1);
	  }
	}
	catch(string errorString){
	  cerr << errorString << endl;
	}
	catch(...){
	  cerr << "Error writing file " << imgWriterEst.getFileName() << endl;
	}
      }
    }
    else{//we have a measurement at end time
      if(verbose_opt[0])
	cout << "we have a measurement at end time" << endl;
      imgReaderObs.open(observation_opt.back());
      imgReaderObs.getGeoTransform(geotransform);
      imgReaderObs.setNoData(obsnodata_opt);
      
      vector< vector<double> > obsLineVector(down_opt[0]);
      vector<double> obsLineBuffer;
      vector<double> obsWindowBuffer;//buffer for observation to calculate average corresponding to model pixel
      vector<double> estReadBuffer;
      vector<double> estWriteBuffer(ncol);
      vector<double> uncertWriteBuffer(ncol);
      vector<double> uncertObsLineBuffer;

      if(verbose_opt[0])
	cout << "initialize obsLineVector" << endl;
      assert(down_opt[0]%2);//window size must be odd 
      for(int iline=-down_opt[0]/2;iline<down_opt[0]/2+1;++iline){
	if(iline<0)//replicate line 0
	  imgReaderObs.readData(obsLineVector[iline+down_opt[0]/2],GDT_Float64,0,0);
	else
	  imgReaderObs.readData(obsLineVector[iline+down_opt[0]/2],GDT_Float64,iline,0);
      }
      for(int jrow=0;jrow<nrow;jrow+=down_opt[0]){
	for(int irow=jrow;irow<jrow+down_opt[0];++irow){
	  imgWriterEst.image2geo(0,irow,geox,geoy);
	  imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
	  assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	  RESAMPLE theResample=BILINEAR;
	  imgReaderModel1.readData(estReadBuffer,GDT_Float64,modRow,0,theResample);
	  int maxRow=(irow+down_opt[0]/2<imgReaderObs.nrOfRow()) ? irow+down_opt[0]/2 : imgReaderObs.nrOfRow()-1;
	  obsLineVector.erase(obsLineVector.begin());
	  imgReaderObs.readData(obsLineBuffer,GDT_Float64,maxRow,0);
	  obsLineVector.push_back(obsLineBuffer);
	  // obsLineBuffer=obsLineVector[down_opt[0]/2];
	  if(imgReaderObs.nrOfBand()>1)
	    imgReaderObs.readData(uncertObsLineBuffer,GDT_Float64,irow,1);

	  for(int jcol=0;jcol<ncol;jcol+=down_opt[0]){
	    for(int icol=jcol;icol<jcol+down_opt[0];++icol){
	      imgWriterEst.image2geo(icol,irow,geox,geoy);
	      imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
	      assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	      lowerCol=modCol-0.5;
	      lowerCol=static_cast<int>(lowerCol);
	      upperCol=modCol+0.5;
	      upperCol=static_cast<int>(upperCol);
	      if(lowerCol<0)
		lowerCol=0;
	      if(upperCol>=imgReaderModel1.nrOfCol())
		upperCol=imgReaderModel1.nrOfCol()-1;
	      double modValue=(modCol-0.5-lowerCol)*estReadBuffer[upperCol]+(1-modCol+0.5+lowerCol)*estReadBuffer[lowerCol];
	      // double modValue=estReadBuffer[modCol];
	      double errMod=uncertModel_opt[0]*stdDev*stdDev;
	      if(imgReaderModel1.isNoData(modValue)){//model is nodata: retain observation 
		if(imgReaderObs.isNoData(obsLineBuffer[icol])){//both model and observation nodata
		  estWriteBuffer[icol]=obsnodata_opt[0];
		  uncertWriteBuffer[icol]=uncertNodata_opt[0];
		}
		else{
		  estWriteBuffer[icol]=obsLineBuffer[icol];
		  if(obsmin_opt.size()){
		    if(estWriteBuffer[icol]<obsmin_opt[0])
		      estWriteBuffer[icol]=obsmin_opt[0];
		  }
		  if(obsmax_opt.size()){
		    if(estWriteBuffer[icol]>obsmax_opt[0])
		      estWriteBuffer[icol]=obsmax_opt[0];
		  }
		  if(uncertObsLineBuffer.size()>icol)
		    uncertWriteBuffer[icol]=uncertObsLineBuffer[icol];
		  else
		    uncertWriteBuffer[icol]=uncertObs_opt[0];
		}
	      }
	      else{//model is valid: calculate estimate from model
		estWriteBuffer[icol]=modValue;
		uncertWriteBuffer[icol]=errMod;//in case observation is not valid
	      }
	      //measurement update
	      if(!imgReaderObs.isNoData(obsLineBuffer[icol])){
		// estWriteBuffer[icol]=estReadBuffer[icol]*modValue1/modValue2
		double kalmanGain=1;
		int minCol=(icol>down_opt[0]/2) ? icol-down_opt[0]/2 : 0;
		int maxCol=(icol+down_opt[0]/2<imgReaderObs.nrOfCol()) ? icol+down_opt[0]/2 : imgReaderObs.nrOfCol()-1;
		int minRow=(irow>down_opt[0]/2) ? irow-down_opt[0]/2 : 0;
		int maxRow=(irow+down_opt[0]/2<imgReaderObs.nrOfRow()) ? irow+down_opt[0]/2 : imgReaderObs.nrOfRow()-1;
		obsWindowBuffer.clear();
		for(int iline=0;iline<obsLineVector.size();++iline){
		  for(int isample=minCol;isample<=maxCol;++isample){
		    assert(isample<obsLineVector[iline].size());
		    obsWindowBuffer.push_back(obsLineVector[iline][isample]);
		  }
		}
		if(!imgReaderModel1.isNoData(modValue)){//model is valid
		  statfactory::StatFactory statobs;
		  statobs.setNoDataValues(obsnodata_opt);
		  double obsMeanValue=statobs.mean(obsWindowBuffer);
		  double difference=0;
		  difference=obsMeanValue-modValue;
		  errObs=uncertObs_opt[0]*difference*difference;//uncertainty of the observation (R in Kalman equations)
		  double errorCovariance=errMod;//assumed initial errorCovariance (P in Kalman equations)
		  if(errorCovariance+errObs>eps_opt[0])
		    kalmanGain=errorCovariance/(errorCovariance+errObs);
		  else 
		    kalmanGain=1;
		  estWriteBuffer[icol]+=kalmanGain*(obsLineBuffer[icol]-estWriteBuffer[icol]);
		  if(obsmin_opt.size()){
		    if(estWriteBuffer[icol]<obsmin_opt[0])
		      estWriteBuffer[icol]=obsmin_opt[0];
		  }
		  if(obsmax_opt.size()){
		    if(estWriteBuffer[icol]>obsmax_opt[0])
		      estWriteBuffer[icol]=obsmax_opt[0];
		  }
		}
		assert(kalmanGain<=1);
	      }
	    }
	  }
	  imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,0);
	  imgWriterEst.writeData(uncertWriteBuffer,GDT_Float64,irow,1);
	}
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
	outputstream << outputbw_opt[0] << "_";
	outputstream << setfill('0') << setw(ndigit) << tmodel_opt[modindex];
	outputstream << ".tif";
	// outputstream << outputbw_opt[0] << "_" << tmodel_opt[modindex] << ".tif";
	output=outputstream.str();
      }

      //two band output band0=estimation, band1=uncertainty
      imgWriterEst.open(output,ncol,nrow,2,GDT_Float64,imageType,option_opt);
      imgWriterEst.setProjectionProj4(projection_opt[0]);
      imgWriterEst.setGeoTransform(geotransform);
      imgWriterEst.GDALSetNoDataValue(obsnodata_opt[0]);

      //calculate regression between two subsequence model inputs
      imgReaderModel1.open(model_opt[modindex+1]);
      imgReaderModel1.setNoData(modnodata_opt);
      imgReaderModel2.open(model_opt[modindex]);
      imgReaderModel2.setNoData(modnodata_opt);
      //calculate regression
      //we could re-use the points from second image from last run, but
      //to keep it general, we must redo it (overlap might have changed)
    
      pfnProgress(progress,pszMessage,pProgressArg);

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
      }
      //prediction (also to fill cloudy pixels in update mode)
      string input;
      if(outputbw_opt.size()==model_opt.size())
	input=outputbw_opt[modindex+1];
      else{
	ostringstream outputstream;
	outputstream << outputbw_opt[0] << "_";
	outputstream << setfill('0') << setw(ndigit) << tmodel_opt[modindex+1];
	outputstream << ".tif";
	// outputstream << outputbw_opt[0] << "_" << tmodel_opt[modindex+1] << ".tif";
	input=outputstream.str();
      }
      if(verbose_opt[0])
	cout << "opening " << input << endl;
      ImgReaderGdal imgReaderEst(input);
      imgReaderEst.setNoData(obsnodata_opt);
      
      vector< vector<double> > obsLineVector(down_opt[0]);
      vector<double> obsLineBuffer;
      vector<double> obsWindowBuffer;//buffer for observation to calculate average corresponding to model pixel
      vector<double> model1LineBuffer;
      vector<double> model2LineBuffer;
      vector<double> model1buffer;//buffer for model 1 to calculate time regression based on window
      vector<double> model2buffer;//buffer for model 2 to calculate time regression based on window
      vector<double> uncertObsLineBuffer;
      vector< vector<double> > estLineVector(down_opt[0]);
      vector<double> estLineBuffer;
      vector<double> estWindowBuffer;//buffer for estimate to calculate average corresponding to model pixel
      vector<double> uncertReadBuffer;
      vector<double> estWriteBuffer(ncol);
      vector<double> uncertWriteBuffer(ncol);

      //initialize obsLineVector
      if(update){
	if(verbose_opt[0])
	  cout << "initialize obsLineVector" << endl;
	assert(down_opt[0]%2);//window size must be odd 
	for(int iline=-down_opt[0]/2;iline<down_opt[0]/2+1;++iline){
	  if(iline<0)//replicate line 0
	    imgReaderObs.readData(obsLineVector[iline+down_opt[0]/2],GDT_Float64,0,0);
	  else
	    imgReaderObs.readData(obsLineVector[iline+down_opt[0]/2],GDT_Float64,iline,0);
	}
      }
      //initialize estLineVector
      if(verbose_opt[0])
	cout << "initialize estLineVector" << endl;
      assert(down_opt[0]%2);//window size must be odd 
      for(int iline=-down_opt[0]/2;iline<down_opt[0]/2+1;++iline){
	if(iline<0)//replicate line 0
	  imgReaderEst.readData(estLineVector[iline+down_opt[0]/2],GDT_Float64,0,0);
	else
	  imgReaderEst.readData(estLineVector[iline+down_opt[0]/2],GDT_Float64,iline,0);
      }
      statfactory::StatFactory statobs;
      statobs.setNoDataValues(obsnodata_opt);

      for(int jrow=0;jrow<nrow;jrow+=down_opt[0]){
	//todo: read entire window for uncertReadBuffer...
	for(int irow=jrow;irow<jrow+down_opt[0];++irow){
	  imgReaderEst.readData(uncertReadBuffer,GDT_Float64,irow,1);
	  imgReaderEst.image2geo(0,irow,geox,geoy);
	  imgReaderModel2.geo2image(geox,geoy,modCol,modRow);
	  assert(modRow>=0&&modRow<imgReaderModel2.nrOfRow());
	  imgReaderModel2.readData(model2LineBuffer,GDT_Float64,modRow,0,theResample);

	  imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
	  assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	  imgReaderModel1.readData(model1LineBuffer,GDT_Float64,modRow,0,theResample);

	  int maxRow=(irow+down_opt[0]/2<imgReaderEst.nrOfRow()) ? irow+down_opt[0]/2 : imgReaderEst.nrOfRow()-1;
	  estLineVector.erase(estLineVector.begin());
	  imgReaderEst.readData(estLineBuffer,GDT_Float64,maxRow,0);
	  estLineVector.push_back(estLineBuffer);
	  estLineBuffer=estLineVector[down_opt[0]/2];

	  if(update){
	    int maxRow=(irow+down_opt[0]/2<imgReaderObs.nrOfRow()) ? irow+down_opt[0]/2 : imgReaderObs.nrOfRow()-1;
	    obsLineVector.erase(obsLineVector.begin());
	    imgReaderObs.readData(obsLineBuffer,GDT_Float64,maxRow,0);
	    obsLineVector.push_back(obsLineBuffer);
	    obsLineBuffer=obsLineVector[down_opt[0]/2];
	    // imgReaderObs.readData(obsLineBuffer,GDT_Float64,irow,0);
	    if(imgReaderObs.nrOfBand()>1)
	      imgReaderObs.readData(uncertObsLineBuffer,GDT_Float64,irow,1);
	  }
	  for(int jcol=0;jcol<ncol;jcol+=down_opt[0]){
	    for(int icol=jcol;icol<jcol+down_opt[0];++icol){
	      imgReaderEst.image2geo(icol,irow,geox,geoy);
	      int minCol=(icol>down_opt[0]/2) ? icol-down_opt[0]/2 : 0;
	      int maxCol=(icol+down_opt[0]/2<imgReaderEst.nrOfCol()) ? icol+down_opt[0]/2 : imgReaderEst.nrOfCol()-1;
	      int minRow=(irow>down_opt[0]/2) ? irow-down_opt[0]/2 : 0;
	      int maxRow=(irow+down_opt[0]/2<imgReaderEst.nrOfRow()) ? irow+down_opt[0]/2 : imgReaderEst.nrOfRow()-1;
	      estWindowBuffer.clear();
	      for(int iline=0;iline<estLineVector.size();++iline){
		for(int isample=minCol;isample<=maxCol;++isample){
		  assert(isample<estLineVector[iline].size());
		  estWindowBuffer.push_back(estLineVector[iline][isample]);
		}
	      }
	      if(update){
		obsWindowBuffer.clear();
		for(int iline=0;iline<obsLineVector.size();++iline){
		  for(int isample=minCol;isample<=maxCol;++isample){
		    assert(isample<obsLineVector[iline].size());
		    obsWindowBuffer.push_back(obsLineVector[iline][isample]);
		  }
		}
	      }
	      double estValue=estLineBuffer[icol];
	      imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
	      assert(modRow>=0&&modRow<imgReaderModel2.nrOfRow());
	      lowerCol=modCol-0.5;
	      lowerCol=static_cast<int>(lowerCol);
	      upperCol=modCol+0.5;
	      upperCol=static_cast<int>(upperCol);
	      if(lowerCol<0)
		lowerCol=0;
	      if(upperCol>=imgReaderModel1.nrOfCol())
		upperCol=imgReaderModel1.nrOfCol()-1;
	      double modValue1=(modCol-0.5-lowerCol)*model1LineBuffer[upperCol]+(1-modCol+0.5+lowerCol)*model1LineBuffer[lowerCol];
	      // double modValue1=model1LineBuffer[modCol];
	      imgReaderModel2.geo2image(geox,geoy,modCol,modRow);
	      assert(modRow>=0&&modRow<imgReaderModel2.nrOfRow());
	      lowerCol=modCol-0.5;
	      lowerCol=static_cast<int>(lowerCol);
	      upperCol=modCol+0.5;
	      upperCol=static_cast<int>(upperCol);
	      if(lowerCol<0)
		lowerCol=0;
	      if(upperCol>=imgReaderModel1.nrOfCol())
		upperCol=imgReaderModel1.nrOfCol()-1;
	      double modValue2=(modCol-0.5-lowerCol)*model2LineBuffer[upperCol]+(1-modCol+0.5+lowerCol)*model2LineBuffer[lowerCol];
	      // double modValue2=model2LineBuffer[modCol];
	      if(imgReaderEst.isNoData(estValue)){
		//we have not found any valid data yet, better here to take the current model value if valid
		if(imgReaderModel2.isNoData(modValue2)){//if both estimate and model are no-data, set obs to nodata
		  estWriteBuffer[icol]=obsnodata_opt[0];
		  uncertWriteBuffer[icol]=uncertNodata_opt[0];
		}
		else{
		  estWriteBuffer[icol]=modValue2;
		  if(obsmin_opt.size()){
		    if(estWriteBuffer[icol]<obsmin_opt[0])
		      estWriteBuffer[icol]=obsmin_opt[0];
		  }
		  if(obsmax_opt.size()){
		    if(estWriteBuffer[icol]>obsmax_opt[0])
		      estWriteBuffer[icol]=obsmax_opt[0];
		  }
		  uncertWriteBuffer[icol]=uncertModel_opt[0]*stdDev*stdDev;
		}
	      }
	      else{//previous estimate is valid
		double estMeanValue=statobs.mean(estWindowBuffer);
		double nvalid=0;
		//time update
		double processNoiseVariance=processNoise_opt[0];
		//todo: estimate process noise variance expressing instabilityof weights over time
		//estimate stability of weight distribution from model (low resolution) data in a window mod1 -> mod2 and assume distribution holds at fine spatial resolution. 

		if(imgReaderModel1.isNoData(modValue1)||imgReaderModel2.isNoData(modValue2)){
		  estWriteBuffer[icol]=estValue;
		  uncertWriteBuffer[icol]=uncertReadBuffer[icol]+processNoiseVariance;
		}
		else{//model is good
		  double modRatio=modValue2/modValue1;//transition matrix A in Kalman equations
		  estWriteBuffer[icol]=estValue*modRatio;
		  uncertWriteBuffer[icol]=uncertReadBuffer[icol]*modRatio*modRatio+processNoiseVariance;
		}
		if(obsmin_opt.size()){
		  if(estWriteBuffer[icol]<obsmin_opt[0])
		    estWriteBuffer[icol]=obsmin_opt[0];
		}
		if(obsmax_opt.size()){
		  if(estWriteBuffer[icol]>obsmax_opt[0])
		    estWriteBuffer[icol]=obsmax_opt[0];
		}
	      }
	      //measurement update
	      if(update&&!imgReaderObs.isNoData(obsLineBuffer[icol])){
		double kalmanGain=1;
		if(!imgReaderModel2.isNoData(modValue2)){//model is valid
		  statfactory::StatFactory statobs;
		  statobs.setNoDataValues(obsnodata_opt);
		  double obsMeanValue=statobs.mean(obsWindowBuffer);
		  double difference=0;
		  difference=obsMeanValue-modValue2;
		  errObs=uncertObs_opt[0]*difference*difference;//uncertainty of the observation (R in Kalman equations)

		  double errorCovariance=uncertWriteBuffer[icol];//P in Kalman equations

		  if(errorCovariance+errObs>eps_opt[0])
		    kalmanGain=errorCovariance/(errorCovariance+errObs);
		  else 
		    kalmanGain=1;
		  estWriteBuffer[icol]+=kalmanGain*(obsLineBuffer[icol]-estWriteBuffer[icol]);
		  if(obsmin_opt.size()){
		    if(estWriteBuffer[icol]<obsmin_opt[0])
		      estWriteBuffer[icol]=obsmin_opt[0];
		  }
		  if(obsmax_opt.size()){
		    if(estWriteBuffer[icol]>obsmax_opt[0])
		      estWriteBuffer[icol]=obsmax_opt[0];
		  }
		  uncertWriteBuffer[icol]*=(1-kalmanGain);
		}
		assert(kalmanGain<=1);
	      }
	    }
	  }
	  imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,0);
	  imgWriterEst.writeData(uncertWriteBuffer,GDT_Float64,irow,1);
	  progress=static_cast<float>((irow+1.0)/imgWriterEst.nrOfRow());
	  pfnProgress(progress,pszMessage,pProgressArg);
	}
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
    cout << "Running smooth model" << endl;
    obsindex=0;
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
	outputstream << outputfb_opt[0] << "_";
	outputstream << setfill('0') << setw(ndigit) << tmodel_opt[modindex];
	outputstream << ".tif";
	// outputstream << outputfb_opt[0] << "_" << tmodel_opt[modindex] << ".tif";
	output=outputstream.str();
      }
    
      //two band output band0=estimation, band1=uncertainty
      imgWriterEst.open(output,ncol,nrow,2,GDT_Float64,imageType,option_opt);
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
	outputstream << outputfw_opt[0] << "_";
	outputstream << setfill('0') << setw(ndigit) << tmodel_opt[modindex];
	outputstream << ".tif";
	// outputstream << outputfw_opt[0] << "_" << tmodel_opt[modindex] << ".tif";
	inputfw=outputstream.str();
      }
      if(outputbw_opt.size()==model_opt.size())
	inputbw=outputbw_opt[modindex];
      else{
	ostringstream outputstream;
	outputstream << outputbw_opt[0] << "_";
	outputstream << setfill('0') << setw(ndigit) << tmodel_opt[modindex];
	outputstream << ".tif";
	// outputstream << outputbw_opt[0] << "_" << tmodel_opt[modindex] << ".tif";
	inputbw=outputstream.str();
      }
      ImgReaderGdal imgReaderForward(inputfw);
      ImgReaderGdal imgReaderBackward(inputbw);
      imgReaderForward.setNoData(obsnodata_opt);
      imgReaderBackward.setNoData(obsnodata_opt);
      
      vector<double> estForwardBuffer;
      vector<double> estBackwardBuffer;
      vector<double> uncertObsLineBuffer;
      vector<double> uncertForwardBuffer;
      vector<double> uncertBackwardBuffer;
      vector<double> uncertReadBuffer;
      vector<double> estWriteBuffer(ncol);
      vector<double> uncertWriteBuffer(ncol);
      // vector<double> lineMask;

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

	// double oldRowMask=-1;//keep track of row mask to optimize number of line readings
	for(int icol=0;icol<imgWriterEst.nrOfCol();++icol){
	  imgWriterEst.image2geo(icol,irow,geox,geoy);
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
	      if(obsmin_opt.size()){
		if(estWriteBuffer[icol]<obsmin_opt[0])
		estWriteBuffer[icol]=obsmin_opt[0];
	      }
	      if(obsmax_opt.size()){
		if(estWriteBuffer[icol]>obsmax_opt[0])
		estWriteBuffer[icol]=obsmax_opt[0];
	      }
	      double P=0;
	      if(C>eps_opt[0])
		P+=1.0/C;
	      if(D>eps_opt[0])
		P+=1.0/D;
	      if(uncertObs*uncertObs>eps_opt[0])
		P-=1.0/uncertObs/uncertObs;
	      if(P>eps_opt[0])
		P=1.0/P;
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
  // if(mask_opt.size())
  //   maskReader.close();
}

