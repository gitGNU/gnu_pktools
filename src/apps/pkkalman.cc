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
#include "imageclasses/ImgUpdaterGdal.h"
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
 | mod    | model                | std::string |       |coarse spatial resolution input datasets(s) used as model. Use either multi-band input (-model multiband_model.tif) or multiple single-band inputs (-mod model1 -mod model2 etc.) | 
 | modmask| modmask              | std::string |       |model mask datasets(s). Must have same dimension as model input. Use either multi-band input or multiple single-band inputs| 
 | obs    | observation          | std::string |       |fine spatial resolution input dataset(s) used as observation. Use either multi-band input (-obs multiband_obs.tif) or multiple single-band inputs (-obs obs1 -obs obs2 etc.) | 
 | obsmask| obsmask              | std::string |       |observation mask dataset(s). Must have same dimension as observation input (use multi-band input or multiple single-band inputs | 
 | tmod   | tmodel               | int  |       |time sequence of model input. Sequence must have exact same length as model input. Leave empty to have default sequence 0,1,2,etc. | 
 | tobs   | tobservation         | int  |       |time sequence of observation input. Sequence must have exact same length as observation input | 
 | a_srs  | a_srs                | std::string |       |Override the projection for the output file (leave blank to copy from input file, use epsg:3035 to use European projection and force to European grid | 
 | ofw    | outputfw             | std::string |       |Output raster dataset for forward model | 
 | u_ofw  | u_outputfw           | std::string |       |Uncertainty output raster dataset for forward model | 
 | obw    | outputbw             | std::string |       |Output raster dataset for backward model | 
 | u_obw  | u_outputbw           | std::string |       |Uncertainty output raster dataset for backward model | 
 | ofb    | outputfb             | std::string |       |Output raster dataset for smooth model | 
 | u_ofb  | u_outputfb           | std::string |       |Uncertainty output raster dataset for smooth model | 
 | modnodata | modnodata         | double | 0     |invalid value for model input | 
 | obsnodata | obsnodata         | double | 0     |invalid value for observation input | 
 | msknodata | msknodata         | float | 0     |Mask value not to consider
 | mskband | mskband             | short | 0     |Mask band to read (0 indexed) | 
 | obsmin | obsmin               | double |      |Minimum value for observation data | 
 | obsmax | obsmax               | double |      |Maximum value for observation data | 
 | eps    | eps                  | double | 1e-05 |epsilon for non zero division | 
 | um     | uncertmodel          | double | 1     |Uncertainty of the model | 
 | uo     | uncertobs            | double | 1     |Uncertainty of valid observations | 
 | unodata | uncertnodata         | double | 100 |Uncertainty in case of no-data values in observation | 
 | q      | q                    | double | 1     |Process noise: expresses instability (variance) of proportions of fine res pixels within a moderate resolution pixel | 
 | down   | down                 | int  |       |Downsampling factor for reading model data to calculate regression (default is ratio between coarse (model) and fine (obs) resolution raster datasets)| 
 | ot     | otype                | std::string |       |Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image | 
 | of     | oformat              | std::string | GTiff |Output image format (see also gdal_translate).| 
 | co     | co                   | std::string |       |Creation option for output file. Multiple options can be specified. | 
 | v      | verbose              | short | 0     |verbose mode when positive | 

**/

using namespace std;
/*------------------
  Main procedure
  ----------------*/
int main(int argc,char **argv) {
  Optionpk<string> direction_opt("dir","direction","direction to run model (forward|backward|smooth)","forward");
  Optionpk<string> model_opt("mod","model","coarse spatial resolution input datasets(s) used as model. Use either multi-band input (-model multiband_model.tif) or multiple single-band inputs (-mod model1 -mod model2 etc.)");
  Optionpk<string> modelmask_opt("modmask","modmask","model mask datasets(s). Must have same dimension as model input. Use either multi-band input or multiple single-band inputs");
  Optionpk<string> observation_opt("obs","observation","fine spatial resolution input dataset(s) used as observation. Use either multi-band input (-obs multiband_obs.tif) or multiple single-band inputs (-obs obs1 -obs obs2 etc.)");
  Optionpk<string> observationmask_opt("obsmask","obsmask","observation mask dataset(s). Must have same dimension as observation input (use multi-band input or multiple single-band inputs");
  Optionpk<int> tmodel_opt("tmod","tmodel","time sequence of model input. Sequence must have exact same length as model input. Leave empty to have default sequence 0,1,2,etc."); 
  Optionpk<int> tobservation_opt("tobs","tobservation","time sequence of observation input. Sequence must have exact same length as observation input)"); 
  Optionpk<string>  projection_opt("a_srs", "a_srs", "Override the projection for the output file (leave blank to copy from input file, use epsg:3035 to use European projection and force to European grid");
  Optionpk<string> outputfw_opt("ofw", "outputfw", "Output raster dataset for forward model");
  Optionpk<string> uncertfw_opt("u_ofw", "u_outputfw", "Uncertainty output raster dataset for forward model");
  Optionpk<string> outputbw_opt("obw", "outputbw", "Output raster dataset for backward model");
  Optionpk<string> uncertbw_opt("u_obw", "u_outputbw", "Uncertainty output raster dataset for backward model");
  Optionpk<string> outputfb_opt("ofb", "outputfb", "Output raster dataset for smooth model");
  Optionpk<string> uncertfb_opt("u_ofb", "u_outputfb", "Uncertainty output raster dataset for smooth model");
  Optionpk<string> gain_opt("gain", "gain", "Output raster dataset for gain");
  Optionpk<double> modnodata_opt("modnodata", "modnodata", "invalid value for model input", 0);
  Optionpk<double> obsnodata_opt("obsnodata", "obsnodata", "invalid value for observation input", 0);
  Optionpk<double> obsmin_opt("obsmin", "obsmin", "Minimum value for observation data");
  Optionpk<double> obsmax_opt("obsmax", "obsmax", "Maximum value for observation data");
  Optionpk<double> msknodata_opt("msknodata", "msknodata", "Mask value not to consider", 0);
  Optionpk<short> mskband_opt("mskband", "mskband", "Mask band to read (0 indexed)", 0);
  Optionpk<double> eps_opt("eps", "eps", "epsilon for non zero division", 0.00001);
  Optionpk<double> uncertModel_opt("um", "uncertmodel", "Uncertainty of model",1);
  Optionpk<double> uncertObs_opt("uo", "uncertobs", "Uncertainty of valid observations",1);
  Optionpk<double> processNoise_opt("q", "q", "Process noise: expresses instability (variance) of proportions of fine res pixels within a moderate resolution pixel",1);
  Optionpk<double> uncertNodata_opt("unodata", "uncertnodata", "Uncertainty in case of no-data values in observation", 100);
  Optionpk<int> down_opt("down", "down", "Downsampling factor for reading model data to calculate regression");
  Optionpk<string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image","");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate).","GTiff",2);
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode when positive", 0);

  observationmask_opt.setHide(1);
  modelmask_opt.setHide(1);
  tmodel_opt.setHide(1);
  tobservation_opt.setHide(1);
  projection_opt.setHide(1);
  uncertfw_opt.setHide(1);
  uncertbw_opt.setHide(1);
  uncertfb_opt.setHide(1);
  obsmin_opt.setHide(1);
  obsmax_opt.setHide(1);
  msknodata_opt.setHide(1);
  mskband_opt.setHide(1);
  eps_opt.setHide(1);
  uncertNodata_opt.setHide(1);
  down_opt.setHide(1);
  otype_opt.setHide(1);
  oformat_opt.setHide(1);
  option_opt.setHide(1);
  verbose_opt.setHide(1);
  gain_opt.setHide(2);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=direction_opt.retrieveOption(argc,argv);
    model_opt.retrieveOption(argc,argv);
    modelmask_opt.retrieveOption(argc,argv);
    observation_opt.retrieveOption(argc,argv);
    observationmask_opt.retrieveOption(argc,argv);
    tmodel_opt.retrieveOption(argc,argv);
    tobservation_opt.retrieveOption(argc,argv);
    projection_opt.retrieveOption(argc,argv);
    outputfw_opt.retrieveOption(argc,argv);
    uncertfw_opt.retrieveOption(argc,argv);
    outputbw_opt.retrieveOption(argc,argv);
    uncertbw_opt.retrieveOption(argc,argv);
    outputfb_opt.retrieveOption(argc,argv);
    uncertfb_opt.retrieveOption(argc,argv);
    gain_opt.retrieveOption(argc,argv);
    modnodata_opt.retrieveOption(argc,argv);
    obsnodata_opt.retrieveOption(argc,argv);
    obsmin_opt.retrieveOption(argc,argv);
    obsmax_opt.retrieveOption(argc,argv);
    msknodata_opt.retrieveOption(argc,argv);
    mskband_opt.retrieveOption(argc,argv);
    eps_opt.retrieveOption(argc,argv);
    uncertModel_opt.retrieveOption(argc,argv);
    uncertObs_opt.retrieveOption(argc,argv);
    processNoise_opt.retrieveOption(argc,argv);
    uncertNodata_opt.retrieveOption(argc,argv);
    down_opt.retrieveOption(argc,argv);
    otype_opt.retrieveOption(argc,argv);
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

  try{
    ostringstream errorStream;
    if(model_opt.size()<1){
      errorStream << "Error: no model dataset selected, use option -mod" << endl;
      throw(errorStream.str());
    }
    if(observation_opt.size()<1){
      errorStream << "Error: no observation dataset selected, use option -obs" << endl;
      throw(errorStream.str());
    }
    if(find(direction_opt.begin(),direction_opt.end(),"forward")!=direction_opt.end()){
      if(outputfw_opt.empty()){
	errorStream << "Error: output forward datasets is not provided, use option -ofw" << endl;
	throw(errorStream.str());
      }
      if(uncertfw_opt.empty()){
	ostringstream uncertStream;
	uncertStream << outputfw_opt[0] << "_uncert";
	uncertfw_opt.push_back(uncertStream.str());
      }
    }
    if(find(direction_opt.begin(),direction_opt.end(),"backward")!=direction_opt.end()){
      if(outputbw_opt.empty()){
	errorStream << "Error: output backward datasets is not provided, use option -obw" << endl;
	throw(errorStream.str());
      }
      if(uncertbw_opt.empty()){
	ostringstream uncertStream;
	uncertStream << outputbw_opt[0] << "_uncert";
	uncertbw_opt.push_back(uncertStream.str());
      }
    }
    // if(model_opt.size()<observation_opt.size()){
    // 	errorStream << "Error: sequence of models should be larger than observations" << endl;
    // 	throw(errorStream.str());
    // }
    if(tmodel_opt.empty()){
      cout << "Warning: model time sequence is not provided, self generating time sequence from model input"  << endl;
    }
    if(tobservation_opt.empty()){
      cout << "Warning: observation time sequence is not provided, self generating time sequence from observation input" << endl;
    }
    if(find(direction_opt.begin(),direction_opt.end(),"smooth")!=direction_opt.end()){
      if(outputfw_opt.empty()){
	errorStream << "Error: output forward dataset is not provided, use option -ofw" << endl;
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
      if(uncertfb_opt.empty()){
	ostringstream uncertStream;
	uncertStream << outputfb_opt[0] << "_uncert";
	uncertfb_opt.push_back(uncertStream.str());
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
  ImgReaderGdal imgReaderModel1Mask;
  ImgReaderGdal imgReaderModel2Mask;
  ImgReaderGdal imgReaderObs;
  ImgReaderGdal imgReaderObsMask;
  //test
  ImgWriterGdal imgWriterGain;

  imgReaderModel1.open(model_opt[0]);
  imgReaderModel1.setNoData(modnodata_opt);
  imgReaderObs.open(observation_opt[0]);
  imgReaderObs.setNoData(obsnodata_opt);
  // if(observationmask_opt.empty())
  //   observationmask_opt=observation_opt;
  if(modelmask_opt.size()){
    imgReaderModel1Mask.open(modelmask_opt[0]);
    imgReaderModel1Mask.setNoData(msknodata_opt);
  }
  if(observationmask_opt.size()){
    imgReaderObsMask.open(observationmask_opt[0]);
    imgReaderObsMask.setNoData(msknodata_opt);
  }

  unsigned int nobs=(observation_opt.size()>1)? observation_opt.size() : imgReaderObs.nrOfBand();
  unsigned int nmodel=(model_opt.size()>1)? model_opt.size() : imgReaderModel1.nrOfBand();

  if(verbose_opt[0]){
    cout << "number of observations: " << nobs << endl;
    cout << "number of models: " << nmodel << endl;
  }

  int ncol=imgReaderObs.nrOfCol();
  int nrow=imgReaderObs.nrOfRow();
  if(projection_opt.empty())
    projection_opt.push_back(imgReaderObs.getProjection());
  double geotransform[6];
  imgReaderObs.getGeoTransform(geotransform);

  GDALDataType theType=GDT_Unknown;
  if(verbose_opt[0])
    cout << "possible output data types: ";
  for(int iType = 0; iType < GDT_TypeCount; ++iType){
    if(verbose_opt[0])
      cout << " " << GDALGetDataTypeName((GDALDataType)iType);
    if( GDALGetDataTypeName((GDALDataType)iType) != NULL
        && EQUAL(GDALGetDataTypeName((GDALDataType)iType),
                 otype_opt[0].c_str()))
      theType=(GDALDataType) iType;
  }
  if(theType==GDT_Unknown)
    theType=imgReaderObs.getDataType();

  string imageType;//=imgReaderObs.getImageType();
  if(oformat_opt.size())//default
    imageType=oformat_opt[0];
  if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
    string theInterleave="INTERLEAVE=";
    theInterleave+=imgReaderObs.getInterleave();
    option_opt.push_back(theInterleave);
  }

  if(down_opt.empty()){
    double resModel=imgReaderModel1.getDeltaX();
    double resObs=imgReaderObs.getDeltaX();
    int down=static_cast<int>(ceil(resModel/resObs));
    if(!(down%2))
      down+=1;
    down_opt.push_back(down);
  }

  int obsindex=0;

  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;

  double errObs=uncertNodata_opt[0];//start with high initial value in case we do not have first observation at time 0

  while(tmodel_opt.size()<nmodel)
    tmodel_opt.push_back(tmodel_opt.size());
  try{
    if(tobservation_opt.size()<nobs){
      if(nobs==nmodel){
	while(tobservation_opt.size()<nobs)
	  tobservation_opt.push_back(tobservation_opt.size());
      }
      else{
	ostringstream errorStream;
	errorStream << "Error: please provide time sequence for observation using option tobs" << endl;
	throw(errorStream.str());
      }
    }
  }
  catch(string errorString){
    std::cout << errorString << std::endl;
    exit(1);
  }
  
  vector<int> relobsindex;

  for(int tindex=0;tindex<tobservation_opt.size();++tindex){
    vector<int>::iterator modit;
    modit=upper_bound(tmodel_opt.begin(),tmodel_opt.end(),tobservation_opt[tindex]);
    int relpos=modit-tmodel_opt.begin()-1;
    assert(relpos>=0);//todo: for now, we assume model is available at time before first measurement
    relobsindex.push_back(relpos);
    if(verbose_opt[0]){
      cout << "observation " << tindex << ": " << "relative position in model time series is " << relpos << ", date of observation is (tobservation_opt[tindex]): " << tobservation_opt[tindex] << ", relobsindex.back(): " << relobsindex.back();
      if(observation_opt.size()>tindex)
	cout << ", filename observation: " << observation_opt[tindex];
      else
	cout << ", observation band index: " << tindex;
      if(model_opt.size()>relpos)
	cout << ", filename of corresponding model: " << model_opt[relpos] << endl;
      else
	cout << ", band index of corresponding model: " << relpos;
    }
  }

  int ndigit=log(1.0*tmodel_opt.back())/log(10.0)+1;

  double geox=0;
  double geoy=0;

  if(model_opt.size()==nmodel)
    imgReaderModel1.close();
  if(modelmask_opt.size()==nmodel)
    imgReaderModel1Mask.close();
  if(observation_opt.size()==nobs)
    imgReaderObs.close();
  if(observationmask_opt.size()==nobs)
    imgReaderObsMask.close();

  try{
    if(find(direction_opt.begin(),direction_opt.end(),"forward")!=direction_opt.end()){
      ///////////////////////////// forward model /////////////////////////
      cout << "Running forward model" << endl;
      obsindex=0;
      if(verbose_opt[0])
	cout << "Opening image " << outputfw_opt[0] << " for writing " << endl << flush;
    
      // imgWriterEst.open(theOutput,ncol,nrow,2,theType,imageType,option_opt);
      ImgWriterGdal imgWriterEst;
      ImgWriterGdal imgWriterUncert;
      imgWriterEst.open(outputfw_opt[0],ncol,nrow,nmodel,theType,imageType,option_opt);
      imgWriterEst.setProjectionProj4(projection_opt[0]);
      imgWriterEst.setGeoTransform(geotransform);
      imgWriterEst.GDALSetNoDataValue(obsnodata_opt[0]);
      imgWriterUncert.open(uncertfw_opt[0],ncol,nrow,nmodel,theType,imageType,option_opt);
      imgWriterUncert.setProjectionProj4(projection_opt[0]);
      imgWriterUncert.setGeoTransform(geotransform);

      try{
	//test
	if(gain_opt.size()){
	  imgWriterGain.open(gain_opt[0],ncol,nrow,nmodel,GDT_Float64,imageType,option_opt);
	  imgWriterGain.setProjectionProj4(projection_opt[0]);
	  imgWriterGain.setGeoTransform(geotransform);
	  imgWriterGain.GDALSetNoDataValue(obsnodata_opt[0]);
	}

	if(verbose_opt[0]){
	  cout << "processing time " << tmodel_opt[0] << endl;
	  if(obsindex<relobsindex.size()){
	    assert(tmodel_opt.size()>relobsindex[obsindex]);
	    cout << "next observation " << tmodel_opt[relobsindex[obsindex]] << endl;
	  }
	  else
	    cout << "There is no next observation" << endl;
	}
	if(model_opt.size()==nmodel){
	  imgReaderModel1.open(model_opt[0]);
	  imgReaderModel1.setNoData(modnodata_opt);
	}
	if(modelmask_opt.size()==nmodel){
	  imgReaderModel1Mask.open(modelmask_opt[0]);
	  imgReaderModel1Mask.setNoData(msknodata_opt);
	}
      }
      catch(string errorString){
	cerr << errorString << endl;
      }
      catch(...){
	cerr << "Error opening file " << model_opt[0] << endl;
      }
      
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
	  vector<double> lineModelMask;
	  vector<double> estWriteBuffer(ncol);
	  vector<double> uncertWriteBuffer(ncol);
	  //test
	  vector<double> gainWriteBuffer(ncol);
	  try{
	    for(int irow=jrow;irow<jrow+down_opt[0]&&irow<nrow;++irow){
	      imgWriterEst.image2geo(0,irow,geox,geoy);
	      imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
	      if(modRow<0||modRow>=imgReaderModel1.nrOfRow()){
		cerr << "Error: geo coordinates (" << geox << "," << geoy << ") not covered in model image " << imgReaderModel1.getFileName() << endl;
		assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	      }
	      // imgReaderModel1.readData(estReadBuffer,GDT_Float64,modRow,0,theResample);

	      int readModelBand=(model_opt.size()==nmodel)? 0:0;
	      int readModelMaskBand=(modelmask_opt.size()==nmodel)? mskband_opt[0]:0;
	      imgReaderModel1.readData(estReadBuffer,GDT_Float64,modRow,readModelBand,theResample);
	      if(modelmask_opt.size())
		imgReaderModel1Mask.readData(lineModelMask,GDT_Float64,modRow,readModelMaskBand,theResample);
	      for(int jcol=0;jcol<ncol;jcol+=down_opt[0]){
		for(int icol=jcol;icol<jcol+down_opt[0]&&icol<ncol;++icol){
		  imgWriterEst.image2geo(icol,irow,geox,geoy);
		  imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
		  if(modelmask_opt.size()){
		    if(imgReaderModel1Mask.isNoData(lineModelMask[modCol])){
		      estWriteBuffer[icol]=obsnodata_opt[0];
		      uncertWriteBuffer[icol]=uncertNodata_opt[0];
		      //test
		      gainWriteBuffer[icol]=obsnodata_opt[0];
		      continue;
		    }
		  }
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
		    //test
		    gainWriteBuffer[icol]=obsnodata_opt[0];
		    continue;
		  }
		  estWriteBuffer[icol]=modValue;
		  if(obsmin_opt.size()){
		    if(estWriteBuffer[icol]<obsmin_opt[0])
		      estWriteBuffer[icol]=obsmin_opt[0];
		  }
		  if(obsmax_opt.size()){
		    if(estWriteBuffer[icol]>obsmax_opt[0])
		      estWriteBuffer[icol]=obsmax_opt[0];
		  }
		  uncertWriteBuffer[icol]=uncertModel_opt[0];
		  //test
		  gainWriteBuffer[icol]=0;
		}
	      }
	      imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,0);
	      imgWriterUncert.writeData(uncertWriteBuffer,GDT_Float64,irow,0);
	      //test
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
	if(observation_opt.size()==nobs){
	  imgReaderObs.open(observation_opt[0]);
	  imgReaderObs.setNoData(obsnodata_opt);
	}
	if(observationmask_opt.size()==nobs){
	  imgReaderObsMask.open(observationmask_opt[0]);
	  imgReaderObsMask.setNoData(msknodata_opt);
	}
	imgReaderObs.getGeoTransform(geotransform);

	vector< vector<double> > obsLineVector(down_opt[0]);
	vector<double> obsLineBuffer;
	vector<double> obsMaskLineBuffer;
	vector<double> modelMaskLineBuffer;
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
	int readObsBand=(observation_opt.size()==nobs)? 0:0;
	int readObsMaskBand=(observationmask_opt.size()==nobs)? mskband_opt[0]:0;
	int readModelBand=(model_opt.size()==nmodel)? 0:0;
	int readModelMaskBand=(modelmask_opt.size()==nmodel)? mskband_opt[0]:0;
	for(int iline=-down_opt[0]/2;iline<down_opt[0]/2+1;++iline){
	  if(iline<0)//replicate line 0
	    imgReaderObs.readData(obsLineVector[iline+down_opt[0]/2],GDT_Float64,0,readObsBand);
	  else
	    imgReaderObs.readData(obsLineVector[iline+down_opt[0]/2],GDT_Float64,iline,readObsBand);
	}
	for(int jrow=0;jrow<nrow;jrow+=down_opt[0]){
	  for(int irow=jrow;irow<jrow+down_opt[0]&&irow<nrow;++irow){
	    imgWriterEst.image2geo(0,irow,geox,geoy);
	    imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
	    assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	    imgReaderModel1.readData(estReadBuffer,GDT_Float64,modRow,readModelBand,theResample);
	    if(modelmask_opt.size())
	      imgReaderModel1Mask.readData(modelMaskLineBuffer,GDT_Float64,modRow,readModelMaskBand);
	    int maxRow=(irow+down_opt[0]/2<imgReaderObs.nrOfRow()) ? irow+down_opt[0]/2 : imgReaderObs.nrOfRow()-1;
	    obsLineVector.erase(obsLineVector.begin());
	    imgReaderObs.readData(obsLineBuffer,GDT_Float64,maxRow,readObsBand);
	    obsLineVector.push_back(obsLineBuffer);

	    if(observationmask_opt.size())
	      imgReaderObsMask.readData(obsMaskLineBuffer,GDT_Float64,irow,readObsMaskBand);

	    for(int jcol=0;jcol<ncol;jcol+=down_opt[0]){
	      for(int icol=jcol;icol<jcol+down_opt[0]&&icol<ncol;++icol){
		imgWriterEst.image2geo(icol,irow,geox,geoy);
		imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
		assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
		bool modelIsNoData=false;
		if(modelmask_opt.size())
		  modelIsNoData=imgReaderModel1Mask.isNoData(modelMaskLineBuffer[modCol]);
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
		double errMod=uncertModel_opt[0];//*stdDev*stdDev;
		modelIsNoData=modelIsNoData||imgReaderModel1.isNoData(modValue);
		bool obsIsNoData=false;
		if(observationmask_opt.size())
		  obsIsNoData=imgReaderObsMask.isNoData(obsMaskLineBuffer[icol]);
		obsIsNoData=obsIsNoData||imgReaderObs.isNoData(obsLineBuffer[icol]);
		if(modelIsNoData){//model is nodata: retain observation 
		  if(obsIsNoData){//both model and observation nodata
		    estWriteBuffer[icol]=obsnodata_opt[0];
		    uncertWriteBuffer[icol]=uncertNodata_opt[0];
		    //test
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
		    uncertWriteBuffer[icol]=uncertObs_opt[0];
		  }
		}
		else{//model is valid: calculate estimate from model
		  estWriteBuffer[icol]=modValue;
		  uncertWriteBuffer[icol]=errMod;//in case observation is not valid
		  //test
		  gainWriteBuffer[icol]=0;
		}
		//measurement update
		if(!obsIsNoData){
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
		  if(!modelIsNoData){//model is valid
		    statfactory::StatFactory statobs;
		    statobs.setNoDataValues(obsnodata_opt);
		    double obsMeanValue=0;
		    double obsVarValue=0;
		    statobs.meanVar(obsWindowBuffer,obsMeanValue,obsVarValue);
		    double difference=0;
		    difference=obsMeanValue-modValue;
		    // errObs=uncertObs_opt[0]*sqrt(difference*difference);
		    errObs=uncertObs_opt[0]*difference*difference;//uncertainty of the observation (R in Kalman equations)
		    // double errorCovariance=errMod;
		    double errorCovariance=processNoise_opt[0]*obsVarValue;//assumed initial errorCovariance (P in Kalman equations)
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
		      if(uncertWriteBuffer[icol]>obsmax_opt[0])
			uncertWriteBuffer[icol]=obsmax_opt[0];
		    }
		  }
		  assert(kalmanGain<=1);
		  //test
		  gainWriteBuffer[icol]=kalmanGain;
		}
	      }
	    }
	    imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,0);
	    imgWriterUncert.writeData(uncertWriteBuffer,GDT_Float64,irow,0);
	    //test
	    if(gain_opt.size())
	      imgWriterGain.writeData(gainWriteBuffer,GDT_Float64,irow,0);
	  }
	}
	if(observation_opt.size()==nobs)
	  imgReaderObs.close();
	if(observationmask_opt.size()==nobs)
	  imgReaderObsMask.close();
	++obsindex;
      }
      if(model_opt.size()==nmodel)
	imgReaderModel1.close();
      if(modelmask_opt.size()==nmodel)
	imgReaderModel1Mask.close();
      imgWriterEst.close();
      imgWriterUncert.close();

      ImgUpdaterGdal imgUpdaterEst;
      ImgUpdaterGdal imgUpdaterUncert;
      for(int modindex=1;modindex<nmodel;++modindex){
	imgUpdaterEst.open(outputfw_opt[0]);
	imgUpdaterEst.setNoData(obsnodata_opt);
	imgUpdaterUncert.open(uncertfw_opt[0]);
	if(verbose_opt[0]){
	  cout << "processing time " << tmodel_opt[modindex] << endl;
	  if(obsindex<relobsindex.size())
	    cout << "next observation " << tmodel_opt[relobsindex[obsindex]] << endl;
	  else
	    cout << "There is no next observation" << endl;
	}

	//calculate regression between two subsequence model inputs
	if(model_opt.size()==nmodel){
	  imgReaderModel1.open(model_opt[modindex-1]);
	  imgReaderModel1.setNoData(modnodata_opt);
	  imgReaderModel2.open(model_opt[modindex]);
	  imgReaderModel2.setNoData(modnodata_opt);
	}
	if(modelmask_opt.size()==nmodel){
	  imgReaderModel1Mask.open(modelmask_opt[modindex-1]);
	  imgReaderModel1Mask.setNoData(msknodata_opt);
	  imgReaderModel2Mask.open(modelmask_opt[modindex]);
	  imgReaderModel2Mask.setNoData(msknodata_opt);
	}
	
	pfnProgress(progress,pszMessage,pProgressArg);

	bool update=false;
	if(obsindex<relobsindex.size()){
	  update=(relobsindex[obsindex]==modindex);
	}
	if(update){
	  if(observation_opt.size()==nobs){
	    if(verbose_opt[0])
	      cout << "***update " << relobsindex[obsindex] << " = " << modindex << " " << observation_opt[obsindex] << " ***" << endl;
	    imgReaderObs.open(observation_opt[obsindex]);
	    imgReaderObs.getGeoTransform(geotransform);
	    imgReaderObs.setNoData(obsnodata_opt);
	  }
	  if(observationmask_opt.size()==nobs){
	    imgReaderObsMask.open(observationmask_opt[obsindex]);
	    imgReaderObsMask.setNoData(msknodata_opt);
	  }
	}
	//prediction (also to fill cloudy pixels in measurement update mode)
	string input;
	input=outputfw_opt[0];

	vector< vector<double> > obsLineVector(down_opt[0]);
	vector<double> obsLineBuffer;
	vector<double> obsMaskLineBuffer;
	vector<double> model1MaskLineBuffer;
	vector<double> model2MaskLineBuffer;
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
	//test
	vector<double> gainWriteBuffer(ncol);

	int readObsBand=(observation_opt.size()==nobs)? 0:obsindex;
	int readObsMaskBand=(observationmask_opt.size()==nobs)? mskband_opt[0]:obsindex;
	int readModel1Band=(model_opt.size()==nmodel)? 0:modindex-1;
	int readModel2Band=(model_opt.size()==nmodel)? 0:modindex;
	int readModel1MaskBand=(modelmask_opt.size()==nmodel)? mskband_opt[0]:modindex-1;
	int readModel2MaskBand=(modelmask_opt.size()==nmodel)? mskband_opt[0]:modindex;

	//initialize obsLineVector if update
	if(update){
	  if(verbose_opt[0])
	    cout << "initialize obsLineVector" << endl;
	  assert(down_opt[0]%2);//window size must be odd 
	  for(int iline=-down_opt[0]/2;iline<down_opt[0]/2+1;++iline){
	    if(iline<0)//replicate line 0
	      imgReaderObs.readData(obsLineVector[iline+down_opt[0]/2],GDT_Float64,0,readObsBand);
	    else
	      imgReaderObs.readData(obsLineVector[iline+down_opt[0]/2],GDT_Float64,iline,readObsBand);
	  }
	}
	//initialize estLineVector
	if(verbose_opt[0])
	  cout << "initialize estLineVector" << endl;
	assert(down_opt[0]%2);//window size must be odd 

	for(int iline=-down_opt[0]/2;iline<down_opt[0]/2+1;++iline){
	  if(iline<0)//replicate line 0
	    imgUpdaterEst.readData(estLineVector[iline+down_opt[0]/2],GDT_Float64,0,modindex-1);
	  else
	    imgUpdaterEst.readData(estLineVector[iline+down_opt[0]/2],GDT_Float64,iline,modindex-1);
	}
	statfactory::StatFactory statobs;
	statobs.setNoDataValues(obsnodata_opt);
	for(int jrow=0;jrow<nrow;jrow+=down_opt[0]){
	  //todo: read entire window for uncertReadBuffer...
	  for(int irow=jrow;irow<jrow+down_opt[0]&&irow<nrow;++irow){
	    imgUpdaterUncert.readData(uncertReadBuffer,GDT_Float64,irow,modindex-1);
	    imgUpdaterUncert.image2geo(0,irow,geox,geoy);
	    if(model_opt.size()==nmodel){
	      imgReaderModel2.geo2image(geox,geoy,modCol,modRow);
	      assert(modRow>=0&&modRow<imgReaderModel2.nrOfRow());
	      imgReaderModel2.readData(model2LineBuffer,GDT_Float64,modRow,readModel2Band,theResample);
	      imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
	    }
	    else{
	      imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
	      imgReaderModel1.readData(model2LineBuffer,GDT_Float64,modRow,readModel2Band,theResample);
	    }
	    assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	    imgReaderModel1.readData(model1LineBuffer,GDT_Float64,modRow,readModel1Band,theResample);

	    if(modelmask_opt.size()){
	      imgReaderModel1Mask.readData(model1MaskLineBuffer,GDT_Float64,modRow,readModel1MaskBand);
	      if(modelmask_opt.size()==nmodel)
		imgReaderModel2Mask.readData(model2MaskLineBuffer,GDT_Float64,modRow,readModel2MaskBand);
	      else
		imgReaderModel1Mask.readData(model2MaskLineBuffer,GDT_Float64,modRow,readModel2MaskBand);
	    }

	    int maxRow=(irow+down_opt[0]/2<imgUpdaterEst.nrOfRow()) ? irow+down_opt[0]/2 : imgUpdaterEst.nrOfRow()-1;
	    estLineVector.erase(estLineVector.begin());
	    imgUpdaterEst.readData(estLineBuffer,GDT_Float64,maxRow,modindex-1);
	    estLineVector.push_back(estLineBuffer);
	    estLineBuffer=estLineVector[down_opt[0]/2];

	    if(update){
	      int maxRow=(irow+down_opt[0]/2<imgReaderObs.nrOfRow()) ? irow+down_opt[0]/2 : imgReaderObs.nrOfRow()-1;
	      obsLineVector.erase(obsLineVector.begin());
	      imgReaderObs.readData(obsLineBuffer,GDT_Float64,maxRow,readObsBand);
	      obsLineVector.push_back(obsLineBuffer);
	      obsLineBuffer=obsLineVector[down_opt[0]/2];

	      if(observationmask_opt.size())
		imgReaderObsMask.readData(obsMaskLineBuffer,GDT_Float64,irow,readObsMaskBand);
	    }

	    for(int jcol=0;jcol<ncol;jcol+=down_opt[0]){
	      for(int icol=jcol;icol<jcol+down_opt[0]&&icol<ncol;++icol){
		imgUpdaterEst.image2geo(icol,irow,geox,geoy);
		int minCol=(icol>down_opt[0]/2) ? icol-down_opt[0]/2 : 0;
		int maxCol=(icol+down_opt[0]/2<imgUpdaterEst.nrOfCol()) ? icol+down_opt[0]/2 : imgUpdaterEst.nrOfCol()-1;
		int minRow=(irow>down_opt[0]/2) ? irow-down_opt[0]/2 : 0;
		int maxRow=(irow+down_opt[0]/2<imgUpdaterEst.nrOfRow()) ? irow+down_opt[0]/2 : imgUpdaterEst.nrOfRow()-1;
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
		bool model1IsNoData=false;
		if(modelmask_opt.size())
		  model1IsNoData=imgReaderModel1Mask.isNoData(model1MaskLineBuffer[modCol]);
		lowerCol=modCol-0.5;
		lowerCol=static_cast<int>(lowerCol);
		upperCol=modCol+0.5;
		upperCol=static_cast<int>(upperCol);
		if(lowerCol<0)
		  lowerCol=0;
		if(upperCol>=imgReaderModel1.nrOfCol())
		  upperCol=imgReaderModel1.nrOfCol()-1;
		double modValue1=(modCol-0.5-lowerCol)*model1LineBuffer[upperCol]+(1-modCol+0.5+lowerCol)*model1LineBuffer[lowerCol];
		model1IsNoData=model1IsNoData||imgReaderModel1.isNoData(modValue1);
		if(model_opt.size()==nmodel)
		  imgReaderModel2.geo2image(geox,geoy,modCol,modRow);
		else
		  imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
		bool model2IsNoData=false;
		if(modelmask_opt.size())
		  model2IsNoData=imgReaderModel1Mask.isNoData(model2MaskLineBuffer[modCol]);
		lowerCol=modCol-0.5;
		lowerCol=static_cast<int>(lowerCol);
		upperCol=modCol+0.5;
		upperCol=static_cast<int>(upperCol);
		if(lowerCol<0)
		  lowerCol=0;
		if(upperCol>=imgReaderModel1.nrOfCol())
		  upperCol=imgReaderModel1.nrOfCol()-1;
		double modValue2=(modCol-0.5-lowerCol)*model2LineBuffer[upperCol]+(1-modCol+0.5+lowerCol)*model2LineBuffer[lowerCol];
		model2IsNoData=model2IsNoData||imgReaderModel1.isNoData(modValue2);
		bool obsIsNoData=false;
		if(observationmask_opt.size())
		  obsIsNoData=imgReaderObsMask.isNoData(obsMaskLineBuffer[icol]);
		obsIsNoData=obsIsNoData||imgReaderObs.isNoData(obsLineBuffer[icol]);

		if(imgUpdaterEst.isNoData(estValue)){
		  //we have not found any valid data yet, better here to take the current model value if valid
		  if(model2IsNoData){//if both estimate and model are no-data, set obs to nodata
		    estWriteBuffer[icol]=obsnodata_opt[0];
		    uncertWriteBuffer[icol]=uncertNodata_opt[0];
		    //test
		    gainWriteBuffer[icol]=0;
		  }
		  else{
		    estWriteBuffer[icol]=modValue2;
		    uncertWriteBuffer[icol]=uncertModel_opt[0];//*stdDev*stdDev;
		    if(obsmin_opt.size()){
		      if(estWriteBuffer[icol]<obsmin_opt[0])
			estWriteBuffer[icol]=obsmin_opt[0];
		    }
		    if(obsmax_opt.size()){
		      if(estWriteBuffer[icol]>obsmax_opt[0])
			estWriteBuffer[icol]=obsmax_opt[0];
		      if(uncertWriteBuffer[icol]>obsmax_opt[0])
			uncertWriteBuffer[icol]=obsmax_opt[0];
		    }
		    //test
		    gainWriteBuffer[icol]=0;
		  }
		}
		else{//previous estimate is valid
		  double estMeanValue=0;
		  double estVarValue=0;
		  statobs.meanVar(estWindowBuffer,estMeanValue,estVarValue);
		  double nvalid=0;
		  //time update
		  double processNoiseVariance=processNoise_opt[0]*estVarValue;
		  //estimate stability of weight distribution from model (low resolution) data in a window mod1 -> mod2 and assume distribution holds at fine spatial resolution. 

		  if(model1IsNoData||model2IsNoData){
		    estWriteBuffer[icol]=estValue;
		    // uncertWriteBuffer[icol]=uncertReadBuffer[icol]+processNoiseVariance;
		    //todo: check following line if makes sense
		    uncertWriteBuffer[icol]=uncertReadBuffer[icol]+uncertObs_opt[0];
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
		    if(uncertWriteBuffer[icol]>obsmax_opt[0])
		      uncertWriteBuffer[icol]=obsmax_opt[0];
		  }
		}
		//measurement update
		if(update&&!obsIsNoData){
		  double kalmanGain=1;
		  if(!model2IsNoData){//model is valid
		    statfactory::StatFactory statobs;
		    statobs.setNoDataValues(obsnodata_opt);
		    double obsMeanValue=0;
		    double obsVarValue=0;
		    double difference=0;
		    statobs.meanVar(obsWindowBuffer,obsMeanValue,obsVarValue);
		    difference=obsMeanValue-modValue2;
		    // errObs=uncertObs_opt[0]*sqrt(difference*difference);
		    errObs=uncertObs_opt[0]*difference*difference;//uncertainty of the observation (R in Kalman equations)

		    if(errObs<eps_opt[0])
		      errObs=eps_opt[0];
		    double errorCovariance=uncertWriteBuffer[icol];//P in Kalman equations

		    if(errorCovariance+errObs>eps_opt[0])
		      kalmanGain=errorCovariance/(errorCovariance+errObs);
		    else 
		      kalmanGain=1;
		    estWriteBuffer[icol]+=kalmanGain*(obsLineBuffer[icol]-estWriteBuffer[icol]);
		    uncertWriteBuffer[icol]*=(1-kalmanGain);
		    if(obsmin_opt.size()){
		      if(estWriteBuffer[icol]<obsmin_opt[0])
			estWriteBuffer[icol]=obsmin_opt[0];
		    }
		    if(obsmax_opt.size()){
		      if(estWriteBuffer[icol]>obsmax_opt[0])
			estWriteBuffer[icol]=obsmax_opt[0];
		      if(uncertWriteBuffer[icol]>obsmax_opt[0])
			uncertWriteBuffer[icol]=obsmax_opt[0];
		    }
		  }
		  assert(kalmanGain<=1);
		  //test
		  gainWriteBuffer[icol]=kalmanGain;
		}
	      }
	    }

	    //test
	    if(gain_opt.size())
	      imgWriterGain.writeData(gainWriteBuffer,GDT_Float64,irow,modindex);
	    imgUpdaterEst.writeData(estWriteBuffer,GDT_Float64,irow,modindex);
	    imgUpdaterUncert.writeData(uncertWriteBuffer,GDT_Float64,irow,modindex);
	    progress=static_cast<float>((irow+1.0)/imgUpdaterEst.nrOfRow());
	    pfnProgress(progress,pszMessage,pProgressArg);
	  }
	}

	//must close writers to ensure flush
	imgUpdaterEst.close();
	imgUpdaterUncert.close();
	// imgWriterEst.close();
	// imgReaderEst.close();

	if(update){
	  if(observation_opt.size()==nobs)
	    imgReaderObs.close();
	  if(observationmask_opt.size()==nobs)
	    imgReaderObsMask.close();
	  ++obsindex;
	}
	if(model_opt.size()==nmodel){
	  imgReaderModel1.close();
	  imgReaderModel2.close();
	}
	if(modelmask_opt.size()==nmodel){
	  imgReaderModel1Mask.close();
	  imgReaderModel2Mask.close();
	}
      }
      //test
      if(gain_opt.size())
	imgWriterGain.close();
    }
  }
  catch(string errorString){
    cerr << errorString << endl;
    exit(1);
  }
  catch(...){
    cerr << "Error in forward direction " << endl;
    exit(2);
  }
  try{
    if(find(direction_opt.begin(),direction_opt.end(),"backward")!=direction_opt.end()){
      ///////////////////////////// backward model /////////////////////////
      cout << "Running backward model" << endl;
      obsindex=relobsindex.size()-1;
      if(verbose_opt[0])
	cout << "Opening image " << outputbw_opt[0] << " for writing " << endl;

      // imgWriterEst.open(theOutput,ncol,nrow,2,theType,imageType,option_opt);
      ImgWriterGdal imgWriterEst;
      ImgWriterGdal imgWriterUncert;
      imgWriterEst.open(outputbw_opt[0],ncol,nrow,nmodel,theType,imageType,option_opt);
      imgWriterEst.setProjectionProj4(projection_opt[0]);
      imgWriterEst.setGeoTransform(geotransform);
      imgWriterEst.GDALSetNoDataValue(obsnodata_opt[0]);
      imgWriterUncert.open(uncertbw_opt[0],ncol,nrow,nmodel,theType,imageType,option_opt);
      imgWriterUncert.setProjectionProj4(projection_opt[0]);
      imgWriterUncert.setGeoTransform(geotransform);

      try{
	// //test
	// if(gain_opt.size()){
	//   imgWriterGain.open(gain_opt[0],ncol,nrow,nmodel,GDT_Float64,imageType,option_opt);
	//   imgWriterGain.setProjectionProj4(projection_opt[0]);
	//   imgWriterGain.setGeoTransform(geotransform);
	//   imgWriterGain.GDALSetNoDataValue(obsnodata_opt[0]);
	// }

	if(verbose_opt[0]){
	  cout << "processing time " << tmodel_opt.back() << endl;
	  if(obsindex<relobsindex.size()){
	    assert(tmodel_opt.size()>relobsindex[obsindex]);
	    cout << "next observation " << tmodel_opt[relobsindex[obsindex]] << endl;
	  }
	  else
	    cout << "There is no next observation" << endl;
	}
	if(model_opt.size()==nmodel){
	  imgReaderModel1.open(model_opt.back());
	  imgReaderModel1.setNoData(modnodata_opt);
	}
	if(modelmask_opt.size()==nmodel){
	  imgReaderModel1Mask.open(modelmask_opt[0]);
	  imgReaderModel1Mask.setNoData(msknodata_opt);
	}
      }
      catch(string errorString){
	cerr << errorString << endl;
      }
      catch(...){
	cerr << "Error opening file " << model_opt[0] << endl;
      }

      double modRow=0;
      double modCol=0;
      double lowerCol=0;
      double upperCol=0;
      RESAMPLE theResample=BILINEAR;

      if(relobsindex.back()<nmodel-1){//initialize output_opt.back() as last model
	//write last model as output
	if(verbose_opt[0])
	  cout << "write last model as output" << endl;
	for(int jrow=0;jrow<nrow;jrow+=down_opt[0]){
	  vector<double> estReadBuffer;
	  vector<double> lineModelMask;
	  vector<double> estWriteBuffer(ncol);
	  vector<double> uncertWriteBuffer(ncol);
	  // //test
	  // vector<double> gainWriteBuffer(ncol);
	  try{
	    for(int irow=jrow;irow<jrow+down_opt[0]&&irow<nrow;++irow){
	      imgWriterEst.image2geo(0,irow,geox,geoy);
	      imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
	      if(modRow<0||modRow>=imgReaderModel1.nrOfRow()){
		cerr << "Error: geo coordinates (" << geox << "," << geoy << ") not covered in model image " << imgReaderModel1.getFileName() << endl;
		assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	      }
	      // imgReaderModel1.readData(estReadBuffer,GDT_Float64,modRow,0,theResample);
	      int readModelBand=(model_opt.size()==nmodel)? 0:nmodel-1;
	      int readModelMaskBand=(modelmask_opt.size()==nmodel)? mskband_opt[0]:0;
	      imgReaderModel1.readData(estReadBuffer,GDT_Float64,modRow,readModelBand,theResample);
	      if(modelmask_opt.size())
		imgReaderModel1Mask.readData(lineModelMask,GDT_Float64,modRow,readModelMaskBand,theResample);
	      for(int jcol=0;jcol<ncol;jcol+=down_opt[0]){
		for(int icol=jcol;icol<jcol+down_opt[0]&&icol<ncol;++icol){
		  imgWriterEst.image2geo(icol,irow,geox,geoy);
		  imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
		  if(lineModelMask.size()>modCol){
		    if(imgReaderModel1Mask.isNoData(lineModelMask[modCol])){
		      estWriteBuffer[icol]=obsnodata_opt[0];
		      uncertWriteBuffer[icol]=uncertNodata_opt[0];
		      //test
		      // gainWriteBuffer[icol]=obsnodata_opt[0];
		      continue;
		    }
		  }
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
		    //test
		    // gainWriteBuffer[icol]=obsnodata_opt[0];
		    continue;
		  }
		  estWriteBuffer[icol]=modValue;
		  if(obsmin_opt.size()){
		    if(estWriteBuffer[icol]<obsmin_opt[0])
		      estWriteBuffer[icol]=obsmin_opt[0];
		  }
		  if(obsmax_opt.size()){
		    if(estWriteBuffer[icol]>obsmax_opt[0])
		      estWriteBuffer[icol]=obsmax_opt[0];
		  }
		  uncertWriteBuffer[icol]=uncertModel_opt[0];//*stdDev*stdDev;
		  //test
		  // gainWriteBuffer[icol]=0;
		}
	      }
	      imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,nmodel-1);
	      imgWriterUncert.writeData(uncertWriteBuffer,GDT_Float64,irow,nmodel-1);
	      // //test
	      // if(gain_opt.size())
	      //   imgWriterGain.writeData(gainWriteBuffer,GDT_Float64,irow,0);
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
	if(observation_opt.size()==nobs){
	  imgReaderObs.open(observation_opt.back());
	  imgReaderObs.setNoData(obsnodata_opt);
	}
	if(observationmask_opt.size()==nobs){
	  imgReaderObsMask.open(observationmask_opt.back());
	  imgReaderObsMask.setNoData(msknodata_opt);
	}
	imgReaderObs.getGeoTransform(geotransform);
      
	vector< vector<double> > obsLineVector(down_opt[0]);
	vector<double> obsLineBuffer;
	vector<double> obsMaskLineBuffer;
	vector<double> modelMaskLineBuffer;
	vector<double> obsWindowBuffer;//buffer for observation to calculate average corresponding to model pixel
	vector<double> estReadBuffer;
	vector<double> estWriteBuffer(ncol);
	vector<double> uncertWriteBuffer(ncol);
	vector<double> uncertObsLineBuffer;
	// //test
	// vector<double> gainWriteBuffer(ncol);

	if(verbose_opt[0])
	  cout << "initialize obsLineVector" << endl;
	assert(down_opt[0]%2);//window size must be odd 
	int readObsBand=(observation_opt.size()==nobs)? 0:nobs-1;
	int readObsMaskBand=(observationmask_opt.size()==nobs)? mskband_opt[0]:nobs-1;
	int readModelBand=(model_opt.size()==nmodel)? 0:nmodel-1;
	int readModelMaskBand=(modelmask_opt.size()==nmodel)? mskband_opt[0]:nmodel-1;
	for(int iline=-down_opt[0]/2;iline<down_opt[0]/2+1;++iline){
	  if(iline<0)//replicate line 0
	    imgReaderObs.readData(obsLineVector[iline+down_opt[0]/2],GDT_Float64,0,readObsBand);
	  else
	    imgReaderObs.readData(obsLineVector[iline+down_opt[0]/2],GDT_Float64,iline,readObsBand);
	}
	for(int jrow=0;jrow<nrow;jrow+=down_opt[0]){
	  for(int irow=jrow;irow<jrow+down_opt[0]&&irow<nrow;++irow){
	    imgWriterEst.image2geo(0,irow,geox,geoy);
	    imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
	    assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	    imgReaderModel1.readData(estReadBuffer,GDT_Float64,modRow,readModelBand,theResample);
	    if(modelmask_opt.size())
	      imgReaderModel1Mask.readData(modelMaskLineBuffer,GDT_Float64,modRow,readModelMaskBand);
	    int maxRow=(irow+down_opt[0]/2<imgReaderObs.nrOfRow()) ? irow+down_opt[0]/2 : imgReaderObs.nrOfRow()-1;
	    obsLineVector.erase(obsLineVector.begin());
	    imgReaderObs.readData(obsLineBuffer,GDT_Float64,maxRow,readObsBand);
	    obsLineVector.push_back(obsLineBuffer);

	    if(observationmask_opt.size())
	      imgReaderObsMask.readData(obsMaskLineBuffer,GDT_Float64,irow,readObsMaskBand);

	    for(int jcol=0;jcol<ncol;jcol+=down_opt[0]){
	      for(int icol=jcol;icol<jcol+down_opt[0]&&icol<ncol;++icol){
		imgWriterEst.image2geo(icol,irow,geox,geoy);
		imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
		assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
		bool modelIsNoData=false;
		if(modelmask_opt.size())
		  modelIsNoData=imgReaderModel1Mask.isNoData(modelMaskLineBuffer[modCol]);
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
		double errMod=uncertModel_opt[0];//*stdDev*stdDev;
		modelIsNoData=modelIsNoData||imgReaderModel1.isNoData(modValue);
		bool obsIsNoData=false;
		if(observationmask_opt.size())
		  obsIsNoData=imgReaderObsMask.isNoData(obsMaskLineBuffer[icol]);
		obsIsNoData=obsIsNoData||imgReaderObs.isNoData(obsLineBuffer[icol]);
		if(modelIsNoData){//model is nodata: retain observation 
		  if(obsIsNoData){//both model and observation nodata
		    estWriteBuffer[icol]=obsnodata_opt[0];
		    uncertWriteBuffer[icol]=uncertNodata_opt[0];
		    //test
		    // gainWriteBuffer[icol]=obsnodata_opt[0];
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
		    uncertWriteBuffer[icol]=uncertObs_opt[0];
		  }
		}
		else{//model is valid: calculate estimate from model
		  estWriteBuffer[icol]=modValue;
		  uncertWriteBuffer[icol]=errMod;//in case observation is not valid
		  //test
		  // gainWriteBuffer[icol]=0;
		}
		//measurement update
		if(!obsIsNoData){
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
		  if(!modelIsNoData){//model is valid
		    statfactory::StatFactory statobs;
		    statobs.setNoDataValues(obsnodata_opt);
		    double obsMeanValue=0;
		    double obsVarValue=0;
		    statobs.meanVar(obsWindowBuffer,obsMeanValue,obsVarValue);
		    double difference=0;
		    difference=obsMeanValue-modValue;
		    // errObs=uncertObs_opt[0]*sqrt(difference*difference);
		    errObs=uncertObs_opt[0]*difference*difference;//uncertainty of the observation (R in Kalman equations)
		    // double errorCovariance=errMod;
		    double errorCovariance=processNoise_opt[0]*obsVarValue;//assumed initial errorCovariance (P in Kalman equations)
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
		      if(uncertWriteBuffer[icol]>obsmax_opt[0])
			uncertWriteBuffer[icol]=obsmax_opt[0];
		    }
		  }
		  assert(kalmanGain<=1);
		  //test
		  // gainWriteBuffer[icol]=kalmanGain;
		}
	      }
	    }
	    imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,nmodel-1);
	    imgWriterUncert.writeData(uncertWriteBuffer,GDT_Float64,irow,nmodel-1);
	    // //test
	    // if(gain_opt.size())
	    //   imgWriterGain.writeData(gainWriteBuffer,GDT_Float64,irow,0);
	  }
	}
	if(observation_opt.size()==nobs)
	  imgReaderObs.close();
	if(observationmask_opt.size()==nobs)
	  imgReaderObsMask.close();
	--obsindex;
      }

      if(model_opt.size()==nmodel)
	imgReaderModel1.close();
      if(modelmask_opt.size()==nmodel)
	imgReaderModel1Mask.close();
      imgWriterEst.close();
      imgWriterUncert.close();

      ImgUpdaterGdal imgUpdaterEst;
      ImgUpdaterGdal imgUpdaterUncert;
      for(int modindex=nmodel-2;modindex>=0;--modindex){
	imgUpdaterEst.open(outputbw_opt[0]);
	imgUpdaterEst.setNoData(obsnodata_opt);
	imgUpdaterUncert.open(uncertbw_opt[0]);
	if(verbose_opt[0]){
	  cout << "processing time " << tmodel_opt[modindex] << endl;
	  if(obsindex<relobsindex.size())
	    cout << "next observation " << tmodel_opt[relobsindex[obsindex]] << endl;
	  else
	    cout << "There is no next observation" << endl;
	}

	//calculate regression between two subsequence model inputs
	if(model_opt.size()==nmodel){
	  imgReaderModel1.open(model_opt[modindex+1]);
	  imgReaderModel1.setNoData(modnodata_opt);
	  imgReaderModel2.open(model_opt[modindex]);
	  imgReaderModel2.setNoData(modnodata_opt);
	}
	if(modelmask_opt.size()==nmodel){
	  imgReaderModel1Mask.open(modelmask_opt[modindex-1]);
	  imgReaderModel1Mask.setNoData(msknodata_opt);
	  imgReaderModel2Mask.open(modelmask_opt[modindex]);
	  imgReaderModel2Mask.setNoData(msknodata_opt);
	}
    
	pfnProgress(progress,pszMessage,pProgressArg);

	bool update=false;
	if(obsindex<relobsindex.size()){
	  update=(relobsindex[obsindex]==modindex);
	}
	if(update){
	  if(observation_opt.size()==nobs){
	    if(verbose_opt[0])
	      cout << "***update " << relobsindex[obsindex] << " = " << modindex << " " << observation_opt[obsindex] << " ***" << endl;
	    imgReaderObs.open(observation_opt[obsindex]);
	    imgReaderObs.getGeoTransform(geotransform);
	    imgReaderObs.setNoData(obsnodata_opt);
	  }
	  if(observationmask_opt.size()==nobs){
	    imgReaderObsMask.open(observationmask_opt[obsindex]);
	    imgReaderObsMask.setNoData(msknodata_opt);
	  }
	}
	//prediction (also to fill cloudy pixels in update mode)
	string input;
	input=outputbw_opt[0];
      
	vector< vector<double> > obsLineVector(down_opt[0]);
	vector<double> obsLineBuffer;
	vector<double> obsMaskLineBuffer;
	vector<double> model1MaskLineBuffer;
	vector<double> model2MaskLineBuffer;
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
	//test
	// vector<double> gainWriteBuffer(ncol);

	int readObsBand=(observation_opt.size()==nobs)? 0:obsindex;
	int readObsMaskBand=(observationmask_opt.size()==nobs)? mskband_opt[0]:obsindex;
	int readModel1Band=(model_opt.size()==nmodel)? 0:modindex+1;
	int readModel2Band=(model_opt.size()==nmodel)? 0:modindex;
	int readModel1MaskBand=(modelmask_opt.size()==nmodel)? mskband_opt[0]:modindex+1;
	int readModel2MaskBand=(modelmask_opt.size()==nmodel)? mskband_opt[0]:modindex;

	//initialize obsLineVector
	if(update){
	  if(verbose_opt[0])
	    cout << "initialize obsLineVector" << endl;
	  assert(down_opt[0]%2);//window size must be odd 
	  for(int iline=-down_opt[0]/2;iline<down_opt[0]/2+1;++iline){
	    if(iline<0)//replicate line 0
	      imgReaderObs.readData(obsLineVector[iline+down_opt[0]/2],GDT_Float64,0,readObsBand);
	    else
	      imgReaderObs.readData(obsLineVector[iline+down_opt[0]/2],GDT_Float64,iline,readObsBand);
	  }
	}
	//initialize estLineVector
	if(verbose_opt[0])
	  cout << "initialize estLineVector" << endl;
	assert(down_opt[0]%2);//window size must be odd 

	for(int iline=-down_opt[0]/2;iline<down_opt[0]/2+1;++iline){
	  if(iline<0)//replicate line 0
	    imgUpdaterEst.readData(estLineVector[iline+down_opt[0]/2],GDT_Float64,0,modindex+1);
	  else
	    imgUpdaterEst.readData(estLineVector[iline+down_opt[0]/2],GDT_Float64,iline,modindex+1);
	}
	statfactory::StatFactory statobs;
	statobs.setNoDataValues(obsnodata_opt);

	for(int jrow=0;jrow<nrow;jrow+=down_opt[0]){
	  //todo: read entire window for uncertReadBuffer...
	  for(int irow=jrow;irow<jrow+down_opt[0]&&irow<nrow;++irow){
	    imgUpdaterUncert.readData(uncertReadBuffer,GDT_Float64,irow,modindex+1);
	    imgUpdaterUncert.image2geo(0,irow,geox,geoy);
	    if(model_opt.size()==nmodel){
	      imgReaderModel2.geo2image(geox,geoy,modCol,modRow);
	      assert(modRow>=0&&modRow<imgReaderModel2.nrOfRow());
	      imgReaderModel2.readData(model2LineBuffer,GDT_Float64,modRow,readModel2Band,theResample);
	      imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
	    }
	    else{
	      imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
	      imgReaderModel1.readData(model2LineBuffer,GDT_Float64,modRow,readModel2Band,theResample);
	    }

	    assert(modRow>=0&&modRow<imgReaderModel1.nrOfRow());
	    imgReaderModel1.readData(model1LineBuffer,GDT_Float64,modRow,readModel1Band,theResample);
	    if(modelmask_opt.size()){
	      imgReaderModel1Mask.readData(model1MaskLineBuffer,GDT_Float64,modRow,readModel1MaskBand);
	      if(modelmask_opt.size()==nmodel)
		imgReaderModel2Mask.readData(model2MaskLineBuffer,GDT_Float64,modRow,readModel2MaskBand);
	      else
		imgReaderModel1Mask.readData(model2MaskLineBuffer,GDT_Float64,modRow,readModel2MaskBand);
	    }
	    int maxRow=(irow+down_opt[0]/2<imgUpdaterEst.nrOfRow()) ? irow+down_opt[0]/2 : imgUpdaterEst.nrOfRow()-1;
	    estLineVector.erase(estLineVector.begin());
	    imgUpdaterEst.readData(estLineBuffer,GDT_Float64,maxRow,modindex+1);
	    estLineVector.push_back(estLineBuffer);
	    estLineBuffer=estLineVector[down_opt[0]/2];

	    if(update){
	      int maxRow=(irow+down_opt[0]/2<imgReaderObs.nrOfRow()) ? irow+down_opt[0]/2 : imgReaderObs.nrOfRow()-1;
	      obsLineVector.erase(obsLineVector.begin());
	      imgReaderObs.readData(obsLineBuffer,GDT_Float64,maxRow,readObsBand);
	      obsLineVector.push_back(obsLineBuffer);
	      obsLineBuffer=obsLineVector[down_opt[0]/2];

	      if(observationmask_opt.size())
		imgReaderObsMask.readData(obsMaskLineBuffer,GDT_Float64,irow,readObsBand);
	    }
	    for(int jcol=0;jcol<ncol;jcol+=down_opt[0]){
	      for(int icol=jcol;icol<jcol+down_opt[0]&&icol<ncol;++icol){
		imgUpdaterEst.image2geo(icol,irow,geox,geoy);
		int minCol=(icol>down_opt[0]/2) ? icol-down_opt[0]/2 : 0;
		int maxCol=(icol+down_opt[0]/2<imgUpdaterEst.nrOfCol()) ? icol+down_opt[0]/2 : imgUpdaterEst.nrOfCol()-1;
		int minRow=(irow>down_opt[0]/2) ? irow-down_opt[0]/2 : 0;
		int maxRow=(irow+down_opt[0]/2<imgUpdaterEst.nrOfRow()) ? irow+down_opt[0]/2 : imgUpdaterEst.nrOfRow()-1;
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
		bool model1IsNoData=false;

		if(modelmask_opt.size())
		  model1IsNoData=imgReaderModel1Mask.isNoData(model1MaskLineBuffer[modCol]);

		lowerCol=modCol-0.5;
		lowerCol=static_cast<int>(lowerCol);
		upperCol=modCol+0.5;
		upperCol=static_cast<int>(upperCol);
		if(lowerCol<0)
		  lowerCol=0;
		if(upperCol>=imgReaderModel1.nrOfCol())
		  upperCol=imgReaderModel1.nrOfCol()-1;
		double modValue1=(modCol-0.5-lowerCol)*model1LineBuffer[upperCol]+(1-modCol+0.5+lowerCol)*model1LineBuffer[lowerCol];
		model1IsNoData=model1IsNoData||imgReaderModel1.isNoData(modValue1);
		if(model_opt.size()==nmodel)
		  imgReaderModel2.geo2image(geox,geoy,modCol,modRow);
		else
		  imgReaderModel1.geo2image(geox,geoy,modCol,modRow);
		bool model2IsNoData=false;

		if(modelmask_opt.size())
		  model2IsNoData=imgReaderModel1Mask.isNoData(model2MaskLineBuffer[modCol]);
		lowerCol=modCol-0.5;
		lowerCol=static_cast<int>(lowerCol);
		upperCol=modCol+0.5;
		upperCol=static_cast<int>(upperCol);
		if(lowerCol<0)
		  lowerCol=0;
		if(upperCol>=imgReaderModel1.nrOfCol())
		  upperCol=imgReaderModel1.nrOfCol()-1;
		double modValue2=(modCol-0.5-lowerCol)*model2LineBuffer[upperCol]+(1-modCol+0.5+lowerCol)*model2LineBuffer[lowerCol];
		model2IsNoData=model2IsNoData||imgReaderModel1.isNoData(modValue2);
		bool obsIsNoData=false;
		if(observationmask_opt.size())
		  obsIsNoData=imgReaderObsMask.isNoData(obsMaskLineBuffer[icol]);
		obsIsNoData=obsIsNoData||imgReaderObs.isNoData(obsLineBuffer[icol]);

		if(imgUpdaterEst.isNoData(estValue)){
		  //we have not found any valid data yet, better here to take the current model value if valid
		  if(model2IsNoData){//if both estimate and model are no-data, set obs to nodata
		    estWriteBuffer[icol]=obsnodata_opt[0];
		    uncertWriteBuffer[icol]=uncertNodata_opt[0];
		    //test
		    // gainWriteBuffer[icol]=0;
		  }
		  else{
		    estWriteBuffer[icol]=modValue2;
		    uncertWriteBuffer[icol]=uncertModel_opt[0];//*stdDev*stdDev;
		    if(obsmin_opt.size()){
		      if(estWriteBuffer[icol]<obsmin_opt[0])
			estWriteBuffer[icol]=obsmin_opt[0];
		    }
		    if(obsmax_opt.size()){
		      if(estWriteBuffer[icol]>obsmax_opt[0])
			estWriteBuffer[icol]=obsmax_opt[0];
		      if(uncertWriteBuffer[icol]>obsmax_opt[0])
			uncertWriteBuffer[icol]=obsmax_opt[0];
		    }
		    //test
		    // gainWriteBuffer[icol]=0;
		  }
		}
		else{//previous estimate is valid
		  double estMeanValue=0;
		  double estVarValue=0;
		  statobs.meanVar(estWindowBuffer,estMeanValue,estVarValue);
		  double nvalid=0;
		  //time update
		  double processNoiseVariance=processNoise_opt[0]*estVarValue;
		  //estimate stability of weight distribution from model (low resolution) data in a window mod1 -> mod2 and assume distribution holds at fine spatial resolution. 

		  if(model1IsNoData||model2IsNoData){
		    estWriteBuffer[icol]=estValue;
		    // uncertWriteBuffer[icol]=uncertReadBuffer[icol]+processNoiseVariance;
		    //todo: check following line if makes sense
		    uncertWriteBuffer[icol]=uncertReadBuffer[icol]+uncertObs_opt[0];
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
		    if(uncertWriteBuffer[icol]>obsmax_opt[0])
		      uncertWriteBuffer[icol]=obsmax_opt[0];
		  }
		}
		//measurement update
		if(update&&!imgReaderObs.isNoData(obsLineBuffer[icol])){
		  double kalmanGain=1;
		  if(!imgReaderModel1.isNoData(modValue2)){//model is valid
		    statfactory::StatFactory statobs;
		    statobs.setNoDataValues(obsnodata_opt);
		    double obsMeanValue=0;
		    double obsVarValue=0;
		    double difference=0;
		    statobs.meanVar(obsWindowBuffer,obsMeanValue,obsVarValue);
		    difference=obsMeanValue-modValue2;
		    // errObs=uncertObs_opt[0]*sqrt(difference*difference);
		    errObs=uncertObs_opt[0]*difference*difference;//uncertainty of the observation (R in Kalman equations)

		    if(errObs<eps_opt[0])
		      errObs=eps_opt[0];
		    double errorCovariance=uncertWriteBuffer[icol];//P in Kalman equations

		    if(errorCovariance+errObs>eps_opt[0])
		      kalmanGain=errorCovariance/(errorCovariance+errObs);
		    else 
		      kalmanGain=1;
		    estWriteBuffer[icol]+=kalmanGain*(obsLineBuffer[icol]-estWriteBuffer[icol]);
		    uncertWriteBuffer[icol]*=(1-kalmanGain);
		    if(obsmin_opt.size()){
		      if(estWriteBuffer[icol]<obsmin_opt[0])
			estWriteBuffer[icol]=obsmin_opt[0];
		    }
		    if(obsmax_opt.size()){
		      if(estWriteBuffer[icol]>obsmax_opt[0])
			estWriteBuffer[icol]=obsmax_opt[0];
		      if(uncertWriteBuffer[icol]>obsmax_opt[0])
			uncertWriteBuffer[icol]=obsmax_opt[0];
		    }
		  }
		  assert(kalmanGain<=1);
		  //test
		  // gainWriteBuffer[icol]=kalmanGain;
		}
	      }
	    }
	    // //test
	    // if(gain_opt.size())
	    //   imgWriterGain.writeData(gainWriteBuffer,GDT_Float64,irow,modindex);
	    imgUpdaterEst.writeData(estWriteBuffer,GDT_Float64,irow,modindex);
	    imgUpdaterUncert.writeData(uncertWriteBuffer,GDT_Float64,irow,modindex);
	    progress=static_cast<float>((irow+1.0)/imgUpdaterEst.nrOfRow());
	    pfnProgress(progress,pszMessage,pProgressArg);
	  }
	}
	//must close writers to ensure flush
	imgUpdaterEst.close();
	imgUpdaterUncert.close();
	// imgWriterEst.close();
	// imgReaderEst.close();

	if(update){
	  if(observation_opt.size()==nobs)
	    imgReaderObs.close();
	  if(observationmask_opt.size()==nobs)
	    imgReaderObsMask.close();
	  --obsindex;
	}
	if(model_opt.size()==nmodel){
	  imgReaderModel1.close();
	  imgReaderModel2.close();
	}
	if(modelmask_opt.size()==nmodel){
	  imgReaderModel1Mask.close();
	  imgReaderModel2Mask.close();
	}
      }
      // //test
      // if(gain_opt.size())
      // 	imgWriterGain.close();
    }
  }
  catch(string errorString){
    cerr << errorString << endl;
    exit(1);
  }
  catch(...){
    cerr << "Error in backward direction " << endl;
    exit(2);
  }
  if(find(direction_opt.begin(),direction_opt.end(),"smooth")!=direction_opt.end()){
    ///////////////////////////// smooth model /////////////////////////
    cout << "Running smooth model" << endl;
    obsindex=0;

    ImgReaderGdal imgReaderForward(outputfw_opt[0]);
    ImgReaderGdal imgReaderBackward(outputbw_opt[0]);
    ImgReaderGdal imgReaderForwardUncert(uncertfw_opt[0]);
    ImgReaderGdal imgReaderBackwardUncert(uncertbw_opt[0]);
    imgReaderForward.setNoData(obsnodata_opt);
    imgReaderBackward.setNoData(obsnodata_opt);
      
    assert(imgReaderForward.nrOfBand()==nmodel);
    assert(imgReaderForwardUncert.nrOfBand()==nmodel);
    assert(imgReaderBackward.nrOfBand()==nmodel);
    assert(imgReaderBackwardUncert.nrOfBand()==nmodel);
    ImgWriterGdal imgWriterEst;
    imgWriterEst.setNoData(obsnodata_opt);
    ImgWriterGdal imgWriterUncert;
    imgWriterEst.open(outputfb_opt[0],ncol,nrow,nmodel,theType,imageType,option_opt);
    imgWriterEst.setProjectionProj4(projection_opt[0]);
    imgWriterEst.setGeoTransform(geotransform);
    imgWriterEst.GDALSetNoDataValue(obsnodata_opt[0]);

    imgWriterUncert.open(uncertfb_opt[0],ncol,nrow,nmodel,theType,imageType,option_opt);
    imgWriterUncert.setProjectionProj4(projection_opt[0]);
    imgWriterUncert.setGeoTransform(geotransform);
    for(int modindex=0;modindex<nmodel;++modindex){
      if(verbose_opt[0]){
	cout << "processing time " << tmodel_opt[modindex] << endl;
	if(obsindex<relobsindex.size())
	  cout << "next observation " << tmodel_opt[relobsindex[obsindex]] << endl;
	else
	  cout << "There is no next observation" << endl;
      }
      // if(outputfb_opt.size()==model_opt.size())
      // 	theOutput=outputfb_opt[modindex];
      // else{
      // 	ostringstream outputstream;
      // 	outputstream << outputfb_opt[0] << "_";
      // 	outputstream << setfill('0') << setw(ndigit) << tmodel_opt[modindex];
      // 	outputstream << ".tif";
      // 	// outputstream << outputfb_opt[0] << "_" << tmodel_opt[modindex] << ".tif";
      // 	theOutput=outputstream.str();
      // }
    
      //two band output band0=estimation, band1=uncertainty

      //open forward and backward estimates
      //we assume forward in model and backward in observation...

      // string inputfw;
      // string inputbw;
      // if(outputfw_opt.size()==model_opt.size())
      // 	inputfw=outputfw_opt[modindex];
      // else{
      // 	ostringstream outputstream;
      // 	outputstream << outputfw_opt[0] << "_";
      // 	outputstream << setfill('0') << setw(ndigit) << tmodel_opt[modindex];
      // 	outputstream << ".tif";
      // 	// outputstream << outputfw_opt[0] << "_" << tmodel_opt[modindex] << ".tif";
      // 	inputfw=outputstream.str();
      // }
      // if(outputbw_opt.size()==model_opt.size())
      // 	inputbw=outputbw_opt[modindex];
      // else{
      // 	ostringstream outputstream;
      // 	outputstream << outputbw_opt[0] << "_";
      // 	outputstream << setfill('0') << setw(ndigit) << tmodel_opt[modindex];
      // 	outputstream << ".tif";
      // 	// outputstream << outputbw_opt[0] << "_" << tmodel_opt[modindex] << ".tif";
      // 	inputbw=outputstream.str();
      // }
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

      int readObsBand=(observation_opt.size()==nobs)? 0:obsindex;

      if(update){
	if(observation_opt.size()==nobs){
	  if(verbose_opt[0])
	    cout << "***update " << relobsindex[obsindex] << " = " << modindex << " " << observation_opt[obsindex] << " ***" << endl;
	  imgReaderObs.open(observation_opt[obsindex]);
	  imgReaderObs.setNoData(obsnodata_opt);
	  imgReaderObs.getGeoTransform(geotransform);
	}
	if(observationmask_opt.size()==nobs){
	  imgReaderObsMask.open(observationmask_opt[obsindex]);
	  imgReaderObsMask.setNoData(msknodata_opt);
	}
	// imgReaderObs.open(observation_opt[obsindex]);
	// imgReaderObs.getGeoTransform(geotransform);
	// imgReaderObs.setNoData(obsnodata_opt);
	//calculate regression between model and observation
      }

      pfnProgress(progress,pszMessage,pProgressArg);

      for(int irow=0;irow<imgWriterEst.nrOfRow();++irow){
	assert(irow<imgReaderForward.nrOfRow());
	assert(irow<imgReaderBackward.nrOfRow());
	imgReaderForward.readData(estForwardBuffer,GDT_Float64,irow,modindex);
	imgReaderBackward.readData(estBackwardBuffer,GDT_Float64,irow,modindex);
	imgReaderForwardUncert.readData(uncertForwardBuffer,GDT_Float64,irow,modindex);
	imgReaderBackwardUncert.readData(uncertBackwardBuffer,GDT_Float64,irow,modindex);
	// imgReaderForward.readData(estForwardBuffer,GDT_Float64,irow,0);
	// imgReaderBackward.readData(estBackwardBuffer,GDT_Float64,irow,0);
	// imgReaderForward.readData(uncertForwardBuffer,GDT_Float64,irow,1);
	// imgReaderBackward.readData(uncertBackwardBuffer,GDT_Float64,irow,1);

	if(update){
	  if(observation_opt.size()==nobs)
	    imgReaderObs.readData(estWriteBuffer,GDT_Float64,irow,readObsBand);
	  if(observationmask_opt.size())
	    imgReaderObsMask.readData(uncertObsLineBuffer,GDT_Float64,irow,readObsBand);
	}

	// double oldRowMask=-1;//keep track of row mask to optimize number of line readings
	for(int icol=0;icol<imgWriterEst.nrOfCol();++icol){
	  imgWriterEst.image2geo(icol,irow,geox,geoy);
	  double A=estForwardBuffer[icol];
	  double B=estBackwardBuffer[icol];
	  double C=uncertForwardBuffer[icol];
	  double D=uncertBackwardBuffer[icol];
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
	      uncertWriteBuffer[icol]=eps_opt[0];
	    }
	    else{
	      estWriteBuffer[icol]=(A*D+B*C)/noemer;
	      uncertWriteBuffer[icol]=C*D/noemer;
	      if(obsmin_opt.size()){
		if(estWriteBuffer[icol]<obsmin_opt[0])
		estWriteBuffer[icol]=obsmin_opt[0];
	      }
	      if(obsmax_opt.size()){
		if(estWriteBuffer[icol]>obsmax_opt[0])
		  estWriteBuffer[icol]=obsmax_opt[0];
		if(uncertWriteBuffer[icol]>obsmax_opt[0])
		  uncertWriteBuffer[icol]=obsmax_opt[0];
	      }
	      // double P=0;
	      // if(C>eps_opt[0])
	      // 	P+=1.0/C;
	      // if(D>eps_opt[0])
	      // 	P+=1.0/D;
	      // if(P>eps_opt[0])
	      // 	P=1.0/P;
	      // else
	      // 	P=0;
	      // uncertWriteBuffer[icol]=P;
	    }
	  }
	}
	imgWriterEst.writeData(estWriteBuffer,GDT_Float64,irow,modindex);
	imgWriterUncert.writeData(uncertWriteBuffer,GDT_Float64,irow,modindex);
	progress=static_cast<float>((irow+1.0)/imgWriterEst.nrOfRow());
	pfnProgress(progress,pszMessage,pProgressArg);
      }
      if(update){
	if(observation_opt.size()==nobs)
	  imgReaderObs.close();
	++obsindex;
      }
    }
    imgReaderForward.close();
    imgReaderBackward.close();
    imgWriterEst.close();
    imgWriterUncert.close();
  }
  if(observation_opt.size()<nobs)
    imgReaderObs.close();
  if(model_opt.size()<nmodel)
    imgReaderModel1.close();
  // if(mask_opt.size())
  //   maskReader.close();
}

