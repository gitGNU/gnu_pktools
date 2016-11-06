/**********************************************************************
pkgetmask_lib.cc: program to create mask image based on values in input raster image
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
#include <assert.h>
#include <vector>
#include "imageclasses/ImgRasterGdal.h"
#include "base/Optionpk.h"
#include "apps/AppFactory.h"

using namespace std;
using namespace app;

/**
 * @param app application specific option arguments
 * @return output image
 **/
shared_ptr<ImgRasterGdal> ImgRasterGdal::getMask(const AppFactory& app){
  shared_ptr<ImgRasterGdal> imgWriter=ImgRasterGdal::createImg();
  getMask(*imgWriter, app);
  return(imgWriter);
}

/**
 * @param imgWriter output raster setmask dataset
 * @return CE_None if successful, CE_Failure if failed
 **/
CPLErr ImgRasterGdal::getMask(ImgRasterGdal& imgWriter, const AppFactory& app){
  Optionpk<string> input_opt("i", "input", "Input image file");
  Optionpk<short> band_opt("b", "band", "band(s) used for mask", 0);
  Optionpk<double> min_opt("min", "min", "Values smaller than min threshold(s) are masked as invalid. Use one threshold for each band");
  Optionpk<double> max_opt("max", "max", "Values greater than max threshold(s) are masked as invalid. Use one threshold for each band");
  Optionpk<string> operator_opt("p", "operator", "Operator: [AND,OR].", "OR");
  Optionpk<unsigned short> data_opt("data", "data", "value(s) for valid pixels: between min and max", 1);
  Optionpk<unsigned short> nodata_opt("nodata", "nodata", "value(s) for invalid pixels: not between min and max", 0);
  Optionpk<string> output_opt("o", "output", "Output mask file");
  Optionpk<string> otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image", "Byte");
  Optionpk<string> oformat_opt("of", "oformat", "Output image format (see also gdal_translate).","GTiff");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<string> colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)");
  Optionpk<short> verbose_opt("v", "verbose", "verbose", 0,2);

  band_opt.setHide(1);
  operator_opt.setHide(1);
  otype_opt.setHide(1);
  oformat_opt.setHide(1);
  option_opt.setHide(1);
  colorTable_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=min_opt.retrieveOption(app.getArgc(),app.getArgv());
    max_opt.retrieveOption(app.getArgc(),app.getArgv());
    data_opt.retrieveOption(app.getArgc(),app.getArgv());
    nodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    band_opt.retrieveOption(app.getArgc(),app.getArgv());
    operator_opt.retrieveOption(app.getArgc(),app.getArgv());
    otype_opt.retrieveOption(app.getArgc(),app.getArgv());
    colorTable_opt.retrieveOption(app.getArgc(),app.getArgv());
    verbose_opt.retrieveOption(app.getArgc(),app.getArgv());
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    return(CE_Failure);
  }
  if(!doProcess){
    cout << endl;
    std::ostringstream helpStream;
    helpStream << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    throw(helpStream.str());//help was invoked, stop processing
  }

  GDALDataType theType=GDT_Unknown;
  if(otype_opt.size()){
    theType=getGDALDataType(otype_opt[0]);
    if(theType==GDT_Unknown)
      std::cout << "Warning: unknown output pixel type: " << otype_opt[0] << ", using input type as default" << std::endl;
  }
  //if output type not set, get type from input image
  if(theType==GDT_Unknown){
    theType=getDataType();
    if(verbose_opt[0])
      cout << "Using data type from input image: " << GDALGetDataTypeName(theType) << endl;
  }
  if(verbose_opt[0])
    cout << "Output pixel type:  " << GDALGetDataTypeName(theType) << endl;
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);

  if(band_opt.empty()){
    std::string errorString="Error: band is empty, use option -b";
    throw(errorString);
  }
  if(band_opt.size()>=nrOfBand()){
    std::string errorString="Error: bands exceed number of band in input image";
    throw(errorString);
  }

  if(min_opt.size()&&max_opt.size()){
    if(min_opt.size()!=max_opt.size()){
      std::string errorString="Error: number of min and max options must correspond if both min and max options are provided";
      throw(errorString);
    }
  }
  if(min_opt.size()){
    while(band_opt.size()>min_opt.size())
      min_opt.push_back(min_opt[0]);
    while(min_opt.size()>data_opt.size())
      data_opt.push_back(data_opt[0]);
  }
  if(max_opt.size()){
    while(band_opt.size()>max_opt.size())
      max_opt.push_back(max_opt[0]);
    while(max_opt.size()>data_opt.size())
      data_opt.push_back(data_opt[0]);
  }

  vector< vector<float> > lineBuffer(band_opt.size());
  for(unsigned int iband=0;iband<band_opt.size();++iband)
    lineBuffer.resize(nrOfCol());
  //if output type not set, get type from input image
  if(theType==GDT_Unknown){
    theType=getDataType();
    if(verbose_opt[0])
      cout << "Using data type from input image: " << GDALGetDataTypeName(theType) << endl;
  }
  string imageType;//=getImageType();
  if(oformat_opt.size())//default
    imageType=oformat_opt[0];
  if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
    string theInterleave="INTERLEAVE=";
    theInterleave+=getInterleave();
    option_opt.push_back(theInterleave);
  }
  imgWriter.open(nrOfCol(),nrOfRow(),1,theType);
  if(colorTable_opt.size()){
    if(colorTable_opt[0]!="none")
      imgWriter.setColorTable(colorTable_opt[0]);
  }
  else if (getColorTable()!=NULL)//copy colorTable from input image
    imgWriter.setColorTable(getColorTable());

  imgWriter.setProjection(getProjection());
  double gt[6];
  getGeoTransform(gt);
  imgWriter.setGeoTransform(gt);//ulx,uly,getDeltaX(),getDeltaY(),0,0);

  if(nodata_opt.size())
    imgWriter.GDALSetNoDataValue(nodata_opt[0]);

  vector<char> writeBuffer(imgWriter.nrOfCol());
  for(unsigned int irow=0;irow<nrOfRow();++irow){
    for(unsigned int iband=0;iband<band_opt.size();++iband)
      readData(lineBuffer[iband],irow,band_opt[iband]);
    for(unsigned int icol=0;icol<nrOfCol();++icol){
      bool valid=(operator_opt[0]=="OR")?false:true;
      unsigned short validValue=data_opt[0];
      if(min_opt.size()&&max_opt.size()){
        assert(max_opt.size()==min_opt.size());
        for(unsigned int ivalid=0;ivalid<min_opt.size();++ivalid){
          bool validBand=false;
          // for(unsigned int iband=0;iband<band_opt.size();++iband){
          unsigned short theBand=(band_opt.size()==min_opt.size())? ivalid:0;
          if(lineBuffer[theBand][icol]>=min_opt[ivalid]&&lineBuffer[theBand][icol]<=max_opt[ivalid]){
            validValue=data_opt[ivalid];
            validBand=true;
          }
          valid=(operator_opt[0]=="OR")?valid||validBand : valid&&validBand;
        }
      }
      else if(min_opt.size()){
        for(int ivalid=0;ivalid<min_opt.size();++ivalid){
          bool validBand=false;
          // for(int iband=0;iband<band_opt.size();++iband){
          unsigned short theBand=(band_opt.size()==min_opt.size())? ivalid:0;
          if(lineBuffer[theBand][icol]>=min_opt[ivalid]){
            validValue=data_opt[ivalid];
            validBand=true;
          }
          valid=(operator_opt[0]=="OR")?valid||validBand : valid&&validBand;
        }
      }
      else if(max_opt.size()){
        for(int ivalid=0;ivalid<max_opt.size();++ivalid){
          bool validBand=false;
          // for(int iband=0;iband<band_opt.size();++iband){
          unsigned short theBand=(band_opt.size()==max_opt.size())? ivalid:0;
          if(lineBuffer[theBand][icol]<=max_opt[ivalid]){
            validValue=data_opt[ivalid];
            validBand=true;
          }
          valid=(operator_opt[0]=="OR")?valid||validBand : valid&&validBand;
        }
      }
      if(valid)
        writeBuffer[icol]=validValue;
      else
        writeBuffer[icol]=nodata_opt[0];
    }
    imgWriter.writeData(writeBuffer,irow);
    progress=(1.0+irow)/imgWriter.nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
  return(CE_None);
}
