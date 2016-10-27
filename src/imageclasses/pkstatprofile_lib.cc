/**********************************************************************
pkstatprofile_lib.cc: program to calculate statistics in temporal or spectral profile
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
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <sys/types.h>
#include <stdio.h>
#include "base/Optionpk.h"
#include "base/Vector2d.h"
#include "algorithms/Filter2d.h"
#include "algorithms/Filter.h"
#include "imageclasses/ImgRasterGdal.h"
#include "imageclasses/ImgCollection.h"
#include "algorithms/StatFactory.h"
#include "apps/AppFactory.h"

using namespace std;
using namespace app;

/**
 * @param app application specific option arguments
 * @return output image
 **/
shared_ptr<ImgRasterGdal> ImgCollection::statProfile(const AppFactory& app){
  shared_ptr<ImgRasterGdal> imgWriter=ImgRasterGdal::createImg();
  statProfile(*imgWriter, app);
  return(imgWriter);
}

/**
 * @param imgWriter output raster profile dataset
 * @param app application specific option arguments
 * @return CE_None if successful, CE_Failure if failed
 **/
CPLErr ImgCollection::statProfile(ImgRasterGdal& imgWriter, const AppFactory& app){
  Optionpk<std::string> function_opt("f", "function", "Statistics function (mean, median, var, stdev, min, max, sum, mode (provide classes), ismin, ismax, proportion (provide classes), percentile, nvalid");
  Optionpk<double> percentile_opt("perc","perc","Percentile value(s) used for rule percentile");
  // Optionpk<short> class_opt("class", "class", "class value(s) to use for mode, proportion");
  Optionpk<double> nodata_opt("nodata", "nodata", "nodata value(s)");
  Optionpk<std::string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image","GDT_Unknown");
  // Optionpk<short> down_opt("d", "down", "down sampling factor. Use value 1 for no downsampling). Use value n>1 for downsampling (aggregation)", 1);
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0,2);

  percentile_opt.setHide(1);
  // class_opt.setHide(1);
  otype_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=function_opt.retrieveOption(app.getArgc(),app.getArgv());
    percentile_opt.retrieveOption(app.getArgc(),app.getArgv());
    nodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    otype_opt.retrieveOption(app.getArgc(),app.getArgv());
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

  if(empty()){
    std::string errorString="Input collection is empty";
    throw(errorString);
  }

  if(function_opt.empty()){
    std::ostringstream errorStream;
    errorStream << "Error: no function selected, use option -f" << endl;
    throw(errorStream.str());
  }

  GDALDataType theType=getGDALDataType(otype_opt[0]);
  if(theType==GDT_Unknown)
    theType=this->front()->getDataType();

  if(verbose_opt[0])
    std::cout << std::endl << "Output pixel type:  " << GDALGetDataTypeName(theType) << endl;

  if(verbose_opt[0])
    cout << "Calculating statistic metrics: " << function_opt.size() << endl;
  //todo: expand for collections that have different image dimensions and geotransforms?
  imgWriter.open(this->front()->nrOfCol(),this->front()->nrOfRow(),function_opt.size(),theType);
  imgWriter.setProjection(this->front()->getProjection());
  double gt[6];
  this->front()->getGeoTransform(gt);
  imgWriter.setGeoTransform(gt);

  // if(colorTable_opt.size()){
  //   if(colorTable_opt[0]!="none"){
  //     if(verbose_opt[0])
  //       cout << "set colortable " << colorTable_opt[0] << endl;
  //     assert(imgWriter.getDataType()==GDT_Byte);
  //     imgWriter.setColorTable(colorTable_opt[0]);
  //   }
  // }
  // else if(this->front()->getColorTable()!=NULL)
  //   imgWriter.setColorTable(this->front()->getColorTable());

  if(nodata_opt.size()){
    for(int iband=0;iband<imgWriter.nrOfBand();++iband)
      imgWriter.setNoData(nodata_opt);
  }

  assert(imgWriter.nrOfBand()==function_opt.size());
  Vector2d<double> lineInput(this->size(),this->front()->nrOfCol());
  assert(imgWriter.nrOfCol()==this->front()->nrOfCol());
  Vector2d<double> lineOutput(function_opt.size(),imgWriter.nrOfCol());
  statfactory::StatFactory stat;
  stat.setNoDataValues(nodata_opt);
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(unsigned int y=0;y<this->front()->nrOfRow();++y){
    std::vector<std::shared_ptr<ImgRasterGdal> >::const_iterator imit=begin();
    for(int ifile=0;ifile<size();++ifile)
      at(ifile)->readData(lineInput[ifile],y);
    vector<double> pixelInput(this->size());
    for(unsigned int x=0;x<this->front()->nrOfCol();++x){
      pixelInput=lineInput.selectCol(x);
      int ithreshold=0;//threshold to use for percentiles
      for(int imethod=0;imethod<function_opt.size();++imethod){
        switch(filter::Filter::getFilterType(function_opt[imethod])){
        case(filter::first):
          lineOutput[imethod][x]=lineInput[0][x];
          break;
        case(filter::last):
          lineOutput[imethod][x]=lineInput.back()[x];
          break;
        case(filter::nvalid):
          lineOutput[imethod][x]=stat.nvalid(pixelInput);
          break;
        case(filter::median):
          lineOutput[imethod][x]=stat.median(pixelInput);
          break;
        case(filter::min):
          lineOutput[imethod][x]=stat.mymin(pixelInput);
          break;
        case(filter::max):
          lineOutput[imethod][x]=stat.mymax(pixelInput);
          break;
        case(filter::sum):
          lineOutput[imethod][x]=stat.sum(pixelInput);
          break;
        case(filter::var):
          lineOutput[imethod][x]=stat.var(pixelInput);
          break;
        case(filter::stdev):
          lineOutput[imethod][x]=sqrt(stat.var(pixelInput));
          break;
        case(filter::mean):
          lineOutput[imethod][x]=stat.mean(pixelInput);
          break;
        case(filter::percentile):{
          double threshold=(ithreshold<percentile_opt.size())? percentile_opt[ithreshold] : percentile_opt[0];
          lineOutput[imethod][x]=stat.percentile(pixelInput,pixelInput.begin(),pixelInput.end(),threshold);
          ++ithreshold;
          break;
        }
        default:
          std::string errorString="method not supported";
          throw(errorString);
          break;
        }
      }
    }
    for(int imethod=0;imethod<function_opt.size();++imethod){
      imgWriter.writeData(lineOutput[imethod],y,imethod);

    }
    progress=(1.0+y)/imgWriter.nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
  // filter1d.stats(input,output,function_opt);

  // this->close();
  // imgWriter.close();
  return(CE_None);
}

/**
 * @param app application specific option arguments
 * @return output image
 **/
shared_ptr<ImgRasterGdal> ImgRasterGdal::statProfile(const AppFactory& app){
  try{
    shared_ptr<ImgRasterGdal> imgWriter;
    imgWriter=this->clone();//create clone to first object, allowing for polymorphism in case of derived ImgRasterGdal objects
    statProfile(*imgWriter, app);
    return(imgWriter);
  }
  catch(string helpString){
    cerr << helpString << endl;
    return(0);
  }
}

/**
 * @param imgWriter output raster profile dataset
 * @param app application specific option arguments
 * @return CE_None if successful, CE_Failure if failed
 **/
CPLErr ImgRasterGdal::statProfile(ImgRasterGdal& imgWriter, const AppFactory& app){
  Optionpk<std::string> function_opt("f", "function", "Statistics function (mean, median, var, stdev, min, max, sum, mode (provide classes), ismin, ismax, proportion (provide classes), percentile, nvalid");
  Optionpk<double> percentile_opt("perc","perc","Percentile value(s) used for rule percentile");
  // Optionpk<short> class_opt("class", "class", "class value(s) to use for mode, proportion");
  Optionpk<double> nodata_opt("nodata", "nodata", "nodata value(s)");
  Optionpk<std::string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image","GDT_Unknown");
  // Optionpk<short> down_opt("d", "down", "down sampling factor. Use value 1 for no downsampling). Use value n>1 for downsampling (aggregation)", 1);
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0,2);

  percentile_opt.setHide(1);
  // class_opt.setHide(1);
  otype_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=function_opt.retrieveOption(app.getArgc(),app.getArgv());
    percentile_opt.retrieveOption(app.getArgc(),app.getArgv());
    nodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    otype_opt.retrieveOption(app.getArgc(),app.getArgv());
    verbose_opt.retrieveOption(app.getArgc(),app.getArgv());
    if(!doProcess){
      cout << endl;
      std::ostringstream helpStream;
      helpStream << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
      throw(helpStream.str());//help was invoked, stop processing
    }

    if(function_opt.empty()){
      std::ostringstream errorStream;
      errorStream << "Error: no function selected, use option -f" << endl;
      throw(errorStream.str());
    }
    while(function_opt.countSubstring("percentile")<percentile_opt.size())
      function_opt.push_back("percentile");

    GDALDataType theType=getGDALDataType(otype_opt[0]);
    if(theType==GDT_Unknown)
      theType=this->getDataType();

    if(verbose_opt[0])
      std::cout << std::endl << "Output pixel type:  " << GDALGetDataTypeName(theType) << endl;

    if(verbose_opt[0])
      cout << "Calculating statistic metrics: " << function_opt.size() << endl;
    imgWriter.open(this->nrOfCol(),this->nrOfRow(),function_opt.size(),theType);
    imgWriter.setProjection(this->getProjection());
    double gt[6];
    this->getGeoTransform(gt);
    imgWriter.setGeoTransform(gt);

    // if(colorTable_opt.size()){
    //   if(colorTable_opt[0]!="none"){
    //     if(verbose_opt[0])
    //       cout << "set colortable " << colorTable_opt[0] << endl;
    //     assert(imgWriter.getDataType()==GDT_Byte);
    //     imgWriter.setColorTable(colorTable_opt[0]);
    //   }
    // }
    // else if(this->getColorTable()!=NULL)
    //   imgWriter.setColorTable(this->getColorTable());

    if(nodata_opt.size()){
      for(int iband=0;iband<imgWriter.nrOfBand();++iband)
        imgWriter.GDALSetNoDataValue(nodata_opt[0],iband);
    }

    assert(imgWriter.nrOfBand()==function_opt.size());
    Vector2d<double> lineInput(this->nrOfBand(),this->nrOfCol());
    assert(imgWriter.nrOfCol()==this->nrOfCol());
    Vector2d<double> lineOutput(function_opt.size(),imgWriter.nrOfCol());
    statfactory::StatFactory stat;
    stat.setNoDataValues(nodata_opt);
    const char* pszMessage;
    void* pProgressArg=NULL;
    GDALProgressFunc pfnProgress=GDALTermProgress;
    double progress=0;
    pfnProgress(progress,pszMessage,pProgressArg);
    for(unsigned int y=0;y<this->nrOfRow();++y){
      for(unsigned int iband=0;iband<this->nrOfBand();++iband)
        this->readData(lineInput[iband],y,iband);
      vector<double> pixelInput(this->nrOfBand());
      for(unsigned int x=0;x<this->nrOfCol();++x){
        pixelInput=lineInput.selectCol(x);
        int ithreshold=0;//threshold to use for percentiles
        for(int imethod=0;imethod<function_opt.size();++imethod){
          switch(filter::Filter::getFilterType(function_opt[imethod])){
          case(filter::first):
            lineOutput[imethod][x]=lineInput[0][x];
            break;
          case(filter::last):
            lineOutput[imethod][x]=lineInput.back()[x];
            break;
          case(filter::nvalid):
            lineOutput[imethod][x]=stat.nvalid(pixelInput);
            break;
          case(filter::median):
            lineOutput[imethod][x]=stat.median(pixelInput);
            break;
          case(filter::min):
            lineOutput[imethod][x]=stat.mymin(pixelInput);
            break;
          case(filter::max):
            lineOutput[imethod][x]=stat.mymax(pixelInput);
            break;
          case(filter::sum):
            lineOutput[imethod][x]=stat.sum(pixelInput);
            break;
          case(filter::var):
            lineOutput[imethod][x]=stat.var(pixelInput);
            break;
          case(filter::stdev):
            lineOutput[imethod][x]=sqrt(stat.var(pixelInput));
            break;
          case(filter::mean):
            lineOutput[imethod][x]=stat.mean(pixelInput);
            break;
          case(filter::percentile):{
            double threshold=(ithreshold<percentile_opt.size())? percentile_opt[ithreshold] : percentile_opt[0];
            lineOutput[imethod][x]=stat.percentile(pixelInput,pixelInput.begin(),pixelInput.end(),threshold);
            ++ithreshold;
            break;
          }
          default:
            std::string errorString="method not supported";
            throw(errorString);
            break;
          }
        }
      }
      for(int imethod=0;imethod<function_opt.size();++imethod){
        imgWriter.writeData(lineOutput[imethod],y,imethod);

      }
      progress=(1.0+y)/imgWriter.nrOfRow();
      pfnProgress(progress,pszMessage,pProgressArg);
    }
    // filter1d.stats(input,output,function_opt);

    // this->close();
    // imgWriter.close();
    return(CE_None);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    return(CE_Failure);
  }
}
