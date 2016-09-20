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
#include "imageclasses/ImgRaster.h"
#include "algorithms/StatFactory.h"
#include "apps/AppFactory.h"

using namespace std;
using namespace app;

/**
 * @param app application specific option arguments
 * @return output image
 **/
shared_ptr<ImgRaster> ImgRaster::statProfile(const AppFactory& app){
  try{
    shared_ptr<ImgRaster> imgWriter;
    imgWriter=this->clone();//create clone to first object, allowing for polymorphism in case of derived ImgRaster objects
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
CPLErr ImgRaster::statProfile(ImgRaster& imgWriter, const AppFactory& app){
  Optionpk<std::string> function_opt("f", "function", "Statistics function (mean, median, var, stdev, min, max, sum, mode (provide classes), ismin, ismax, proportion (provide classes), percentile, nvalid");
  Optionpk<double> percentile_opt("perc","perc","Percentile value(s) used for rule percentile");
  // Optionpk<short> class_opt("class", "class", "class value(s) to use for mode, proportion");
  Optionpk<double> nodata_opt("nodata", "nodata", "nodata value(s)");
  Optionpk<std::string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image","");
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

    // Filter filter1d;
    // filter1d.setNoDataValues(nodata_opt);
    // for(int iclass=0;iclass<class_opt.size();++iclass)
    //   filter1d.pushClass(class_opt[iclass]);

    // filter1d.setThresholds(percentile_opt);

    // ImgRaster input;
    // ImgRaster output;
    // try{
    //   input.open(input_opt[0],memory_opt[0]);
    // }
    // catch(string errorstring){
    //   cout << errorstring << endl;
    //   exit(1);
    // }

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
      theType=this->getDataType();

    if(verbose_opt[0])
      std::cout << std::endl << "Output pixel type:  " << GDALGetDataTypeName(theType) << endl;

    if(verbose_opt[0])
      cout << "Calculating statistic metrics: " << function_opt.size() << endl;
    // imgWriter.open(output_opt[0],this->nrOfCol(),this->nrOfRow(),function_opt.size(),theType,imageType,memory_opt[0],option_opt);
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
