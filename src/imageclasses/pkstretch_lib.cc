/**********************************************************************
pkstretch_lib.cc: program to stretch raster images: histogram stretching
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
#include <string>
#include <iostream>
#include <memory>
#include "base/Optionpk.h"
#include "imageclasses/ImgRaster.h"
#include "apps/AppFactory.h"

using namespace std;
using namespace app;

/**
 * @param app application specific option arguments
 * @return output stretched raster dataset
 **/
shared_ptr<ImgRaster> ImgRaster::stretch(const AppFactory& app){
  try{
    shared_ptr<ImgRaster> imgWriter;
    // imgWriter=this->clone();//create clone to first object, allowing for polymorphism in case of derived ImgRaster objects
    imgWriter=ImgRaster::createImg();
    stretch(*imgWriter, app);
    return(imgWriter);
  }
  catch(string helpString){
    cerr << helpString << endl;
    return(0);
  }
}

/**
 * @param app application specific option arguments
 * @return output stretched raster dataset
 **/
shared_ptr<ImgRaster> ImgRaster::equalize(const AppFactory& app){
  try{
    shared_ptr<ImgRaster> imgWriter;
    // imgWriter=this->clone();//create clone to first object, allowing for polymorphism in case of derived ImgRaster objects
    imgWriter=ImgRaster::createImg();
    equalize(*imgWriter, app);
    return(imgWriter);
  }
  catch(string helpString){
    cerr << helpString << endl;
    return(0);
  }
}

/**
 * @param imgWriter output stretched raster dataset
 * @param app application specific option arguments
 * @return CE_None if successful, CE_Failure if failed
 **/
CPLErr ImgRaster::equalize(ImgRaster& imgWriter, const AppFactory& app){
  Optionpk<double> nodata_opt("nodata", "nodata", "nodata value(s) (first value will be put in output image)");
  Optionpk<double> src_min_opt("src_min","src_min","clip source below this minimum value");
  Optionpk<double> src_max_opt("src_max","src_max","clip source above this maximum value");
  Optionpk<double> fromValue_opt("dst_min", "dst_min", "mininum value in output image", 0);
  Optionpk<double> toValue_opt("dst_max", "dst_max", "maximum value in output image", 255);
  Optionpk<double> minp_opt("cc_min", "cc_min", "cumulative count cut from");
  Optionpk<double> maxp_opt("cc_max", "cc_max", "cumulative count cut to");
  Optionpk<int> band_opt("b", "b", "band");
  Optionpk<std::string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image");
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0,2);

  src_min_opt.setHide(1);
  src_max_opt.setHide(1);
  otype_opt.setHide(1);
  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=minp_opt.retrieveOption(app.getArgc(),app.getArgv());
    nodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    src_min_opt.retrieveOption(app.getArgc(),app.getArgv());
    src_max_opt.retrieveOption(app.getArgc(),app.getArgv());
    maxp_opt.retrieveOption(app.getArgc(),app.getArgv());
    fromValue_opt.retrieveOption(app.getArgc(),app.getArgv());
    toValue_opt.retrieveOption(app.getArgc(),app.getArgv());
    band_opt.retrieveOption(app.getArgc(),app.getArgv());
    otype_opt.retrieveOption(app.getArgc(),app.getArgv());
    verbose_opt.retrieveOption(app.getArgc(),app.getArgv());
    if(!doProcess){
      cout << endl;
      std::ostringstream helpStream;
      helpStream << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
      throw(helpStream.str());//help was invoked, stop processing
    }
    if(src_min_opt.size()){
      while(src_min_opt.size()<band_opt.size())
        src_min_opt.push_back(src_min_opt[0]);
    }
    if(src_max_opt.size()){
      while(src_max_opt.size()<band_opt.size())
        src_max_opt.push_back(src_max_opt[0]);
    }
    if(minp_opt.empty()&&maxp_opt.empty())
      minp_opt.push_back(0);
    if(minp_opt.size()&&maxp_opt.empty())
      maxp_opt.push_back(100.0-minp_opt[0]);
    if(maxp_opt.size()&&minp_opt.empty())
      minp_opt.push_back(100.0-maxp_opt[0]);
    if(maxp_opt[0]<minp_opt[0]){
      minp_opt.push_back(maxp_opt[0]);
      maxp_opt[0]=minp_opt[0];
      minp_opt[0]=minp_opt.back();
      minp_opt.erase(minp_opt.begin()+1);
    }
    if(band_opt.empty()){
      while(band_opt.size()<nrOfBand())
        band_opt.push_back(band_opt.size());
    }
    GDALDataType theType=getDataType();
    if(otype_opt.size())
      theType=getGDALDataType(otype_opt[0]);
    imgWriter.open(this->nrOfCol(),this->nrOfRow(),band_opt.size(),theType);
    imgWriter.setProjection(this->getProjection());
    imgWriter.copyGeoTransform(*this);

    if(nodata_opt.size()){
      setNoData(nodata_opt);
      imgWriter.setNoData(nodata_opt);
      for(unsigned int iband=0;iband<imgWriter.nrOfBand();++iband)
        imgWriter.GDALSetNoDataValue(nodata_opt[0],iband);
    }
    for(int iband=0;iband<band_opt.size();++iband){
      int theBand=band_opt[iband];
      if(theBand<0||theBand>nrOfBand()){
        string ErrorString="Error: band exceeds number of bands";
        throw(ErrorString);
      }
      double minValue=(src_min_opt.size())? src_min_opt[0] : 0;
      double maxValue=(src_max_opt.size())? src_max_opt[0] : 0;
      getMinMax(minValue,maxValue,theBand);
      if(src_min_opt.size())
        minValue=src_min_opt[0];
      if(src_max_opt.size())
        maxValue=src_max_opt[0];
      if(minValue>=maxValue)
       getMinMax(minValue,maxValue,theBand);
      int nbin=0;//calculate automatically from min and max
      getMinMax(minValue,maxValue,theBand);
      //nbin is calculated automatically in case of Byte/short integer types, else set to 100
      if((GDALGetDataTypeSize(getDataType())>>3)>2)
        nbin=255;
      std::vector<double> histogramValues;
      std::vector<double> histogramOutput;
      double nvalid=getHistogram(histogramOutput, minValue, maxValue, nbin, theBand);
      double scale=static_cast<double>(nbin-1)/(maxValue-minValue);
      //create cumulative histogram
      std::vector<double>::iterator histit=histogramOutput.begin();
      double value=0;
      int ibin=0;
      double minpv=minValue;//value corresponding to minp_opt[0]
      double maxpv=maxValue;//value corresponding to maxp_opt[0]
      while(histit!=histogramOutput.end()){
          histogramValues.push_back(ibin/scale+minValue);
          value+=*histit/nvalid;
          *histit=value;
          if(100*(*histit)<=minp_opt[0])
            minpv=histogramValues.back();
          if(100*(*histit)<=maxp_opt[0])
            maxpv=histogramValues.back();
          ++histit;
          ++ibin;
        // }
        // else
        //   histogramOutput.erase(histit);
      }
      nbin=histogramOutput.size();
      if(verbose_opt[0]){
        std::cout << "minValue: " << minValue << std::endl;
        std::cout << "maxValue: " << maxValue << std::endl;
        std::cout << "nbin: " << nbin << std::endl;
        std::cout << "total count cumulative histogram: " << value << std::endl;
        std::cout << "nvalid: " << nvalid << std::endl;
        std::cout << "minpv: " << minpv << std::endl;
        std::cout << "maxpv: " << maxpv << std::endl;
      }
      std::vector<double> lineBuffer(nrOfCol());
      for(int irow=0;irow<nrOfRow();++irow){
        readData(lineBuffer,irow,theBand);
        for(int icol=0;icol<nrOfCol();++icol){
          if(isNoData(lineBuffer[icol]))
            continue;
          else if(lineBuffer[icol]>maxValue)
            lineBuffer[icol]=toValue_opt[0];
          else if(lineBuffer[icol]<minValue)
            lineBuffer[icol]=fromValue_opt[0];
          else{
            int theBin=static_cast<unsigned long int>(scale*(lineBuffer[icol]-minValue));
            double cdf_min=100*histogramOutput[0];
            if(cdf_min<minp_opt[0])
              cdf_min=minp_opt[0];
            double cdf_max=maxp_opt[0];
            double cdf=100*histogramOutput[theBin];
            double value=0;
            if(lineBuffer[icol]<minpv)
              value=fromValue_opt[0];
            else if(lineBuffer[icol]>maxpv)
              value=toValue_opt[0];
            else
              value=(cdf-cdf_min)/(cdf_max-cdf_min)*(toValue_opt[0]-fromValue_opt[0])+fromValue_opt[0];
            //round values if output data type is not float or double
            if(getGDALDataType(otype_opt[0])==GDT_Float32||getGDALDataType(otype_opt[0])==GDT_Float64)
              lineBuffer[icol]=value;
            else
              lineBuffer[icol]=static_cast<int>(value+0.5);
          }
        }
        imgWriter.writeData(lineBuffer,irow,theBand);
      }
    }
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    return(CE_Failure);
  }
}

/**
 * @param imgWriter output stretched raster dataset
 * @param app application specific option arguments
 * @return CE_None if successful, CE_Failure if failed
 **/
CPLErr ImgRaster::stretch(ImgRaster& imgWriter, const AppFactory& app){
  Optionpk<double> nodata_opt("nodata", "nodata", "nodata value(s) (first value will be put in output image)");
  Optionpk<double> src_min_opt("src_min","src_min","clip source below this minimum value");
  Optionpk<double> src_max_opt("src_max","src_max","clip source above this maximum value");
  Optionpk<double> fromValue_opt("dst_min", "dst_min", "mininum value in output image", 0);
  Optionpk<double> toValue_opt("dst_max", "dst_max", "maximum value in output image", 255);
  Optionpk<double> minp_opt("cc_min", "cc_min", "cumulative count cut from");
  Optionpk<double> maxp_opt("cc_max", "cc_max", "cumulative count cut to");
  Optionpk<int> band_opt("b", "b", "band");
  Optionpk<bool> equalize_opt("eq", "eq", "Histogram equalization",false);
  Optionpk<std::string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image");
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0,2);

  src_min_opt.setHide(1);
  src_max_opt.setHide(1);
  otype_opt.setHide(1);
  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=minp_opt.retrieveOption(app.getArgc(),app.getArgv());
    nodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    src_min_opt.retrieveOption(app.getArgc(),app.getArgv());
    src_max_opt.retrieveOption(app.getArgc(),app.getArgv());
    maxp_opt.retrieveOption(app.getArgc(),app.getArgv());
    fromValue_opt.retrieveOption(app.getArgc(),app.getArgv());
    toValue_opt.retrieveOption(app.getArgc(),app.getArgv());
    band_opt.retrieveOption(app.getArgc(),app.getArgv());
    equalize_opt.retrieveOption(app.getArgc(),app.getArgv());
    otype_opt.retrieveOption(app.getArgc(),app.getArgv());
    verbose_opt.retrieveOption(app.getArgc(),app.getArgv());
    if(!doProcess){
      cout << endl;
      std::ostringstream helpStream;
      helpStream << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
      throw(helpStream.str());//help was invoked, stop processing
    }
    if(src_min_opt.size()){
      while(src_min_opt.size()<band_opt.size())
        src_min_opt.push_back(src_min_opt[0]);
    }
    if(src_max_opt.size()){
      while(src_max_opt.size()<band_opt.size())
        src_max_opt.push_back(src_max_opt[0]);
    }
    if(minp_opt.empty()&&maxp_opt.empty())
      minp_opt.push_back(0);
    if(minp_opt.size()&&maxp_opt.empty())
      maxp_opt.push_back(100.0-minp_opt[0]);
    if(maxp_opt.size()&&minp_opt.empty())
      minp_opt.push_back(100.0-maxp_opt[0]);
    if(maxp_opt[0]<minp_opt[0]){
      minp_opt.push_back(maxp_opt[0]);
      maxp_opt[0]=minp_opt[0];
      minp_opt[0]=minp_opt.back();
      minp_opt.erase(minp_opt.begin()+1);
    }
    if(band_opt.empty()){
      while(band_opt.size()<nrOfBand())
        band_opt.push_back(band_opt.size());
    }
    GDALDataType theType=getDataType();
    if(otype_opt.size())
      theType=getGDALDataType(otype_opt[0]);
    imgWriter.open(this->nrOfCol(),this->nrOfRow(),band_opt.size(),theType);
    imgWriter.setProjection(this->getProjection());
    imgWriter.copyGeoTransform(*this);

    if(nodata_opt.size()){
      setNoData(nodata_opt);
      imgWriter.setNoData(nodata_opt);
      for(unsigned int iband=0;iband<imgWriter.nrOfBand();++iband)
        imgWriter.GDALSetNoDataValue(nodata_opt[0],iband);
    }
    for(int iband=0;iband<band_opt.size();++iband){
      int theBand=band_opt[iband];
      if(theBand<0||theBand>nrOfBand()){
        string ErrorString="Error: band exceeds number of bands";
        throw(ErrorString);
      }
      double minValue=(src_min_opt.size())? src_min_opt[0] : 0;
      double maxValue=(src_max_opt.size())? src_max_opt[0] : 0;
      getMinMax(minValue,maxValue,theBand);
      if(src_min_opt.size())
        minValue=src_min_opt[0];
      if(src_max_opt.size())
        maxValue=src_max_opt[0];
      if(minValue>=maxValue)
       getMinMax(minValue,maxValue,theBand);
      int nbin=0;//calculate automatically from min and max
      getMinMax(minValue,maxValue,theBand);
      //nbin is calculated automatically in case of Byte/short integer types, else set to 100
      if((GDALGetDataTypeSize(getDataType())>>3)>2)
        nbin=255;
      std::vector<double> histogramValues;
      std::vector<double> histogramOutput;
      double nvalid=getHistogram(histogramOutput, minValue, maxValue, nbin, theBand);
      double scale=static_cast<double>(nbin-1)/(maxValue-minValue);
      //create cumulative histogram
      std::vector<double>::iterator histit=histogramOutput.begin();
      double value=0;
      int ibin=0;
      double minpv=minValue;//value corresponding to minp_opt[0]
      double maxpv=maxValue;//value corresponding to maxp_opt[0]
      while(histit!=histogramOutput.end()){
          histogramValues.push_back(ibin/scale+minValue);
          value+=*histit/nvalid;
          *histit=value;
          if(100*(*histit)<=minp_opt[0])
            minpv=histogramValues.back();
          if(100*(*histit)<=maxp_opt[0])
            maxpv=histogramValues.back();
          ++histit;
          ++ibin;
        // }
        // else
        //   histogramOutput.erase(histit);
      }
      nbin=histogramOutput.size();
      if(verbose_opt[0]){
        std::cout << "minValue: " << minValue << std::endl;
        std::cout << "maxValue: " << maxValue << std::endl;
        std::cout << "nbin: " << nbin << std::endl;
        std::cout << "total count cumulative histogram: " << value << std::endl;
        std::cout << "nvalid: " << nvalid << std::endl;
        std::cout << "minpv: " << minpv << std::endl;
        std::cout << "maxpv: " << maxpv << std::endl;
      }
      std::vector<double> lineBuffer(nrOfCol());
      for(int irow=0;irow<nrOfRow();++irow){
        readData(lineBuffer,irow,theBand);
        for(int icol=0;icol<nrOfCol();++icol){
          if(isNoData(lineBuffer[icol]))
            continue;
          else if(lineBuffer[icol]>maxValue)
            lineBuffer[icol]=toValue_opt[0];
          else if(lineBuffer[icol]<minValue)
            lineBuffer[icol]=fromValue_opt[0];
          else{
            int theBin=static_cast<unsigned long int>(scale*(lineBuffer[icol]-minValue));
            double cdf_min=100*histogramOutput[0];
            if(cdf_min<minp_opt[0])
              cdf_min=minp_opt[0];
            double cdf_max=maxp_opt[0];
            double cdf=100*histogramOutput[theBin];
            double value=0;
            if(lineBuffer[icol]<minpv)
              value=fromValue_opt[0];
            else if(lineBuffer[icol]>maxpv)
              value=toValue_opt[0];
            else if(equalize_opt[0])
              value=(cdf-cdf_min)/(cdf_max-cdf_min)*(toValue_opt[0]-fromValue_opt[0])+fromValue_opt[0];
            else
              value=(lineBuffer[icol]-minpv)/(maxpv-minpv)*(toValue_opt[0]-fromValue_opt[0])+fromValue_opt[0];
            //round values if output data type is not float or double
            if(getGDALDataType(otype_opt[0])==GDT_Float32||getGDALDataType(otype_opt[0])==GDT_Float64)
              lineBuffer[icol]=value;
            else
              lineBuffer[icol]=static_cast<int>(value+0.5);
          }
        }
        imgWriter.writeData(lineBuffer,irow,theBand);
      }
    }
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    return(CE_Failure);
  }
}
