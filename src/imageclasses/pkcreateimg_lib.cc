
/**********************************************************************
pkcreateimg_lib.cc: program to mosaic and composite geo-referenced images
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
#include <iostream>
#include <string>
#include <memory>
#include "ImgRaster.h"
#include "base/Optionpk.h"
#include "algorithms/StatFactory.h"
#include "apps/AppFactory.h"

using namespace std;
using namespace app;
using namespace statfactory;

/**
 * @param app application specific option arguments
 * @return output image
 **/
shared_ptr<ImgRaster> ImgRaster::createImg(const AppFactory& app){
  shared_ptr<ImgRaster> pRaster=createImg();
  createImg(*pRaster, app);
  return(pRaster);
}

/**
 * @param imgWriter output raster dataset
 * @param app application specific option arguments
 * @return CE_None if successful, CE_Failure if failed
 **/
CPLErr ImgRaster::createImg(ImgRaster& imgWriter, const AppFactory& app){

  Optionpk<unsigned int> nsample_opt("ns", "nsample", "Number of samples");
  Optionpk<unsigned int> nline_opt("nl", "nline", "Number of lines");
  Optionpk<unsigned int> nband_opt("b", "nband", "Number of bands",1);
  Optionpk<double> ulx_opt("ulx", "ulx", "Upper left x value bounding box", 0.0);
  Optionpk<double> uly_opt("uly", "uly", "Upper left y value bounding box", 0.0);
  Optionpk<double> lrx_opt("lrx", "lrx", "Lower right x value bounding box", 0.0);
  Optionpk<double> lry_opt("lry", "lry", "Lower right y value bounding box", 0.0);
  Optionpk<double> dx_opt("dx", "dx", "Resolution in x");
  Optionpk<double> dy_opt("dy", "dy", "Resolution in y");
  Optionpk<string> otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image","Byte");
  Optionpk<double> nodata_opt("nodata", "nodata", "Nodata value to put in image if out of bounds.");
  Optionpk<unsigned long int> seed_opt("seed", "seed", "seed value for random generator",0);
  Optionpk<double> mean_opt("mean", "mean", "Mean value for random generator",0);
  Optionpk<double> sigma_opt("sigma", "sigma", "Sigma value for random generator",0);
  Optionpk<string> description_opt("d", "description", "Set image description");
  Optionpk<string> projection_opt("a_srs", "a_srs", "Override the spatial reference for the output file (leave blank to copy from input file, use epsg:3035 to use European projection and force to European grid");

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=nsample_opt.retrieveOption(app.getArgc(),app.getArgv());
    nline_opt.retrieveOption(app.getArgc(),app.getArgv());
    nband_opt.retrieveOption(app.getArgc(),app.getArgv());
    ulx_opt.retrieveOption(app.getArgc(),app.getArgv());
    uly_opt.retrieveOption(app.getArgc(),app.getArgv());
    lrx_opt.retrieveOption(app.getArgc(),app.getArgv());
    lry_opt.retrieveOption(app.getArgc(),app.getArgv());
    dx_opt.retrieveOption(app.getArgc(),app.getArgv());
    dy_opt.retrieveOption(app.getArgc(),app.getArgv());
    otype_opt.retrieveOption(app.getArgc(),app.getArgv());
    nodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    seed_opt.retrieveOption(app.getArgc(),app.getArgv());
    mean_opt.retrieveOption(app.getArgc(),app.getArgv());
    sigma_opt.retrieveOption(app.getArgc(),app.getArgv());
    description_opt.retrieveOption(app.getArgc(),app.getArgv());
    projection_opt.retrieveOption(app.getArgc(),app.getArgv());
    if(!doProcess){
      cout << endl;
      std::ostringstream helpStream;
      helpStream << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
      throw(helpStream.str());//help was invoked, stop processing
    }
    else if(nsample_opt.empty()){
      std::ostringstream errorStream;
      errorStream << "Warning: no number of samples (use option -ns). Returning empty image" << std::endl;
      throw(errorStream.str());
    }
    else if(nline_opt.empty()){
      std::ostringstream errorStream;
      errorStream << "Warning: no number of lines (use option -nl). Returning empty image" << std::endl;
      throw(errorStream.str());
    }
    GDALDataType theType=GDT_Unknown;
    for(int iType = 0; iType < GDT_TypeCount; ++iType){
      if( GDALGetDataTypeName((GDALDataType)iType) != NULL
          && EQUAL(GDALGetDataTypeName((GDALDataType)iType),
                   otype_opt[0].c_str()))
        theType=(GDALDataType) iType;
    }
    imgWriter.open(nsample_opt[0],nline_opt[0],nband_opt[0],theType);
    imgWriter.setNoData(nodata_opt);
    if(description_opt.size())
      imgWriter.setImageDescription(description_opt[0]);
    double gt[6];
    if(ulx_opt[0]<lrx_opt[0])
      gt[0]=ulx_opt[0];
    else
      gt[0]=0;
    if(dx_opt.size())
      gt[1]=dx_opt[0];
    else if(lrx_opt[0]>0){
      gt[1]=lrx_opt[0]-ulx_opt[0];
      gt[1]/=imgWriter.nrOfCol();
    }
    else
      gt[1]=1;
    gt[2]=0;
    if(uly_opt[0]>lry_opt[0])
      gt[3]=uly_opt[0];
    else
      gt[3]=0;
    gt[4]=0;
    if(dy_opt.size())
      gt[5]=-dy_opt[0];
    else if(lry_opt[0]>0){
      gt[5]=lry_opt[0]-uly_opt[0];
      gt[5]/=imgWriter.nrOfRow();
    }
    else
      gt[5]=1;
    imgWriter.setGeoTransform(gt);
    if(projection_opt.size())
      imgWriter.setProjectionProj4(projection_opt[0]);
    StatFactory stat;
    gsl_rng* rndgen=stat.getRandomGenerator(seed_opt[0]);
    vector<double> lineBuffer(imgWriter.nrOfCol());
    double value=stat.getRandomValue(rndgen,"gaussian",mean_opt[0],sigma_opt[0]);
    for(unsigned int iband=0;iband<imgWriter.nrOfBand();++iband){
      for(unsigned int irow=0;irow<imgWriter.nrOfRow();++irow){
        for(unsigned int icol=0;icol<imgWriter.nrOfCol();++icol){
          if(sigma_opt[0]>0||(!irow&&!iband)){
            value=stat.getRandomValue(rndgen,"gaussian",mean_opt[0],sigma_opt[0]);
            lineBuffer[icol]=value;
          }
        }
        imgWriter.writeData(lineBuffer,irow,iband);
      }
    }
    return(CE_None);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    return(CE_Failure);
  }
}
