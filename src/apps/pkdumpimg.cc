/**********************************************************************
pkdumpimg.cc: program to dump image content to ascii or std out
Copyright (C) 2008-2012 Pieter Kempeneers

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
#include <fstream>
#include <vector>
#include <iostream>
#include <assert.h>
#include "Optionpk.h"
#include "imageclasses/ImgWriterGdal.h"
#include "imageclasses/ImgWriterOgr.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

int main(int argc, char *argv[])
{
  Optionpk<std::string> input_opt("i","input","input image file","");
  Optionpk<string> output_opt("o", "output", "Output ascii file (Default is empty: use stdout", "");
  Optionpk<string> otype_opt("ot", "otype", "Data type for output ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image", "");
  Optionpk<string> oformat_opt("of", "oformat", "Output format (matrix form or list (x,y,z) form. Default is matrix form", "matrix");
  Optionpk<int> band_opt("b", "band", "band index to crop (-1: crop all bands)", 0);
  Optionpk<string> extent_opt("e", "extent", "get boundary from extent from polygons in vector file", "");
  Optionpk<double> ulx_opt("ulx", "ulx", "Upper left x value bounding box (in geocoordinates if georef is true)",0.0);
  Optionpk<double> uly_opt("uly", "uly", "Upper left y value bounding box (in geocoordinates if georef is true)",0.0);
  Optionpk<double> lrx_opt("lrx", "lrx", "Lower left x value bounding box (in geocoordinates if georef is true)",0.0);
  Optionpk<double> lry_opt("lry", "lry", "Lower left y value bounding box (in geocoordinates if georef is true)",0.0);
  Optionpk<double> dx_opt("dx", "dx", "Output resolution in x (in meter) (0.0: keep original resolution)",0.0);
  Optionpk<double> dy_opt("dy", "dy", "Output resolution in y (in meter) (0.0: keep original resolution)",0.0);
  Optionpk<string> resample_opt("r", "resampling-method", "Resampling method (near: nearest neighbour, bilinear: bi-linear interpolation).", "near");
  Optionpk<short> flag_opt("f", "flag", "Flag value to put in image if out of bounds.", 0);
  Optionpk<double> nodata_opt("nodata", "nodata", "set no data value(s) for calculations (flags in input image)");
  Optionpk<short> verbose_opt("v", "verbose", "verbose (Default: 0)", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    otype_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    extent_opt.retrieveOption(argc,argv);
    ulx_opt.retrieveOption(argc,argv);
    uly_opt.retrieveOption(argc,argv);
    lrx_opt.retrieveOption(argc,argv);
    lry_opt.retrieveOption(argc,argv);
    dx_opt.retrieveOption(argc,argv);
    dy_opt.retrieveOption(argc,argv);
    resample_opt.retrieveOption(argc,argv);
    flag_opt.retrieveOption(argc,argv);
    nodata_opt.retrieveOption(argc,argv);
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

  ofstream outputStream;
  if(!output_opt[0].empty())
    outputStream.open(output_opt[0].c_str());
  
  RESAMPLE theResample;
  if(resample_opt[0]=="near"){
    theResample=NEAR;
    if(verbose_opt[0])
      cout << "resampling: nearest neighbour" << endl;
  }
  else if(resample_opt[0]=="bilinear"){
    theResample=BILINEAR;
    if(verbose_opt[0])
      cout << "resampling: bilinear interpolation" << endl;
  }
  else{
    std::cout << "Error: resampling method " << resample_opt[0] << " not supported" << std::endl;
    exit(1);
  }

  ImgWriterGdal imgWriter;
  GDALDataType theType;

  ImgReaderGdal imgReader(input_opt[0]);
  for(int inodata=0;inodata<nodata_opt.size();++inodata)
    imgReader.pushNoDataValue(nodata_opt[inodata]);

  ImgWriterGdal virtualWriter;//only for coordinate conversion (no output file defined)
  
  int nband=imgReader.nrOfBand();
  //get number of lines
  int nrow=imgReader.nrOfRow();
  int ncol=imgReader.nrOfCol();
  int ncropcol=0;
  int ncroprow=0;
  double dx=dx_opt[0];
  double dy=dy_opt[0];
  if(!dx||!dy){
    dx=imgReader.getDeltaX();
    dy=imgReader.getDeltaY();
  }
  //bounding box of cropped image
  double cropulx=ulx_opt[0];
  double cropuly=uly_opt[0];
  double croplrx=lrx_opt[0];
  double croplry=lry_opt[0];
  //get bounding box from extentReader if defined
  ImgReaderOgr extentReader;
  if(extent_opt[0]!=""){
    for(int iextent=0;iextent<extent_opt.size();++iextent){
      extentReader.open(extent_opt[iextent]);
      if(!(extentReader.getExtent(ulx_opt[0],uly_opt[0],lrx_opt[0],lry_opt[0]))){
        cerr << "Error: could not get extent from " << extent_opt[0] << endl;
        exit(1);
      }
      if(!iextent){
        cropulx=ulx_opt[0];
        cropuly=uly_opt[0];
        croplry=lry_opt[0];
        croplrx=lrx_opt[0];
      }
      else{
        if(ulx_opt[0]<cropulx)
          cropulx=ulx_opt[0];
        if(uly_opt[0]>cropuly)
          cropuly=uly_opt[0];
        if(lry_opt[0]<croplry)
          croplry=lry_opt[0];
        if(lrx_opt[0]>croplrx)
          croplrx=lrx_opt[0];
      }
      extentReader.close();
    }
  }
     double uli,ulj,lri,lrj;//image coordinates
     double magicX=1,magicY=1;
  if(ulx_opt[0]>=lrx_opt[0]){//default bounding box: no cropping
    uli=0;
    lri=imgReader.nrOfCol()-1;
    ulj=0;
    lrj=imgReader.nrOfRow()-1;
    ncropcol=imgReader.nrOfCol();
    ncroprow=imgReader.nrOfRow();
    imgReader.getBoundingBox(cropulx,cropuly,croplrx,croplry);
    imgReader.geo2image(cropulx+(magicX-1.0)*imgReader.getDeltaX(),cropuly-(magicY-1.0)*imgReader.getDeltaY(),uli,ulj);
    imgReader.geo2image(croplrx+(magicX-2.0)*imgReader.getDeltaX(),croplry-(magicY-2.0)*imgReader.getDeltaY(),lri,lrj);
    ncropcol=ceil((croplrx-cropulx)/dx);
    ncroprow=ceil((cropuly-croplry)/dy);
    ncropcol=ceil((croplrx-cropulx)/dx);
    ncroprow=ceil((cropuly-croplry)/dy);
  }
  else{
    imgReader.geo2image(cropulx+(magicX-1.0)*imgReader.getDeltaX(),cropuly-(magicY-1.0)*imgReader.getDeltaY(),uli,ulj);
    imgReader.geo2image(croplrx+(magicX-2.0)*imgReader.getDeltaX(),croplry-(magicY-2.0)*imgReader.getDeltaY(),lri,lrj);
    
    ncropcol=ceil((croplrx-cropulx)/dx);
    ncroprow=ceil((cropuly-croplry)/dy);
    uli=floor(uli);
    ulj=floor(ulj);
    lri=floor(lri);
    lrj=floor(lrj);
  }
  double startCol=uli;
  double endCol=lri;
  double dcropcol=0;
  double dcroprow=0;
  double deltaX=imgReader.getDeltaX();
  double deltaY=imgReader.getDeltaY();
  dcropcol=(lri-uli+1)/(dx/deltaX);
  dcroprow=(lrj-ulj+1)/(dy/deltaY);
  if(verbose_opt[0]){
    cout << "cropulx: " << cropulx << endl;
    cout << "cropuly: " << cropuly << endl;
    cout << "croplrx: " << croplrx << endl;
    cout << "croplry: " << croplry << endl;
    cout << "ncropcol: " << ncropcol << endl;
    cout << "ncroprow: " << ncroprow << endl;
    cout << "cropulx+ncropcol*dx: " << cropulx+ncropcol*dx << endl;
    cout << "cropuly-ncroprow*dy: " << cropuly-ncroprow*dy << endl;
    cout << "selected upper left column in input image: " << uli << endl;
    cout << "selected upper left row of input image: " << ulj << endl;
    cout << "selected lower right column in input image: " << lri << endl;
    cout << "selected lower right row in input image: " << lrj << endl;
  }

  imgWriter.setGeoTransform(cropulx,cropuly,dx,dy,0,0);
  // imgWriter.setProjection(imgReader.getProjection());

  double readRow=0;
  double readCol=0;
  double lowerCol=0;
  double upperCol=0;
  int readncol=endCol-startCol+1;
  //test
  // vector<double> readBuffer(readncol);
  vector<double> readBuffer(readncol+1);
  assert(band_opt[0]>=0);
  assert(band_opt[0]<imgReader.nrOfBand());
  assert(imgWriter.isGeoRef());
  for(int irow=0;irow<ncroprow;++irow){
    if(verbose_opt[0])
      std::cout << irow << std::endl;
    double x=0;
    double y=0;
    //convert irow to geo
    imgWriter.image2geo(0,irow,x,y);
    //lookup corresponding row for irow in this file
    imgReader.geo2image(x,y,readCol,readRow);
    if(readRow<0||readRow>=imgReader.nrOfRow()){
      if(verbose_opt[0]>1)
        std::cout << "skipping row " << readRow << std::endl;
      continue;
    }
    try{
      //test
      // imgReader.readData(readBuffer,GDT_Float64,startCol,endCol,readRow,band_opt[0],theResample);
      if(endCol<imgReader.nrOfCol()-1)
        imgReader.readData(readBuffer,GDT_Float64,startCol,endCol+1,readRow,band_opt[0],theResample);
      else
        imgReader.readData(readBuffer,GDT_Float64,startCol,endCol,readRow,band_opt[0],theResample);
      for(int ib=0;ib<ncropcol;++ib){
        assert(imgWriter.image2geo(ib,irow,x,y));
        //lookup corresponding row for irow in this file
        imgReader.geo2image(x,y,readCol,readRow);
        if(readCol<0||readCol>=imgReader.nrOfCol()){
          if(oformat_opt[0]=="matrix"){
            if(output_opt[0].empty())
              std::cout << flag_opt[0] << " ";
            else
              outputStream << flag_opt[0] << " ";
          }
          else{
            if(output_opt[0].empty())
              std::cout << x << " " << y << " " << flag_opt[0] << endl;
            else
              outputStream << x << " " << y << " " << flag_opt[0] << endl;
          }
        }
        else{
          switch(theResample){
          case(BILINEAR):
            lowerCol=readCol-0.5;
            lowerCol=static_cast<int>(lowerCol);
            upperCol=readCol+0.5;
            upperCol=static_cast<int>(upperCol);
            if(lowerCol<0)
              lowerCol=0;
            if(upperCol>=imgReader.nrOfCol())
              upperCol=imgReader.nrOfCol()-1;
            if(oformat_opt[0]=="matrix"){
              if(output_opt[0].empty())
                std::cout << (readCol-0.5-lowerCol)*readBuffer[upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[lowerCol-startCol] << " ";
              else
                outputStream << (readCol-0.5-lowerCol)*readBuffer[upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[lowerCol-startCol] << " ";
            }
            else if(!imgReader.isNoData(readBuffer[upperCol-startCol])&&!imgReader.isNoData(readBuffer[lowerCol-startCol])){
                if(output_opt[0].empty())
                  std::cout << x << " " << y << " " << (readCol-0.5-lowerCol)*readBuffer[upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[lowerCol-startCol] << " " << endl;
                else
                  outputStream << x << " " << y << " " << (readCol-0.5-lowerCol)*readBuffer[upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[lowerCol-startCol] << " " << endl;
            }
            break;
          default:
            readCol=static_cast<int>(readCol);
            readCol-=startCol;//we only start reading from startCol
            assert(readCol>=0);
            if(readCol>=readBuffer.size())
              std::cout << "Error: " << readCol << " >= " << readBuffer.size() << std::endl;
            assert(readCol<readBuffer.size());
            if(oformat_opt[0]=="matrix"){
              if(output_opt[0].empty())
                std::cout << readBuffer[readCol] << " ";
              else
                outputStream << readBuffer[readCol] << " ";
            }
            else if(!imgReader.isNoData(readBuffer[readCol])){
              if(output_opt[0].empty())
                std::cout << x << " " << y << " " << readBuffer[readCol] << std::endl;
              else
                outputStream << x << " " << y << " " << readBuffer[readCol] << std::endl;
            }
            break;
          }
        }
      }
    }
    catch(string errorstring){
      cout << errorstring << endl;
      exit(1);
    }
    if(oformat_opt[0]=="matrix"){
      if(output_opt[0].empty())
        std::cout << std::endl;
      else
        outputStream << std::endl;
    }
  }
  if(extent_opt[0]!=""){
    extentReader.close();
  }
  imgReader.close();
  if(!output_opt[0].empty())
    outputStream.close();
}
