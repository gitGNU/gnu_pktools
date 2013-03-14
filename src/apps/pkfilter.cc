/**********************************************************************
pkfilter.cc: program to filter raster images: median, min/max, morphological, filtering
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
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"


/*------------------
  Main procedure
  ----------------*/
int main(int argc,char **argv) {
  Optionpk<std::string> input_opt("i","input","input image file");
  Optionpk<std::string> output_opt("o", "output", "Output image file");
  Optionpk<bool> disc_opt("c", "circular", "circular disc kernel for dilation and erosion", false);
  Optionpk<double> angle_opt("a", "angle", "angle used for directional filtering in dilation.");
  Optionpk<std::string> method_opt("f", "filter", "filter function (median,variance,min,max,sum,mean,minmax,dilation,erosion,closing,opening,spatially homogeneous (central pixel must be identical to all other pixels within window),SobelX edge detection in X,SobelY edge detection in Y,SobelXY,SobelYX,smooth,density,majority voting (only for classes),forest aggregation (mixed),smooth no data (mask) values,threshold local filtering,ismin,ismax,heterogeneous (central pixel must be different than all other pixels within window),order,stdev)", "median");
  Optionpk<int> dimX_opt("dx", "dx", "filter kernel size in x, must be odd", 3);
  Optionpk<int> dimY_opt("dy", "dy", "filter kernel size in y, must be odd", 3);
  Optionpk<int> dimZ_opt("dz", "dz", "filter kernel size in z (band or spectral dimension), must be odd (example: 3).. Set dz>0 if 1-D filter must be used in band domain");
  Optionpk<short> class_opt("class", "class", "class value(s) to use for density, erosion, dilation, openening and closing, thresholding");
  Optionpk<double> threshold_opt("t", "threshold", "threshold value(s) to use for threshold filter (one for each class)", 0);
  Optionpk<short> mask_opt("m", "mask", "mask value(s) ");
  Optionpk<std::string> tap_opt("tap", "tap", "text file containing taps used for spatial filtering (from ul to lr). Use dimX and dimY to specify tap dimensions in x and y. Leave empty for not using taps");
  Optionpk<double> tapz_opt("tapz", "tapz", "taps used for spectral filtering");
  Optionpk<double> fwhm_opt("fwhm", "fwhm", "list of full width half to apply spectral filtering (-fwhm band1 -fwhm band2 ...)");
  Optionpk<std::string> srf_opt("srf", "srf", "list of ASCII files containing spectral response functions (two columns: wavelength response)");
  Optionpk<double> wavelengthIn_opt("win", "wavelengthIn", "list of wavelengths in input spectrum (-w band1 -w band2 ...)");
  Optionpk<double> wavelengthOut_opt("wout", "wavelengthOut", "list of wavelengths in output spectrum (-w band1 -w band2 ...)");
  Optionpk<std::string> interpolationType_opt("interp", "interp", "type of interpolation for spectral filtering (see http://www.gnu.org/software/gsl/manual/html_node/Interpolation-Types.html)","akima");
  Optionpk<std::string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image","");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image");
  Optionpk<string>  colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)");
  Optionpk<std::string> option_opt("co", "co", "options: NAME=VALUE [-co COMPRESS=LZW] [-co INTERLEAVE=BAND]");
  Optionpk<short> down_opt("d", "down", "down sampling factor. Use value 1 for no downsampling)", 1);
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    disc_opt.retrieveOption(argc,argv);
    angle_opt.retrieveOption(argc,argv);
    method_opt.retrieveOption(argc,argv);
    dimX_opt.retrieveOption(argc,argv);
    dimY_opt.retrieveOption(argc,argv);
    dimZ_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    class_opt.retrieveOption(argc,argv);
    threshold_opt.retrieveOption(argc,argv);
    mask_opt.retrieveOption(argc,argv);
    tap_opt.retrieveOption(argc,argv);
    tapz_opt.retrieveOption(argc,argv);
    fwhm_opt.retrieveOption(argc,argv);
    srf_opt.retrieveOption(argc,argv);
    wavelengthIn_opt.retrieveOption(argc,argv);
    wavelengthOut_opt.retrieveOption(argc,argv);
    down_opt.retrieveOption(argc,argv);
    interpolationType_opt.retrieveOption(argc,argv);
    otype_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
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

  ImgReaderGdal input;
  ImgWriterGdal output;
  assert(input_opt.size());
  input.open(input_opt[0]);
  // output.open(output_opt[0],input);
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
    theType=input.getDataType();

  if(verbose_opt[0])
    std::cout << std::endl << "Output pixel type:  " << GDALGetDataTypeName(theType) << endl;

  string imageType=input.getImageType();
  if(oformat_opt.size())
    imageType=oformat_opt[0];

  if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
    string theInterleave="INTERLEAVE=";
    theInterleave+=input.getInterleave();
    option_opt.push_back(theInterleave);
  }
  try{
    assert(output_opt.size());
    if(fwhm_opt.size()||srf_opt.size()){
      //todo: support down and offset
      int nband=fwhm_opt.size()? fwhm_opt.size():srf_opt.size();
      output.open(output_opt[0],input.nrOfCol(),input.nrOfRow(),nband,theType,imageType,option_opt);
    }
    else
      output.open(output_opt[0],(input.nrOfCol()+down_opt[0]-1)/down_opt[0],(input.nrOfRow()+down_opt[0]-1)/down_opt[0],input.nrOfBand(),theType,imageType,option_opt);
  }
  catch(string errorstring){
    cout << errorstring << endl;
    exit(4);
  }
  if(input.isGeoRef()){
    output.setProjection(input.getProjection());
    double ulx,uly,deltaX,deltaY,rot1,rot2;
    input.getGeoTransform(ulx,uly,deltaX,deltaY,rot1,rot2);
    output.setGeoTransform(ulx,uly,deltaX*down_opt[0],deltaY*down_opt[0],rot1,rot2);
  }
  if(colorTable_opt.size()){
    if(verbose_opt[0])
      cout << "set colortable " << colorTable_opt[0] << endl;
    assert(output.getDataType()==GDT_Byte);
    output.setColorTable(colorTable_opt[0]);
  }
  else if(input.getColorTable()!=NULL)
    output.setColorTable(input.getColorTable());

  if(input.isGeoRef()){
    output.setProjection(input.getProjection());
    output.copyGeoTransform(input);
  }

  filter2d::Filter2d filter2d;
  filter::Filter filter1d;
  if(class_opt.size()&&verbose_opt[0]){
    std::cout<< "class values: ";
    for(int iclass=0;iclass<class_opt.size();++iclass){
      if(!dimZ_opt.size())
        filter2d.pushClass(class_opt[iclass]);
      else
        filter1d.pushClass(class_opt[iclass]);
      if(verbose_opt[0])
        std::cout<< class_opt[iclass] << " ";
    }
    if(verbose_opt[0])
      std::cout<< std::endl;
  }
  if(mask_opt.size()&&verbose_opt[0]){
    std::cout<< "mask values: ";
    for(int imask=0;imask<mask_opt.size();++imask){
      if(verbose_opt[0])
        std::cout<< mask_opt[imask] << " ";
      filter2d.pushMask(mask_opt[imask]);
    }
    if(verbose_opt[0])
      std::cout<< std::endl;
  }
  if(tap_opt.size()){
    ifstream tapfile(tap_opt[0].c_str());
    assert(tapfile);
    Vector2d<double> taps(dimY_opt[0],dimX_opt[0]);

    for(int j=0;j<dimY_opt[0];++j){
      for(int i=0;i<dimX_opt[0];++i){
        tapfile >> taps[j][i];
      }
    }
    if(verbose_opt[0]){
      std::cout << "taps: ";
      for(int j=0;j<dimY_opt[0];++j){
        for(int i=0;i<dimX_opt[0];++i){
          std::cout<< taps[j][i] << " ";
        }
        std::cout<< std::endl;
      }
    }
    filter2d.setTaps(taps);    
    filter2d.filter(input,output);
    tapfile.close();
  }
  else if(tapz_opt.size()){
    filter1d.setTaps(tapz_opt);    
    filter1d.doit(input,output,down_opt[0]);
  }
  else if(fwhm_opt.size()){
    if(verbose_opt[0])
      std::cout << "spectral filtering to " << fwhm_opt.size() << " bands with provided fwhm " << std::endl;
    assert(wavelengthOut_opt.size()==fwhm_opt.size());
    assert(wavelengthIn_opt.size());

    Vector2d<double> lineInput(input.nrOfBand(),input.nrOfCol());
    Vector2d<double> lineOutput(wavelengthOut_opt.size(),input.nrOfCol());
    const char* pszMessage;
    void* pProgressArg=NULL;
    GDALProgressFunc pfnProgress=GDALTermProgress;
    double progress=0;
    pfnProgress(progress,pszMessage,pProgressArg);
    for(int y=0;y<input.nrOfRow();++y){
      for(int iband=0;iband<input.nrOfBand();++iband)
        input.readData(lineInput[iband],GDT_Float64,y,iband);
      filter1d.applyFwhm<double>(wavelengthIn_opt,lineInput,wavelengthOut_opt,fwhm_opt, interpolationType_opt[0], lineOutput, verbose_opt[0]);
      for(int iband=0;iband<output.nrOfBand();++iband){
        try{
          output.writeData(lineOutput[iband],GDT_Float64,y,iband);
        }
        catch(string errorstring){
          cerr << errorstring << "in band " << iband << ", line " << y << endl;
        }
      }
      progress=(1.0+y)/output.nrOfRow();
      pfnProgress(progress,pszMessage,pProgressArg);
    }
    // filter1d.applyFwhm(wavelengthIn_opt,input,wavelengthOut_opt,fwhm_opt,interpolationType_opt[0],output,verbose_opt[0]);
  }
  else if(srf_opt.size()){
    if(verbose_opt[0])
      std::cout << "spectral filtering to " << srf_opt.size() << " bands with provided SRF " << std::endl;
    assert(wavelengthIn_opt.size());
    vector< Vector2d<double> > srf(srf_opt.size());//[0] srf_nr, [1]: wavelength, [2]: response
    ifstream srfFile;
    for(int isrf=0;isrf<srf_opt.size();++isrf){
      srf[isrf].resize(2);
      srfFile.open(srf_opt[isrf].c_str());
      double v;
      //add 0 to make sure srf is 0 at boundaries after interpolation step
      srf[isrf][0].push_back(0);
      srf[isrf][1].push_back(0);
      srf[isrf][0].push_back(1);
      srf[isrf][1].push_back(0);
      while(srfFile >> v){
        srf[isrf][0].push_back(v);
        srfFile >> v;
        srf[isrf][1].push_back(v);
      }
      srfFile.close();
      //add 0 to make sure srf[isrf] is 0 at boundaries after interpolation step
      srf[isrf][0].push_back(srf[isrf][0].back()+1);
      srf[isrf][1].push_back(0);    
      srf[isrf][0].push_back(srf[isrf][0].back()+1);
      srf[isrf][1].push_back(0);
      if(verbose_opt[0])
        cout << "srf file details: " << srf[isrf][0].size() << " wavelengths defined" << endl;    
    }
    assert(output.nrOfBand()==srf.size());
    double centreWavelength=0;
    Vector2d<double> lineInput(input.nrOfBand(),input.nrOfCol());
    const char* pszMessage;
    void* pProgressArg=NULL;
    GDALProgressFunc pfnProgress=GDALTermProgress;
    double progress=0;
    pfnProgress(progress,pszMessage,pProgressArg);
    for(int y=0;y<input.nrOfRow();++y){
      for(int iband=0;iband<input.nrOfBand();++iband)
        input.readData(lineInput[iband],GDT_Float64,y,iband);
      for(int isrf=0;isrf<srf.size();++isrf){
        vector<double> lineOutput(input.nrOfCol());
        double delta=1.0;
        bool normalize=true;
        centreWavelength=filter1d.applySrf<double>(wavelengthIn_opt,lineInput,srf[isrf], interpolationType_opt[0], lineOutput, delta, normalize, verbose_opt[0]);
        if(verbose_opt[0])
          std::cout << "centre wavelength srf " << isrf << ": " << centreWavelength << std::endl;
        try{
          output.writeData(lineOutput,GDT_Float64,y,isrf);
        }
        catch(string errorstring){
          cerr << errorstring << "in srf " << srf_opt[isrf] << ", line " << y << endl;
        }

      }
      progress=(1.0+y)/output.nrOfRow();
      pfnProgress(progress,pszMessage,pProgressArg);
    }

    // filter1d.applySrf(wavelengthIn_opt,input,srf,interpolationType_opt[0],output,verbose_opt[0]);
  }
  else{
    if(colorTable_opt.size())
      output.setColorTable(colorTable_opt[0]);
    switch(filter2d::Filter2d::getFilterType(method_opt[0])){
    case(filter2d::dilate):
      if(dimZ_opt.size()){
        if(verbose_opt[0])
          std::cout<< "1-D filtering: dilate" << std::endl;
        filter1d.morphology(input,output,"dilate",dimZ_opt[0]);
      }
      else{
        if(angle_opt.size())
          filter2d.morphology(input,output,"dilate",dimX_opt[0],dimY_opt[0],disc_opt[0],angle_opt[0]);
        else
          filter2d.morphology(input,output,"dilate",dimX_opt[0],dimY_opt[0],disc_opt[0]);
      }
      break;
    case(filter2d::erode):
      if(dimZ_opt[0]>0){
        if(verbose_opt[0])
          std::cout<< "1-D filtering: dilate" << std::endl;
        filter1d.morphology(input,output,"erode",dimZ_opt[0]);
      }
      else{
        if(angle_opt.size())
          filter2d.morphology(input,output,"erode",dimX_opt[0],dimY_opt[0],disc_opt[0],angle_opt[0]);
        else
          filter2d.morphology(input,output,"erode",dimX_opt[0],dimY_opt[0],disc_opt[0]);
      }
      break;
    case(filter2d::close):{//closing
      ostringstream tmps;
      tmps << "/tmp/dilation_" << getpid() << ".tif";
      ImgWriterGdal tmpout;
      tmpout.open(tmps.str(),input);
      try{
        if(dimZ_opt.size()){
          filter1d.morphology(input,tmpout,"dilate",dimZ_opt[0]);
        }
        else{
          if(angle_opt.size())
            filter2d.morphology(input,tmpout,"dilate",dimX_opt[0],dimY_opt[0],disc_opt[0],angle_opt[0]);
          else
            filter2d.morphology(input,tmpout,"dilate",dimX_opt[0],dimY_opt[0],disc_opt[0]);
        }
      }
      catch(std::string errorString){
	std::cout<< errorString;
	exit(1);
      }
      tmpout.close();
      ImgReaderGdal tmpin;
      tmpin.open(tmps.str());
      if(dimZ_opt.size()){
        filter1d.morphology(tmpin,output,"erode",dimZ_opt[0]);
      }
      else{
        if(angle_opt.size())
          filter2d.morphology(tmpin,output,"erode",dimX_opt[0],dimY_opt[0],disc_opt[0],angle_opt[0]);
        else
          filter2d.morphology(tmpin,output,"erode",dimX_opt[0],dimY_opt[0],disc_opt[0]);
      }
      tmpin.close();
      if(remove(tmps.str().c_str( )) !=0){
        cerr << "could not remove " << tmps.str() << std::endl;
      }
      break;
    }
    case(filter2d::open):{//opening
      ostringstream tmps;
      tmps << "/tmp/erosion_" << getpid() << ".tif";
      ImgWriterGdal tmpout;
      tmpout.open(tmps.str(),input);
      if(dimZ_opt.size()){
        filter1d.morphology(input,tmpout,"erode",dimZ_opt[0]);
      }
      else{
        if(angle_opt.size())
          filter2d.morphology(input,tmpout,"erode",dimX_opt[0],dimY_opt[0],disc_opt[0],angle_opt[0]);
        else
          filter2d.morphology(input,tmpout,"erode",dimX_opt[0],dimY_opt[0],disc_opt[0]);
      }
      tmpout.close();
      ImgReaderGdal tmpin;
      tmpin.open(tmps.str());
      if(dimZ_opt.size()){
        filter1d.morphology(tmpin,output,"dilate",dimZ_opt[0]);
      }
      else{
        if(angle_opt.size())
          filter2d.morphology(tmpin,output,"dilate",dimX_opt[0],dimY_opt[0],disc_opt[0],angle_opt[0]);
        else
          filter2d.morphology(tmpin,output,"dilate",dimX_opt[0],dimY_opt[0],disc_opt[0]);
      }
      tmpin.close();
      if(remove(tmps.str().c_str( )) !=0){
        cerr << "could not remove " << tmps.str() << std::endl;
      }
      break;
    }
    case(filter2d::homog):{//spatially homogeneous
      filter2d.doit(input,output,"homog",dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      // filter2d.homogeneousSpatial(input,output,dimX_opt[0],disc_opt[0]);
      break;
    }
    case(filter2d::heterog):{//spatially heterogeneous
      filter2d.doit(input,output,"heterog",dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    }
    case(filter2d::sobelx):{//Sobel edge detection in X
      Vector2d<double> theTaps(3,3);
      theTaps[0][0]=-1.0;
      theTaps[0][1]=0.0;
      theTaps[0][2]=1.0;
      theTaps[1][0]=-2.0;
      theTaps[1][1]=0.0;
      theTaps[1][2]=2.0;
      theTaps[2][0]=-1.0;
      theTaps[2][1]=0.0;
      theTaps[2][2]=1.0;
      filter2d.setTaps(theTaps);
      filter2d.filter(input,output,true);
      break;
    }
    case(filter2d::sobely):{//Sobel edge detection in Y
      Vector2d<double> theTaps(3,3);
      theTaps[0][0]=1.0;
      theTaps[0][1]=2.0;
      theTaps[0][2]=1.0;
      theTaps[1][0]=0.0;
      theTaps[1][1]=0.0;
      theTaps[1][2]=0.0;
      theTaps[2][0]=-1.0;
      theTaps[2][1]=-2.0;
      theTaps[2][2]=-1.0;
      filter2d.setTaps(theTaps);
      filter2d.filter(input,output,true);
      break;
    }
    case(filter2d::sobelxy):{//Sobel edge detection in XY
      Vector2d<double> theTaps(3,3);
      theTaps[0][0]=0.0;
      theTaps[0][1]=1.0;
      theTaps[0][2]=2.0;
      theTaps[1][0]=-1.0;
      theTaps[1][1]=0.0;
      theTaps[1][2]=1.0;
      theTaps[2][0]=-2.0;
      theTaps[2][1]=-1.0;
      theTaps[2][2]=0.0;
      filter2d.setTaps(theTaps);
      filter2d.filter(input,output,true);
      break;
    }
    case(filter2d::sobelyx):{//Sobel edge detection in XY
      Vector2d<double> theTaps(3,3);
      theTaps[0][0]=2.0;
      theTaps[0][1]=1.0;
      theTaps[0][2]=0.0;
      theTaps[1][0]=1.0;
      theTaps[1][1]=0.0;
      theTaps[1][2]=-1.0;
      theTaps[2][0]=0.0;
      theTaps[2][1]=-1.0;
      theTaps[2][2]=-2.0;
      filter2d.setTaps(theTaps);
      filter2d.filter(input,output,true);
      break;
    }
    case(filter2d::smooth):{//Smoothing filter
      filter2d.smooth(input,output,dimX_opt[0],dimY_opt[0]);
      break;
    }
    case(filter2d::smoothnodata):{//Smoothing filter
      filter2d.smoothNoData(input,output,dimX_opt[0],dimY_opt[0]);
      break;
    }
    case(filter2d::threshold):
      filter2d.setThresholds(threshold_opt);
      filter2d.setClasses(class_opt);//deliberate fall through
    default:
      filter2d.doit(input,output,method_opt[0],dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    }
  }
  input.close();
  output.close();
  return 0;
}
