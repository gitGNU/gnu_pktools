/**********************************************************************
pkfilter.cc: program to filter raster images (e.g., median, min/max, morphological filtering)
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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*------------------
  Main procedure
  ----------------*/
int main(int argc,char **argv) {
  std::string versionString="version ";
  versionString+=VERSION;
  versionString+=", Copyright (C) 2008-2012 Pieter Kempeneers.\n\
   This program comes with ABSOLUTELY NO WARRANTY; for details type use option -h.\n\
   This is free software, and you are welcome to redistribute it\n\
   under certain conditions; use option --license for details.";
  Optionpk<bool> version_opt("\0","version",versionString,false);
  Optionpk<bool> license_opt("lic","license","show license information",false);
  Optionpk<bool> help_opt("h","help","shows this help info",false);
  Optionpk<bool> todo_opt("\0","todo","",false);
  Optionpk<std::string> input_opt("i","input","input image file","");
  Optionpk<std::string> output_opt("o", "output", "Output image file", "");
  Optionpk<bool> disc_opt("c", "circular", "circular disc kernel for dilation and erosion", false);
  Optionpk<double> angle_opt("a", "angle", "angle used for directional filtering in dilation.");
  Optionpk<int> function_opt("f", "filter", "filter function (0: median, 1: variance, 2: min, 3: max, 4: sum, 5: mean, 6: minmax, 7: dilation, 8: erosion, 9: closing, 10: opening, 11: spatially homogeneous (central pixel must be identical to all other pixels within window), 12: SobelX edge detection in X, 13: SobelY edge detection in Y, 14: SobelXY, -14: SobelYX, 15: smooth, 16: density, 17: majority voting (only for classes), 18: forest aggregation (mixed), 19: smooth no data (mask) values), 20: threshold local filtering, 21: ismin, 22: ismax, 23: heterogeneous (central pixel must be different than all other pixels within window), 24: order, 25: stdev", 0);
  Optionpk<int> dimX_opt("dx", "dx", "filter kernel size in x, must be odd", 3);
  Optionpk<int> dimY_opt("dy", "dy", "filter kernel size in y, must be odd", 3);
  Optionpk<int> dimZ_opt("dz", "dz", "filter kernel size in z (band or spectral dimension), must be odd (example: 3).. Set dz>0 if 1-D filter must be used in band domain");
  Optionpk<short> class_opt("class", "class", "class value(s) to use for density, erosion, dilation, openening and closing, thresholding");
  Optionpk<double> threshold_opt("t", "threshold", "threshold value(s) to use for threshold filter (one for each class)", 0);
  Optionpk<short> mask_opt("\0", "mask", "mask value(s) ");
  Optionpk<std::string> tap_opt("tap", "tap", "text file containing taps used for spatial filtering (from ul to lr). Use dimX and dimY to specify tap dimensions in x and y. Leave empty for not using taps");
  Optionpk<double> tapz_opt("tapz", "tapz", "taps used for spectral filtering");
  Optionpk<std::string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image", "");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image");
  Optionpk<string>  colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)", "");
  Optionpk<std::string> option_opt("co", "co", "options: NAME=VALUE [-co COMPRESS=LZW] [-co INTERLEAVE=BAND]");
  Optionpk<short> down_opt("d", "down", "down sampling factor. Use value 1 for no downsampling)", 1);
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0);

  version_opt.retrieveOption(argc,argv);
  license_opt.retrieveOption(argc,argv);
  help_opt.retrieveOption(argc,argv);
  todo_opt.retrieveOption(argc,argv);

  input_opt.retrieveOption(argc,argv);
  output_opt.retrieveOption(argc,argv);
  disc_opt.retrieveOption(argc,argv);
  angle_opt.retrieveOption(argc,argv);
  function_opt.retrieveOption(argc,argv);
  dimX_opt.retrieveOption(argc,argv);
  dimY_opt.retrieveOption(argc,argv);
  dimZ_opt.retrieveOption(argc,argv);
  option_opt.retrieveOption(argc,argv);
  class_opt.retrieveOption(argc,argv);
  threshold_opt.retrieveOption(argc,argv);
  mask_opt.retrieveOption(argc,argv);
  tap_opt.retrieveOption(argc,argv);
  tapz_opt.retrieveOption(argc,argv);
  down_opt.retrieveOption(argc,argv);
  otype_opt.retrieveOption(argc,argv);
  oformat_opt.retrieveOption(argc,argv);
  colorTable_opt.retrieveOption(argc,argv);
  verbose_opt.retrieveOption(argc,argv);

  if(version_opt[0]||todo_opt[0]){
    cout << version_opt.getHelp() << endl;
    cout << "todo: " << todo_opt.getHelp() << endl;
    exit(0);
  }
  if(license_opt[0]){
    cout << Optionpk<bool>::getGPLv3License() << endl;
    exit(0);
  }
  if(help_opt[0]){
    cout << "usage: pkfilter -i inputimage -o outputimage [OPTIONS]" << endl;
    exit(0);
  }

  ImgReaderGdal input;
  ImgWriterGdal output;
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
  if(input.getColorTable()!=NULL)
    output.setColorTable(input.getColorTable());

  // if(input.isGeoRef()){
  //   output.setProjection(input.getProjection());
  //   output.copyGeoTransform(input);
  // }

  Filter2d::Filter2d filter2d;
  Filter filter1d;
  if(verbose_opt[0])
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
  if(verbose_opt[0])
    std::cout<< "mask values: ";
  for(int imask=0;imask<mask_opt.size();++imask){
    if(verbose_opt[0])
      std::cout<< mask_opt[imask] << " ";
    filter2d.pushMask(mask_opt[imask]);
  }
  if(verbose_opt[0])
    std::cout<< std::endl;
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
  else{
    if(colorTable_opt[0]!="")
      output.setColorTable(colorTable_opt[0]);
    switch(function_opt[0]){
    case(Filter2d::MEDIAN):
      filter2d.doit(input,output,Filter2d::MEDIAN,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    case(Filter2d::VAR):
      filter2d.doit(input,output,Filter2d::VAR,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    case(Filter2d::STDEV):
      filter2d.doit(input,output,Filter2d::STDEV,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    case(Filter2d::MIN):
      filter2d.doit(input,output,Filter2d::MIN,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    case(Filter2d::ISMIN):
      filter2d.doit(input,output,Filter2d::ISMIN,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    case(Filter2d::MAX):
      filter2d.doit(input,output,Filter2d::MAX,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    case(Filter2d::ISMAX):
      filter2d.doit(input,output,Filter2d::ISMAX,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    case(Filter2d::MINMAX):
      filter2d.doit(input,output,Filter2d::MINMAX,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    case(Filter2d::SUM):
      filter2d.doit(input,output,Filter2d::SUM,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    case(Filter2d::MEAN):
      filter2d.doit(input,output,Filter2d::MEAN,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    case(Filter2d::MAJORITY):
      filter2d.doit(input,output,Filter2d::MAJORITY,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    case(Filter2d::THRESHOLD):
      filter2d.setThresholds(threshold_opt);
      filter2d.setClasses(class_opt);
      filter2d.doit(input,output,Filter2d::THRESHOLD,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    case(Filter2d::MIXED):
      filter2d.doit(input,output,Filter2d::MIXED,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    case(Filter2d::DILATE):
      if(dimZ_opt.size()){
        if(verbose_opt[0])
          std::cout<< "1-D filtering: dilate" << std::endl;
        filter1d.morphology(input,output,Filter::DILATE,dimZ_opt[0]);
      }
      else{
        if(angle_opt.size())
          filter2d.morphology(input,output,Filter2d::DILATE,dimX_opt[0],dimY_opt[0],disc_opt[0],angle_opt[0]);
        else
          filter2d.morphology(input,output,Filter2d::DILATE,dimX_opt[0],dimY_opt[0],disc_opt[0]);
      }
      break;
    case(Filter2d::ERODE):
      if(dimZ_opt[0]>0){
        if(verbose_opt[0])
          std::cout<< "1-D filtering: dilate" << std::endl;
        filter1d.morphology(input,output,Filter::ERODE,dimZ_opt[0]);
      }
      else{
        if(angle_opt.size())
          filter2d.morphology(input,output,Filter2d::ERODE,dimX_opt[0],dimY_opt[0],disc_opt[0],angle_opt[0]);
        else
          filter2d.morphology(input,output,Filter2d::ERODE,dimX_opt[0],dimY_opt[0],disc_opt[0]);
      }
      break;
    case(Filter2d::CLOSE):{//closing
      ostringstream tmps;
      tmps << "/tmp/dilation_" << getpid() << ".tif";
      ImgWriterGdal tmpout;
      tmpout.open(tmps.str(),input);
      try{
        if(dimZ_opt.size()){
          filter1d.morphology(input,tmpout,Filter::DILATE,dimZ_opt[0]);
        }
        else{
          if(angle_opt.size())
            filter2d.morphology(input,tmpout,Filter2d::DILATE,dimX_opt[0],dimY_opt[0],disc_opt[0],angle_opt[0]);
          else
            filter2d.morphology(input,tmpout,Filter2d::DILATE,dimX_opt[0],dimY_opt[0],disc_opt[0]);
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
        filter1d.morphology(tmpin,output,Filter::ERODE,dimZ_opt[0]);
      }
      else{
        if(angle_opt.size())
          filter2d.morphology(tmpin,output,Filter2d::ERODE,dimX_opt[0],dimY_opt[0],disc_opt[0],angle_opt[0]);
        else
          filter2d.morphology(tmpin,output,Filter2d::ERODE,dimX_opt[0],dimY_opt[0],disc_opt[0]);
      }
      tmpin.close();
      if(remove(tmps.str().c_str( )) !=0){
        cerr << "could not remove " << tmps.str() << std::endl;
      }
      break;
    }
    case(Filter2d::OPEN):{//opening
      ostringstream tmps;
      tmps << "/tmp/erosion_" << getpid() << ".tif";
      ImgWriterGdal tmpout;
      tmpout.open(tmps.str(),input);
      if(dimZ_opt.size()){
        filter1d.morphology(input,tmpout,Filter::ERODE,dimZ_opt[0]);
      }
      else{
        if(angle_opt.size())
          filter2d.morphology(input,tmpout,Filter2d::ERODE,dimX_opt[0],dimY_opt[0],disc_opt[0],angle_opt[0]);
        else
          filter2d.morphology(input,tmpout,Filter2d::ERODE,dimX_opt[0],dimY_opt[0],disc_opt[0]);
      }
      tmpout.close();
      ImgReaderGdal tmpin;
      tmpin.open(tmps.str());
      if(dimZ_opt.size()){
        filter1d.morphology(tmpin,output,Filter::DILATE,dimZ_opt[0]);
      }
      else{
        if(angle_opt.size())
          filter2d.morphology(tmpin,output,Filter2d::DILATE,dimX_opt[0],dimY_opt[0],disc_opt[0],angle_opt[0]);
        else
          filter2d.morphology(tmpin,output,Filter2d::DILATE,dimX_opt[0],dimY_opt[0],disc_opt[0]);
      }
      tmpin.close();
      if(remove(tmps.str().c_str( )) !=0){
        cerr << "could not remove " << tmps.str() << std::endl;
      }
      break;
    }
    case(Filter2d::HOMOG):{//spatially homogeneous
      filter2d.doit(input,output,Filter2d::HOMOG,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      // filter2d.homogeneousSpatial(input,output,dimX_opt[0],disc_opt[0]);
      break;
    }
    case(Filter2d::HETEROG):{//spatially heterogeneous
      filter2d.doit(input,output,Filter2d::HETEROG,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    }
    case(Filter2d::SOBELX):{//Sobel edge detection in X
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
    case(Filter2d::SOBELY):{//Sobel edge detection in Y
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
    case(Filter2d::SOBELXY):{//Sobel edge detection in XY
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
    case(Filter2d::SOBELYX):{//Sobel edge detection in XY
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
    case(Filter2d::SMOOTH):{//Smoothing filter
      filter2d.smooth(input,output,dimX_opt[0],dimY_opt[0]);
      break;
    }
    case(Filter2d::SMOOTHNODATA):{//Smoothing filter
      filter2d.smoothNoData(input,output,dimX_opt[0],dimY_opt[0]);
      break;
    }
    case(Filter2d::DENSITY):{//estimation of forest density
      filter2d.doit(input,output,Filter2d::DENSITY,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    }
    case(Filter2d::ORDER):{//order centre pixel with respect to values in window
      assert(dimX_opt[0]);
      filter2d.doit(input,output,Filter2d::ORDER,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    }
    default:
      cerr << "filter function " << function_opt[0] << " not supported" << std::endl;
      return 1;
    }
  }
  input.close();
  output.close();
  return 0;
}
