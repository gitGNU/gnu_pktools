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

/*------------------
  Main procedure
  ----------------*/
int main(int argc,char **argv) {
  Optionpk<bool> version_opt("\0","version","version 20120625, Copyright (C) 2008-2012 Pieter Kempeneers.\n\
   This program comes with ABSOLUTELY NO WARRANTY; for details type use option -h.\n\
   This is free software, and you are welcome to redistribute it\n\
   under certain conditions; use option --license for details.",false);
  Optionpk<bool> license_opt("lic","license","show license information",false);
  Optionpk<bool> help_opt("h","help","shows this help info",false);
  Optionpk<bool> todo_opt("\0","todo","",false);
  Optionpk<std::string> input_opt("i","input","input image file","");
  Optionpk<std::string> output_opt("o", "output", "Output image file", "");
  Optionpk<bool> disc_opt("c", "circular", "circular disc kernel for dilation and erosion (default is false)", false);
  Optionpk<double> angle_opt("a", "angle", "angle used for directional filtering in dilation. Default is less than -180 (no directional effect)", -190);
  Optionpk<int> function_opt("f", "filter", "filter function (0: median, 1: variance, 2: min, 3: max, 4: sum, 5: mean, 6: minmax, 7: dilation, 8: erosion, 9: closing, 10: opening, 11: spatially homogeneous, 12: SobelX edge detection in X, 13: SobelY edge detection in Y, 14: SobelXY (not supported), 15: smooth, 16: density, 17: majority voting (only for classes), 18: forest aggregation (mixed), 19: smooth no data (mask) values), 20: threshold local filtering", 0);
  Optionpk<int> dimX_opt("dx", "dx", "filter kernel size in x, must be odd (example: 3). Default is 0", 0);
  Optionpk<int> dimY_opt("dy", "dy", "filter kernel size in y, must be odd (example: 3). Default is 0. Set dy=0 if 1-D filter must be used in band domain", 0);
  Optionpk<string> option_opt("co", "co", "options: NAME=VALUE [-co COMPRESS=LZW] [-co INTERLEAVE=BAND]", "INTERLEAVE=BAND");
  Optionpk<short> class_opt("class", "class", "class value(s) to use for density, erosion, dilation, openening and closing, thresholding (default is 1)", 1);
  Optionpk<double> threshold_opt("t", "threshold", "threshold value(s) to use for threshold filter (one for each class)", 0);
  Optionpk<short> mask_opt("\0", "mask", "mask value(s) ", -1);
  Optionpk<std::string> tap_opt("tap", "tap", "text file conttaining taps used for filtering (from ul to lr). Default is empty: do not use taps. Use dimX and dimY to specify tap dimensions in x and y", "");
  Optionpk<std::string> colorTable_opt("\0", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)", "");
  Optionpk<short> down_opt("d", "down", "down sampling factor. Default is 1: no downsample)", 1);
  Optionpk<short> verbose_opt("v", "verbose", "verbose (default is 0)", 0);

  version_opt.retrieveOption(argc,argv);
  license_opt.retrieveOption(argc,argv);
  help_opt.retrieveOption(argc,argv);
  todo_opt.retrieveOption(argc,argv);

  input_opt.retrieveOption(argc,argv);
  output_opt.retrieveOption(argc,argv);
  verbose_opt.retrieveOption(argc,argv);
  disc_opt.retrieveOption(argc,argv);
  angle_opt.retrieveOption(argc,argv);
  function_opt.retrieveOption(argc,argv);
  dimX_opt.retrieveOption(argc,argv);
  dimY_opt.retrieveOption(argc,argv);
  option_opt.retrieveOption(argc,argv);
  class_opt.retrieveOption(argc,argv);
  threshold_opt.retrieveOption(argc,argv);
  mask_opt.retrieveOption(argc,argv);
  tap_opt.retrieveOption(argc,argv);
  down_opt.retrieveOption(argc,argv);
  colorTable_opt.retrieveOption(argc,argv);

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

  output.open(output_opt[0],(input.nrOfCol()+down_opt[0]-1)/down_opt[0],(input.nrOfRow()+down_opt[0]-1)/down_opt[0],input.nrOfBand(),input.getDataType(),input.getImageType(),option_opt);
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
    if(dimY_opt[0]>0)
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
  if(tap_opt[0]!=""){
    ifstream tapfile(tap_opt[0].c_str());
    assert(tapfile);
    Vector2d<double> taps(dimY_opt[0],dimX_opt[0]);

    for(int j=0;j<dimY_opt[0];++j){
      for(int i=0;i<dimX_opt[0];++i){
        tapfile >> taps[j][i];
      }
    }
    if(verbose_opt[0]){
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
    case(Filter2d::MIN):
      filter2d.doit(input,output,Filter2d::MIN,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    case(Filter2d::MAX):
      filter2d.doit(input,output,Filter2d::MAX,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
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
      if(dimY_opt[0]>0)
        filter2d.morphology(input,output,Filter2d::DILATE,dimX_opt[0],dimY_opt[0],disc_opt[0],angle_opt[0]);
      else{
        if(verbose_opt[0])
          std::cout<< "1-D filtering: dilate" << std::endl;
        filter1d.morphology(input,output,Filter::DILATE,dimX_opt[0]);
      }
      break;
    case(Filter2d::ERODE):
      if(dimY_opt[0]>0)
        filter2d.morphology(input,output,Filter2d::ERODE,dimX_opt[0],dimY_opt[0],disc_opt[0],angle_opt[0]);
      else{
        if(verbose_opt[0])
          std::cout<< "1-D filtering: dilate" << std::endl;
        filter1d.morphology(input,output,Filter::ERODE,dimX_opt[0]);
      }
      break;
    case(Filter2d::CLOSE):{//closing
      ostringstream tmps;
      tmps << "/tmp/dilation_" << getpid() << ".tif";
      ImgWriterGdal tmpout;
      tmpout.open(tmps.str(),input);
      try{
        if(dimY_opt[0]>0)
          filter2d.morphology(input,tmpout,Filter2d::DILATE,dimX_opt[0],dimY_opt[0],disc_opt[0],angle_opt[0]);
        else
          filter1d.morphology(input,tmpout,Filter::DILATE,dimX_opt[0]);
      }
      catch(std::string errorString){
	std::cout<< errorString;
	exit(1);
      }
      tmpout.close();
      ImgReaderGdal tmpin;
      tmpin.open(tmps.str());
      if(dimY_opt[0]>0)
        filter2d.morphology(tmpin,output,Filter2d::ERODE,dimX_opt[0],dimY_opt[0],disc_opt[0],angle_opt[0]);
      else
        filter1d.morphology(tmpin,output,Filter::ERODE,dimX_opt[0]);
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
      if(dimY_opt[0]>0)
        filter2d.morphology(input,tmpout,Filter2d::ERODE,dimX_opt[0],dimY_opt[0],disc_opt[0],angle_opt[0]);
      else
        filter1d.morphology(input,tmpout,Filter::ERODE,dimX_opt[0]);
      tmpout.close();
      ImgReaderGdal tmpin;
      tmpin.open(tmps.str());
      if(dimY_opt[0]>0)
        filter2d.morphology(tmpin,output,Filter2d::DILATE,dimX_opt[0],dimY_opt[0],disc_opt[0],angle_opt[0]);
      else
        filter1d.morphology(tmpin,output,Filter::DILATE,dimX_opt[0]);
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
    default:
      cerr << "filter function " << function_opt[0] << " not supported" << std::endl;
      return 1;
    }
  }
  input.close();
  output.close();
  return 0;
}
