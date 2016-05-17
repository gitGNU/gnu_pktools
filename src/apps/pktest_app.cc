/**********************************************************************
pktest_app.cc: program to mosaic and composite geo-referenced images
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
#include <algorithm>
#include <vector>
#include <iostream>
#include <string>
#include "imageclasses/ImgRasterGdal.h"
#include "imageclasses/ImgReaderOgr.h"
#include "base/Vector2d.h"
#include "base/Optionpk.h"
#include "algorithms/StatFactory.h"
#include "algorithms/Filter2d.h"
#include "AppFactory.h"
using namespace std;

int main(int argc, char *argv[])
{
  Optionpk<string>  input_opt("i", "input", "Input image file(s). If input contains multiple images, a multi-band output is created");
  Optionpk<string>  output_opt("o", "output", "Output image file");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate).","GTiff");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<string>  projection_opt("a_srs", "a_srs", "Override the spatial reference for the output file (leave blank to copy from input file, use epsg:3035 to use European projection and force to European grid");
  Optionpk<double> scale_opt("scale", "scale", "output=scale*input+offset");
  Optionpk<double> offset_opt("offset", "offset", "output=scale*input+offset");
  Optionpk<unsigned long int>  memory_opt("mem", "mem", "Buffer size (in MB) to read image data blocks in memory",0,1);

  option_opt.setHide(1);
  scale_opt.setHide(1);
  offset_opt.setHide(1);
  memory_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    projection_opt.retrieveOption(argc,argv);
    scale_opt.retrieveOption(argc,argv);
    offset_opt.retrieveOption(argc,argv);
    memory_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    cout << endl;
    cout << "Usage: pktest -i input -o output" << endl;
    cout << endl;
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }

  if(input_opt.empty()){
    std::cerr << "No input file provided (use option -i). Use --help for help information" << std::endl;
    exit(0);
  }

  if(output_opt.empty()){
    std::cerr << "No output file provided (use option -o). Use --help for help information" << std::endl;
    exit(0);
  }

  appfactory::AppFactory app;
  filter2d::Filter2d filter;

  app.setOptions(argc,argv);

  vector<ImgRasterGdal> imgReader(1);
  try{
    imgReader[0].open(input_opt[0],GA_ReadOnly,memory_opt[0]);
    for(int iband=0;iband<scale_opt.size();++iband)
      imgReader[0].setScale(scale_opt[iband],iband);
    for(int iband=0;iband<offset_opt.size();++iband)
      imgReader[0].setOffset(offset_opt[iband],iband);
  }
  catch(string errorstring){
    cerr << errorstring << " " << input_opt[0] << endl;
  }

  ImgRasterGdal imgRaster1;

  app.pkcrop(imgReader,imgRaster1);
  //todo: try inplace?
  filter.smooth(imgRaster1,imgRaster1,5);
  filter.morphology(imgRaster1,imgRaster1,"erode",3,3);
  filter.morphology(imgRaster1,imgRaster1,"dilate",3,3);
  string imageType;
  if(oformat_opt.size())//default
    imageType=oformat_opt[0];
  else
    imageType=imgReader[0].getImageType();

  if(projection_opt.size())
    imgRaster1.setProjectionProj4(projection_opt[0]);
  else if(imgReader[0].getProjection()!="")
    imgRaster1.setProjection(imgReader[0].getProjection());

  imgRaster1.setFile(output_opt[0],imageType);

  imgReader[0].close();

  imgRaster1.close();
  return 0;
}
  
