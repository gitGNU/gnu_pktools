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
#include "imageclasses/ImgRaster.h"
#include "imageclasses/ImgCollection.h"
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

  if(doProcess){
    if(input_opt.empty()){
      std::cerr << "Error: no input file provided (use option -i). Use --help for help information" << std::endl;
      exit(1);
    }
    else if(output_opt.empty()){
      std::cerr << "Error: no output file provided (use option -o). Use --help for help information" << std::endl;
      exit(1);
    }
  }

  app::AppFactory app(argc,argv);
  app.setOption("dx","100");
  app.setOption("dy","100");
  
  filter2d::Filter2d filter;
  
  try{
    ImgCollection imgCollection;
    shared_ptr<ImgRaster> inputImg(new ImgRaster);
    inputImg->open(input_opt[0]);
    imgCollection.pushImage(inputImg);

    double theMin=0;
    double theMax=0;
    imgCollection[0]->getMinMax(theMin,theMax,0);
    cout << "min, max: " << theMin << ", " << theMax << endl;
    cout << "file1: " << imgCollection[0]->getFileName() << endl;
    shared_ptr<ImgRaster> imgRaster(new ImgRaster());
    app.showOptions();
    cout << "file2: " << imgCollection[0]->getFileName() << endl;
    imgCollection.crop(imgRaster,app);
    filter.smooth(*imgRaster,*imgRaster,5);
    // filter.morphology(imgRaster,imgRaster,"erode",3,3);
    // filter.morphology(imgRaster,imgRaster,"dilate",3,3);

    string imageType;
    if(oformat_opt.size())//default
      imageType=oformat_opt[0];
    else
      imageType=imgCollection[0]->getImageType();

    // if(projection_opt.size())
    //   imgRaster.setProjectionProj4(projection_opt[0]);
    // else if(imgCollection[0].getProjection()!="")
    //   imgRaster.setProjection(imgCollection[0].getProjection());

    cout << "imgRaster->nrOfCol(): " << imgRaster->nrOfCol() << endl;
    imgRaster->setFile(output_opt[0],imageType);
    imgRaster->close();
    inputImg->close();
  }
  catch(string helpString){//help was invoked
    std::cout << helpString << std::endl;
    return(1);
  }
  return(0);
}
  
