/**********************************************************************
pkegcs.cc: Utility for raster files in European Grid Coordinate System
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
#include "base/Optionpk.h"
#include "imageclasses/ImgReaderGdal.h"
#include "algorithms/Egcs.h"

int main(int argc, char *argv[])
{
  Optionpk<std::string> image_opt("i","image","input image to analyse","");
  Optionpk<unsigned short>  band_opt("b", "band", "Band specific information", 0);
  Optionpk<std::string> cell2bb_opt("c2b","cell2bb","convert cell code to geo coordinates of boundingbox (e.g. 32-AB)","");
  Optionpk<std::string> cell2mid_opt("c2m","cell2mid","convert cell code to centre in geo coordinates (e.g. 32-AB)","");
  Optionpk<bool> refpixel_opt("\0", "ref", "get reference pixel (lower left corner of centre of gravity pixel)", false);
  Optionpk<double> maskValue_opt("m", "mask", "mask value(s) for no data to calculate reference pixel in image",0);
  Optionpk<int> dx_opt("dx","dx","resolution",250);
  Optionpk<bool> geo2cell_opt("g2c", "geo2cell", "get cell code for coordinates in x_opt and y_opt given the resolution in dx_opt", false);
  Optionpk<double> x_opt("x","x","x coordinate in epsg:3035",0);
  Optionpk<double> y_opt("y","y","y coordinate in epsg:3035",0);


  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=image_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    cell2bb_opt.retrieveOption(argc,argv);
    cell2mid_opt.retrieveOption(argc,argv);
    geo2cell_opt.retrieveOption(argc,argv);
    refpixel_opt.retrieveOption(argc,argv);
    maskValue_opt.retrieveOption(argc,argv);
    dx_opt.retrieveOption(argc,argv);
    x_opt.retrieveOption(argc,argv);
    y_opt.retrieveOption(argc,argv);
  }
  catch(std::string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }
  
  Egcs egcs;
  if(cell2bb_opt[0]!=""){
    int theULX, theULY, theLRX, theLRY;
    egcs.setLevel(egcs.cell2level(cell2bb_opt[0]));
    egcs.cell2bb(cell2bb_opt[0],theULX,theULY,theLRX,theLRY);
    std::cout << std::setprecision(12) << "--ulx=" << theULX << " --uly=" << theULY << " --lrx=" << theLRX << " --lry=" << theLRY << std::endl;
  }
  if(cell2mid_opt[0]!=""){
    double midX, midY;
    egcs.setLevel(egcs.cell2level(cell2mid_opt[0]));
    egcs.cell2mid(cell2mid_opt[0],midX,midY);
    std::cout << std::setprecision(12) << "-x=" << midX << " -y=" << midY << std::endl;
  }
  if(geo2cell_opt[0]){
    egcs.setLevel(egcs.res2level(dx_opt[0]));
    std::cout << egcs.geo2cell(x_opt[0],y_opt[0]) << std::endl;
  }
  if(image_opt[0]!=""){
    ImgReaderGdal imgReader;
    imgReader.open(image_opt[0]);
    if(refpixel_opt[0]){
      assert(band_opt[0]<imgReader.nrOfBand());
      for(int inodata=0;inodata<maskValue_opt.size();++inodata)
        imgReader.pushNoDataValue(maskValue_opt[inodata]);
      // if(verbose_opt[0]){
      //   vector<double> noData;
      //   imgReader.getNoDataValues(noData,band_opt[0]);
      //   std::cout << "number of no data values: " << noData.size() << std::endl;
      // }
      double refX,refY;
      //get centre of reference (centre of gravity) pixel in image
      imgReader.getRefPix(refX,refY,band_opt[0]);
      std::cout << std::setprecision(12) << "--x " << refX << " --y " << refY << std::endl;
      egcs.setLevel(egcs.res2level(imgReader.getDeltaX()));
      // unsigned short theLevel=egcs.getLevel(imgReader.getDeltaX());
      // egcs.setLevel(theLevel);
      std::cout << "cell code at level " << egcs.getLevel() << " (resolution is " << egcs.getResolution() << "): " << egcs.geo2cell(refX,refY) << std::endl;
    }
  }
}
