/**********************************************************************
pkdsm2shadow.cc: program to calculate sun shadow based on digital surface model and sun angles
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
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"

using namespace std;

/*------------------
  Main procedure
  ----------------*/
int main(int argc,char **argv) {
  Optionpk<std::string> input_opt("i","input","input image file");
  Optionpk<std::string> output_opt("o", "output", "Output image file");
  Optionpk<double> sza_opt("sza", "sza", "Sun zenith angle.");
  Optionpk<double> saa_opt("saa", "saa", "Sun azimuth angle (N=0 E=90 S=180 W=270).");
  Optionpk<int> flag_opt("f", "flag", "Flag to put in image if pixel shadow", 0);
  Optionpk<std::string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image", "");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image");
  Optionpk<string>  colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    sza_opt.retrieveOption(argc,argv);
    saa_opt.retrieveOption(argc,argv);
    flag_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
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
  assert(output_opt.size());
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
    output.open(output_opt[0],input.nrOfCol(),input.nrOfRow(),input.nrOfBand(),theType,imageType,option_opt);
  }
  catch(string errorstring){
    cout << errorstring << endl;
    exit(4);
  }
  if(input.isGeoRef()){
    output.setProjection(input.getProjection());
    double gt[6];
    input.getGeoTransform(gt);
    output.setGeoTransform(gt);
  }
  if(input.getColorTable()!=NULL)
    output.setColorTable(input.getColorTable());

  filter2d::Filter2d filter2d;
  if(verbose_opt[0])
    std::cout<< "class values: ";
  if(colorTable_opt.size())
    output.setColorTable(colorTable_opt[0]);
  filter2d.shadowDsm(input,output,sza_opt[0],saa_opt[0],input.getDeltaX(),flag_opt[0]);
  input.close();
  output.close();
  return 0;
}
