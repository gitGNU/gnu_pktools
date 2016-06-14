/**********************************************************************
pkcreatect.cc: program to create and import colour table to GTiff image
Copyright (C) 2008-2014 Pieter Kempeneers

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
#include "imageclasses/ImgRaster.h"
#include "base/Optionpk.h"

/******************************************************************************/
/*! \page pkcreatect pkcreatect
 program to create and import colour table to GTiff image
## SYNOPSIS

<code>
  Usage: pkcreatect -i input.txt -o output [-ct colortable | -min value -max value]
</code>

<code>

  
  Options: [--grey] [-of GDALformat] [-co option]* [-d description]

  Advanced options: [--legend filename [--dim cols --dim rows]] 
</code>

\section pkascii2ogr_description Description

Utility to include a color table to a raster dataset. You can either define an existing color table (ASCII text file) with the option -ct or define a minimum (-min) and maximum (-max) value.

\section pkcreatect_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |Input image file | 
 | o      | output               | std::string |       |Output image file | 
 | l      | legend               | std::string |       |Create legend as png file | 
 | dim    | dim                  | short | 100   |number of columns and rows in legend. | 
 | min    | min                  | double | 0     |minimum value | 
 | max    | max                  | double | 100   |maximum value | 
 | g      | grey                 | bool | false |grey scale | 
 | ct     | ct                   | std::string |       |color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid) | 
 | d      | description          | std::string |       |Set image description | 
 | of     | oformat              | std::string | GTiff |Output image format (see also gdal_translate).| 
 | co     | co                   | std::string |       |Creation option for output file. Multiple options can be specified. | 
 | mem    | mem                  | unsigned long int | 0 |Buffer size (in MB) to read image data blocks in memory | 

Usage: pkcreatect -i input.txt -o output [-ct colortable | -min value -max value]


Examples
========
Some examples how to use pkcreatect can be found \ref examples_pkcreatect "here"
**/

using namespace std;

int main(int argc,char **argv) {

  short red=-1;
  short green=-1;
  short blue=-1;

  Optionpk<string>  input_opt("i", "input", "Input image file");
  Optionpk<string>  output_opt("o", "output", "Output image file");
  Optionpk<string>  legend_opt("l", "legend", "Create legend as png file");
  Optionpk<short>  dim_opt("dim", "dim", "number of columns and rows in legend.", 100);
  Optionpk<double>  min_opt("min", "min", "minimum value", 0);
  Optionpk<double>  max_opt("max", "max", "maximum value", 100);
  Optionpk<bool>  grey_opt("g", "grey", "grey scale", false);
  Optionpk<string>  colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)");
  Optionpk<string> oformat_opt("of", "oformat", "Output image format (see also gdal_translate).", "GTiff");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<string>  description_opt("d", "description", "Set image description");
  Optionpk<unsigned long int>  memory_opt("mem", "mem", "Buffer size (in MB) to read image data blocks in memory",0,1);
  Optionpk<bool>  verbose_opt("v", "verbose", "verbose", false,2);
  
  legend_opt.setHide(1);
  dim_opt.setHide(1);
  memory_opt.setHide(1);
  
  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    legend_opt.retrieveOption(argc,argv);
    dim_opt.retrieveOption(argc,argv);
    min_opt.retrieveOption(argc,argv);
    max_opt.retrieveOption(argc,argv);
    grey_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
    description_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    memory_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    cout << endl;
    cout << "Usage: pkcreatect -i input.txt -o output [-ct colortable | -min value -max value]" << endl;
    cout << endl;
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }

  GDALColorTable colorTable;
  GDALColorEntry sEntry;
  if(colorTable_opt.empty()){
    sEntry.c4=255;
    for(int i=min_opt[0];i<=max_opt[0];++i){
      if(grey_opt[0]){
        sEntry.c1=255*(i-min_opt[0])/(max_opt[0]-min_opt[0]);
        sEntry.c2=255*(i-min_opt[0])/(max_opt[0]-min_opt[0]);
        sEntry.c3=255*(i-min_opt[0])/(max_opt[0]-min_opt[0]);
      }
      else{//hot to cold colour ramp
        sEntry.c1=255;
        sEntry.c2=255;
        sEntry.c3=255;
        double delta=max_opt[0]-min_opt[0];
        if(i<(min_opt[0]+0.25*delta)){
          sEntry.c1=0;
          sEntry.c2=255*4*(i-min_opt[0])/delta;
        }
        else if(i<(min_opt[0]+0.5*delta)){
          sEntry.c1=0;
          sEntry.c3=255*(1+4*(min_opt[0]+0.25*delta-i)/delta);
        }
        else if(i<(min_opt[0]+0.75*delta)){
          sEntry.c1=255*4*(i-min_opt[0]-0.5*delta)/delta;
          sEntry.c3=0;
        }
        else{
          sEntry.c2=255*(1+4*(min_opt[0]+0.75*delta-i)/delta);
          sEntry.c3=0;
        }
      }
      colorTable.SetColorEntry(i,&sEntry);
      if(output_opt.empty())
        cout << i << " " << sEntry.c1 << " " << sEntry.c2 << " " << sEntry.c3 << " " << sEntry.c4 << endl;
    }
  }
  ImgRaster legendWriter;
  short ncol=dim_opt[0];
  short nrow;
  if(legend_opt.size()){
    if(dim_opt.size()>1)
      nrow=dim_opt[1];
    else{
      nrow=max_opt[0]-min_opt[0]+1;
      ncol=dim_opt[0];
    }
    vector<string> pngOption;
    // pngOption.push_back("-co worldfile=no");
    pngOption.push_back("");
    legendWriter.open(legend_opt[0],ncol,nrow,1,GDT_Byte,oformat_opt[0],memory_opt[0],option_opt);
    if(colorTable_opt.size()){
      if(colorTable_opt[0]!="none")
        legendWriter.setColorTable(colorTable_opt[0]);
    }
    else
      legendWriter.setColorTable(&colorTable);
    if(legend_opt.size()){
      for(unsigned int irow=0;irow<legendWriter.nrOfRow();++irow){
        vector<char> buffer(legendWriter.nrOfCol());
        for(unsigned int icol=0;icol<legendWriter.nrOfCol();++icol)
          buffer[icol]=min_opt[0]+irow*static_cast<short>(max_opt[0]-min_opt[0]+1)/legendWriter.nrOfRow();
        legendWriter.writeData(buffer,legendWriter.nrOfRow()-1-irow);
      }
    }
  }

  // const char* pszMessage;
  // void* pProgressArg=NULL;
  // GDALProgressFunc pfnProgress=GDALTermProgress;
  // double progress=0;
  // pfnProgress(progress,pszMessage,pProgressArg);
  if(input_opt.size()&&output_opt.size()){
    ImgRaster imgReader(input_opt[0],memory_opt[0]);
    ImgRaster imgWriter;
    if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
      string theInterleave="INTERLEAVE=";
      theInterleave+=imgReader.getInterleave();
      option_opt.push_back(theInterleave);
    }

    imgWriter.open(output_opt[0],imgReader.nrOfCol(),imgReader.nrOfRow(),1,GDT_Byte,oformat_opt[0],memory_opt[0],option_opt);

    imgWriter.copyGeoTransform(imgReader);
    if(colorTable_opt.size()){
      if(colorTable_opt[0]!="none")
        imgWriter.setColorTable(colorTable_opt[0]);
    }
    else
      imgWriter.setColorTable(&colorTable);
    if(description_opt.size())
      imgWriter.setImageDescription(description_opt[0]);
    switch(imgReader.getDataType()){
    case(GDT_Byte):{
      vector<char> buffer;
      for(unsigned int irow=0;irow<imgReader.nrOfRow();++irow){
        imgReader.readData(buffer,irow);
        imgWriter.writeData(buffer,irow);
      }
      break;
    }
    case(GDT_Int16):{
      vector<short> buffer;
      cout << "Warning: copying short to unsigned short without conversion, use gdal_translate -scale if needed..." << endl;
      for(unsigned int irow=0;irow<imgReader.nrOfRow();++irow){
        imgReader.readData(buffer,irow);
        imgWriter.writeData(buffer,irow);
      }
      break;
    }
    case(GDT_UInt16):{
      vector<unsigned short> buffer;
      for(unsigned int irow=0;irow<imgReader.nrOfRow();++irow){
        imgReader.readData(buffer,irow);
        imgWriter.writeData(buffer,irow);
      }
      break;
    }
    default:
      cerr << "data type " << imgReader.getDataType() << " not supported for adding a colortable" << endl;
      break;
    }
    imgReader.close();
    imgWriter.close();
  }
  if(legend_opt.size())
    legendWriter.close();
}

