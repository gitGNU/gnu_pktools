/**********************************************************************
pkascii2img.cc: program to create raster image based on ascii file
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
#include <string>
#include <fstream>
#include "base/Optionpk.h"
#include <assert.h>
#include "imageclasses/ImgWriterGdal.h"

using namespace std;

int main(int argc, char *argv[])
{
  Optionpk<std::string> input_opt("i","input","input ASCII file");
  Optionpk<string> output_opt("o", "output", "Output image file");
  Optionpk<string> dataType_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image","Byte");
  Optionpk<string> imageType_opt("of", "oformat", "image type string (see also gdal_translate)", "GTiff");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<double> ulx_opt("ulx", "ulx", "Upper left x value bounding box (in geocoordinates if georef is true)", 0.0);
  Optionpk<double> uly_opt("uly", "uly", "Upper left y value bounding box (in geocoordinates if georef is true)", 0.0);
  Optionpk<double> dx_opt("dx", "dx", "Output resolution in x (in meter)");
  Optionpk<double> dy_opt("dy", "dy", "Output resolution in y (in meter)");
  Optionpk<string> projection_opt("a_srs", "a_srs", "Override the projection for the output file");
  Optionpk<string> colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)");
  Optionpk<string> description_opt("d", "description", "Set image description");
  Optionpk<bool> verbose_opt("v", "verbose", "verbose", false);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    dataType_opt.retrieveOption(argc,argv);
    imageType_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    ulx_opt.retrieveOption(argc,argv);
    uly_opt.retrieveOption(argc,argv);
    dx_opt.retrieveOption(argc,argv);
    dy_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
    projection_opt.retrieveOption(argc,argv);
    description_opt.retrieveOption(argc,argv);
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

  assert(input_opt.size());
  assert(output_opt.size());
  ImgWriterGdal imgWriter;
  ifstream ifile(input_opt[0].c_str(),ios::in);
  //get number of lines
  string line;
  int nrow=0;
  int ncol=0;
  int irow=0;
  string interleave="BAND";
  vector< vector<double> > data;
  vector<double> row;
  double value;
  try{
    while(getline(ifile,line)){
      row.clear();
      //read data from ascii file
      istringstream ist(line);
      while(ist>>value)
        row.push_back(value);
      if(!ncol){
        ncol=row.size();
        data.push_back(row);
      }
      else
        data.push_back(row);
      if(verbose_opt[0]){
        for(int icol=0;icol<row.size();++icol)
          cout << row[icol] << " ";
        cout << endl;
      }
      ++irow;
    }
    nrow=irow;
    assert(ncol);
    assert(nrow);
    if(verbose_opt[0]){
      cout << "nrow: " << nrow << endl;
      cout << "ncol: " << ncol << endl;
    }
  }
  catch(string theError){
    cout << theError << endl;
  }

  GDALDataType dataType=GDT_Unknown;
  if(verbose_opt[0])
    cout << "possible output data types: ";
  for(int iType = 0; iType < GDT_TypeCount; ++iType){
    if(verbose_opt[0])
      cout << " " << GDALGetDataTypeName((GDALDataType)iType);
    if( GDALGetDataTypeName((GDALDataType)iType) != NULL
        && EQUAL(GDALGetDataTypeName((GDALDataType)iType),
                 dataType_opt[0].c_str()))
      dataType=(GDALDataType) iType;
  }
  if(verbose_opt[0])
    cout << endl;
  if(verbose_opt[0]){
    if(dataType==GDT_Unknown)
      cout << "Unknown output pixel type: " << dataType_opt[0] << endl;
    else
      cout << "Output pixel type:  " << GDALGetDataTypeName(dataType) << endl;
  }

  imgWriter.open(output_opt[0],ncol,nrow,1,dataType,imageType_opt[0],option_opt);
  if(description_opt.size())
    imgWriter.setImageDescription(description_opt[0]);
  if(projection_opt.size()){
    assert(dx_opt.size());
    assert(dy_opt.size());
    if(verbose_opt[0])
      cout << output_opt[0] << " is georeferenced." << endl;
    double gt[6];
    gt[0]=ulx_opt[0];
    gt[1]=dx_opt[0];
    gt[2]=0;
    gt[3]=uly_opt[0];
    gt[4]=0;
    gt[5]=-dy_opt[0];
    imgWriter.setGeoTransform(gt);
    imgWriter.setProjectionProj4(projection_opt[0]);
  }
  else{
    if(verbose_opt[0])
      cout << output_opt[0] << " is not georeferenced." << endl;
    assert(!imgWriter.isGeoRef());
  }
  if(colorTable_opt.size()){
    assert(imgWriter.getDataType()==GDT_Byte);
    imgWriter.setColorTable(colorTable_opt[0]);
  }
  assert(data.size()==nrow);
  for(irow=0;irow<nrow;++irow)
    imgWriter.writeData(data[irow],GDT_Float64,irow);
  imgWriter.close();
}

