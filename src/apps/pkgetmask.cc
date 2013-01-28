/**********************************************************************
pkgetmask.cc: program to create mask image based on values in input raster image
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
#include <vector>
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"
#include "base/Optionpk.h"

using namespace std;
int main(int argc,char **argv) {
  Optionpk<string>  input_opt("i", "input", "Input image file", "");
  Optionpk<short>  band_opt("b", "band", "band(s) used for mask", 0);
  Optionpk<double> min_opt("min", "min", "Values smaller than min threshold(s) are masked as invalid. Use one threshold for each band", 0);
  Optionpk<double> max_opt("max", "max", "Values greater than max threshold(s) are masked as invalid. Use one threshold for each band", 0);
  Optionpk<string>  operator_opt("p", "operator", "Operator: [AND,OR].", "OR");
  Optionpk<unsigned short> valid_opt("t", "tvalue", "value(s) for valid pixels: between min and max", 1);
  Optionpk<unsigned short> invalid_opt("f", "fvalue", "value(s) for invalid pixels: not between min and max", 0);
  Optionpk<string> output_opt("o", "output", "Output mask file", "");
  Optionpk<string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image", "Byte");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image", "");
  Optionpk<string> option_opt("co", "co", "options: NAME=VALUE [-co COMPRESS=LZW] [-co INTERLEAVE=BAND]", "INTERLEAVE=BAND");
  Optionpk<string> colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)", "");
  Optionpk<short> verbose_opt("v", "verbose", "verbose", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    min_opt.retrieveOption(argc,argv);
    max_opt.retrieveOption(argc,argv);
    operator_opt.retrieveOption(argc,argv);
    valid_opt.retrieveOption(argc,argv);
    invalid_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    otype_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
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

  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
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
  if(verbose_opt[0]){
    cout << endl;
    if(theType==GDT_Unknown)
      cout << "Unknown output pixel type: " << otype_opt[0] << endl;
    else
      cout << "Output pixel type:  " << GDALGetDataTypeName(theType) << endl;
  }

  ImgReaderGdal imgReader(input_opt[0]);
  assert(band_opt.size()>=0);
  assert(band_opt.size()<=imgReader.nrOfBand());

  while(band_opt.size()>min_opt.size())
    min_opt.push_back(min_opt[0]);
  while(band_opt.size()>max_opt.size())
    max_opt.push_back(max_opt[0]);
  while(min_opt.size()>valid_opt.size())
    valid_opt.push_back(valid_opt[0]);
  assert(min_opt.size()==max_opt.size());
  if(verbose_opt[0]){
    cout << "min,max values: ";
    for(int imin=0;imin<min_opt.size();++imin){
      cout << min_opt[imin] << "," << max_opt[imin];
      if(imin<min_opt.size()-1)
	cout << " ";
      else
	cout << endl;
    }
  }
  
  vector< vector<float> > lineBuffer(band_opt.size());
  for(int iband=0;iband<band_opt.size();++iband)
    lineBuffer.resize(imgReader.nrOfCol());
  ImgWriterGdal imgWriter;
  //if output type not set, get type from input image
  if(theType==GDT_Unknown){
    theType=imgReader.getDataType();
    if(verbose_opt[0])
      cout << "Using data type from input image: " << GDALGetDataTypeName(theType) << endl;
  }
  string imageType=imgReader.getImageType();
  if(oformat_opt[0]!="")//default
    imageType=oformat_opt[0];
  if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
    string theInterleave="INTERLEAVE=";
    theInterleave+=imgReader.getInterleave();
    option_opt.push_back(theInterleave);
  }
  imgWriter.open(output_opt[0],imgReader.nrOfCol(),imgReader.nrOfRow(),1,theType,imageType,option_opt);
  if(colorTable_opt[0]!=""){
    if(colorTable_opt[0]!="none")
      imgWriter.setColorTable(colorTable_opt[0]);
  }
  else if (imgReader.getColorTable()!=NULL)//copy colorTable from input image
    imgWriter.setColorTable(imgReader.getColorTable());
  if(imgReader.isGeoRef()){
    imgWriter.setProjection(imgReader.getProjection());
    double ulx,uly,lrx,lry;
    imgReader.getBoundingBox(ulx,uly,lrx,lry);
    imgWriter.setGeoTransform(ulx,uly,imgReader.getDeltaX(),imgReader.getDeltaY(),0,0);
  }
  vector<char> writeBuffer(imgWriter.nrOfCol());
  for(int irow=0;irow<imgReader.nrOfRow();++irow){
    for(int iband=0;iband<band_opt.size();++iband)
      imgReader.readData(lineBuffer[iband],GDT_Float32,irow,band_opt[iband]);
    for(int icol=0;icol<imgReader.nrOfCol();++icol){
      bool valid=(operator_opt[0]=="OR")?false:true;
      unsigned short validValue=valid_opt[0];
      for(int ivalid=0;ivalid<min_opt.size();++ivalid){
        bool validBand=false;
      // for(int iband=0;iband<band_opt.size();++iband){
        unsigned short theBand=(band_opt.size()==min_opt.size())? ivalid:0;
        if(lineBuffer[theBand][icol]>=min_opt[ivalid]&&lineBuffer[theBand][icol]<=max_opt[ivalid]){
          validValue=valid_opt[ivalid];
          validBand=true;
        }
        valid=(operator_opt[0]=="OR")?valid||validBand : valid&&validBand;
      }
      if(valid)
	writeBuffer[icol]=validValue;
      else
	writeBuffer[icol]=invalid_opt[0];
    }
    imgWriter.writeData(writeBuffer,GDT_Byte,irow);
    progress=(1.0+irow)/imgWriter.nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }

  imgReader.close();
  imgWriter.close();
}
