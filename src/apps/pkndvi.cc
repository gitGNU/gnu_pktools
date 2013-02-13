/**********************************************************************
pkndvi.cc: program to calculate vegetation index image
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

int main(int argc, char *argv[])
{
  //command line options
  Optionpk<string> input_opt("i","input","input image file","");
  Optionpk<string> output_opt("o","output","output image file containing ndvi","");
  Optionpk<short> band_opt("b", "band", "Bands to be used for vegetation index (see rule option)", 0);
  Optionpk<string> rule_opt("r", "rule", "Rule for index. ndvi (b1-b0)/(b1+b0), ndvi2 (b1-b0)/(b2+b3), gvmi (b0+0.1)-(b1+0.02))/((b0+0.1)+(b1+0.02))), vari (b1-b2)/(b1+b2-b0), osavi, mcari, tcari, diff (b1-b0), scale, ratio.", "ndvi");
  Optionpk<double> invalid_opt("t", "invalid", "Mask value where image is invalid.", 0);
  Optionpk<int> flag_opt("f", "flag", "Flag value to put in image if not valid (0)", 0);
  Optionpk<string> colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)", "");
  Optionpk<string> description_opt("d", "description", "Set image description", "");
  Optionpk<double> minmax_opt("m", "minmax", "minimum and maximum values for ndvi (limit all values smaller/larger to min/max", 0);
  Optionpk<double> eps_opt("e", "eps", "epsilon, contraint division by zero", 0);
  Optionpk<double> scale_opt("s", "scale", "scale[0] is used for input, scale[1] is used for output: DN=scale[1]*ndvi+offset[1]", 1);
  Optionpk<double> offset_opt("off", "offset", "offset[0] is used for input, offset[1] is used for output (see also scale option", 0);
  Optionpk<string> otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image", "Byte");
  Optionpk<string> oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image", "GTiff");
  Optionpk<string> option_opt("co", "co", "options: NAME=VALUE [-co COMPRESS=LZW] [-co INTERLEAVE=BAND]");
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    rule_opt.retrieveOption(argc,argv);
    invalid_opt.retrieveOption(argc,argv);
    flag_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
    description_opt.retrieveOption(argc,argv);
    minmax_opt.retrieveOption(argc,argv);
    eps_opt.retrieveOption(argc,argv);
    scale_opt.retrieveOption(argc,argv);
    offset_opt.retrieveOption(argc,argv);
    otype_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
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

  if(scale_opt.size()<2){
    if(input_opt.size()<2)
      scale_opt.push_back(1);
    else
      scale_opt.push_back(scale_opt[0]);
  }
  if(verbose_opt[0])
    std::cout << scale_opt;
  if(offset_opt.size()<2){
    if(input_opt.size()<2)
      offset_opt.push_back(0);
    else
      offset_opt.push_back(offset_opt[0]);
  }

  if(verbose_opt[0])
    std::cout << offset_opt;
  int reqBand=0;
  if(rule_opt[0]=="scale")
    reqBand=1;
  else if(rule_opt[0]=="vari"||rule_opt[0]=="mcari"||rule_opt[0]=="tcari")
    reqBand=3;
  else if(rule_opt[0]=="ndvi2")
    reqBand=4;
  else
    reqBand=2;
  while(band_opt.size()<reqBand)//bands can be explicitly provided by user or
    band_opt.push_back(band_opt[0]);//default is to use band 0 for each input
  if(verbose_opt[0])
    std::cout << band_opt;

  //todo: a bit stupid to duplicate input reader, but it works
  while(input_opt.size()<reqBand)
    input_opt.push_back(input_opt[0]);
  if(verbose_opt[0])
    std::cout << input_opt;

  vector<ImgReaderGdal> inputReader(reqBand);
  for(int ifile=0;ifile<reqBand;++ifile){
    inputReader[ifile].open(input_opt[ifile]);
    assert(inputReader[ifile].nrOfBand()>band_opt[ifile]);
  }

  if(verbose_opt[0]){
    cout << "opening output image file " << output_opt[0] << endl;
    cout << "data type: " << otype_opt[0] << endl;
  }
  //create output image with user defined data type 
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
    theType=inputReader[0].getDataType();
  if(verbose_opt[0])
    cout << endl << "Output pixel type:  " << GDALGetDataTypeName(theType) << endl;

  ImgWriterGdal outputWriter;
  if(verbose_opt[0])
    cout << "opening output image file " << output_opt[0] << endl;

  if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
    string theInterleave="INTERLEAVE=";
    theInterleave+=inputReader[0].getInterleave();
    option_opt.push_back(theInterleave);
  }
  outputWriter.open(output_opt[0],inputReader[0].nrOfCol(),inputReader[0].nrOfRow(),1,theType,oformat_opt[0],option_opt);

  if(description_opt[0]!="")
      outputWriter.setImageDescription(description_opt[0]);
  //if input image is georeferenced, copy projection info to output image
  if(inputReader[0].isGeoRef()){
    outputWriter.setProjection(inputReader[0].getProjection());
    double ulx,uly,lrx,lry;
    inputReader[0].getBoundingBox(ulx,uly,lrx,lry);
    outputWriter.copyGeoTransform(inputReader[0]);
  }
  if(colorTable_opt[0]!=""){
    if(colorTable_opt[0]!="none")
      outputWriter.setColorTable(colorTable_opt[0]);
  }
  else if (inputReader[0].getColorTable()!=NULL)//copy colorTable from first input image
    outputWriter.setColorTable(inputReader[0].getColorTable());
  
  Vector2d<double> lineInput(reqBand,inputReader[0].nrOfCol());
  vector<double> lineOutput(outputWriter.nrOfCol());

  int irow=0;
  int icol=0;
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  float progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(irow=0;irow<inputReader[0].nrOfRow();++irow){
    //read line in lineInput buffer
    try{
      if(rule_opt[0]=="scale")
        inputReader[0].readData(lineInput[0],GDT_Float64,irow,band_opt[0]);
      else if(rule_opt[0]=="vari"||rule_opt[0]=="tcari"){
        inputReader[0].readData(lineInput[0],GDT_Float64,irow,band_opt[0]);
        inputReader[1].readData(lineInput[1],GDT_Float64,irow,band_opt[1]);
        inputReader[2].readData(lineInput[2],GDT_Float64,irow,band_opt[2]);
      }
      else if(rule_opt[0]=="ndvi2"){
        inputReader[0].readData(lineInput[0],GDT_Float64,irow,band_opt[0]);
        inputReader[1].readData(lineInput[1],GDT_Float64,irow,band_opt[1]);
        inputReader[2].readData(lineInput[2],GDT_Float64,irow,band_opt[2]);
        inputReader[3].readData(lineInput[3],GDT_Float64,irow,band_opt[3]);
      }
      else{
        inputReader[0].readData(lineInput[0],GDT_Float64,irow,band_opt[0]);
        inputReader[1].readData(lineInput[1],GDT_Float64,irow,band_opt[1]);
      }
    }
    catch(string errorstring){
      cerr << errorstring << endl;
      exit(1);
    }
    assert(invalid_opt.size()==flag_opt.size());
    for(icol=0;icol<inputReader[0].nrOfCol();++icol){
      double ndvi=minmax_opt[0];
      double flagValue=flag_opt[0];
      bool valid=true;
      for(int iflag=0;valid&&iflag<invalid_opt.size();++iflag){
        for(int iband=0;iband<lineInput.size();++iband){
          if(lineInput[iband][icol]==invalid_opt[iflag]){
            flagValue=flag_opt[iflag];
            valid=false;
            break;
          }
        }
      }
      double denom;
      double nom;
      if(valid){
        if(rule_opt[0]=="ndvi"){
          //Example of indices addressed by ndvi:
          //structural indices
          //NDVI (Rouse1974): b0=b_680, b1=b_800
          //Chlorophyll indices:
          //Normalized Phaeophytinization index (NPQI Barnes1992): b0=R_435, b1=R_415
          //Photochemical Reflectance index (PRI1 Gamon1992): b0=R_567, b1=R_528
          //Photochemical Reflectance index (PRI2 Gamon1992): b0=R_570, b1=R_531
          //Normalized Phaeophytinization index (NPQI Barnes1992): b0=R_435, b1=R_415
          //Normalized Pigment Chlorophyll index (NPCI Penuelas1994): b0=R_430, b1=R_680
          //Structure Intensive Pigment index (SIPI Penuelas 1995): b0=R_450, b1=R_800
          //Lichtenthaler index 1 (Lic1 Lichtenthaler1996): b0=R_680, b2=R_800
          denom=(lineInput[1][icol]-offset_opt[0])/scale_opt[0]-(lineInput[0][icol]-offset_opt[0])/scale_opt[0];
          nom=(lineInput[1][icol]-offset_opt[0])/scale_opt[0]+(lineInput[0][icol]-offset_opt[0])/scale_opt[0];
        }
        else if(rule_opt[0]=="ndvi2"){//normalized difference with different wavelengths used in denom and nom
          //Example of indices addressed by ndvi2
          //Structure Intensive Pigment index (SIPI Penuelas 1995): b0=R_450, b1=R_800, b2=R_650, b=R_800
          //Vogelmann index 2 (Vog2 Vogelmann1993): b0=R_747, b1=R_735, b2=R_715, b3=R_726
          //Vogelmann index 3 (Vog3 Vogelmann1993): b0=R_747, b1=R_734, b2=R_715, b3=R_720
          denom=(lineInput[1][icol]-offset_opt[0])/scale_opt[0]-(lineInput[0][icol]-offset_opt[0])/scale_opt[0];
          nom=(lineInput[2][icol]-offset_opt[0])/scale_opt[0]+(lineInput[3][icol]-offset_opt[0])/scale_opt[0];
        }
        else if(rule_opt[0]=="gvmi"){
          denom=((lineInput[0][icol]-offset_opt[0])/scale_opt[0]+0.1)-((lineInput[1][icol]-offset_opt[0])/scale_opt[0]+0.02);
          nom=((lineInput[0][icol]-offset_opt[0])/scale_opt[0]+0.1)+((lineInput[1][icol]-offset_opt[0])/scale_opt[0]+0.02);
        }
        else if(rule_opt[0]=="vari"){
          denom=(lineInput[1][icol]-offset_opt[0])/scale_opt[0]-(lineInput[2][icol]-offset_opt[0])/scale_opt[0];
          nom=(lineInput[1][icol]-offset_opt[0])/scale_opt[0]+(lineInput[2][icol]-offset_opt[0])/scale_opt[0]-(lineInput[0][icol]-offset_opt[0])/scale_opt[0];
        }
        else if(rule_opt[0]=="osavi"){//structural index (Rondeaux1996): //b0=R_670, b1=R_800
          denom=(1.0+0.16)*(lineInput[1][icol]-offset_opt[0])/scale_opt[0]-(lineInput[0][icol]-offset_opt[0])/scale_opt[0];
          nom=(lineInput[1][icol]-offset_opt[0])/scale_opt[0]+(lineInput[0][icol]-offset_opt[0])/scale_opt[0]+0.16;
        }
        else if(rule_opt[0]=="mcari"){//chlorophyll index (Daughtry2000): b0=R_550, b1=R_670, b2=R_700
          denom=((lineInput[2][icol]-offset_opt[0])/scale_opt[0]-(lineInput[1][icol]-offset_opt[0])/scale_opt[0]-0.2*((lineInput[2][icol]-offset_opt[0])/scale_opt[0]-(lineInput[0][icol]-offset_opt[0])/scale_opt[0]))*(lineInput[2][icol]-offset_opt[0])/scale_opt[0];
          nom=(lineInput[1][icol]-offset_opt[0])/scale_opt[0];
        }
        else if(rule_opt[0]=="tcari"){//chlorophyll index (Haboudane2002): b0=R_550, b1=R_670, B2=R_700
          denom=3*((lineInput[1][icol]-offset_opt[0])/scale_opt[0]*(lineInput[2][icol]-offset_opt[0])/scale_opt[0]-(lineInput[1][icol]-offset_opt[0])/scale_opt[0]-0.2*((lineInput[2][icol]-offset_opt[0])/scale_opt[0]-(lineInput[0][icol]-offset_opt[0])/scale_opt[0])*(lineInput[2][icol]-offset_opt[0])/scale_opt[0]);
          nom=(lineInput[1][icol]-offset_opt[0])/scale_opt[0];
        }
        else if(rule_opt[0]=="diff"){
          denom=(lineInput[1][icol]-offset_opt[0])/scale_opt[0]-(lineInput[0][icol]-offset_opt[0])/scale_opt[0];
          nom=1.0;
        }
        else if(rule_opt[0]=="scale"){
          denom=(lineInput[0][icol]-offset_opt[0])/scale_opt[0];
          nom=1.0;
        }
        else if(rule_opt[0]=="ratio"){
          //Examples of indices addressed by ratio:
          //structural indices:
          //Simple Ratio Index (SR Jordan1969, Rouse1974): b0=R_NIR/R_RED
          //chlorophyll indices:
          //Greenness Index: b0=R_554, b1=R_677; 
          //Zarco-Tejada&Miller (Zarco2001): b0=R_750,b1=R_710
          //Simple Red Pigment Index (SRPI Penuelas1995): b0=R_430, b1=R_680
          //Carter index 1 (Ctr1 Carter1994): b0=R_695, b1=R_420
          //Carter index 2 (Ctr2 Carter1994): b0=R_695, b1=R_760
          //Lichtenthaler index 2 (Lic2 Lichtenthaler1996): b0=R_440, b2=R_690
          //Vogelmann index 1 (Vog1 Vogelmann1993): b0=R_740, b1=R_720
          //Gitelson and Merzlyak 1 (GM1 Gitelson1997): b0=R_750 b1=R_550
          //Gitelson and Merzlyak (GM2 Gitelson1997) b0=R_750 b1=R_700
          denom=(lineInput[0][icol]-offset_opt[0])/scale_opt[0];
          nom=(lineInput[1][icol]-offset_opt[0])/scale_opt[0];
        }
        else{
          std::cout << "Error: rule " << rule_opt[0] << " not supported" << std::endl;
          exit(1);
        }
        if(nom>eps_opt[0]||nom<-eps_opt[0])
        ndvi=denom/nom;
        if(ndvi<minmax_opt[0])
          ndvi=minmax_opt[0];
        else if(minmax_opt.size()>1){
          if(ndvi>minmax_opt[1])
            ndvi=minmax_opt[1];
        }
        switch(theType){
        case(GDT_Byte):
        case(GDT_Int16):
        case(GDT_UInt16):
        case(GDT_UInt32):
        case(GDT_Int32):
          lineOutput[icol]=static_cast<int>(0.5+ndvi*scale_opt[1]+offset_opt[1]);
          break;
        default:
          lineOutput[icol]=ndvi*scale_opt[1]+offset_opt[1];
        break;
        }
      }
      else
        lineOutput[icol]=flagValue;
    }
    //write buffer lineOutput to output file
    try{
      outputWriter.writeData(lineOutput,GDT_Float64,irow);
    }
    catch(string errorstring){
      cerr << errorstring << endl;
      exit(1);
    }
    //progress bar
    progress=static_cast<float>(irow+1.0)/outputWriter.nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
  for(int ifile=0;ifile<inputReader.size();++ifile)
    inputReader[ifile].close();
  outputWriter.close();
}
