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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

using namespace std;

int main(int argc, char *argv[])
{
  //command line options
  std::string versionString="version ";
  versionString+=VERSION;
  versionString+=", Copyright (C) 2008-2012 Pieter Kempeneers.\n\
   This program comes with ABSOLUTELY NO WARRANTY; for details type use option -h.\n\
   This is free software, and you are welcome to redistribute it\n\
   under certain conditions; use option --license for details.";
  Optionpk<bool> version_opt("\0","version",versionString,false);
  Optionpk<bool> license_opt("lic","license","show license information",false);
  Optionpk<bool> help_opt("h","help","shows this help info",false);
  Optionpk<string> input_opt("i","input","input image file","");
  Optionpk<string> output_opt("o","output","output image file containing ndvi","");
  Optionpk<short> band_opt("b", "band", "Bands to be used for vegetation index (see rule option)", 0);
  Optionpk<string> rule_opt("r", "rule", "Rule for index. [ndvi (b1-b0)/(b1+b0)|gvmi (b0+0.1)-(b1+0.02))/((b0+0.1)+(b1+0.02)))|vari (b1-b2)/(b1+b2-b0)|diff (b1-b0)|scale|ratio.", "ndvi");
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
  Optionpk<string> option_opt("co", "co", "options: NAME=VALUE [-co COMPRESS=LZW] [-co INTERLEAVE=BAND]", "INTERLEAVE=BAND");
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0);

  version_opt.retrieveOption(argc,argv);
  license_opt.retrieveOption(argc,argv);
  help_opt.retrieveOption(argc,argv);
  input_opt.retrieveOption(argc,argv);
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

  if(version_opt[0]){
    cout << version_opt.getHelp() << endl;
    exit(0);
  }
  if(license_opt[0]){
    cout << Optionpk<bool>::getGPLv3License() << endl;
    exit(0);
  }
  if(help_opt[0]){
    cout << "usage: pkinfo -i inputimage -o outputimage [OPTIONS]" << endl;
    exit(0);
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
  else if(rule_opt[0]=="vari")
    reqBand=3;
  else
    reqBand=2;
  while(band_opt.size()<reqBand)
    band_opt.push_back(band_opt[0]);
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
      else if(rule_opt[0]=="vari"){
        inputReader[0].readData(lineInput[0],GDT_Float64,irow,band_opt[0]);
        inputReader[1].readData(lineInput[1],GDT_Float64,irow,band_opt[1]);
        inputReader[2].readData(lineInput[2],GDT_Float64,irow,band_opt[2]);
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
          denom=(lineInput[1][icol]-offset_opt[0])/scale_opt[0]-(lineInput[0][icol]-offset_opt[0])/scale_opt[0];
          nom=(lineInput[1][icol]-offset_opt[0])/scale_opt[0]+(lineInput[0][icol]-offset_opt[0])/scale_opt[0];
        }
        else if(rule_opt[0]=="gvmi"){
          denom=((lineInput[0][icol]-offset_opt[0])/scale_opt[0]+0.1)-((lineInput[1][icol]-offset_opt[0])/scale_opt[0]+0.02);
          nom=((lineInput[0][icol]-offset_opt[0])/scale_opt[0]+0.1)+((lineInput[1][icol]-offset_opt[0])/scale_opt[0]+0.02);
        }
        else if(rule_opt[0]=="vari"){
          denom=(lineInput[1][icol]-offset_opt[0])/scale_opt[0]-(lineInput[2][icol]-offset_opt[0])/scale_opt[0];
          nom=(lineInput[1][icol]-offset_opt[0])/scale_opt[0]+(lineInput[2][icol]-offset_opt[0])/scale_opt[0]-(lineInput[0][icol]-offset_opt[0])/scale_opt[0];
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
