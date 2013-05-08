/**********************************************************************
pksetchandelier.cc: program to apply model parameters for brdf correction
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
#include <string>
#include <list>
#include <iostream>
#include "base/PointData.h"
#include "imageclasses/ImgReaderOgr.h"
#include "imageclasses/ImgWriterOgr.h"
#include "Optionpk.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

int main(int argc, char *argv[])
{
  std::string versionString="version ";
  versionString+=VERSION;
  versionString+=", Copyright (C) 2008-2012 Pieter Kempeneers.\n\
   This program comes with ABSOLUTELY NO WARRANTY; for details type use option -h.\n\
   This is free software, and you are welcome to redistribute it\n\
   under certain conditions; use option --license for details.";
  Optionpk<bool> version_opt("\0","version",versionString,false);
  Optionpk<bool> license_opt("lic","license","show license information",false);
  Optionpk<bool> help_opt("h","help","shows this help info",false);
  Optionpk<bool> todo_opt("\0","todo","",false);
  Optionpk<string> input_opt("i", "input", "Reflectance input","");
  Optionpk<string> geom_opt("g", "geom", "Geometry (ogr vector) file","");
  Optionpk<string> output_opt("o", "output", "corrected image output","");
  Optionpk<string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image", "Float32");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image", "GTiff");
  Optionpk<string> option_opt("co", "co", "options: NAME=VALUE [-co COMPRESS=LZW] [-co INTERLEAVE=BAND]");
  Optionpk<unsigned short> band_opt("b", "band", "field name of Reflectance band in reflectance file",0);
  Optionpk<unsigned short> sza_opt("sza", "sza", "band number (starting from 0) for Sun Zenith Angle in geometry image file",3);
  // Optionpk<unsigned short> vza_opt("vza", "vza", "band number (starting from 0) for Sun Zenith Angle in geometry image file",5);
  Optionpk<unsigned short> saa_opt("saa", "saa", "band number (starting from 0) for Sun Azimuth Angle in geometry image file",6);
  Optionpk<unsigned short> vaa_opt("vaa", "vaa", "band number (starting from 0) for View Azimuth Angle in geometry image file",7);
  Optionpk<double> tau_opt("tau", "tau", "Optical depth",0.2);
  Optionpk<double> deltaT_opt("deltaT", "deltaT", "Exposure time",1.0);
  Optionpk<double> state_opt("s","state","state vector: -s k -s e -s a -s haze",0);
  Optionpk<bool> model_opt("m","model","forward model reflectance",false);
  Optionpk<double> flag_opt("f", "flag", "No value flag",1000);
  Optionpk<short> verbose_opt("v", "verbose", "verbose", 0);

  version_opt.retrieveOption(argc,argv);
  license_opt.retrieveOption(argc,argv);
  help_opt.retrieveOption(argc,argv);
  todo_opt.retrieveOption(argc,argv);

  input_opt.retrieveOption(argc,argv);
  geom_opt.retrieveOption(argc,argv);
  output_opt.retrieveOption(argc,argv);
  otype_opt.retrieveOption(argc,argv);
  oformat_opt.retrieveOption(argc,argv);
  option_opt.retrieveOption(argc,argv);
  band_opt.retrieveOption(argc,argv);
  sza_opt.retrieveOption(argc,argv);
  // vza_opt.retrieveOption(argc,argv);
  saa_opt.retrieveOption(argc,argv);
  vaa_opt.retrieveOption(argc,argv);
  tau_opt.retrieveOption(argc,argv);
  deltaT_opt.retrieveOption(argc,argv);
  state_opt.retrieveOption(argc,argv);
  model_opt.retrieveOption(argc,argv);
  flag_opt.retrieveOption(argc,argv);
  verbose_opt.retrieveOption(argc,argv);

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
    cout << "usage: pkbrdf -i input -o output -g geom [OPTIONS]" << endl;
    exit(0);
  }

  if(verbose_opt[0])
    std::cout << "state_opt.size(): " << state_opt.size() << std::endl;
  assert(state_opt.size()==4);
  ImgReaderGdal inputReader(input_opt[0]);
  ImgReaderGdal geomReader(geom_opt[0]);
  vector<double> inputBuffer(inputReader.nrOfCol());
  vector<double> szaBuffer(geomReader.nrOfCol());
  vector<double> saaBuffer(geomReader.nrOfCol());
  // vector<double> vzaBuffer(geomReader.nrOfCol());
  vector<double> vaaBuffer(geomReader.nrOfCol());
  vector<double> fiBuffer(geomReader.nrOfCol());
  vector<double> xBuffer(inputReader.nrOfCol());//x pos
  vector<double> yBuffer(inputReader.nrOfCol());//y pos
  ImgWriterGdal outputWriter;
  string imageType=inputReader.getImageType();
  if(oformat_opt[0]!="")//default
    imageType=oformat_opt[0];
  GDALDataType theType=GDT_Unknown;
  if(verbose_opt[0]){
    std::cout << "Image type: " << imageType << std::endl;
    std::cout << "possible output data types: ";
  }
  for(int iType = 0; iType < GDT_TypeCount; ++iType){
    if(verbose_opt[0])
      cout << " " << GDALGetDataTypeName((GDALDataType)iType);
    if( GDALGetDataTypeName((GDALDataType)iType) != NULL
        && EQUAL(GDALGetDataTypeName((GDALDataType)iType),
                 otype_opt[0].c_str()))
      theType=(GDALDataType) iType;
  }
  if(theType==GDT_Unknown)
    theType=inputReader.getDataType();

  if(verbose_opt[0]){
    std::cout << std::endl << "Output data type:  " << GDALGetDataTypeName(theType) << std::endl;
    std::cout << "opening output image for writing: " << output_opt[0] << std::endl;
  }
  if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
    string theInterleave="INTERLEAVE=";
    theInterleave+=inputReader.getInterleave();
    option_opt.push_back(theInterleave);
  }
  try{
    outputWriter.open(output_opt[0],inputReader.nrOfCol(),inputReader.nrOfRow(),1,theType,imageType,option_opt);
    outputWriter.setProjection(inputReader.getProjection());
    outputWriter.copyGeoTransform(inputReader);
  }
  catch(string errorstring){
    cout << errorstring << endl;
    exit(1);
  }
  vector<double> writeBuffer(outputWriter.nrOfCol());

  PointData::m_tau=tau_opt[0];
  PointData::m_deltaT=deltaT_opt[0];

  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(int irow=0;irow<inputReader.nrOfRow();++irow){
    try{
      geomReader.readData(szaBuffer,GDT_Float64,irow,sza_opt[0]);
      geomReader.readData(saaBuffer,GDT_Float64,irow,saa_opt[0]);
      // geomReader.readData(vzaBuffer,GDT_Float64,irow,vza_opt[0]);
      geomReader.readData(vaaBuffer,GDT_Float64,irow,vaa_opt[0]);
    }
    catch(string errorstring){
      cout << errorstring << endl;
      exit(1);
    }
    //convert angles to radiance values
    for(int icol=0;icol<inputReader.nrOfCol();++icol){
      inputReader.image2geo(icol,irow,xBuffer[icol],yBuffer[icol]);
      double theRelativeAzimuth=saaBuffer[icol]-vaaBuffer[icol];
      fiBuffer[icol]=theRelativeAzimuth;
    }
    assert(band_opt[0]<inputReader.nrOfBand());
    try{
      inputReader.readData(inputBuffer,GDT_Float64,irow,band_opt[0]);
    }
    catch(string errorstring){
      cout << errorstring << endl;
      exit(1);
    }
    PointData pd;
    State theState;
    for(int index=0;index<inputBuffer.size();++index){
      pd.setReflectance(inputBuffer[index]);
      pd.setCosSZA(szaBuffer[index]);
      pd.setCosFi(fiBuffer[index]);
      double centreX,centreY;
      inputReader.getCentrePos(centreX,centreY);
      double ulx=inputReader.getUlx();
      double uly=inputReader.getUly();
      double r0=((ulx-centreX)*(ulx-centreX)+(uly-centreY)*(uly-centreY));
      double r=((xBuffer[index]-centreX)*(xBuffer[index]-centreX)+(yBuffer[index]-centreY)*(yBuffer[index]-centreY));
      pd.setR(r/r0);
      theState.k=state_opt[0];
      theState.e=state_opt[1];
      theState.a=state_opt[2];
      theState.haze=state_opt[3];
      if(model_opt[0])
        writeBuffer[index]=pd.modelReflectance(theState);
      else
        writeBuffer[index]=pd.correctReflectance(theState);
    }
    try{
      outputWriter.writeData(writeBuffer,GDT_Float64,irow,0);
    }
    catch(string errorstring){
      cout << errorstring << endl;
      exit(1);
    }
    progress=(1.0+irow)/inputReader.nrOfRow();
    assert(progress>=0);
    assert(progress<=1);
    pfnProgress(progress,pszMessage,pProgressArg);
  }
  inputReader.close();
  geomReader.close();
  if(output_opt[0]!="")
    outputWriter.close();
}
