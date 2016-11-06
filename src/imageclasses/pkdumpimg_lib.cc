/**********************************************************************
pkdumpimg_lib.cc: dump image on screen or ASCII file
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
#include <string>
#include <iostream>
#include <memory>
#include "imageclasses/ImgRasterGdal.h"
#include "imageclasses/ImgReaderOgr.h"
#include "base/Optionpk.h"
#include "apps/AppFactory.h"

using namespace std;
using namespace app;

/**
 * @param imgWriter output raster dumpimg dataset
 * @return CE_None if successful, CE_Failure if failed
 **/
CPLErr ImgRasterGdal::dumpimg(const AppFactory& app){

  Optionpk<string> output_opt("o", "output", "Output ascii file (Default is empty: use stdout");
  Optionpk<string> oformat_opt("of", "oformat", "Output format (matrix form or list (x,y,z) form. Default is matrix form", "matrix");
  Optionpk<int> band_opt("b", "band", "band index to crop");
  Optionpk<short> dstnodata_opt("dstnodata", "dstnodata", "nodata value for output if out of bounds.", 0);
  Optionpk<double> srcnodata_opt("srcnodata", "srcnodata", "set no data value(s) for input image");
  Optionpk<short> verbose_opt("v", "verbose", "verbose (Default: 0)", 0,2);

  srcnodata_opt.setHide(1);
  dstnodata_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=output_opt.retrieveOption(app.getArgc(),app.getArgv());
    oformat_opt.retrieveOption(app.getArgc(),app.getArgv());
    band_opt.retrieveOption(app.getArgc(),app.getArgv());
    srcnodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    dstnodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    verbose_opt.retrieveOption(app.getArgc(),app.getArgv());
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    return(CE_Failure);
  }
  if(!doProcess){
    cout << endl;
    std::ostringstream helpStream;
    helpStream << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    //test
    helpStream << "throwing exception" << std::endl;
    throw(helpStream.str());//help was invoked, stop processing
  }
  ofstream outputStream;
  if(output_opt.size())
    outputStream.open(output_opt[0].c_str());

  for(int inodata=0;inodata<srcnodata_opt.size();++inodata)
    pushNoDataValue(srcnodata_opt[inodata]);

  if(band_opt.empty()){
    for(int iband=0;iband<nrOfBand();++iband)
      band_opt.push_back(iband);
  }
  std::vector<double> readBuffer(nrOfCol());
  for(int iband=0;iband<band_opt.size();++iband){
    assert(band_opt[iband]>=0);
    assert(band_opt[iband]<nrOfBand());
    for(int irow=0;irow<nrOfRow();++irow){
      readData(readBuffer,irow,band_opt[iband]);
      for(int icol=0;icol<nrOfCol();++icol){
        if(oformat_opt[0]=="matrix"){
          if(output_opt.empty())
            std::cout << readBuffer[icol] << " ";
          else
            outputStream << readBuffer[icol] << " ";
        }
        else if(!isNoData(readBuffer[icol])){
          if(output_opt.empty())
            std::cout << icol << " " << irow << " " << readBuffer[icol] << std::endl;
          else
            outputStream << icol << " " << irow << " " << readBuffer[icol] << std::endl;
        }
      }
      std::cout << std::endl;
    }
  }
  if(oformat_opt[0]=="matrix"){
    if(output_opt.empty())
      std::cout << std::endl;
    else
      outputStream << std::endl;
  }
  if(!output_opt.empty())
    outputStream.close();
  return(CE_None);
}
