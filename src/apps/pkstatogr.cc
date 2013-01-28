/**********************************************************************
pkstatogr.cc: program to calculate basic statistics from vector file
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
#include <math.h>
#include <string>
#include <fstream>
#include <assert.h>
#include "base/Optionpk.h"
#include "imageclasses/ImgReaderOgr.h"
#include "algorithms/Histogram.h"

int main(int argc, char *argv[])
{
  Optionpk<string> input_opt("i", "input", "Input shape file", "");
  Optionpk<string> fieldname_opt("n", "fname", "fields on which to calculate statistics", "");
  Optionpk<short> nbin_opt("nbin", "nbin", "number of bins", 0);
  Optionpk<bool> min_opt("min","min","calculate minimum value",false);
  Optionpk<bool> max_opt("max","max","calculate maximum value",false);
  Optionpk<bool> mean_opt("mean","mean","calculate mean value",false);
  Optionpk<bool> median_opt("median","median","calculate median value",false);
  Optionpk<bool> stdev_opt("stdev","stdev","calculate standard deviation",false);
  Optionpk<bool> size_opt("s","size","sample size (number of points)",false);
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    fieldname_opt.retrieveOption(argc,argv);
    nbin_opt.retrieveOption(argc,argv);
    min_opt.retrieveOption(argc,argv);
    max_opt.retrieveOption(argc,argv);
    mean_opt.retrieveOption(argc,argv);
    median_opt.retrieveOption(argc,argv);
    stdev_opt.retrieveOption(argc,argv);
    size_opt.retrieveOption(argc,argv);
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

  ImgReaderOgr imgReader;
  try{
    imgReader.open(input_opt[0]);
  }
  catch(string errorstring){
    std::cerr << errorstring << std::endl;
  }

  ImgReaderOgr inputReader(input_opt[0]);
  vector<double> theData;
  Histogram hist;
  //todo: implement ALL

  for(int ifield=0;ifield<fieldname_opt.size();++ifield){
    if(verbose_opt[0])
      std::cout << "field: " << ifield << std::endl;
    theData.clear();
    inputReader.readData(theData,OFTReal,fieldname_opt[ifield],0,verbose_opt[0]);
    vector<int> binData;
    double minimum=0;
    double maximum=0;
    int nbin=(nbin_opt[0]>1)? nbin_opt[0] : 2;
    std::cout << " --fname " << fieldname_opt[ifield];
    try{
      hist.distribution(theData,theData.begin(),theData.end(),binData,nbin,minimum,maximum);
      double theMean=0;
      double theVar=0;
      hist.meanVar(theData,theMean,theVar);
      if(mean_opt[0])
        std::cout << " --mean " << theMean;
      if(stdev_opt[0])
        std::cout << " --stdev " << sqrt(theVar);
      if(min_opt[0])
        cout << " -m " << minimum;
      if(max_opt[0])
        cout << " -M " << maximum;
      if(median_opt[0])
        std::cout << " -median " << hist.median(theData);
      if(size_opt[0])
        std::cout << " -size " << theData.size();
      std::cout << std::endl;
      if(nbin_opt[0]>1){
        for(int bin=0;bin<nbin;++bin)
          std::cout << (maximum-minimum)*bin/(nbin-1)+minimum << " " << static_cast<double>(binData[bin])/theData.size() << std::endl;
      }
    }
    catch(string theError){
      if(mean_opt[0])
        std::cout << " --mean " << theData.back();
      if(stdev_opt[0])
        std::cout << " --stdev " << "0";
      if(min_opt[0])
        cout << " -m " << theData.back();
      if(max_opt[0])
        cout << " -M " << theData.back();
      if(median_opt[0])
        std::cout << " -median " << theData.back();
      if(size_opt[0])
        std::cout << " -size " << theData.size();
      std::cout << std::endl;
      std::cerr << "Warning: all identical values in data" << std::endl;
    }
  }
  imgReader.close();
}

