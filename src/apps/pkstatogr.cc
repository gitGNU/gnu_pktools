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
#include "algorithms/StatFactory.h"

int main(int argc, char *argv[])
{
  Optionpk<string> input_opt("i", "input", "Input shape file", "");
  Optionpk<string> fieldname_opt("n", "fname", "fields on which to calculate statistics", "");
  Optionpk<bool> minmax_opt("mm","minmax","calculate minimum and maximum value",false);
  Optionpk<double> min_opt("min","min","set minimum value",0);
  Optionpk<double> max_opt("max","max","set maximum value",0);
  Optionpk<bool> histogram_opt("hist","hist","calculate histogram",false);
  Optionpk<short> nbin_opt("nbin", "nbin", "number of bins", 0);
  Optionpk<bool> relative_opt("rel","relative","use percentiles for histogram to calculate histogram",false);
  Optionpk<bool> mean_opt("mean","mean","calculate mean value",false);
  Optionpk<bool> median_opt("median","median","calculate median value",false);
  Optionpk<bool> stdev_opt("stdev","stdev","calculate standard deviation",false);
  Optionpk<bool> size_opt("s","size","sample size (number of points)",false);
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    fieldname_opt.retrieveOption(argc,argv);
    minmax_opt.retrieveOption(argc,argv);
    min_opt.retrieveOption(argc,argv);
    max_opt.retrieveOption(argc,argv);
    histogram_opt.retrieveOption(argc,argv);
    nbin_opt.retrieveOption(argc,argv);
    relative_opt.retrieveOption(argc,argv);
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
  statfactory::StatFactory stat;
  //todo: implement ALL

  for(int ifield=0;ifield<fieldname_opt.size();++ifield){
    if(verbose_opt[0])
      std::cout << "field: " << ifield << std::endl;
    theData.clear();
    inputReader.readData(theData,OFTReal,fieldname_opt[ifield],0,verbose_opt[0]);
    vector<int> binData;
    double minValue=min_opt[0];
    double maxValue=max_opt[0];
    if(histogram_opt[0]){
      if(nbin_opt[0]<1){
        if(maxValue<=minValue)
          stat.minmax(theData,theData.begin(),theData.end(),minValue,maxValue);
        nbin_opt[0]=maxValue-minValue+1;
      }
      assert(nbin_opt[0]);
      try{
        stat.distribution(theData,theData.begin(),theData.end(),binData,nbin_opt[0],minValue,maxValue);
      }
      catch(string theError){
        std::cerr << "Warning: all identical values in data" << std::endl;
        exit(1);
      }
    }
    // int nbin=(nbin_opt[0]>1)? nbin_opt[0] : 2;
    std::cout << " --fname " << fieldname_opt[ifield];
    try{
      double theMean=0;
      double theVar=0;
      stat.meanVar(theData,theMean,theVar);
      if(mean_opt[0])
        std::cout << " --mean " << theMean;
      if(stdev_opt[0])
        std::cout << " --stdev " << sqrt(theVar);
      if(minmax_opt[0]){
        cout << " -min " << stat.min(theData);
        cout << " -max " << stat.max(theData);
      }
      if(median_opt[0])
        std::cout << " -median " << stat.median(theData);
      if(size_opt[0])
        std::cout << " -size " << theData.size();
      std::cout << std::endl;
      if(histogram_opt[0]){
        for(int bin=0;bin<nbin_opt[0];++bin){
          if(relative_opt[0])
            std::cout << (maxValue-minValue)*bin/(nbin_opt[0]-1)+minValue << " " << 100.0*static_cast<double>(binData[bin])/theData.size() << std::endl;
          else
            std::cout << (maxValue-minValue)*bin/(nbin_opt[0]-1)+minValue << " " << binData[bin] << std::endl;
        }
      }
    }
    catch(string theError){
      if(mean_opt[0])
        std::cout << " --mean " << theData.back();
      if(stdev_opt[0])
        std::cout << " --stdev " << "0";
      if(min_opt[0])
        cout << " -min " << theData.back();
      if(max_opt[0])
        cout << " -max " << theData.back();
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

