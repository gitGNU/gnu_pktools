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
  Optionpk<std::string> input_opt("i", "input", "Input shape file", "");
  Optionpk<std::string> fieldname_opt("n", "fname", "fields on which to calculate statistics", "");
  Optionpk<bool> minmax_opt("mm","minmax","calculate minimum and maximum value",false);
  Optionpk<bool> min_opt("min","min","calculate minimum value",0);
  Optionpk<bool> max_opt("max","max","calculate maximum value",0);
  Optionpk<double> src_min_opt("src_min","src_min","set minimum value for histogram");
  Optionpk<double> src_max_opt("src_max","src_max","set maximum value for histogram");
  Optionpk<double> nodata_opt("nodata","nodata","set nodata value(s)");
  Optionpk<bool> histogram_opt("hist","hist","calculate histogram",false);
  Optionpk<unsigned int> nbin_opt("nbin", "nbin", "number of bins");
  Optionpk<bool> relative_opt("rel","relative","use percentiles for histogram to calculate histogram",false);
  Optionpk<double> kde_opt("kde","kde","bandwith of kernel density when producing histogram, use 0 for practical estimation based on Silverman's rule of thumb. Leave empty if no kernel density is required");
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
    src_min_opt.retrieveOption(argc,argv);
    src_max_opt.retrieveOption(argc,argv);
    nodata_opt.retrieveOption(argc,argv);
    histogram_opt.retrieveOption(argc,argv);
    nbin_opt.retrieveOption(argc,argv);
    relative_opt.retrieveOption(argc,argv);
    kde_opt.retrieveOption(argc,argv);
    mean_opt.retrieveOption(argc,argv);
    median_opt.retrieveOption(argc,argv);
    stdev_opt.retrieveOption(argc,argv);
    size_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(std::string predefinedString){
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
  catch(std::string errorstring){
    std::cerr << errorstring << std::endl;
  }

  ImgReaderOgr inputReader(input_opt[0]);
  std::vector<double> theData;
  statfactory::StatFactory stat;
  //todo: implement ALL

  stat.setNoDataValues(nodata_opt);
  for(int ifield=0;ifield<fieldname_opt.size();++ifield){
    if(verbose_opt[0])
      std::cout << "field: " << ifield << std::endl;
    theData.clear();
    inputReader.readData(theData,OFTReal,fieldname_opt[ifield],0,verbose_opt[0]);
    std::vector<double> binData;
    double minValue=0;
    double maxValue=0;
    stat.minmax(theData,theData.begin(),theData.end(),minValue,maxValue);
    if(src_min_opt.size())
      minValue=src_min_opt[0];
    if(src_max_opt.size())
      maxValue=src_max_opt[0];
    unsigned int nbin=(nbin_opt.size())? nbin_opt[0]:0;

    if(histogram_opt[0]){
      double sigma=0;
      if(kde_opt.size()){
        if(kde_opt[0]>0)
          sigma=kde_opt[0];
        else
          sigma=1.06*sqrt(stat.var(theData))*pow(theData.size(),-0.2);
      }
      if(nbin<1)
        nbin=(maxValue-minValue+1);
      try{
        stat.distribution(theData,theData.begin(),theData.end(),binData,nbin,minValue,maxValue,sigma);
      }
      catch(std::string theError){
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
      if(minmax_opt[0]||min_opt[0]||max_opt[0]){
	if(minmax_opt[0])
	  std::cout << " --min " << minValue << " --max " << maxValue << " ";
	else{
	  if(min_opt[0])
	    std::cout << " --min " << minValue << " ";
	  if(max_opt[0])
	    std::cout << " --max " << maxValue << " ";
	}
      }
      if(median_opt[0])
        std::cout << " -median " << stat.median(theData);
      if(size_opt[0])
        std::cout << " -size " << theData.size();
      std::cout << std::endl;
      if(histogram_opt[0]){
        for(int ibin=0;ibin<nbin;++ibin){
	  double binValue=0;
	  if(nbin==maxValue-minValue+1)
	    binValue=minValue+ibin;
	  else
	    binValue=minValue+static_cast<double>(maxValue-minValue)*(ibin+0.5)/nbin;
	  std::cout << binValue << " ";
          if(relative_opt[0])
            std::cout << 100.0*static_cast<double>(binData[ibin])/theData.size() << std::endl;
          else
            std::cout << binData[ibin] << std::endl;
        }
      }
    }
    catch(std::string theError){
      if(mean_opt[0])
        std::cout << " --mean " << theData.back();
      if(stdev_opt[0])
        std::cout << " --stdev " << "0";
      if(min_opt[0])
        std::cout << " -min " << theData.back();
      if(max_opt[0])
        std::cout << " -max " << theData.back();
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

