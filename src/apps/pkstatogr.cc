/**********************************************************************
pkstatogr.cc: program to calculate basic statistics from vector file
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
#include <math.h>
#include <string>
#include <fstream>
#include <assert.h>
#include "base/Optionpk.h"
#include "imageclasses/ImgReaderOgr.h"
#include "algorithms/StatFactory.h"
/******************************************************************************/
/*! \page pkstatogr pkstatogr
 program to calculate basic statistics from vector file
## SYNOPSIS

<code>
  Usage: pkstatogr -i input [-n attribute]*
</code>

<code>

  Options: [-ln layer]* [-n attribute]* [srcnodata]* [src_min] [src_max] [-s] [-mm] [-min] [-max] [-mean] [-median] [-stdev] [-hist] [-nbin] [-rel] [-kde]
</code>

\section pkstatogr_description Description

The utility pkstatogr calculates basic statistics on attributes of a vector file. Examples of the basic statistics include: minimum, maximum, median, mean and standard deviation. Histograms (in percentage or absolute values) can also be calculated. The attribute of interest can be selected using the option -n|--fname. Values defined by the -nodata option, or not withing the limits set by the options -src_min and -src_max are ignored for the statistics.
\section pkstatogr_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |Input OGR vector file | 
 | n      | fname                | std::string |       |Fields on which to calculate statistics | 
 | ln     | lname                | std::string |       |Layer name(s) in sample (leave empty to select all) | 
 | nodata | nodata               | double |       |Set nodata value(s) | 
 | src_min | src_min              | double |       |Set minimum value for histogram | 
 | src_max | src_max              | double |       |Set maximum value for histogram | 
 | s      | size                 | bool | false |Sample size (number of points) | 
 | mm     | minmax               | bool | false |Calculate minimum and maximum value | 
 | min    | min                  | bool | false |Calculate minimum value | 
 | max    | max                  | bool | false |Calculate maximum value | 
 | mean   | mean                 | bool | false |Calculate mean value | 
 | median | median               | bool | false |Calculate median value | 
 | stdev  | stdev                | bool | false |Calculate standard deviation | 
 | hist   | hist                 | bool | false |Calculate histogram | 
 | nbin   | nbin                 | unsigned int |       |Number of bins | 
 | rel    | relative             | bool | false |Use percentiles for histogram to calculate histogram | 
 | kde    | kde                  | bool | false |Use Kernel density estimation when producing histogram. The standard deviation is estimated based on Silverman's rule of thumb | 

Usage: pkstatogr -i input [-n attribute]*


Examples
========
Some examples how to use pkstatogr can be found \ref examples_pkstatogr "here"
**/

using namespace std;

int main(int argc, char *argv[])
{
  Optionpk<string> input_opt("i", "input", "Input OGR vector file", "");
  Optionpk<string> layer_opt("ln", "lname", "Layer name(s) in sample (leave empty to select all)");
  Optionpk<string> fieldname_opt("n", "fname", "Fields on which to calculate statistics", "");
  Optionpk<double> nodata_opt("nodata","nodata","Set nodata value(s)");
  Optionpk<double> src_min_opt("src_min","src_min","Set minimum value for histogram");
  Optionpk<double> src_max_opt("src_max","src_max","Set maximum value for histogram");
  Optionpk<bool> size_opt("s","size","Sample size (number of points)",false);
  Optionpk<bool> minmax_opt("mm","minmax","Calculate minimum and maximum value",false);
  Optionpk<bool> min_opt("min","min","Calculate minimum value",0);
  Optionpk<bool> max_opt("max","max","Calculate maximum value",0);
  Optionpk<bool> mean_opt("mean","mean","Calculate mean value",false);
  Optionpk<bool> median_opt("median","median","Calculate median value",false);
  Optionpk<bool> stdev_opt("stdev","stdev","Calculate standard deviation",false);
  Optionpk<bool> histogram_opt("hist","hist","Calculate histogram",false);
  Optionpk<unsigned int> nbin_opt("nbin", "nbin", "Number of bins");
  Optionpk<bool> relative_opt("rel","relative","Use percentiles for histogram to calculate histogram",false);
  Optionpk<bool> kde_opt("kde","kde","Use Kernel density estimation when producing histogram. The standard deviation is estimated based on Silverman's rule of thumb",false);
  Optionpk<short> verbose_opt("v", "verbose", "Verbose level", 0,2);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    fieldname_opt.retrieveOption(argc,argv);
    layer_opt.retrieveOption(argc,argv);
    nodata_opt.retrieveOption(argc,argv);
    src_min_opt.retrieveOption(argc,argv);
    src_max_opt.retrieveOption(argc,argv);
    size_opt.retrieveOption(argc,argv);
    minmax_opt.retrieveOption(argc,argv);
    min_opt.retrieveOption(argc,argv);
    max_opt.retrieveOption(argc,argv);
    mean_opt.retrieveOption(argc,argv);
    median_opt.retrieveOption(argc,argv);
    stdev_opt.retrieveOption(argc,argv);
    histogram_opt.retrieveOption(argc,argv);
    nbin_opt.retrieveOption(argc,argv);
    relative_opt.retrieveOption(argc,argv);
    kde_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    cout << predefinedString << endl;
    exit(0);
  }
  if(!doProcess){
    cout << endl;
    cout << "Usage: pkstatogr -i input [-n attribute]*" << endl;
    cout << endl;
    cout << "short option -h shows basic options only, use long option --help to show all options" << endl;
    exit(0);//help was invoked, stop processing
  }

  ImgReaderOgr imgReader;
  try{
    imgReader.open(input_opt[0]);
  }
  catch(string errorstring){
    cerr << errorstring << endl;
  }

  ImgReaderOgr inputReader(input_opt[0]);
  vector<double> theData;
  statfactory::StatFactory stat;
  //todo: implement ALL

  stat.setNoDataValues(nodata_opt);

  //support multiple layers
  int nlayerRead=inputReader.getDataSource()->GetLayerCount();
  if(verbose_opt[0])
    cout << "number of layers: " << nlayerRead << endl;
      
  for(int ilayer=0;ilayer<nlayerRead;++ilayer){
    OGRLayer *readLayer=inputReader.getLayer(ilayer);
    string currentLayername=readLayer->GetName();
    if(layer_opt.size())
      if(find(layer_opt.begin(),layer_opt.end(),currentLayername)==layer_opt.end())
	continue;
    if(verbose_opt[0])
      cout << "processing layer " << currentLayername << endl;
    if(layer_opt.size())
      cout << " --lname " << currentLayername;
      
    for(int ifield=0;ifield<fieldname_opt.size();++ifield){
      if(verbose_opt[0])
	cout << "field: " << ifield << endl;
      theData.clear();
      inputReader.readData(theData,OFTReal,fieldname_opt[ifield],ilayer,verbose_opt[0]);
      vector<double> binData;
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
	if(kde_opt[0]){
	  // if(kde_opt[0]>0)
	  //   sigma=kde_opt[0];
	  // else
	    sigma=1.06*sqrt(stat.var(theData))*pow(theData.size(),-0.2);
	}
	if(nbin<1)
	  nbin=(maxValue-minValue+1);
	try{
	  stat.distribution(theData,theData.begin(),theData.end(),binData,nbin,minValue,maxValue,sigma);
	}
	catch(string theError){
	  cerr << "Warning: all identical values in data" << endl;
	  exit(1);
	}
      }
      // int nbin=(nbin_opt[0]>1)? nbin_opt[0] : 2;
      cout << " --fname " << fieldname_opt[ifield];
      try{
	double theMean=0;
	double theVar=0;
	stat.meanVar(theData,theMean,theVar);
	if(mean_opt[0])
	  cout << " --mean " << theMean;
	if(stdev_opt[0])
	  cout << " --stdev " << sqrt(theVar);
	if(minmax_opt[0]||min_opt[0]||max_opt[0]){
	  if(minmax_opt[0])
	    cout << " --min " << minValue << " --max " << maxValue << " ";
	  else{
	    if(min_opt[0])
	      cout << " --min " << minValue << " ";
	    if(max_opt[0])
	      cout << " --max " << maxValue << " ";
	  }
	}
	if(median_opt[0])
	  cout << " -median " << stat.median(theData);
	if(size_opt[0])
	  cout << " -size " << theData.size();
	cout << endl;
	if(histogram_opt[0]){
	  for(int ibin=0;ibin<nbin;++ibin){
	    double binValue=0;
	    if(nbin==maxValue-minValue+1)
	      binValue=minValue+ibin;
	    else
	      binValue=minValue+static_cast<double>(maxValue-minValue)*(ibin+0.5)/nbin;
	    cout << binValue << " ";
	    if(relative_opt[0]||kde_opt[0])
	      cout << 100.0*static_cast<double>(binData[ibin])/theData.size() << endl;
	    else
	      cout << binData[ibin] << endl;
	  }
	}
      }
      catch(string theError){
	if(mean_opt[0])
	  cout << " --mean " << theData.back();
	if(stdev_opt[0])
	  cout << " --stdev " << "0";
	if(min_opt[0])
	  cout << " -min " << theData.back();
	if(max_opt[0])
	  cout << " -max " << theData.back();
	if(median_opt[0])
	  cout << " -median " << theData.back();
	if(size_opt[0])
	  cout << " -size " << theData.size();
	cout << endl;
	cerr << "Warning: all identical values in data" << endl;
      }
    }
  }
  imgReader.close();
}

