/**********************************************************************
pkstatascii.cc: program to calculate basic statistics from text file
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
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "base/Optionpk.h"
#include "fileclasses/FileReaderAscii.h"
#include "algorithms/StatFactory.h"
/******************************************************************************/
/*! \page pkstatascii pkstatascii
 program to calculate basic statistics from text file
## SYNOPSIS

<code>
  Usage: pkstatascii -i input [-c column]*
</code>

<code>
  
  Options: [-size] [-rnd number [-dist function] [-rnda value -rndb value]] [-mean] [-median] [-var] [-skew] [-stdev] [-sum] [-mm] [-min] [-max] [-hist [-nbin value] [-rel] [-kde]] [-hist2d [-nbin value] [-rel] [-kde]] [-cor] [-rmse] [-reg] [-regerr]

  Advanced options: [-srcmin value] [-srcmax value] [-fs separator] [-r startrow [-r endrow]] [-o [-t]] [--comment character]

</code>

\section pkstatascii_description Description

The utility pkstatascii calculates basic statistics of a data series in a text file.\section pkstatascii_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |name of the input text file | 
 | c      | column               | int  | 0     |column nr, starting from 0 | 
 | size   | size                 | bool | false |sample size | 
 | rnd    | rnd                  | unsigned int | 0     |generate random numbers | 
 | dist   | dist                 | std::string | gaussian |distribution for generating random numbers, see http://www.gn/software/gsl/manual/gsl-ref_toc.html#TOC320 (only uniform and Gaussian supported yet) | 
 | rnda   | rnda                 | double | 0     |first parameter for random distribution (mean value in case of Gaussian) | 
 | rndb   | rndb                 | double | 1     |second parameter for random distribution (standard deviation in case of Gaussian) | 
 | mean   | mean                 | bool | false |calculate median | 
 | median | median               | bool | false |calculate median | 
 | var    | var                  | bool | false |calculate variance | 
 | stdev  | stdev                | bool | false |calculate standard deviation | 
 | skew   | skewness             | bool | false |calculate skewness | 
 | kurt   | kurtosis             | bool | false |calculate kurtosis | 
 | sum    | sum                  | bool | false |calculate sum of column | 
 | mm     | minmax               | bool | false |calculate minimum and maximum value | 
 | min    | min                  | bool | false |calculate minimum value | 
 | max    | max                  | bool | false |calculate maximum value | 
 | hist   | hist                 | bool | false |calculate histogram | 
 | nbin   | nbin                 | short |       |number of bins to calculate histogram | 
 | rel    | relative             | bool | false |use percentiles for histogram to calculate histogram | 
 | kde    | kde                  | bool | false |Use Kernel density estimation when producing histogram. The standard deviation is estimated based on Silverman's rule of thumb | 
 | hist2d | hist2d               | bool | false |calculate 2-dimensional histogram based on two columns | 
 | cor    | correlation          | bool | false |calculate Pearson produc-moment correlation coefficient between two columns (defined by -c <col1> -c <col2> | 
 | rmse   | rmse                 | bool | false |calculate root mean square error between two columns (defined by -c <col1> -c <col2> | 
 | reg    | regression           | bool | false |calculate linear regression between two columns and get correlation coefficient (defined by -c <col1> -c <col2> | 
 | regerr | regerr               | bool | false |calculate linear regression between two columns and get root mean square error (defined by -c <col1> -c <col2> | 
 | src_min | src_min              | double |       |start reading source from this minimum value | 
 | src_max | src_max              | double |       |stop reading source from this maximum value | 
 | fs     | fs                   | char |       |field separator. | 
 | r      | range                | int  | 0     |rows to start/end reading. Use -r 1 -r 10 to read first 10 rows where first row is header. Use 0 to read all rows with no header. | 
 | o      | output               | bool | false |output the selected columns | 
 | t      | transpose            | bool | false |transpose input ascii vector (use in combination with --output) | 
 | comment | comment              | char | #     |comment character | 

Usage: pkstatascii -i input [-c column]*


**/

using namespace std;

int main(int argc, char *argv[])
{
  Optionpk<string> input_opt("i","input","name of the input text file");
  Optionpk<char> fs_opt("fs","fs","field separator.",' ');
  Optionpk<char> comment_opt("comment","comment","comment character",'#');
  Optionpk<bool> output_opt("o","output","output the selected columns",false);
  Optionpk<bool> transpose_opt("t","transpose","transpose input ascii vector (use in combination with --output)",false);
  Optionpk<int> col_opt("c", "column", "column nr, starting from 0", 0);
  Optionpk<int> range_opt("r", "range", "rows to start/end reading. Use -r 1 -r 10 to read first 10 rows where first row is header. Use 0 to read all rows with no header.", 0);
  Optionpk<bool> size_opt("size","size","sample size",false);
  Optionpk<unsigned int> rand_opt("rnd", "rnd", "generate random numbers", 0);
  Optionpk<std::string> randdist_opt("dist", "dist", "distribution for generating random numbers, see http://www.gn/software/gsl/manual/gsl-ref_toc.html#TOC320 (only uniform and Gaussian supported yet)", "gaussian");
  Optionpk<double> randa_opt("rnda", "rnda", "first parameter for random distribution (mean value in case of Gaussian)", 0);
  Optionpk<double> randb_opt("rndb", "rndb", "second parameter for random distribution (standard deviation in case of Gaussian)", 1);
  Optionpk<bool> mean_opt("mean","mean","calculate median",false);
  Optionpk<bool> median_opt("median","median","calculate median",false);
  Optionpk<bool> var_opt("var","var","calculate variance",false);
  Optionpk<bool> skewness_opt("skew","skewness","calculate skewness",false);
  Optionpk<bool> kurtosis_opt("kurt","kurtosis","calculate kurtosis",false);
  Optionpk<bool> stdev_opt("stdev","stdev","calculate standard deviation",false);
  Optionpk<bool> sum_opt("sum","sum","calculate sum of column",false);
  Optionpk<bool> minmax_opt("mm","minmax","calculate minimum and maximum value",false);
  Optionpk<bool> min_opt("min","min","calculate minimum value",false);
  Optionpk<bool> max_opt("max","max","calculate maximum value",false);
  Optionpk<double> src_min_opt("src_min","src_min","start reading source from this minimum value");
  Optionpk<double> src_max_opt("src_max","src_max","stop reading source from this maximum value");
  Optionpk<bool> histogram_opt("hist","hist","calculate histogram",false);
  Optionpk<bool> histogram2d_opt("hist2d","hist2d","calculate 2-dimensional histogram based on two columns",false);
  Optionpk<short> nbin_opt("nbin","nbin","number of bins to calculate histogram");
  Optionpk<bool> relative_opt("rel","relative","use percentiles for histogram to calculate histogram",false);
  Optionpk<bool> kde_opt("kde","kde","Use Kernel density estimation when producing histogram. The standard deviation is estimated based on Silverman's rule of thumb",false);
  Optionpk<bool> correlation_opt("cor","correlation","calculate Pearson produc-moment correlation coefficient between two columns (defined by -c <col1> -c <col2>",false);
  Optionpk<bool> rmse_opt("rmse","rmse","calculate root mean square error between two columns (defined by -c <col1> -c <col2>",false);
  Optionpk<bool> reg_opt("reg","regression","calculate linear regression between two columns and get correlation coefficient (defined by -c <col1> -c <col2>",false);
  Optionpk<bool> regerr_opt("regerr","regerr","calculate linear regression between two columns and get root mean square error (defined by -c <col1> -c <col2>",false);
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode when positive", 0,2);

  src_min_opt.setHide(1);
  src_max_opt.setHide(1);
  fs_opt.setHide(1);
  range_opt.setHide(1);
  output_opt.setHide(1);
  transpose_opt.setHide(1);
  comment_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    //mandatory options
    doProcess=input_opt.retrieveOption(argc,argv);
    col_opt.retrieveOption(argc,argv);
    //optional options
    size_opt.retrieveOption(argc,argv);
    rand_opt.retrieveOption(argc,argv);
    randdist_opt.retrieveOption(argc,argv);
    randa_opt.retrieveOption(argc,argv);
    randb_opt.retrieveOption(argc,argv);
    mean_opt.retrieveOption(argc,argv);
    median_opt.retrieveOption(argc,argv);
    var_opt.retrieveOption(argc,argv);
    stdev_opt.retrieveOption(argc,argv);
    skewness_opt.retrieveOption(argc,argv);
    kurtosis_opt.retrieveOption(argc,argv);
    sum_opt.retrieveOption(argc,argv);
    minmax_opt.retrieveOption(argc,argv);
    min_opt.retrieveOption(argc,argv);
    max_opt.retrieveOption(argc,argv);
    histogram_opt.retrieveOption(argc,argv);
    nbin_opt.retrieveOption(argc,argv);
    relative_opt.retrieveOption(argc,argv);
    kde_opt.retrieveOption(argc,argv);
    histogram2d_opt.retrieveOption(argc,argv);
    correlation_opt.retrieveOption(argc,argv);
    rmse_opt.retrieveOption(argc,argv);
    reg_opt.retrieveOption(argc,argv);
    regerr_opt.retrieveOption(argc,argv);
    //advanced options
    src_min_opt.retrieveOption(argc,argv);
    src_max_opt.retrieveOption(argc,argv);
    fs_opt.retrieveOption(argc,argv);
    range_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    transpose_opt.retrieveOption(argc,argv);
    comment_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    cout << endl;
    cout << "Usage: pkstatascii -i input [-c column]*" << endl;
    cout << endl;
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }

  if(src_min_opt.size()){
    while(src_min_opt.size()<col_opt.size())
      src_min_opt.push_back(src_min_opt[0]);
  }
  if(src_max_opt.size()){
    while(src_max_opt.size()<col_opt.size())
      src_max_opt.push_back(src_max_opt[0]);
  }
  statfactory::StatFactory stat;
  if(rand_opt[0]>0){
    gsl_rng* r=stat.getRandomGenerator(time(NULL));
    //todo: init random number generator using time...
    if(verbose_opt[0])
      std::cout << "generating " << rand_opt[0] << " random numbers: " << std::endl;
    for(unsigned int i=0;i<rand_opt[0];++i)
      std::cout << i << " " << stat.getRandomValue(r,randdist_opt[0],randa_opt[0],randb_opt[0]) << std::endl;
  }
  vector< vector<double> > dataVector(col_opt.size());
  vector< vector<double> > statVector(col_opt.size());

  if(!input_opt.size())
    exit(0);
  FileReaderAscii asciiReader(input_opt[0]);
  asciiReader.setFieldSeparator(fs_opt[0]);
  asciiReader.setComment(comment_opt[0]);
  asciiReader.setMinRow(range_opt[0]);
  if(range_opt.size()>1)
    asciiReader.setMaxRow(range_opt[1]);
  asciiReader.readData(dataVector,col_opt);
  assert(dataVector.size());
  double minValue=0;
  double maxValue=0;
  unsigned int nbin=0;
  if(nbin_opt.size())
    nbin=nbin_opt[0];
  if(histogram_opt[0]){
    stat.minmax(dataVector[0],dataVector[0].begin(),dataVector[0].end(),minValue,maxValue);
    if(src_min_opt.size())
      minValue=src_min_opt[0];
    if(src_max_opt.size())
      maxValue=src_max_opt[0];
    if(nbin<1){
      std::cerr << "Warning: number of bins not defined, calculating bins from min and max value" << std::endl;
      nbin=maxValue-minValue+1;
    }
  }
  double minX=0;
  double minY=0;
  double maxX=0;
  double maxY=0;
  if(histogram2d_opt[0]){
    assert(col_opt.size()==2);
    if(nbin<1){
      std::cerr << "Warning: number of bins not defined, calculating bins from min and max value" << std::endl;
      stat.minmax(dataVector[0],dataVector[0].begin(),dataVector[0].end(),minX,maxX);
      stat.minmax(dataVector[1],dataVector[1].begin(),dataVector[1].end(),minY,maxY);
      if(src_min_opt.size())
	minX=src_min_opt[0];
      if(src_min_opt.size()>1)
	minY=src_min_opt[1];
      if(src_max_opt.size())
	maxX=src_max_opt[0];
      if(src_max_opt.size()>1)
	maxY=src_max_opt[1];
      minValue=(minX<minY)? minX:minY;
      maxValue=(maxX>maxY)? maxX:maxY;
      if(verbose_opt[0])
        std::cout << "min and max values: " << minValue << ", " << maxValue << std::endl;
      nbin=maxValue-minValue+1;
    }
  }
  for(int icol=0;icol<col_opt.size();++icol){
    if(!dataVector[icol].size()){
      std::cerr << "Warning: dataVector[" << icol << "] is empty" << std::endl;
      continue;
    }
    if(size_opt[0])
      cout << "sample size column " << col_opt[icol] << ": " << dataVector[icol].size() << endl;
    if(mean_opt[0])
      cout << "mean value column " << col_opt[icol] << ": " << stat.mean(dataVector[icol]) << endl;
    if(var_opt[0])
      cout << "variance value column " << col_opt[icol] << ": " << stat.var(dataVector[icol]) << endl;
    if(stdev_opt[0])
      cout << "standard deviation column " << col_opt[icol] << ": " << sqrt(stat.var(dataVector[icol])) << endl;
    if(skewness_opt[0])
      cout << "skewness value column " << col_opt[icol] << ": " << stat.skewness(dataVector[icol]) << endl;
    if(kurtosis_opt[0])
      cout << "kurtosis value column " << col_opt[icol] << ": " << stat.kurtosis(dataVector[icol]) << endl;
    if(sum_opt[0]){
      cout << setprecision(2);
      cout << fixed << "sum column " << col_opt[icol] << ": " << (stat.sum(dataVector[icol])) << endl;
    }
    if(median_opt[0])
      cout << "median value column " << col_opt[icol] << ": " << stat.median(dataVector[icol]) << endl;
    if(minmax_opt[0]){
      cout << "min value column " << col_opt[icol] << ": " << stat.mymin(dataVector[icol]) << endl;
      cout << "max value column " << col_opt[icol] << ": " << stat.mymax(dataVector[icol]) << endl;
    }
    if(min_opt[0])
      cout << "min value column " << col_opt[icol] << ": " << stat.mymin(dataVector[icol]) << endl;
    if(max_opt[0])
      cout << "max value column " << col_opt[icol] << ": " << stat.mymax(dataVector[icol]) << endl;
    if(histogram_opt[0]){
      //todo: support kernel density function and estimate sigma as in practical estimate of the bandwith in http://en.wikipedia.org/wiki/Kernel_density_estimation
      double sigma=0;
      if(kde_opt[0]){//.size()){
        // if(kde_opt[0]>0)
        //   sigma=kde_opt[0];
        // else
          sigma=1.06*sqrt(stat.var(dataVector[icol]))*pow(dataVector[icol].size(),-0.2);
      }
      assert(nbin);
      if(verbose_opt[0]){
        if(sigma>0)
          std::cout << "calculating kernel density estimate with sigma " << sigma << " for col " << icol << std::endl;
        else
          std::cout << "calculating histogram for col " << icol << std::endl;
      }
      //test
      // cout << "debug0" << endl;
      // cout << "dataVector.size(): " << dataVector.size() << endl;
      // cout << "statVector.size(): " << statVector.size() << endl;

      // double theMinValue=0;
      // double theMaxValue=0;
      
      // stat.minmax(dataVector[icol],dataVector[icol].begin(),dataVector[icol].end(),theMinValue,theMaxValue);
      // if(minValue<maxValue&&minValue>theMinValue)
      // 	theMinValue=minValue;
      // if(minValue<maxValue&&maxValue<theMaxValue)
      // 	theMaxValue=maxValue;

      // //todo: check...
      // minValue=theMinValue;
      // maxValue=theMaxValue;

      // if(maxValue<=minValue){
      // 	std::ostringstream s;
      // 	s<<"Error: could not calculate distribution (min>=max)";
      // 	throw(s.str());
      // }
      // assert(nbin);
      // assert(dataVector[icol].size());
      // if(statVector[icol].size()!=nbin){
      // 	statVector[icol].resize(nbin);
      // 	for(int i=0;i<nbin;statVector[icol][i++]=0);
      // }
      // typename std::vector<double>::const_iterator it;
      // for(it=dataVector[icol].begin();it!=dataVector[icol].end();++it){
      // 	if(*it<minValue)
      // 	  continue;
      // 	if(*it>maxValue)
      // 	  continue;
      // 	if(stat.isNoData(*it))
      // 	  continue;
      // 	int theBin=0;
      // 	if(*it==maxValue)
      // 	  theBin=nbin-1;
      // 	else if(*it>minValue && *it<maxValue)
      // 	  theBin=static_cast<int>(static_cast<double>((nbin-1)*(*it)-minValue)/(maxValue-minValue));
      // 	assert(theBin<statVector[icol].size());
      // 	++statVector[icol][theBin];
      // 	// if(*it==maxValue)
      // 	//   ++statVector[icol][nbin-1];
      // 	// else if(*it>=minValue && *it<maxValue)
      // 	//   ++statVector[icol][static_cast<int>(static_cast<double>((*it)-minValue)/(maxValue-minValue)*nbin)];
      // }

      // exit(0);
      //end test
      
      stat.distribution(dataVector[icol],dataVector[icol].begin(),dataVector[icol].end(),statVector[icol],nbin,minValue,maxValue,sigma);
      if(verbose_opt[0])
        std::cout << "min and max values: " << minValue << ", " << maxValue << std::endl;
    }
  }
  if(correlation_opt[0]){
    assert(dataVector.size()==2);
    cout << "correlation between columns " << col_opt[0] << " and " << col_opt[1] << ": " << stat.correlation(dataVector[0],dataVector[1]) << endl;
  }
  if(rmse_opt[0]){
    assert(dataVector.size()==2);
    cout << "root mean square error between columns " << col_opt[0] << " and " << col_opt[1] << ": " << stat.rmse(dataVector[0],dataVector[1]) << endl;
  }
  if(reg_opt[0]){
    assert(dataVector.size()==2);
    double c0=0;
    double c1=0;
    double r2=stat.linear_regression(dataVector[0],dataVector[1],c0,c1);
    cout << "linear regression between columns: " << col_opt[0] << " and " << col_opt[1] << ": " << c0 << "+" << c1 << " * x " << " with R^2 (square correlation coefficient): " << r2 << endl;
  }
  if(regerr_opt[0]){
    assert(dataVector.size()==2);
    double c0=0;
    double c1=0;
    double err=stat.linear_regression_err(dataVector[0],dataVector[1],c0,c1);
    if(verbose_opt[0])
      cout << "linear regression between columns: " << col_opt[0] << " and " << col_opt[1] << ": " << c0 << "+" << c1 << " * x " << " with rmse: " << err << endl;
    else
      cout << c0 << " " << c1 << " " << err << endl;
  }
  if(histogram_opt[0]){
    for(int irow=0;irow<statVector.begin()->size();++irow){
      double binValue=0;
      if(nbin==maxValue-minValue+1)
	binValue=minValue+irow;
      else
	binValue=minValue+static_cast<double>(maxValue-minValue)*(irow+0.5)/nbin;
      std::cout << binValue << " ";
      // std::cout << minValue+static_cast<double>(maxValue-minValue)*(irow+0.5)/nbin << " ";
      for(int icol=0;icol<col_opt.size();++icol){
        if(relative_opt[0])
          std::cout << 100.0*static_cast<double>(statVector[icol][irow])/static_cast<double>(dataVector[icol].size());
        else
          std::cout << statVector[icol][irow];
        if(icol<col_opt.size()-1)
          cout << " ";
      }
      cout << endl;
    }
  }
  if(histogram2d_opt[0]){
    assert(nbin);
    assert(dataVector.size()==2);
    assert(dataVector[0].size()==dataVector[1].size());
    double sigma=0;
    //kernel density estimation as in http://en.wikipedia.org/wiki/Kernel_density_estimation
    if(kde_opt[0]){
      // if(kde_opt[0]>0)
      //   sigma=kde_opt[0];
      // else
        sigma=1.06*sqrt(sqrt(stat.var(dataVector[0]))*sqrt(stat.var(dataVector[1])))*pow(dataVector[0].size(),-0.2);
    }
    assert(nbin);
    if(verbose_opt[0]){
      if(sigma>0)
        std::cout << "calculating 2d kernel density estimate with sigma " << sigma << " for cols " << col_opt[0] << " and " << col_opt[1] << std::endl;
      else
        std::cout << "calculating 2d histogram for cols " << col_opt[0] << " and " << col_opt[1] << std::endl;
      std::cout << "nbin: " << nbin << std::endl;
    }
    std::vector< std::vector<double> > histVector;
    stat.distribution2d(dataVector[0],dataVector[1],histVector,nbin,minX,maxX,minY,maxY,sigma);
    for(int binX=0;binX<nbin;++binX){
      std::cout << std::endl;
      for(int binY=0;binY<nbin;++binY){
	double binValueX=0;
	if(nbin==maxX-minX+1)
	  binValueX=minX+binX;
	else
	  binValueX=minX+static_cast<double>(maxX-minX)*(binX+0.5)/nbin;
	double binValueY=0;
	if(nbin==maxY-minY+1)
	  binValueY=minY+binY;
	else
	  binValueY=minY+static_cast<double>(maxY-minY)*(binY+0.5)/nbin;
        double value=0;
        value=static_cast<double>(histVector[binX][binY])/dataVector[0].size();
	std::cout << binValueX << " " << binValueY << " " << value << std::endl;
	// std::cout << minX+static_cast<double>(maxX-minX)*(binX+0.5)/nbin << " " << minY+static_cast<double>(maxY-minY)*(binY+0.5)/nbin << " " << value << std::endl;
      }
    }
  }
  
  if(output_opt[0]){
    if(transpose_opt[0]){
      for(int icol=0;icol<col_opt.size();++icol){
        for(int irow=0;irow<dataVector.begin()->size();++irow){
          cout << dataVector[icol][irow];
          if(irow<dataVector.begin()->size()-1)
            cout << " ";
        }
        cout << endl;
      }
    }
    else{
      for(int irow=0;irow<dataVector.begin()->size();++irow){
        for(int icol=0;icol<col_opt.size();++icol){
          cout << dataVector[icol][irow];
          if(icol<col_opt.size()-1)
            cout << " ";
        }
        cout << endl;
      }
    }
  }
}
