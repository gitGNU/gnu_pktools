/**********************************************************************
pkstat.cc: program to calculate basic statistics from raster image
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
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "base/Optionpk.h"
#include "algorithms/Histogram.h"

using namespace std;

int main(int argc, char *argv[])
{
  Optionpk<string> input_opt("i","input","name of the input text file","");
  Optionpk<char> fs_opt("fs","fs","field separator.",' ');
  Optionpk<bool> output_opt("o","output","output the selected columns",false);
  Optionpk<short> col_opt("c", "column", "column nr, starting from 0", 0);
  Optionpk<int> range_opt("r", "range", "rows to start/end reading. Use -r 1 -r 10 to read first 10 rows where first row is header. Use 0 to read all rows with no header.", 0);
  Optionpk<bool> size_opt("size","size","sample size",false);
  Optionpk<bool> mean_opt("m","mean","calculate mean value",false);
  Optionpk<bool> median_opt("med","median","calculate median",false);
  Optionpk<bool> var_opt("var","var","calculate variance",false);
  Optionpk<bool> skewness_opt("skew","skewness","calculate skewness",false);
  Optionpk<bool> kurtosis_opt("kurt","kurtosis","calculate kurtosis",false);
  Optionpk<bool> stdev_opt("stdev","stdev","calculate standard deviation",false);
  Optionpk<bool> sum_opt("s","sum","calculate sum of column",false);
  Optionpk<bool> minmax_opt("mm","minmax","calculate minimum and maximum value",false);
  Optionpk<double> min_opt("min","min","calculate minimum value",0);
  Optionpk<double> max_opt("max","max","calculate maximum value",0);
  Optionpk<bool> histogram_opt("hist","hist","calculate histogram",false);
  Optionpk<short> nbin_opt("bin","bin","number of bins to calculate histogram",10);
  Optionpk<bool> relative_opt("rel","relative","use percentiles for histogram to calculate histogram",false);
  Optionpk<bool> correlation_opt("cor","correlation","calculate Pearson produc-moment correlation coefficient between two columns (defined by -c <col1> -c <col2>",false);
  Optionpk<bool> rmse_opt("e","rmse","calculate root mean square error between two columns (defined by -c <col1> -c <col2>",false);
  Optionpk<bool> reg_opt("reg","regression","calculate linear regression error between two columns (defined by -c <col1> -c <col2>",false);
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode when > 0", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    fs_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    col_opt.retrieveOption(argc,argv);
    range_opt.retrieveOption(argc,argv);
    size_opt.retrieveOption(argc,argv);
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
    correlation_opt.retrieveOption(argc,argv);
    rmse_opt.retrieveOption(argc,argv);
    reg_opt.retrieveOption(argc,argv);
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

  vector< vector<double> > dataVector(col_opt.size());
  vector< vector<int> > histVector(col_opt.size());
  ifstream dataFile;
  if(verbose_opt[0])
    cout << "opening file " << input_opt[0] << endl;
  dataFile.open(input_opt[0].c_str());

  int nrow=0;
  bool withinRange=true;
  
  if(fs_opt[0]>' '&&fs_opt[0]<='~'){//field separator is a regular character (minimum ASCII code is space, maximum ASCII code is tilde)
  // if(input_opt[0].find(".csv")!=string::npos){
    if(verbose_opt[0])
      cout << "reading csv file " << input_opt[0] << endl;
    string csvRecord;
    while(getline(dataFile,csvRecord)){//read a line
      withinRange=true;
      if(nrow<range_opt[0])
        withinRange=false;
      if(range_opt.size()>1)
        if(nrow>range_opt[1])
          withinRange=false;
      if(withinRange){
        istringstream csvstream(csvRecord);
        string item;
        int ncol=0;
        while(getline(csvstream,item,fs_opt[0])){//read a column
          if(verbose_opt[0])
            cout << item << " ";
          for(int icol=0;icol<col_opt.size();++icol){
            if(ncol==col_opt[icol]){
              double value=atof(item.c_str());
              if((value>=min_opt[0]&&value<=max_opt[0])||max_opt[0]<=min_opt[0])
                dataVector[icol].push_back(value);
            }
          }
          ++ncol;
        }
        if(verbose_opt[0])
          cout << endl;
        assert(ncol>=col_opt[0]);
      }
      ++nrow;
    }
    assert(dataVector.size());
  }
  else{//space or tab delimited fields
    string spaceRecord;
    while(!getline(dataFile, spaceRecord).eof()){
      withinRange=true;
      if(nrow<range_opt[0])
        withinRange=false;
      if(range_opt.size()>1)
        if(nrow>range_opt[1])
          withinRange=false;
      if(withinRange){
        if(verbose_opt[0]>1)
          cout << spaceRecord << endl;
        istringstream lineStream(spaceRecord);
        string item;
        int ncol=0;
        while(lineStream >> item){
          if(verbose_opt[0]>1)
            cout << item << " ";
          istringstream itemStream(item);
          double value;
          itemStream >> value;
          for(int icol=0;icol<col_opt.size();++icol){
            if(ncol==col_opt[icol]){
              if((value>=min_opt[0]&&value<=max_opt[0])||max_opt[0]<=min_opt[0])
                dataVector[icol].push_back(value);
            }
          }
          ++ncol;
        }
        if(verbose_opt[0]>1)
          cout << endl;
        if(verbose_opt[0])
          cout << "number of columns: " << ncol << endl;
        assert(ncol>=col_opt[0]);
      }
      ++nrow;
    }
  }
  assert(dataVector.size());
  dataFile.close();
  double minValue=min_opt[0];
  double maxValue=max_opt[0];
  Histogram hist;
  for(int icol=0;icol<col_opt.size();++icol){
    assert(dataVector[icol].size());
    if(size_opt[0])
      cout << "sample size column " << col_opt[icol] << ": " << dataVector[icol].size() << endl;
    if(mean_opt[0])
      cout << "mean value column " << col_opt[icol] << ": " << hist.mean(dataVector[icol]) << endl;
    if(var_opt[0])
      cout << "variance value column " << col_opt[icol] << ": " << hist.var(dataVector[icol]) << endl;
    if(stdev_opt[0])
      cout << "standard deviation column " << col_opt[icol] << ": " << sqrt(hist.var(dataVector[icol])) << endl;
    if(skewness_opt[0])
      cout << "skewness value column " << col_opt[icol] << ": " << hist.skewness(dataVector[icol]) << endl;
    if(kurtosis_opt[0])
      cout << "kurtosis value column " << col_opt[icol] << ": " << hist.kurtosis(dataVector[icol]) << endl;
    if(sum_opt[0]){
      cout << setprecision(2);
      cout << fixed << "sum column " << col_opt[icol] << ": " << (hist.sum(dataVector[icol])) << endl;
    }
    if(median_opt[0])
      cout << "median value column " << col_opt[icol] << ": " << hist.median(dataVector[icol]) << endl;
    if(minmax_opt[0]){
      cout << "min value  column " << col_opt[icol] << ": " << hist.min(dataVector[icol]) << endl;
      cout << "max value column " << col_opt[icol] << ": " << hist.max(dataVector[icol]) << endl;
    }
    if(histogram_opt[0]){
      if(verbose_opt[0])
        std::cout << "calculating histogram for col " << icol << std::endl;
      hist.distribution(dataVector[icol],dataVector[icol].begin(),dataVector[icol].end(),histVector[icol],nbin_opt[0],minValue,maxValue);
      if(verbose_opt[0])
        std::cout << "min and max values: " << minValue << ", " << maxValue << std::endl;
    }
  }
  if(correlation_opt[0]){
    assert(dataVector.size()==2);
    cout << "correlation between columns " << col_opt[0] << " and " << col_opt[1] << ": " << hist.correlation(dataVector[0],dataVector[1]) << endl;
  }
  if(rmse_opt[0]){
    assert(dataVector.size()==2);
    cout << "root mean square error between columns " << col_opt[0] << " and " << col_opt[1] << ": " << hist.rmse(dataVector[0],dataVector[1]) << endl;
  }
  if(reg_opt[0]){
    assert(dataVector.size()==2);
    double c0=0;
    double c1=0;
    double r2=hist.linear_regression(dataVector[0],dataVector[1],c0,c1);
    cout << "linear regression between columns: " << col_opt[0] << " and " << col_opt[1] << ": " << c0 << "+" << c1 << " * x " << " with R^2 (square correlation coefficient): " << r2 << endl;
  }
  if(histogram_opt[0]){
    for(int irow=0;irow<histVector.begin()->size();++irow){
      std::cout << (maxValue-minValue)*irow/(nbin_opt[0]-1)+minValue << " ";
      for(int icol=0;icol<col_opt.size();++icol){
        if(relative_opt[0])
          std::cout << 100.0*static_cast<double>(histVector[icol][irow])/static_cast<double>(dataVector[icol].size());
        else
          std::cout << histVector[icol][irow];
        if(icol<col_opt.size()-1)
          cout << " ";
      }
      cout << endl;
    }
  }
  if(output_opt[0]){
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
