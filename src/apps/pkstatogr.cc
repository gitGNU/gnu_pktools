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
  Optionpk<bool> version_opt("\0","version","version 20120625, Copyright (C) 2008-2012 Pieter Kempeneers.\n\
   This program comes with ABSOLUTELY NO WARRANTY; for details type use option -h.\n\
   This is free software, and you are welcome to redistribute it\n\
   under certain conditions; use option --license for details.",false);
  Optionpk<bool> license_opt("lic","license","show license information",false);
  Optionpk<bool> help_opt("h","help","shows this help info",false);
  Optionpk<bool> todo_opt("\0","todo","",false);
  Optionpk<string> input_opt("i", "input", "Input shape file", "");
  Optionpk<string> field_opt("f", "fields", "fields on which to calculate statistics", "");
  Optionpk<short> nbin_opt("n", "nbin", "number of bins", 0);
  Optionpk<bool> min_opt("m","min","calculate minimum value",false);
  Optionpk<bool> max_opt("M","max","calculate maximum value",false);
  Optionpk<bool> mean_opt("mean","mean","calculate mean value",false);
  Optionpk<bool> stdev_opt("stdev","stdev","calculate standard deviation",false);
  Optionpk<short> verbose_opt("v", "verbose", "verbose (Default: 0)", 0);

  version_opt.retrieveOption(argc,argv);
  license_opt.retrieveOption(argc,argv);
  help_opt.retrieveOption(argc,argv);
  todo_opt.retrieveOption(argc,argv);

  input_opt.retrieveOption(argc,argv);
  field_opt.retrieveOption(argc,argv);
  nbin_opt.retrieveOption(argc,argv);
  min_opt.retrieveOption(argc,argv);
  max_opt.retrieveOption(argc,argv);
  mean_opt.retrieveOption(argc,argv);
  stdev_opt.retrieveOption(argc,argv);
  verbose_opt.retrieveOption(argc,argv);

  if(version_opt[0]||todo_opt[0]){
    std::cout << version_opt.getHelp() << std::endl;
    std::cout << "todo: " << todo_opt.getHelp() << std::endl;
    exit(0);
  }
  if(license_opt[0]){
    cout << Optionpk<bool>::getGPLv3License() << endl;
    exit(0);
  }
  if(help_opt[0]){
    std::cout << "usage: pkstatogr -i inputimage -n nbin -f field [-f field]" << std::endl;
    exit(0);
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
  Histogram hist;
  //todo: implement ALL

  for(int ifield=0;ifield<field_opt.size();++ifield){
    if(verbose_opt[0])
      cout << "field: " << ifield << endl;
    theData.clear();
    inputReader.readData(theData,OFTReal,field_opt[ifield],0,verbose_opt[0]);
    vector<int> binData;
    double minimum=0;
    double maximum=0;
    int nbin=(nbin_opt[0]>1)? nbin_opt[0] : 2;
    hist.distribution(theData,theData.begin(),theData.end(),binData,nbin,minimum,maximum);
    double theMean=0;
    double theVar=0;
    hist.meanVar(theData,theMean,theVar);
    std::cout << " -f " << field_opt[ifield];
    if(mean_opt[0])
      std::cout << " --mean " << theMean;
    if(stdev_opt[0])
      std::cout << " --stdev " << sqrt(theVar);
    if(min_opt[0])
      cout << " -m " << minimum;
    if(max_opt[0])
      cout << " -M " << maximum;
    std::cout << std::endl;
    if(nbin_opt[0]>1){
      std::cout << std::endl;
      for(int bin=0;bin<nbin;++bin)
        cout << (maximum-minimum)*bin/(nbin-1)+minimum << " " << static_cast<double>(binData[bin])/theData.size() << endl;
    }
  }
  imgReader.close();
}

