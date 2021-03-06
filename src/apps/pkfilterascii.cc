/**********************************************************************
pkfilterascii.cc: program to filter data in an ASCII file
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
#include <assert.h>
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <sys/types.h>
#include <stdio.h>
#include "base/Optionpk.h"
#include "base/Vector2d.h"
#include "algorithms/Filter.h"
#include "fileclasses/FileReaderAscii.h"

extern "C" {
#include <gsl/gsl_sort.h>
}

/******************************************************************************/
/*! \page pkfilterascii pkfilterascii
 program to filter data in an ASCII file
## SYNOPSIS

<code>
  Usage: pkfilterascii -i input.txt [-ic column]*
</code>

<code>
  
  Options: [-f filter] [-dz value] [-t]

  Advanced options: [-tapz value]* [-fwhm value]* [-srf filename]* [-win col] [-wout value]* [-interp type] [-wt type] [-wf family] [-cut threshold] 

</code>

\section pkfilterascii_description Description

The utility pkfilterascii filters the columns defined by the option -ic. A varietey of filters can be selected from with the option -f. The kernel size is defined with the option -dz. Alternatively, you can define your own filter tap values (use the option -tapz for each tap). In case of spectral filtering, define the full width half max values (-fwhm value) or spectral response functions in ASCII files (-srf filename).\section pkfilterascii_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |input ASCII file | 
 | o      | output               | std::string |       |Output ASCII file | 
 | ic     | inputCols            | int  |       |input columns (e.g., for three dimensional input data in first three columns use: -ic 0 -ic 1 -ic 2 | 
 | f      | filter               | std::string |       |filter function (to be implemented: dwt, dwti,dwt_cut) | 
 | dz     | dz                   | int  |       |filter kernel size in z (band or spectral dimension), must be odd (example: 3).. Set dz>0 if 1-D filter must be used in band domain | 
 | tapz   | tapz                 | double |       |taps used for spectral filtering | 
 | fwhm   | fwhm                 | double |       |list of full width half to apply spectral filtering (-fwhm band1 -fwhm band2 ...) | 
 | srf    | srf                  | std::string |       |list of ASCII files containing spectral response functions (two columns: wavelength response) | 
 | win    | wavelengthIn         | int  |       |column number of input ASCII file containing wavelengths | 
 | wout   | wavelengthOut        | double |       |list of wavelengths in output spectrum (-wout band1 -wout band2 ...) | 
 | interp | interp               | std::string | akima |type of interpolation for spectral filtering (see http://www.gnu.org/software/gsl/manual/html_node/Interpolation-Types.html) | 
 | t      | transpose            | bool | false |transpose output with samples in rows and wavelengths in cols | 
 | wt     | wavelet              | std::string | daubechies |wavelet type: daubechies,daubechies_centered, haar, haar_centered, bspline, bspline_centered | 
 | wf     | family               | int  | 4     |wavelet family (vanishing moment, see also http://www.gnu.org/software/gsl/manual/html_node/DWT-Initialization.html) | 
 | cut    | cut                  | double | 0     |threshold to cut dwt coefficients. Use 0 to keep all. | 

Usage: pkfilterascii -i input.txt [-ic column]*


**/

using namespace std;

/*------------------
  Main procedure
  ----------------*/
int main(int argc,char **argv) {
  Optionpk<std::string> input_opt("i","input","input ASCII file");
  Optionpk<std::string> output_opt("o", "output", "Output ASCII file");
  Optionpk<int> inputCols_opt("ic", "inputCols", "input columns (e.g., for three dimensional input data in first three columns use: -ic 0 -ic 1 -ic 2"); 
  Optionpk<std::string> method_opt("f", "filter", "filter function (to be implemented: dwt, dwti,dwt_cut)");
  Optionpk<std::string> wavelet_type_opt("wt", "wavelet", "wavelet type: daubechies,daubechies_centered, haar, haar_centered, bspline, bspline_centered", "daubechies");
  Optionpk<int> family_opt("wf", "family", "wavelet family (vanishing moment, see also http://www.gnu.org/software/gsl/manual/html_node/DWT-Initialization.html)", 4);
  Optionpk<double> threshold_opt("cut", "cut", "threshold to cut dwt coefficients. Use 0 to keep all.", 0);
  Optionpk<int> dimZ_opt("dz", "dz", "filter kernel size in z (band or spectral dimension), must be odd (example: 3).. Set dz>0 if 1-D filter must be used in band domain");
  Optionpk<double> tapZ_opt("tapz", "tapz", "taps used for spectral filtering");
  Optionpk<double> fwhm_opt("fwhm", "fwhm", "list of full width half to apply spectral filtering (-fwhm band1 -fwhm band2 ...)");
  Optionpk<std::string> srf_opt("srf", "srf", "list of ASCII files containing spectral response functions (two columns: wavelength response)");
  Optionpk<int> wavelengthIn_opt("win", "wavelengthIn", "column number of input ASCII file containing wavelengths");
  Optionpk<double> wavelengthOut_opt("wout", "wavelengthOut", "list of wavelengths in output spectrum (-wout band1 -wout band2 ...)");
  Optionpk<std::string> interpolationType_opt("interp", "interp", "type of interpolation for spectral filtering (see http://www.gnu.org/software/gsl/manual/html_node/Interpolation-Types.html)","akima");
  Optionpk<bool> transpose_opt("t", "transpose", "transpose output with samples in rows and wavelengths in cols", false);
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0,2);

  tapZ_opt.setHide(1);
  fwhm_opt.setHide(1);
  srf_opt.setHide(1);
  wavelengthIn_opt.setHide(1);
  wavelengthOut_opt.setHide(1);
  interpolationType_opt.setHide(1);
  transpose_opt.setHide(1);
  wavelet_type_opt.setHide(1);
  family_opt.setHide(1);
  threshold_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    inputCols_opt.retrieveOption(argc,argv);
    method_opt.retrieveOption(argc,argv);
    dimZ_opt.retrieveOption(argc,argv);
    tapZ_opt.retrieveOption(argc,argv);
    fwhm_opt.retrieveOption(argc,argv);
    srf_opt.retrieveOption(argc,argv);
    wavelengthIn_opt.retrieveOption(argc,argv);
    wavelengthOut_opt.retrieveOption(argc,argv);
    interpolationType_opt.retrieveOption(argc,argv);
    transpose_opt.retrieveOption(argc,argv);
    wavelet_type_opt.retrieveOption(argc,argv);
    family_opt.retrieveOption(argc,argv);
    threshold_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    cout << endl;
    cout << "Usage: pkfilterascii -i input.txt [-ic column]*" << endl;
    cout << endl;
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }
  if(doProcess&&input_opt.empty()){
    std::cerr << "Error: no input file provided (use option -i). Use --help for help information" << std::endl;
    exit(1);
  }

  Vector2d<double> inputData(inputCols_opt.size());
  Vector2d<double> filteredData(inputCols_opt.size());
  std::vector<double> wavelengthIn;
  std::vector<double> wavelengthOut;
  assert(input_opt.size());
  FileReaderAscii asciiReader(input_opt[0]);
  if(wavelengthIn_opt.size())
    asciiReader.readData(wavelengthIn,wavelengthIn_opt[0]);
  assert(inputCols_opt.size());
  asciiReader.readData(inputData,inputCols_opt);
  if(verbose_opt[0]){
    std::cout << "wavelengthIn.size(): " << wavelengthIn.size() << std::endl;
    std::cout << "inputData[0].size(): " << inputData[0].size() << std::endl;
  }
  if(wavelengthIn.size())
    assert(wavelengthIn.size()==inputData[0].size());
  asciiReader.close();
  filter::Filter filter1d;
  if(fwhm_opt.size()){
    filteredData.resize(inputCols_opt.size(),wavelengthOut_opt.size());
    assert(wavelengthIn_opt.size());
    if(verbose_opt[0])
      std::cout << "spectral filtering to " << fwhm_opt.size() << " bands with provided fwhm " << std::endl;
    assert(wavelengthOut_opt.size()==fwhm_opt.size());
    std::vector<double> fwhmData(wavelengthOut_opt.size());
    for(int icol=0;icol<inputCols_opt.size();++icol)
      filter1d.applyFwhm<double>(wavelengthIn,inputData[icol], wavelengthOut_opt,fwhm_opt, interpolationType_opt[0], filteredData[icol],verbose_opt[0]);
    if(verbose_opt[0])
      std::cout << "spectra filtered to " << wavelengthOut_opt.size() << " bands" << std::endl;
    wavelengthOut=wavelengthOut_opt;
  }
  else if(srf_opt.size()){
    wavelengthOut.resize(srf_opt.size());
    filteredData.resize(inputCols_opt.size(),srf_opt.size());
    Vector2d<double> srfData(srf_opt.size(),inputCols_opt.size());//transposed output
    if(verbose_opt[0])
      std::cout << "spectral filtering to " << srf_opt.size() << " bands with provided SRF " << std::endl;
    std::vector< Vector2d<double> > srf(srf_opt.size());//[0] srf_nr, [1]: wavelength, [2]: response
    ifstream srfFile;
    for(int isrf=0;isrf<srf_opt.size();++isrf){
      srf[isrf].resize(2);
      srfFile.open(srf_opt[isrf].c_str());
      double v;
      //add 0 to make sure srf is 0 at boundaries after interpolation step
      srf[isrf][0].push_back(0);
      srf[isrf][1].push_back(0);
      srf[isrf][0].push_back(1);
      srf[isrf][1].push_back(0);
      while(srfFile >> v){
        srf[isrf][0].push_back(v);
        srfFile >> v;
        srf[isrf][1].push_back(v);
      }
      srfFile.close();
      //add 0 to make sure srf[isrf] is 0 at boundaries after interpolation step
      srf[isrf][0].push_back(srf[isrf][0].back()+1);
      srf[isrf][1].push_back(0);    
      srf[isrf][0].push_back(srf[isrf][0].back()+1);
      srf[isrf][1].push_back(0);
      if(verbose_opt[0])
        cout << "srf file details: " << srf[isrf][0].size() << " wavelengths defined" << endl;
      if(verbose_opt[0]>1){
        for(int iw=0;iw<srf[isrf][0].size();++iw)
          std::cout << srf[isrf][0][iw] << " " << srf[isrf][1][iw] << std::endl;
      }
    }
    double centreWavelength=0;
    for(int icol=0;icol<inputCols_opt.size();++icol)
      filteredData[icol].resize(srf.size());

    for(int isrf=0;isrf<srf.size();++isrf){
      double delta=1.0;
      bool normalize=true;
      centreWavelength=filter1d.applySrf<double>(wavelengthIn,inputData,srf[isrf], interpolationType_opt[0], srfData[isrf], delta, normalize,1,true, verbose_opt[0]);
      if(verbose_opt[0])
        std::cout << "centre wavelength srf " << isrf << ": " << centreWavelength << std::endl;
      wavelengthOut[isrf]=static_cast<int>(centreWavelength+0.5);
    }
    srfData.transpose(filteredData);
    if(verbose_opt[0])
      std::cout << "spectra filtered to " << srf.size() << " bands" << std::endl;
  }
  else{//no filtering
    if(verbose_opt[0])
      std::cout << "no filtering selected" << std::endl;
    for(int icol=0;icol<inputCols_opt.size();++icol)
      filteredData[icol]=inputData[icol];
  }
  
  if(method_opt.size()){
    wavelengthOut=wavelengthIn;
    for(int icol=0;icol<inputCols_opt.size();++icol){
      switch(filter::Filter::getFilterType(method_opt[0])){
      case(filter::dwt):
        filter1d.dwtForward(filteredData[icol],wavelet_type_opt[0],family_opt[0]);
        break;
      case(filter::dwti):
        filter1d.dwtInverse(filteredData[icol],wavelet_type_opt[0],family_opt[0]);
        break;
      case(filter::dwt_cut):
	filter1d.dwtCut(filteredData[icol],wavelet_type_opt[0],family_opt[0],threshold_opt[0]);
        break;
      case(filter::smooth):
	if(tapZ_opt.size()){
	  filter1d.setTaps(tapZ_opt);
	  filter1d.filter(inputData[icol],filteredData[icol]);
	}
	else{
	  assert(dimZ_opt.size());
	  filter1d.smooth(inputData[icol],filteredData[icol],dimZ_opt[0]);
	}
	break;
      default:
	assert(tapZ_opt.size());
        filter1d.filter(inputData[icol],filteredData[icol]);
        // if(verbose_opt[0])
        //   std::cout << "method to be implemented" << std::endl;
	// exit(1);
        break;
      }
    }
  }
  ofstream outputStream;
  if(!output_opt.empty())
    outputStream.open(output_opt[0].c_str(),ios::out);

  if(verbose_opt[0])
    std::cout << "stream to output" << std::endl;
  if(transpose_opt[0]){
    for(int icol=0;icol<inputCols_opt.size();++icol){
      for(int iband=0;iband<filteredData[icol].size();++iband){
        if(!output_opt.empty()){
          outputStream << filteredData[icol][iband];
          if(iband<filteredData[icol].size()-1)
            outputStream << " ";
          else
            outputStream << std::endl;
        }
        else{
          std::cout << filteredData[icol][iband];
          if(iband<filteredData[icol].size()-1)
            std::cout << " ";
          else
            std::cout << std::endl;
        }
      }
    }    
  }
  else{
    // int nband=wavelengthOut.size()? wavelengthOut.size() : filteredData[0].size();
    int nband=0;
    if(method_opt.size()){
      switch(filter::Filter::getFilterType(method_opt[0])){
      case(filter::dwt):
        nband=filteredData[0].size();
        break;
      case(filter::dwti):
        nband=filteredData[0].size();
        break;
      case(filter::dwt_cut):
        nband=filteredData[0].size();
        break;
      default:
        nband=wavelengthOut.size();
        break;
      }
    }
    else
      nband=wavelengthOut.size();
    if(verbose_opt[0]){
      std::cout << "number of bands: " << nband << std::endl;
      std::cout << "wavelengthOut.size(): " << wavelengthOut.size() << std::endl;
    }
    for(int iband=0;iband<nband;++iband){
      if(!output_opt.empty()){
        if(wavelengthOut.size())
          outputStream << wavelengthOut[iband] << " ";
        else if(wavelengthIn_opt.size())
          outputStream << iband << " ";
      }
      else{
        if(wavelengthOut.size())
          std::cout << wavelengthOut[iband] << " ";
        else if(wavelengthIn_opt.size())
          std::cout << iband << " ";
      }
      for(int icol=0;icol<inputCols_opt.size();++icol){
        if(!output_opt.empty()){
          outputStream << filteredData[icol][iband];
          if(icol<inputCols_opt.size()-1)
            outputStream << " ";
          else
            outputStream << std::endl;
        }
        else{
          std::cout << filteredData[icol][iband];
          if(icol<inputCols_opt.size()-1)
            std::cout << " ";
          else
            std::cout << std::endl;
        }
      }
    }    
  }
  if(!output_opt.empty())
    outputStream.close();
  return 0;
}
