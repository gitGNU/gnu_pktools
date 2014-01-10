/**********************************************************************
pkenhance.cc: program to enhance raster images: histogram matching
Copyright (C) 2008-2013 Pieter Kempeneers

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
#include "base/Optionpk.h"
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"

using namespace std;

int main(int argc,char **argv) {
  Optionpk<std::string> input_opt("i","input","input image file");
  Optionpk<string> reference_opt("r", "reference", "Reference image file");
  Optionpk<std::string> output_opt("o", "output", "Output image file");
  Optionpk<double> minRef_opt("minref", "minref", "Sets minimum for histogram of reference image");
  Optionpk<double> maxRef_opt("maxref", "maxref", "Sets maximum for histogram of reference image");
  Optionpk<double> minInput_opt("mininput", "mininput", "Sets minimum for histogram of input image");
  Optionpk<double> maxInput_opt("maxinput", "maxinput", "Sets maximum for histogram of input image");
  Optionpk<double> nodata_opt("nodata", "nodata", "Sets no data value(s) for calculations (nodata values in input image)");
  Optionpk<std::string> method_opt("m", "method", "enhancement method (histmatch)", "histmatch");
  // Optionpk<std::string> wavelet_type_opt("wt", "wavelet", "wavelet type: daubechies,daubechies_centered, haar, haar_centered, bspline, bspline_centered", "daubechies");
  // Optionpk<int> family_opt("wf", "family", "wavelet family (vanishing moment, see also http://www.gnu.org/software/gsl/manual/html_node/DWT-Initialization.html)", 4);
  // Optionpk<double> quantize_opt("q", "quantize", "Quantize threshold",0);
  Optionpk<short>  nbin_opt("nbin", "nbin", "Number of bins used in histogram. Use 0 for all input values as integers",0);
  Optionpk<string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image","");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    reference_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    minRef_opt.retrieveOption(argc,argv);
    maxRef_opt.retrieveOption(argc,argv);
    minInput_opt.retrieveOption(argc,argv);
    maxInput_opt.retrieveOption(argc,argv);
    nodata_opt.retrieveOption(argc,argv);
    method_opt.retrieveOption(argc,argv);
    nbin_opt.retrieveOption(argc,argv);
    otype_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
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

  ImgReaderGdal inputImg;
  ImgReaderGdal referenceImg;
  ImgWriterGdal outputImg;
  assert(input_opt.size());
  inputImg.open(input_opt[0]);
  for(int inodata=0;inodata<nodata_opt.size();++inodata){
    if(!inodata)
      inputImg.GDALSetNoDataValue(nodata_opt[0],0);//only single no data can be set in GDALRasterBand (used for ComputeStatistics)
    inputImg.pushNoDataValue(nodata_opt[inodata]);
  }

  int nband=inputImg.nrOfBand();
  GDALDataType theType=GDT_Unknown;
  if(verbose_opt[0])
    cout << "possible output data types: ";
  for(int iType = 0; iType < GDT_TypeCount; ++iType){
    if(verbose_opt[0])
      cout << " " << GDALGetDataTypeName((GDALDataType)iType);
    if( GDALGetDataTypeName((GDALDataType)iType) != NULL
        && EQUAL(GDALGetDataTypeName((GDALDataType)iType),
                 otype_opt[0].c_str()))
      theType=(GDALDataType) iType;
  }
  if(theType==GDT_Unknown)
    theType=inputImg.getDataType();

  if(verbose_opt[0])
    std::cout << std::endl << "Output pixel type:  " << GDALGetDataTypeName(theType) << endl;

  string imageType=inputImg.getImageType();
  if(oformat_opt.size())
    imageType=oformat_opt[0];

  if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
    string theInterleave="INTERLEAVE=";
    theInterleave+=inputImg.getInterleave();
    option_opt.push_back(theInterleave);
  }
  try{
    assert(output_opt.size());
    if(verbose_opt[0])
      std::cout << "opening output image " << output_opt[0] << std::endl;
    outputImg.open(output_opt[0],inputImg.nrOfCol(),inputImg.nrOfRow(),inputImg.nrOfBand(),theType,imageType,option_opt);
  }
  catch(string errorstring){
    cout << errorstring << endl;
    exit(1);
  }

  if(method_opt[0]=="histmatch"){
    assert(reference_opt.size());
    referenceImg.open(reference_opt[0]);
    assert(nband==referenceImg.nrOfBand());
    for(int inodata=0;inodata<nodata_opt.size();++inodata){
      if(!inodata)
        referenceImg.GDALSetNoDataValue(nodata_opt[0],0);//only single no data can be set in GDALRasterBand (used for ComputeStatistics)
      referenceImg.pushNoDataValue(nodata_opt[inodata]);
    }
    const char* pszMessage;
    void* pProgressArg=NULL;
    GDALProgressFunc pfnProgress=GDALTermProgress;
    double progress=0;
    pfnProgress(progress,pszMessage,pProgressArg);
    for(int iband=0;iband<nband;++iband){
      //calculate histograms
      unsigned int nbinRef=nbin_opt[0];
      unsigned int nbinInput=nbin_opt[0];
      std::vector<unsigned long int> histRef(nbinRef);
      std::vector<unsigned long int> histInput(nbinInput);
      double minValueRef=0;
      double maxValueRef=0;
      double minValueInput=0;
      double maxValueInput=0;
      if(minRef_opt.size())
        minValueRef=minRef_opt[0];
      if(maxRef_opt.size())
        maxValueRef=maxRef_opt[0];
      if(minInput_opt.size())
        minValueInput=minInput_opt[0];
      if(maxInput_opt.size())
        maxValueInput=maxInput_opt[0];
      unsigned long int nsampleRef=referenceImg.getHistogram(histRef,minValueRef,maxValueRef,nbinRef,iband);
      unsigned long int nsampleInput=inputImg.getHistogram(histInput,minValueInput,maxValueInput,nbinInput,iband);
      //create cumulative historgrams
      for(unsigned int bin=0;bin<nbinRef;++bin){
        histRef[bin]+=100.0*static_cast<double>(histRef[bin])/static_cast<double>(nsampleRef);
      }
      for(unsigned int bin=0;bin<nbinInput;++bin)
        histInput[bin]+=100.0*static_cast<double>(histInput[bin])/static_cast<double>(nsampleInput);
      //match histograms
      vector<double> lineBuffer(inputImg.nrOfCol());
      for(int irow=0;irow<inputImg.nrOfRow();++irow){
        inputImg.readData(lineBuffer,GDT_Float64, irow, iband);
        for(int icol=0;icol<inputImg.nrOfCol();++icol){
          //find bin in input image histogram
          int inputBin=static_cast<int>(static_cast<double>(lineBuffer[icol]-minValueInput)/(maxValueInput-minValueInput)*(histInput.size()-1));
          //find corresponding bin in reference image histogram
          //todo: optimize with lower_bound?
          // std::vector<unsigned long int>::const_iterator hit=histRef.begin();
          int ibin=0;
          for(ibin=0;ibin<histRef.size();++ibin){
            if(histRef[ibin]>histInput[inputBin])
              break;
          }
          if(ibin)
            --ibin;
          lineBuffer[icol]=(maxValueRef-minValueRef)/(histRef.size()-1)*(ibin)+minValueRef;
          // std::vector<unsigned long int>::const_iterator it=std::lower_bound(histRef.begin(),histRef.end(),culInput);
        }
        outputImg.writeData(lineBuffer,GDT_Float64,irow,iband);
        progress=(1.0+irow);
        progress+=(outputImg.nrOfRow()*iband);
        progress/=outputImg.nrOfBand()*outputImg.nrOfRow();
        assert(progress>=0);
        assert(progress<=1);
        pfnProgress(progress,pszMessage,pProgressArg);
      }
    }
  }

  inputImg.close();
  if(method_opt[0]=="histmatch")
    referenceImg.close();
  outputImg.close();
  return 0;
}
