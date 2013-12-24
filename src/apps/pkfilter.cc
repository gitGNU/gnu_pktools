/**********************************************************************
pkfilter.cc: program to filter raster images: median, min/max, morphological, filtering
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
#include <assert.h>
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <sys/types.h>
#include <stdio.h>
#include "base/Optionpk.h"
#include "base/Vector2d.h"
#include "algorithms/Filter2d.h"
#include "algorithms/Filter.h"
#include "fileclasses/FileReaderAscii.h"
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"

using namespace std;
/*------------------
  Main procedure
  ----------------*/
int main(int argc,char **argv) {
  Optionpk<std::string> input_opt("i","input","input image file");
  Optionpk<std::string> output_opt("o", "output", "Output image file");
  Optionpk<bool> disc_opt("c", "circular", "circular disc kernel for dilation and erosion", false);
  // Optionpk<double> angle_opt("a", "angle", "angle used for directional filtering in dilation (North=0, East=90, South=180, West=270).");
  Optionpk<std::string> method_opt("f", "filter", "filter function (median, var, min, max, sum, mean, dilate, erode, close, open, homog (central pixel must be identical to all other pixels within window), heterog, sobelx (horizontal edge detection), sobely (vertical edge detection), sobelxy (diagonal edge detection NE-SW),sobelyx (diagonal edge detection NW-SE), smooth, density, majority voting (only for classes), smoothnodata (smooth nodata values only) values, threshold local filtering, ismin, ismax, heterogeneous (central pixel must be different than all other pixels within window), order (rank pixels in order), stdev, mrf, dwt, dwti, dwt_cut, scramble, shift, linearfeature)", "median");
  Optionpk<std::string> resample_opt("r", "resampling-method", "Resampling method for shifting operation (near: nearest neighbour, bilinear: bi-linear interpolation).", "near");
  Optionpk<double> dimX_opt("dx", "dx", "filter kernel size in x, better use odd value to avoid image shift", 3);
  Optionpk<double> dimY_opt("dy", "dy", "filter kernel size in y, better use odd value to avoid image shift", 3);
  Optionpk<int> dimZ_opt("dz", "dz", "filter kernel size in z (band or spectral dimension), must be odd (example: 3).. Set dz>0 if 1-D filter must be used in band domain");
  Optionpk<std::string> wavelet_type_opt("wt", "wavelet", "wavelet type: daubechies,daubechies_centered, haar, haar_centered, bspline, bspline_centered", "daubechies");
  Optionpk<int> family_opt("wf", "family", "wavelet family (vanishing moment, see also http://www.gnu.org/software/gsl/manual/html_node/DWT-Initialization.html)", 4);
  Optionpk<short> class_opt("class", "class", "class value(s) to use for density, erosion, dilation, openening and closing, thresholding");
  Optionpk<double> threshold_opt("t", "threshold", "threshold value(s) to use for threshold filter (one for each class), or threshold to cut for dwt_cut (use 0 to keep all), or sigma for shift", 0);
  Optionpk<short> nodata_opt("nodata", "nodata", "nodata value(s) for smoothnodata filter");
  Optionpk<std::string> tap_opt("tap", "tap", "text file containing taps used for spatial filtering (from ul to lr). Use dimX and dimY to specify tap dimensions in x and y. Leave empty for not using taps");
  Optionpk<double> tapz_opt("tapz", "tapz", "taps used for spectral filtering");
  Optionpk<double> fwhm_opt("fwhm", "fwhm", "list of full width half to apply spectral filtering (-fwhm band1 -fwhm band2 ...)");
  Optionpk<std::string> srf_opt("srf", "srf", "list of ASCII files containing spectral response functions (two columns: wavelength response)");
  Optionpk<double> wavelengthIn_opt("win", "wavelengthIn", "list of wavelengths in input spectrum (-win band1 -win band2 ...)");
  Optionpk<double> wavelengthOut_opt("wout", "wavelengthOut", "list of wavelengths in output spectrum (-wout band1 -wout band2 ...)");
  Optionpk<std::string> interpolationType_opt("interp", "interp", "type of interpolation for spectral filtering (see http://www.gnu.org/software/gsl/manual/html_node/Interpolation-Types.html)","akima",1);
  Optionpk<std::string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image","");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image");
  Optionpk<string>  colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid). Use none to ommit color table");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<short> down_opt("d", "down", "down sampling factor. Use value 1 for no downsampling). Use value n>1 for downsampling (aggregation)", 1);
  Optionpk<string> beta_opt("beta", "beta", "ASCII file with beta for each class transition in Markov Random Field");
  Optionpk<double> eps_opt("eps","eps", "error marging for linear feature",0);
  Optionpk<bool> l1_opt("l1","l1", "obtain longest object length for linear feature",false);
  Optionpk<bool> l2_opt("l2","l2", "obtain shortest object length for linear feature",false,2);
  Optionpk<bool> a1_opt("a1","a1", "obtain angle found for longest object length for linear feature",false);
  Optionpk<bool> a2_opt("a2","a2", "obtain angle found for shortest object length for linear feature",false);
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    disc_opt.retrieveOption(argc,argv);
    // angle_opt.retrieveOption(argc,argv);
    method_opt.retrieveOption(argc,argv);
    resample_opt.retrieveOption(argc,argv);
    dimX_opt.retrieveOption(argc,argv);
    dimY_opt.retrieveOption(argc,argv);
    dimZ_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    wavelet_type_opt.retrieveOption(argc,argv);
    family_opt.retrieveOption(argc,argv);
    class_opt.retrieveOption(argc,argv);
    threshold_opt.retrieveOption(argc,argv);
    nodata_opt.retrieveOption(argc,argv);
    tap_opt.retrieveOption(argc,argv);
    tapz_opt.retrieveOption(argc,argv);
    fwhm_opt.retrieveOption(argc,argv);
    srf_opt.retrieveOption(argc,argv);
    wavelengthIn_opt.retrieveOption(argc,argv);
    wavelengthOut_opt.retrieveOption(argc,argv);
    down_opt.retrieveOption(argc,argv);
    beta_opt.retrieveOption(argc,argv);
    eps_opt.retrieveOption(argc,argv);
    l1_opt.retrieveOption(argc,argv);
    l2_opt.retrieveOption(argc,argv);
    a1_opt.retrieveOption(argc,argv);
    a2_opt.retrieveOption(argc,argv);
    interpolationType_opt.retrieveOption(argc,argv);
    otype_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
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

  //not implemented yet, must debug first...
  vector<double> angle_opt;

  ImgReaderGdal input;
  ImgWriterGdal output;
  if(input_opt.empty()){
    cerr << "Error: no input file selected, use option -i" << endl;
    exit(1);
  }
  if(output_opt.empty()){
    cerr << "Error: no output file selected, use option -o" << endl;
    exit(1);
  }
  if(method_opt.empty()){
    cerr << "Error: no filter selected, use option -f" << endl;
    exit(1);
  }
  input.open(input_opt[0]);
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
    theType=input.getDataType();

  if(verbose_opt[0])
    std::cout << std::endl << "Output pixel type:  " << GDALGetDataTypeName(theType) << endl;

  string imageType=input.getImageType();
  if(oformat_opt.size())
    imageType=oformat_opt[0];

  if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
    string theInterleave="INTERLEAVE=";
    theInterleave+=input.getInterleave();
    option_opt.push_back(theInterleave);
  }
  try{
    if(filter2d::Filter2d::getFilterType(method_opt[0])==filter2d::mrf){
      assert(class_opt.size()>1);
      if(verbose_opt[0])
	std::cout << "opening output image " << output_opt[0] << std::endl;
      output.open(output_opt[0],(input.nrOfCol()+down_opt[0]-1)/down_opt[0],(input.nrOfRow()+down_opt[0]-1)/down_opt[0],class_opt.size(),theType,imageType,option_opt);
    }
    else if(filter2d::Filter2d::getFilterType(method_opt[0])==filter2d::linearfeature){
      if(verbose_opt[0])
	std::cout << "opening output image " << output_opt[0] << std::endl;
      int nband=0;
      if(l1_opt[0])
	++nband;
      if(a1_opt[0])
	++nband;
      if(l2_opt[0])
	++nband;
      if(a2_opt[0])
	++nband;
      output.open(output_opt[0],input.nrOfCol(),input.nrOfRow(),nband,theType,imageType,option_opt);
    }
    else if(fwhm_opt.size()||srf_opt.size()){
      //todo: support down and offset
      int nband=fwhm_opt.size()? fwhm_opt.size():srf_opt.size();
      output.open(output_opt[0],(input.nrOfCol()+down_opt[0]-1)/down_opt[0],(input.nrOfRow()+down_opt[0]-1)/down_opt[0],nband,theType,imageType,option_opt);
    }
    else
      output.open(output_opt[0],(input.nrOfCol()+down_opt[0]-1)/down_opt[0],(input.nrOfRow()+down_opt[0]-1)/down_opt[0],input.nrOfBand(),theType,imageType,option_opt);
  }
  catch(string errorstring){
    cout << errorstring << endl;
    exit(4);
  }
  if(input.isGeoRef()){
    output.setProjection(input.getProjection());
    double gt[6];
    input.getGeoTransform(gt);
    output.setGeoTransform(gt);
  }
  if(colorTable_opt.size()){
    if(colorTable_opt[0]!="none"){
      if(verbose_opt[0])
	cout << "set colortable " << colorTable_opt[0] << endl;
      assert(output.getDataType()==GDT_Byte);
      output.setColorTable(colorTable_opt[0]);
    }
  }
  else if(input.getColorTable()!=NULL)
    output.setColorTable(input.getColorTable());

  filter2d::Filter2d filter2d;
  filter::Filter filter1d;
  if(class_opt.size()){
    if(verbose_opt[0])
      std::cout<< "class values: ";
    for(int iclass=0;iclass<class_opt.size();++iclass){
      if(!dimZ_opt.size())
        filter2d.pushClass(class_opt[iclass]);
      else
        filter1d.pushClass(class_opt[iclass]);
      if(verbose_opt[0])
        std::cout<< class_opt[iclass] << " ";
    }
    if(verbose_opt[0])
      std::cout<< std::endl;
  }
  if(nodata_opt.size()){
    if(verbose_opt[0])
      std::cout<< "mask values: ";
    for(int imask=0;imask<nodata_opt.size();++imask){
      if(verbose_opt[0])
        std::cout<< nodata_opt[imask] << " ";
      filter2d.pushNoDataValue(nodata_opt[imask]);
    }
    if(verbose_opt[0])
      std::cout<< std::endl;
  }
  if(tap_opt.size()){
    ifstream tapfile(tap_opt[0].c_str());
    assert(tapfile);
    Vector2d<double> taps(dimY_opt[0],dimX_opt[0]);

    for(int j=0;j<dimY_opt[0];++j){
      for(int i=0;i<dimX_opt[0];++i){
        tapfile >> taps[j][i];
      }
    }
    if(verbose_opt[0]){
      std::cout << "taps: ";
      for(int j=0;j<dimY_opt[0];++j){
        for(int i=0;i<dimX_opt[0];++i){
          std::cout<< taps[j][i] << " ";
        }
        std::cout<< std::endl;
      }
    }
    filter2d.setTaps(taps);    
    filter2d.filter(input,output);
    tapfile.close();
  }
  else if(tapz_opt.size()){
    if(verbose_opt[0]){
      std::cout << "taps: ";
      for(int itap=0;itap<tapz_opt.size();++itap)
	std::cout<< tapz_opt[itap] << " ";
      std::cout<< std::endl;
    }
    filter1d.setTaps(tapz_opt);    
    filter1d.filter(input,output,down_opt[0]);
  }
  else if(fwhm_opt.size()){
    if(verbose_opt[0])
      std::cout << "spectral filtering to " << fwhm_opt.size() << " bands with provided fwhm " << std::endl;
    assert(wavelengthOut_opt.size()==fwhm_opt.size());
    assert(wavelengthIn_opt.size());

    Vector2d<double> lineInput(input.nrOfBand(),input.nrOfCol());
    Vector2d<double> lineOutput(wavelengthOut_opt.size(),input.nrOfCol());
    const char* pszMessage;
    void* pProgressArg=NULL;
    GDALProgressFunc pfnProgress=GDALTermProgress;
    double progress=0;
    pfnProgress(progress,pszMessage,pProgressArg);
    for(int y=0;y<input.nrOfRow();++y){
      if((y+1+down_opt[0]/2)%down_opt[0])
        continue;
      for(int iband=0;iband<input.nrOfBand();++iband)
        input.readData(lineInput[iband],GDT_Float64,y,iband);
      filter1d.applyFwhm<double>(wavelengthIn_opt,lineInput,wavelengthOut_opt,fwhm_opt, interpolationType_opt[0], lineOutput, down_opt[0], verbose_opt[0]);
      for(int iband=0;iband<output.nrOfBand();++iband){
        try{
          output.writeData(lineOutput[iband],GDT_Float64,y/down_opt[0],iband);
        }
        catch(string errorstring){
          cerr << errorstring << "in band " << iband << ", line " << y << endl;
        }
      }
      progress=(1.0+y)/output.nrOfRow();
      pfnProgress(progress,pszMessage,pProgressArg);
    }
    // filter1d.applyFwhm(wavelengthIn_opt,input,wavelengthOut_opt,fwhm_opt,interpolationType_opt[0],output,verbose_opt[0]);
  }
  else if(srf_opt.size()){
    if(verbose_opt[0])
      std::cout << "spectral filtering to " << srf_opt.size() << " bands with provided SRF " << std::endl;
    assert(wavelengthIn_opt.size());
    vector< Vector2d<double> > srf(srf_opt.size());//[0] srf_nr, [1]: wavelength, [2]: response
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
    }
    assert(output.nrOfBand()==srf.size());
    double centreWavelength=0;
    Vector2d<double> lineInput(input.nrOfBand(),input.nrOfCol());
    const char* pszMessage;
    void* pProgressArg=NULL;
    GDALProgressFunc pfnProgress=GDALTermProgress;
    double progress=0;
    pfnProgress(progress,pszMessage,pProgressArg);
    for(int y=0;y<input.nrOfRow();++y){
      if((y+1+down_opt[0]/2)%down_opt[0])
        continue;
      for(int iband=0;iband<input.nrOfBand();++iband)
        input.readData(lineInput[iband],GDT_Float64,y,iband);
      for(int isrf=0;isrf<srf.size();++isrf){
        vector<double> lineOutput(input.nrOfCol());
        double delta=1.0;
        bool normalize=true;
        centreWavelength=filter1d.applySrf<double>(wavelengthIn_opt,lineInput,srf[isrf], interpolationType_opt[0], lineOutput, delta, normalize, verbose_opt[0]);
        if(verbose_opt[0])
          std::cout << "centre wavelength srf " << isrf << ": " << centreWavelength << std::endl;
        try{
          output.writeData(lineOutput,GDT_Float64,y/down_opt[0],isrf);
        }
        catch(string errorstring){
          cerr << errorstring << "in srf " << srf_opt[isrf] << ", line " << y << endl;
        }

      }
      progress=(1.0+y)/output.nrOfRow();
      pfnProgress(progress,pszMessage,pProgressArg);
    }

    // filter1d.applySrf(wavelengthIn_opt,input,srf,interpolationType_opt[0],output,verbose_opt[0]);
  }
  else{
    switch(filter2d::Filter2d::getFilterType(method_opt[0])){
    case(filter2d::dilate):
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for morphological operator" << std::endl;
	exit(1);
      }
      if(dimZ_opt.size()){
        if(verbose_opt[0])
          std::cout<< "1-D filtering: dilate" << std::endl;
        filter1d.morphology(input,output,"dilate",dimZ_opt[0],1,0,verbose_opt[0]);
      }
      else{
	filter2d.morphology(input,output,"dilate",dimX_opt[0],dimY_opt[0],angle_opt,disc_opt[0]);
      }
      break;
    case(filter2d::erode):
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for morphological operator" << std::endl;
	exit(1);
      }
      if(dimZ_opt.size()>0){
        if(verbose_opt[0])
          std::cout<< "1-D filtering: dilate" << std::endl;
        filter1d.morphology(input,output,"erode",dimZ_opt[0]);
      }
      else{
	filter2d.morphology(input,output,"erode",dimX_opt[0],dimY_opt[0],angle_opt,disc_opt[0]);
      }
      break;
    case(filter2d::close):{//closing
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for morphological operator" << std::endl;
	exit(1);
      }
      ostringstream tmps;
      tmps << "/tmp/dilation_" << getpid() << ".tif";
      ImgWriterGdal tmpout;
      tmpout.open(tmps.str(),input);
      try{
        if(dimZ_opt.size()){
          filter1d.morphology(input,tmpout,"dilate",dimZ_opt[0]);
        }
        else{
	  filter2d.morphology(input,tmpout,"dilate",dimX_opt[0],dimY_opt[0],angle_opt,disc_opt[0]);
        }
      }
      catch(std::string errorString){
	std::cout<< errorString;
	exit(1);
      }
      tmpout.close();
      ImgReaderGdal tmpin;
      tmpin.open(tmps.str());
      if(dimZ_opt.size()){
        filter1d.morphology(tmpin,output,"erode",dimZ_opt[0]);
      }
      else{
          filter2d.morphology(tmpin,output,"erode",dimX_opt[0],dimY_opt[0],angle_opt,disc_opt[0]);
      }
      tmpin.close();
      if(remove(tmps.str().c_str( )) !=0){
        cerr << "could not remove " << tmps.str() << std::endl;
      }
      break;
    }
    case(filter2d::open):{//opening
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for morphological operator" << std::endl;
	exit(1);
      }
      ostringstream tmps;
      tmps << "/tmp/erosion_" << getpid() << ".tif";
      ImgWriterGdal tmpout;
      tmpout.open(tmps.str(),input);
      if(dimZ_opt.size()){
        filter1d.morphology(input,tmpout,"erode",dimZ_opt[0]);
      }
      else{
	filter2d.morphology(input,tmpout,"erode",dimX_opt[0],dimY_opt[0],angle_opt,disc_opt[0]);
      }
      tmpout.close();
      ImgReaderGdal tmpin;
      tmpin.open(tmps.str());
      if(dimZ_opt.size()){
        filter1d.morphology(tmpin,output,"dilate",dimZ_opt[0]);
      }
      else{
	filter2d.morphology(tmpin,output,"dilate",dimX_opt[0],dimY_opt[0],angle_opt,disc_opt[0]);
      }
      tmpin.close();
      if(remove(tmps.str().c_str( )) !=0){
        cerr << "could not remove " << tmps.str() << std::endl;
      }
      break;
    }
    case(filter2d::homog):{//spatially homogeneous
      filter2d.doit(input,output,"homog",dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      // filter2d.homogeneousSpatial(input,output,dimX_opt[0],disc_opt[0]);
      break;
    }
    case(filter2d::heterog):{//spatially heterogeneous
      filter2d.doit(input,output,"heterog",dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    }
    case(filter2d::shift):{//shift
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for shift operator" << std::endl;
	exit(1);
      }
      assert(input.nrOfBand());
      assert(input.nrOfCol());
      assert(input.nrOfRow());
      try{
        filter2d.shift(input,output,dimX_opt[0],dimY_opt[0],threshold_opt[0],filter2d::Filter2d::getResampleType(resample_opt[0]));
      }
      catch(string errorstring){
        cerr << errorstring << endl;
      }
      break;
    }
    case(filter2d::linearfeature):{
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for linear feature" << std::endl;
	exit(1);
      }
      assert(input.nrOfBand());
      assert(input.nrOfCol());
      assert(input.nrOfRow());
      float theAngle=361;
      if(angle_opt.size())
	theAngle=angle_opt[0];
      if(verbose_opt[0])
	std::cout << "using angle " << theAngle << std::endl;
      try{
	//using an angle step of 5 degrees and no maximum distance
        filter2d.linearFeature(input,output,theAngle,5,0,eps_opt[0],l1_opt[0],a1_opt[0],l2_opt[0],a2_opt[0],0,verbose_opt[0]);
      }
      catch(string errorstring){
        cerr << errorstring << endl;
      }
      break;
    }
    case(filter2d::mrf):{//Markov Random Field
      if(verbose_opt[0])
	std::cout << "Markov Random Field filtering" << std::endl;
      if(beta_opt.size()){
	//in file: classFrom classTo
	//in variable: beta[classTo][classFrom]
	FileReaderAscii betaReader(beta_opt[0]);
	Vector2d<double> beta(class_opt.size(),class_opt.size());
	vector<int> cols(class_opt.size());
	for(int iclass=0;iclass<class_opt.size();++iclass)
	  cols[iclass]=iclass;
	betaReader.readData(beta,cols);
	if(verbose_opt[0]){
	  std::cout << "using values for beta:" << std::endl;
	  for(int iclass1=0;iclass1<class_opt.size();++iclass1)
	    std::cout << "      " << iclass1 << " (" << class_opt[iclass1] << ")";
	  std::cout << std::endl;
	  for(int iclass1=0;iclass1<class_opt.size();++iclass1){
	    std::cout << iclass1 << " (" << class_opt[iclass1] << ")";
	    for(int iclass2=0;iclass2<class_opt.size();++iclass2)
	      std::cout << " " << beta[iclass2][iclass1] << " (" << class_opt[iclass2] << ")";
	    std::cout << std::endl;
	  }
	}
	filter2d.mrf(input, output, dimX_opt[0], dimY_opt[0], beta, true, down_opt[0], verbose_opt[0]);
      }
      else
	filter2d.mrf(input, output, dimX_opt[0], dimY_opt[0], 1, true, down_opt[0], verbose_opt[0]);
      break;
    }
    case(filter2d::sobelx):{//Sobel edge detection in X
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for sobel edge detection" << std::endl;
	exit(1);
      }
      Vector2d<double> theTaps(3,3);
      theTaps[0][0]=-1.0;
      theTaps[0][1]=0.0;
      theTaps[0][2]=1.0;
      theTaps[1][0]=-2.0;
      theTaps[1][1]=0.0;
      theTaps[1][2]=2.0;
      theTaps[2][0]=-1.0;
      theTaps[2][1]=0.0;
      theTaps[2][2]=1.0;
      filter2d.setTaps(theTaps);
      filter2d.filter(input,output,true);//absolute and normalize
      break;
    }
    case(filter2d::sobely):{//Sobel edge detection in Y
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for sobel edge detection" << std::endl;
	exit(1);
      }
      Vector2d<double> theTaps(3,3);
      theTaps[0][0]=1.0;
      theTaps[0][1]=2.0;
      theTaps[0][2]=1.0;
      theTaps[1][0]=0.0;
      theTaps[1][1]=0.0;
      theTaps[1][2]=0.0;
      theTaps[2][0]=-1.0;
      theTaps[2][1]=-2.0;
      theTaps[2][2]=-1.0;
      filter2d.setTaps(theTaps);
      filter2d.filter(input,output,true);//absolute and normalize
      break;
    }
    case(filter2d::sobelxy):{//Sobel edge detection in XY
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for sobel edge detection" << std::endl;
	exit(1);
      }
      Vector2d<double> theTaps(3,3);
      theTaps[0][0]=0.0;
      theTaps[0][1]=1.0;
      theTaps[0][2]=2.0;
      theTaps[1][0]=-1.0;
      theTaps[1][1]=0.0;
      theTaps[1][2]=1.0;
      theTaps[2][0]=-2.0;
      theTaps[2][1]=-1.0;
      theTaps[2][2]=0.0;
      filter2d.setTaps(theTaps);
      filter2d.filter(input,output,true);//absolute and normalize
      break;
    }
    case(filter2d::sobelyx):{//Sobel edge detection in XY
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for sobel edge detection" << std::endl;
	exit(1);
      }
      Vector2d<double> theTaps(3,3);
      theTaps[0][0]=2.0;
      theTaps[0][1]=1.0;
      theTaps[0][2]=0.0;
      theTaps[1][0]=1.0;
      theTaps[1][1]=0.0;
      theTaps[1][2]=-1.0;
      theTaps[2][0]=0.0;
      theTaps[2][1]=-1.0;
      theTaps[2][2]=-2.0;
      filter2d.setTaps(theTaps);
      filter2d.filter(input,output,true);//absolute and normalize
      break;
    }
    case(filter2d::smooth):{//Smoothing filter
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for this filter" << std::endl;
	exit(1);
      }
      if(dimZ_opt.size()){
        if(verbose_opt[0])
          std::cout<< "1-D filtering: smooth" << std::endl;
        filter1d.smooth(input,output,dimZ_opt[0]);
      }
      else{
	filter2d.smooth(input,output,dimX_opt[0],dimY_opt[0]);
      }
      break;
    }
    case(filter2d::smoothnodata):{//Smoothing filter
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for this filter" << std::endl;
	exit(1);
      }
      filter2d.smoothNoData(input,output,dimX_opt[0],dimY_opt[0]);
      break;
    }
    case(filter2d::dwt):
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for this filter" << std::endl;
	exit(1);
      }
      if(dimZ_opt.size())
	filter1d.dwtForward(input, output, wavelet_type_opt[0], family_opt[0]);
      else
	filter2d.dwtForward(input, output, wavelet_type_opt[0], family_opt[0]);
      break;
    case(filter2d::dwti):
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for this filter" << std::endl;
	exit(1);
      }
      if(dimZ_opt.size())
	filter1d.dwtInverse(input, output, wavelet_type_opt[0], family_opt[0]);
      else
	filter2d.dwtInverse(input, output, wavelet_type_opt[0], family_opt[0]);
      break;
    case(filter2d::dwt_cut):
      if(down_opt[0]!=1){
	std::cerr << "Error: down option not supported for this filter" << std::endl;
	exit(1);
      }
      if(dimZ_opt.size())
	filter1d.dwtCut(input, output, wavelet_type_opt[0], family_opt[0], threshold_opt[0]);
      else
	filter2d.dwtCut(input, output, wavelet_type_opt[0], family_opt[0], threshold_opt[0]);
      break;
    case(filter2d::threshold):
      filter2d.setThresholds(threshold_opt);
      filter2d.setClasses(class_opt);//deliberate fall through
    default:
      filter2d.doit(input,output,method_opt[0],dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);
      break;
    }
  }
  input.close();
  output.close();
  return 0;
}
