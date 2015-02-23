/**********************************************************************
pkstat.cc: program to calculate basic statistics from raster dataset
Copyright (C) 2008-2015 Pieter Kempeneers

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
#include <math.h>
#include "base/Optionpk.h"
#include "algorithms/StatFactory.h"
#include "algorithms/ImgRegression.h"
using namespace std;

int main(int argc, char *argv[])
{
  Optionpk<string> input_opt("i","input","name of the input raster dataset");
  Optionpk<unsigned short> band_opt("b","band","band(s) on which to calculate statistics",0);
  Optionpk<bool>  filename_opt("f", "filename", "Shows image filename ", false);
  Optionpk<bool>  stat_opt("stats", "statistics", "Shows basic statistics (min,max, mean and stdDev of the raster datasets)", false);
  Optionpk<double>  ulx_opt("ulx", "ulx", "Upper left x value bounding box");
  Optionpk<double>  uly_opt("uly", "uly", "Upper left y value bounding box");
  Optionpk<double>  lrx_opt("lrx", "lrx", "Lower right x value bounding box");
  Optionpk<double>  lry_opt("lry", "lry", "Lower right y value bounding box");
  Optionpk<double> nodata_opt("nodata","nodata","Set nodata value(s)");
  Optionpk<short> down_opt("down", "down", "Down sampling factor (for raster sample datasets only). Can be used to create grid points", 1);
  Optionpk<unsigned int> random_opt("rnd", "rnd", "generate random numbers", 0);

  // Optionpk<bool> transpose_opt("t","transpose","transpose output",false);
  // Optionpk<std::string> randdist_opt("dist", "dist", "distribution for generating random numbers, see http://www.gn/software/gsl/manual/gsl-ref_toc.html#TOC320 (only uniform and Gaussian supported yet)", "gaussian");
  // Optionpk<double> randa_opt("rnda", "rnda", "first parameter for random distribution (mean value in case of Gaussian)", 0);
  // Optionpk<double> randb_opt("rndb", "rndb", "second parameter for random distribution (standard deviation in case of Gaussian)", 1);
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
  Optionpk<bool> histogram2d_opt("hist2d","hist2d","calculate 2-dimensional histogram based on two images",false);
  Optionpk<short> nbin_opt("nbin","nbin","number of bins to calculate histogram");
  Optionpk<bool> relative_opt("rel","relative","use percentiles for histogram to calculate histogram",false);
  Optionpk<bool> kde_opt("kde","kde","Use Kernel density estimation when producing histogram. The standard deviation is estimated based on Silverman's rule of thumb",false);
  Optionpk<bool> correlation_opt("cor","correlation","calculate Pearson produc-moment correlation coefficient between two raster datasets (defined by -c <col1> -c <col2>",false);
  Optionpk<bool> rmse_opt("rmse","rmse","calculate root mean square error between two raster datasets",false);
  Optionpk<bool> reg_opt("reg","regression","calculate linear regression between two raster datasets and get correlation coefficient",false);
  Optionpk<bool> regerr_opt("regerr","regerr","calculate linear regression between two raster datasets and get root mean square error",false);
  Optionpk<bool> preg_opt("preg","preg","calculate perpendicular regression between two raster datasets and get correlation coefficient",false);
  Optionpk<bool> pregerr_opt("pregerr","pregerr","calculate perpendicular regression between two raster datasets and get root mean square error",false);
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode when positive", 0,2);

  src_min_opt.setHide(1);
  src_max_opt.setHide(1);
  // range_opt.setHide(1);
  // transpose_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    //mandatory options
    doProcess=input_opt.retrieveOption(argc,argv);
    //optional options
    band_opt.retrieveOption(argc,argv);
    filename_opt.retrieveOption(argc,argv);
    stat_opt.retrieveOption(argc,argv);
    ulx_opt.retrieveOption(argc,argv);
    uly_opt.retrieveOption(argc,argv);
    lrx_opt.retrieveOption(argc,argv);
    lry_opt.retrieveOption(argc,argv);
    nodata_opt.retrieveOption(argc,argv);
    down_opt.retrieveOption(argc,argv);
    random_opt.retrieveOption(argc,argv);
    // randdist_opt.retrieveOption(argc,argv);
    // randa_opt.retrieveOption(argc,argv);
    // randb_opt.retrieveOption(argc,argv);
    mean_opt.retrieveOption(argc,argv);
    median_opt.retrieveOption(argc,argv);
    var_opt.retrieveOption(argc,argv);
    stdev_opt.retrieveOption(argc,argv);
    // skewness_opt.retrieveOption(argc,argv);
    // kurtosis_opt.retrieveOption(argc,argv);
    // sum_opt.retrieveOption(argc,argv);
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
    preg_opt.retrieveOption(argc,argv);
    pregerr_opt.retrieveOption(argc,argv);
    //advanced options
    src_min_opt.retrieveOption(argc,argv);
    src_max_opt.retrieveOption(argc,argv);
    // range_opt.retrieveOption(argc,argv);
    // transpose_opt.retrieveOption(argc,argv);
    // comment_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    cout << endl;
    cout << "Usage: pkstat -i input [-c column]*" << endl;
    cout << endl;
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }

  if(src_min_opt.size()){
    while(src_min_opt.size()<band_opt.size())
      src_min_opt.push_back(src_min_opt[0]);
  }
  if(src_max_opt.size()){
    while(src_max_opt.size()<band_opt.size())
      src_max_opt.push_back(src_max_opt[0]);
  }

  unsigned int nbin=0;
  double minX=0;
  double minY=0;
  double maxX=0;
  double maxY=0;
  double minValue=0;
  double maxValue=0;
  double meanValue=0;
  double stdDev=0;

  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  srand(time(NULL));

  statfactory::StatFactory stat;
  imgregression::ImgRegression imgreg;

  ImgReaderGdal imgReader;
  if(input_opt.empty()){
    std::cerr << "No image dataset provided (use option -i). Use --help for help information";
      exit(0);
  }
  for(int ifile=0;ifile<input_opt.size();++ifile){
    try{
      imgReader.open(input_opt[ifile]);
    }
    catch(std::string errorstring){
      std::cout << errorstring << std::endl;
      exit(0);
    }
    for(int inodata=0;inodata<nodata_opt.size();++inodata){
      if(!inodata)
        imgReader.GDALSetNoDataValue(nodata_opt[0],band_opt[0]);//only single no data can be set in GDALRasterBand (used for ComputeStatistics)
      imgReader.pushNoDataValue(nodata_opt[inodata]);
    }
    if(filename_opt[0])
      std::cout << " --input " << input_opt[ifile] << " ";
    int nband=band_opt.size();
    for(int iband=0;iband<nband;++iband){
      if(stat_opt[0]||mean_opt[0]||var_opt[0]||stdev_opt[0]){
	assert(band_opt[iband]<imgReader.nrOfBand());
	GDALProgressFunc pfnProgress;
	void* pProgressData;
	GDALRasterBand* rasterBand;
	
	rasterBand=imgReader.getRasterBand(band_opt[iband]);
	rasterBand->ComputeStatistics(0,&minValue,&maxValue,&meanValue,&stdDev,pfnProgress,pProgressData);
	if(mean_opt[0])
	  std::cout << "--mean " << meanValue << " ";
	if(stdev_opt[0])
	  std::cout << "--stdDev " << stdDev << " ";
	if(var_opt[0])
	  std::cout << "--var " << stdDev*stdDev << " ";
	if(stat_opt[0])
	  std::cout << "-min " << minValue << " -max " << maxValue << " --mean " << meanValue << " --stdDev " << stdDev << " ";
      }

      if(minmax_opt[0]||min_opt[0]||max_opt[0]){
	assert(band_opt[iband]<imgReader.nrOfBand());
	if((ulx_opt.size()||uly_opt.size()||lrx_opt.size()||lry_opt.size())&&(imgReader.covers(ulx_opt[0],uly_opt[0],lrx_opt[0],lry_opt[0]))){
	  double uli,ulj,lri,lrj;
	  imgReader.geo2image(ulx_opt[0],uly_opt[0],uli,ulj);
	  imgReader.geo2image(lrx_opt[0],lry_opt[0],lri,lrj);
	  imgReader.getMinMax(static_cast<int>(uli),static_cast<int>(lri),static_cast<int>(ulj),static_cast<int>(lrj),band_opt[iband],minValue,maxValue);
	}
	else
	  imgReader.getMinMax(minValue,maxValue,band_opt[iband],true);
	if(minmax_opt[0])
	  std::cout << "-min " << minValue << " -max " << maxValue << " ";
	else{
	  if(min_opt[0])
	    std::cout << "-min " << minValue << " ";
	  if(max_opt[0])
	    std::cout << "-max " << maxValue << " ";
	}
      }
    }
    if(histogram_opt[0]){//only for first selected band
      assert(band_opt[0]<imgReader.nrOfBand());
      nbin=(nbin_opt.size())? nbin_opt[0]:0;
      
      imgReader.getMinMax(minValue,maxValue,band_opt[0]);
      if(src_min_opt.size())
        minValue=src_min_opt[0];
      if(src_max_opt.size())
        maxValue=src_max_opt[0];
      if(minValue>=maxValue)
	imgReader.getMinMax(minValue,maxValue,band_opt[0]);

      if(verbose_opt[0])
	cout << "number of valid pixels in image: " << imgReader.getNvalid(band_opt[0]) << endl;
      std::vector<unsigned long int> output;
      unsigned long int nsample=imgReader.getHistogram(output,minValue,maxValue,nbin,band_opt[0],kde_opt[0]);
      std::cout.precision(10);
      for(int bin=0;bin<nbin;++bin){
	double binValue=0;
	if(nbin==maxValue-minValue+1)
	  binValue=minValue+bin;
	else
	  binValue=minValue+static_cast<double>(maxValue-minValue)*(bin+0.5)/nbin;
	std::cout << binValue << " ";
	if(relative_opt[0])
	  std::cout << 100.0*static_cast<double>(output[bin])/static_cast<double>(nsample) << std::endl;
	else
	  std::cout << static_cast<double>(output[bin]) << std::endl;
      }
    }
    if(histogram2d_opt[0]){
      if(band_opt.size()<2)
	continue;
      imgReader.getMinMax(minX,maxX,band_opt[0]);
      imgReader.getMinMax(minY,maxY,band_opt[1]);
      if(src_min_opt.size()){
	minX=src_min_opt[0];
	minY=src_min_opt[1];
      }
      if(src_max_opt.size()){
	maxX=src_max_opt[0];
	maxY=src_max_opt[1];
      }
      nbin=(nbin_opt.size())? nbin_opt[0]:0;
      if(nbin<=1){
	std::cerr << "Warning: number of bins not defined, calculating bins from min and max value" << std::endl;
	if(minX>=maxX)
	  imgReader.getMinMax(minX,maxX,band_opt[0]);
	if(minY>=maxY)
	  imgReader.getMinMax(minY,maxY,band_opt[1]);

	minValue=(minX<minY)? minX:minY;
	maxValue=(maxX>maxY)? maxX:maxY;
	if(verbose_opt[0])
	  std::cout << "min and max values: " << minValue << ", " << maxValue << std::endl;
	nbin=maxValue-minValue+1;
      }
      assert(nbin>1);
      double sigma=0;
      //kernel density estimation as in http://en.wikipedia.org/wiki/Kernel_density_estimation
      if(kde_opt[0]){
	assert(band_opt[0]<imgReader.nrOfBand());
	assert(band_opt[1]<imgReader.nrOfBand());
	GDALProgressFunc pfnProgress;
	void* pProgressData;
	GDALRasterBand* rasterBand;
	double stdDev1=0;
	double stdDev2=0;
	rasterBand=imgReader.getRasterBand(band_opt[0]);
	rasterBand->ComputeStatistics(0,&minValue,&maxValue,&meanValue,&stdDev1,pfnProgress,pProgressData);
	rasterBand=imgReader.getRasterBand(band_opt[1]);
	rasterBand->ComputeStatistics(0,&minValue,&maxValue,&meanValue,&stdDev2,pfnProgress,pProgressData);

	double estimatedSize=1.0*imgReader.getNvalid(band_opt[0])/down_opt[0]/down_opt[0];
	if(random_opt[0]>0)
	  estimatedSize*=random_opt[0]/100.0;
        sigma=1.06*sqrt(stdDev1*stdDev2)*pow(estimatedSize,-0.2);
      }
      assert(nbin);
      if(verbose_opt[0]){
	if(sigma>0)
	  std::cout << "calculating 2d kernel density estimate with sigma " << sigma << " for bands " << band_opt[0] << " and " << band_opt[1] << std::endl;
	else
	  std::cout << "calculating 2d histogram for bands " << band_opt[0] << " and " << band_opt[1] << std::endl;
	std::cout << "nbin: " << nbin << std::endl;
      }


      vector< vector<double> > output;

      if(maxX<=minX)
	imgReader.getMinMax(minX,maxX,band_opt[0]);
      if(maxY<=minY)
	imgReader.getMinMax(minY,maxY,band_opt[1]);

      if(maxX<=minX){
	std::ostringstream s;
	s<<"Error: could not calculate distribution (minX>=maxX)";
	throw(s.str());
      }
      if(maxY<=minY){
	std::ostringstream s;
	s<<"Error: could not calculate distribution (minY>=maxY)";
	throw(s.str());
      }
      output.resize(nbin);
      for(int i=0;i<nbin;++i){
	output[i].resize(nbin);
	for(int j=0;j<nbin;++j)
	  output[i][j]=0;
      }
      int binX=0;
      int binY=0;
      vector<double> inputX(imgReader.nrOfCol());
      vector<double> inputY(imgReader.nrOfCol());
      unsigned long int nvalid=0;
      for(int irow=0;irow<imgReader.nrOfRow();++irow){
        if(irow%down_opt[0])
          continue;
	imgReader.readData(inputX,GDT_Float64,irow,band_opt[0]);
	imgReader.readData(inputY,GDT_Float64,irow,band_opt[1]);
	for(int icol=0;icol<imgReader.nrOfCol();++icol){
          if(icol%down_opt[0])
            continue;
	  if(random_opt[0]>0){
	    double p=static_cast<double>(rand())/(RAND_MAX);
	    p*=100.0;
	    if(p>random_opt[0])
	      continue;//do not select for now, go to next column
	  }
	  if(imgReader.isNoData(inputX[icol]))
	    continue;
	  if(imgReader.isNoData(inputY[icol]))
	    continue;
	  ++nvalid;
	  if(inputX[icol]>=maxX)
	    binX=nbin-1;
	  else if(inputX[icol]<=minX)
	    binX=0;
	  else
	    binX=static_cast<int>(static_cast<double>(inputX[icol]-minX)/(maxX-minX)*nbin);
	  if(inputY[icol]>=maxY)
	    binY=nbin-1;
	  else if(inputY[icol]<=minX)
	    binY=0;
	  else
	    binY=static_cast<int>(static_cast<double>(inputY[icol]-minY)/(maxY-minY)*nbin);
	  assert(binX>=0);
	  assert(binX<output.size());
	  assert(binY>=0);
	  assert(binY<output[binX].size());
	  if(sigma>0){
	    //create kde for Gaussian basis function
	    //todo: speed up by calculating first and last bin with non-zero contriubtion...
	    for(int ibinX=0;ibinX<nbin;++ibinX){
	      double centerX=minX+static_cast<double>(maxX-minX)*ibinX/nbin;
	      double pdfX=gsl_ran_gaussian_pdf(inputX[icol]-centerX, sigma);
	      for(int ibinY=0;ibinY<nbin;++ibinY){
		//calculate  \integral_ibinX^(ibinX+1)
		double centerY=minY+static_cast<double>(maxY-minY)*ibinY/nbin;
		double pdfY=gsl_ran_gaussian_pdf(inputY[icol]-centerY, sigma);
		output[ibinX][binY]+=pdfX*pdfY;
	      }
	    }
	  }
	  else
	    ++output[binX][binY];
	}
      }
      if(verbose_opt[0])
	cout << "number of valid pixels: " << nvalid << endl;

      for(int binX=0;binX<nbin;++binX){
	cout << endl;
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

	  double value=static_cast<double>(output[binX][binY]);
	  
	  if(relative_opt[0])
	    value*=100.0/nvalid;

	  cout << binValueX << " " << binValueY << " " << value << std::endl;
	  // double value=static_cast<double>(output[binX][binY])/nvalid;
	  // cout << (maxX-minX)*bin/(nbin-1)+minX << " " << (maxY-minY)*bin/(nbin-1)+minY << " " << value << std::endl;
	}
      }
    }
    if(reg_opt[0]){
      if(band_opt.size()<2)
	continue;
      imgreg.setDown(down_opt[0]);
      imgreg.setThreshold(random_opt[0]);
      double c0=0;//offset
      double c1=1;//scale
      double r2=imgreg.getR2(imgReader,band_opt[0],band_opt[1],c0,c1,verbose_opt[0]);
      std::cout << "-c0 " << c0 << " -c1 " << c1 << " -r2 " << r2 << std::endl;
    }
    if(regerr_opt[0]){
      if(band_opt.size()<2)
	continue;
      imgreg.setDown(down_opt[0]);
      imgreg.setThreshold(random_opt[0]);
      double c0=0;//offset
      double c1=1;//scale
      double err=imgreg.getRMSE(imgReader,band_opt[0],band_opt[1],c0,c1,verbose_opt[0]);
      std::cout << "-c0 " << c0 << " -c1 " << c1 << " -rmse " << err << std::endl;
    }
    if(rmse_opt[0]){
      if(band_opt.size()<2)
	continue;
      imgreg.setDown(down_opt[0]);
      imgreg.setThreshold(random_opt[0]);
      double c0=0;//offset
      double c1=1;//scale
      double err=imgreg.getRMSE(imgReader,band_opt[0],band_opt[1],c0,c1,verbose_opt[0]);
      std::cout << " -rmse " << err << std::endl;
    }
    if(preg_opt[0]){
      if(band_opt.size()<2)
	continue;
      imgreg.setDown(down_opt[0]);
      imgreg.setThreshold(random_opt[0]);
      double c0=0;//offset
      double c1=1;//scale
      double r2=imgreg.pgetR2(imgReader,band_opt[0],band_opt[1],c0,c1,verbose_opt[0]);
      std::cout << "-c0 " << c0 << " -c1 " << c1 << " -r2 " << r2 << std::endl;
    }
    imgReader.close();
  }
  if(reg_opt[0]&&(input_opt.size()>1)){
    while(band_opt.size()<input_opt.size())
      band_opt.push_back(band_opt[0]);
    imgreg.setDown(down_opt[0]);
    imgreg.setThreshold(random_opt[0]);
    double c0=0;//offset
    double c1=1;//scale
    ImgReaderGdal imgReader1(input_opt[0]);
    ImgReaderGdal imgReader2(input_opt[1]);
    double r2=imgreg.getR2(imgReader1,imgReader2,c0,c1,band_opt[0],band_opt[1],verbose_opt[0]);
    std::cout << "-c0 " << c0 << " -c1 " << c1 << " -r2 " << r2 << std::endl;
    imgReader1.close();
    imgReader2.close();
  }
  if(preg_opt[0]&&(input_opt.size()>1)){
    while(band_opt.size()<input_opt.size())
      band_opt.push_back(band_opt[0]);
    imgreg.setDown(down_opt[0]);
    imgreg.setThreshold(random_opt[0]);
    double c0=0;//offset
    double c1=1;//scale
    ImgReaderGdal imgReader1(input_opt[0]);
    ImgReaderGdal imgReader2(input_opt[1]);
    double r2=imgreg.pgetR2(imgReader1,imgReader2,c0,c1,band_opt[0],band_opt[1],verbose_opt[0]);
    std::cout << "-c0 " << c0 << " -c1 " << c1 << " -r2 " << r2 << std::endl;
    imgReader1.close();
    imgReader2.close();
  }
  if(regerr_opt[0]&&(input_opt.size()>1)){
    while(band_opt.size()<input_opt.size())
      band_opt.push_back(band_opt[0]);
    imgreg.setDown(down_opt[0]);
    imgreg.setThreshold(random_opt[0]);
    double c0=0;//offset
    double c1=1;//scale
    ImgReaderGdal imgReader1(input_opt[0]);
    ImgReaderGdal imgReader2(input_opt[1]);
    double err=imgreg.getRMSE(imgReader1,imgReader2,c0,c1,band_opt[0],band_opt[1],verbose_opt[0]);
    std::cout << "-c0 " << c0 << " -c1 " << c1 << " -rmse " << err << std::endl;
    imgReader1.close();
    imgReader2.close();
  }
  if(rmse_opt[0]&&(input_opt.size()>1)){
    while(band_opt.size()<input_opt.size())
      band_opt.push_back(band_opt[0]);
    imgreg.setDown(down_opt[0]);
    imgreg.setThreshold(random_opt[0]);
    double c0=0;//offset
    double c1=1;//scale
    ImgReaderGdal imgReader1(input_opt[0]);
    ImgReaderGdal imgReader2(input_opt[1]);
    double err=imgreg.getRMSE(imgReader1,imgReader2,c0,c1,band_opt[0],band_opt[1],verbose_opt[0]);
    std::cout << "-rmse " << err << std::endl;
    imgReader1.close();
    imgReader2.close();
  }
  if(histogram2d_opt[0]&&(input_opt.size()>1)){
    imgReader.getMinMax(minX,maxX,band_opt[0]);
    imgReader.getMinMax(minY,maxY,band_opt[1]);
    if(src_min_opt.size()){
      while(src_min_opt.size()<input_opt.size())
	src_min_opt.push_back(src_min_opt[0]);
    }
    if(src_max_opt.size()){
      while(src_max_opt.size()<input_opt.size())
	src_max_opt.push_back(src_max_opt[0]);
    }
    if(src_min_opt.size()){
      minX=src_min_opt[0];
      minY=src_min_opt[1];
    }
    if(src_max_opt.size()){
      maxX=src_max_opt[0];
      maxY=src_max_opt[1];
    }
    nbin=(nbin_opt.size())? nbin_opt[0]:0;
    ImgReaderGdal imgReader1(input_opt[0]);
    ImgReaderGdal imgReader2(input_opt[1]);
    for(int inodata=0;inodata<nodata_opt.size();++inodata){
      if(!inodata){
        imgReader1.GDALSetNoDataValue(nodata_opt[0],band_opt[0]);//only single no data can be set in GDALRasterBand (used for ComputeStatistics)
        imgReader2.GDALSetNoDataValue(nodata_opt[0],band_opt[0]);//only single no data can be set in GDALRasterBand (used for ComputeStatistics)
      }
      imgReader1.pushNoDataValue(nodata_opt[inodata]);
      imgReader2.pushNoDataValue(nodata_opt[inodata]);
    }
    if(nbin<=1){
      std::cerr << "Warning: number of bins not defined, calculating bins from min and max value" << std::endl;
      // imgReader1.getMinMax(minX,maxX,band_opt[0]);
      // imgReader2.getMinMax(minY,maxY,band_opt[0]);
      if(minX>=maxX)
	imgReader1.getMinMax(minX,maxX,band_opt[0]);
      if(minY>=maxY)
	imgReader2.getMinMax(minY,maxY,band_opt[1]);
      
      minValue=(minX<minY)? minX:minY;
      maxValue=(maxX>maxY)? maxX:maxY;
      if(verbose_opt[0])
        std::cout << "min and max values: " << minValue << ", " << maxValue << std::endl;
      nbin=maxValue-minValue+1;
    }
    assert(nbin>1);
    double sigma=0;
    //kernel density estimation as in http://en.wikipedia.org/wiki/Kernel_density_estimation
    if(kde_opt[0]){
      GDALProgressFunc pfnProgress;
      void* pProgressData;
      GDALRasterBand* rasterBand;
      double stdDev1=0;
      double stdDev2=0;
      rasterBand=imgReader1.getRasterBand(band_opt[0]);
      rasterBand->ComputeStatistics(0,&minValue,&maxValue,&meanValue,&stdDev1,pfnProgress,pProgressData);
      rasterBand=imgReader2.getRasterBand(band_opt[0]);
      rasterBand->ComputeStatistics(0,&minValue,&maxValue,&meanValue,&stdDev2,pfnProgress,pProgressData);

      //todo: think of smarter way how to estimate size (nodata!)
      double estimatedSize=1.0*imgReader.getNvalid(band_opt[0])/down_opt[0]/down_opt[0];
      if(random_opt[0]>0)
	estimatedSize*=random_opt[0]/100.0;
      sigma=1.06*sqrt(stdDev1*stdDev2)*pow(estimatedSize,-0.2);
    }
    assert(nbin);
    if(verbose_opt[0]){
      if(sigma>0)
	std::cout << "calculating 2d kernel density estimate with sigma " << sigma << " for datasets " << input_opt[0] << " and " << input_opt[1] << std::endl;
      else
	std::cout << "calculating 2d histogram for datasets " << input_opt[0] << " and " << input_opt[1] << std::endl;
      std::cout << "nbin: " << nbin << std::endl;
    }

    vector< vector<double> > output;

    if(maxX<=minX)
      imgReader1.getMinMax(minX,maxX,band_opt[0]);
    if(maxY<=minY)
      imgReader2.getMinMax(minY,maxY,band_opt[0]);

    if(maxX<=minX){
      std::ostringstream s;
      s<<"Error: could not calculate distribution (minX>=maxX)";
      throw(s.str());
    }
    if(maxY<=minY){
      std::ostringstream s;
      s<<"Error: could not calculate distribution (minY>=maxY)";
      throw(s.str());
    }
    if(verbose_opt[0]){
      cout << "minX: " << minX << endl;
      cout << "maxX: " << maxX << endl;
      cout << "minY: " << minY << endl;
      cout << "maxY: " << maxY << endl;
    }
    output.resize(nbin);
    for(int i=0;i<nbin;++i){
      output[i].resize(nbin);
      for(int j=0;j<nbin;++j)
	output[i][j]=0;
    }
    int binX=0;
    int binY=0;
    vector<double> inputX(imgReader1.nrOfCol());
    vector<double> inputY(imgReader2.nrOfCol());
    unsigned long int nvalid=0;
    double geoX=0;
    double geoY=0;
    double icol1=0;
    double irow1=0;
    double icol2=0;
    double irow2=0;
    for(int irow=0;irow<imgReader1.nrOfRow();++irow){
      if(irow%down_opt[0])
	continue;
      irow1=irow;
      imgReader1.image2geo(icol1,irow1,geoX,geoY);
      imgReader2.geo2image(geoX,geoY,icol2,irow2);
      irow2=static_cast<int>(irow2);
      imgReader1.readData(inputX,GDT_Float64,irow1,band_opt[0]);
      imgReader2.readData(inputY,GDT_Float64,irow2,band_opt[0]);
      for(int icol=0;icol<imgReader.nrOfCol();++icol){
	if(icol%down_opt[0])
	  continue;
	icol1=icol;
	if(random_opt[0]>0){
	  double p=static_cast<double>(rand())/(RAND_MAX);
	  p*=100.0;
	  if(p>random_opt[0])
	    continue;//do not select for now, go to next column
	}
	if(imgReader1.isNoData(inputX[icol]))
	  continue;
	imgReader1.image2geo(icol1,irow1,geoX,geoY);
	imgReader2.geo2image(geoX,geoY,icol2,irow2);
	icol2=static_cast<int>(icol2);
	if(imgReader2.isNoData(inputY[icol2]))
	  continue;
	++nvalid;
	if(inputX[icol1]>=maxX)
	  binX=nbin-1;
	else if(inputX[icol]<=minX)
	  binX=0;
	else
	  binX=static_cast<int>(static_cast<double>(inputX[icol1]-minX)/(maxX-minX)*nbin);
	if(inputY[icol2]>=maxY)
	  binY=nbin-1;
	else if(inputY[icol2]<=minY)
	  binY=0;
	else
	  binY=static_cast<int>(static_cast<double>(inputY[icol2]-minY)/(maxY-minY)*nbin);
	assert(binX>=0);
	assert(binX<output.size());
	assert(binY>=0);
	assert(binY<output[binX].size());
	if(sigma>0){
	  //create kde for Gaussian basis function
	  //todo: speed up by calculating first and last bin with non-zero contriubtion...
	  for(int ibinX=0;ibinX<nbin;++ibinX){
	    double centerX=minX+static_cast<double>(maxX-minX)*ibinX/nbin;
	    double pdfX=gsl_ran_gaussian_pdf(inputX[icol1]-centerX, sigma);
	    for(int ibinY=0;ibinY<nbin;++ibinY){
	      //calculate  \integral_ibinX^(ibinX+1)
	      double centerY=minY+static_cast<double>(maxY-minY)*ibinY/nbin;
	      double pdfY=gsl_ran_gaussian_pdf(inputY[icol2]-centerY, sigma);
	      output[ibinX][binY]+=pdfX*pdfY;
	    }
	  }
	}
	else
	  ++output[binX][binY];
      }
    }
    if(verbose_opt[0])
      cout << "number of valid pixels: " << nvalid << endl;
    for(int binX=0;binX<nbin;++binX){
      cout << endl;
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
	double value=static_cast<double>(output[binX][binY]);
	  
	if(relative_opt[0])
	  value*=100.0/nvalid;

	cout << binValueX << " " << binValueY << " " << value << std::endl;
	// double value=static_cast<double>(output[binX][binY])/nvalid;
	// cout << (maxX-minX)*bin/(nbin-1)+minX << " " << (maxY-minY)*bin/(nbin-1)+minY << " " << value << std::endl;
      }
    }
    imgReader1.close();
    imgReader2.close();
  }

  if(!histogram_opt[0]||histogram2d_opt[0])
    std::cout << std::endl;
}
  
// int nband=(band_opt.size()) ? band_opt.size() : imgReader.nrOfBand();

// const char* pszMessage;
// void* pProgressArg=NULL;
// GDALProgressFunc pfnProgress=GDALTermProgress;
// double progress=0;
// srand(time(NULL));


// statfactory::StatFactory stat;
// imgregression::ImgRegression imgreg;

// pfnProgress(progress,pszMessage,pProgressArg);
// for(irow=0;irow<classReader.nrOfRow();++irow){
//   if(irow%down_opt[0])
//     continue;
//   // classReader.readData(classBuffer,GDT_Int32,irow);
//   classReader.readData(classBuffer,GDT_Float64,irow);
//   double x,y;//geo coordinates
//   double iimg,jimg;//image coordinates in img image
//   for(icol=0;icol<classReader.nrOfCol();++icol){
//     if(icol%down_opt[0])
  // 	continue;


  // if(rand_opt[0]>0){
  //   gsl_rng* r=stat.getRandomGenerator(time(NULL));
  //   //todo: init random number generator using time...
  //   if(verbose_opt[0])
  //     std::cout << "generating " << rand_opt[0] << " random numbers: " << std::endl;
  //   for(unsigned int i=0;i<rand_opt[0];++i)
  //     std::cout << i << " " << stat.getRandomValue(r,randdist_opt[0],randa_opt[0],randb_opt[0]) << std::endl;
  // }

  // imgreg.setDown(down_opt[0]);
  // imgreg.setThreshold(threshold_opt[0]);
  // double c0=0;//offset
  // double c1=1;//scale
  // double err=uncertNodata_opt[0];//start with high initial value in case we do not have first ob	err=imgreg.getRMSE(imgReaderModel1,imgReader,c0,c1,verbose_opt[0]);

  //   int nband=band_opt.size();
  //   if(band_opt[0]<0)
  //     nband=imgReader.nrOfBand();
  //   for(int iband=0;iband<nband;++iband){
  //     unsigned short band_opt[iband]=(band_opt[0]<0)? iband : band_opt[iband];

  //     if(minmax_opt[0]||min_opt[0]||max_opt[0]){
  // 	assert(band_opt[iband]<imgReader.nrOfBand());
  // 	if((ulx_opt.size()||uly_opt.size()||lrx_opt.size()||lry_opt.size())&&(imgReader.covers(ulx_opt[0],uly_opt[0],lrx_opt[0],lry_opt[0]))){
  // 	  double uli,ulj,lri,lrj;
  // 	  imgReader.geo2image(ulx_opt[0],uly_opt[0],uli,ulj);
  // 	  imgReader.geo2image(lrx_opt[0],lry_opt[0],lri,lrj);
  // 	  imgReader.getMinMax(static_cast<int>(uli),static_cast<int>(lri),static_cast<int>(ulj),static_cast<int>(lrj),band_opt[iband],minValue,maxValue);
  // 	}
  // 	else
  // 	  imgReader.getMinMax(minValue,maxValue,band_opt[iband],true);
  // 	if(minmax_opt[0])
  // 	  std::cout << "-min " << minValue << " -max " << maxValue << " ";
  // 	else{
  // 	  if(min_opt[0])
  // 	    std::cout << "-min " << minValue << " ";
  // 	  if(max_opt[0])
  // 	    std::cout << "-max " << maxValue << " ";
  // 	}
  //     }
  //   }
  //   if(relative_opt[0])
  //     hist_opt[0]=true;
  //   if(hist_opt[0]){
  //     assert(band_opt[0]<imgReader.nrOfBand());
  //     unsigned int nbin=(nbin_opt.size())? nbin_opt[0]:0;
  //     std::vector<unsigned long int> output;
  //     minValue=0;
  //     maxValue=0;
  //     //todo: optimize such that getMinMax is only called once...
  //     imgReader.getMinMax(minValue,maxValue,band_opt[0]);
      
  //     if(src_min_opt.size())
  //       minValue=src_min_opt[0];
  //     if(src_max_opt.size())
  //       maxValue=src_max_opt[0];
  //     unsigned long int nsample=imgReader.getHistogram(output,minValue,maxValue,nbin,band_opt[0]);
  //     std::cout.precision(10);
  //     for(int bin=0;bin<nbin;++bin){
  // 	double binValue=0;
  // 	if(nbin==maxValue-minValue+1)
  // 	  binValue=minValue+bin;
  // 	else
  // 	  binValue=minValue+static_cast<double>(maxValue-minValue)*(bin+0.5)/nbin;
  // 	std::cout << binValue << " ";
  // 	if(relative_opt[0])
  // 	  std::cout << 100.0*static_cast<double>(output[bin])/static_cast<double>(nsample) << std::endl;
  // 	else
  // 	  std::cout << static_cast<double>(output[bin]) << std::endl;
  //     }
  //   }
