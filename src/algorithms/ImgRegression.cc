/**********************************************************************
ImgRegression.cc: class to calculate regression between two raster datasets
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
#include "ImgRegression.h"
#include <iostream>

using namespace imgregression;

ImgRegression::ImgRegression(void)
: m_threshold(0), m_down(1)
{}

ImgRegression::~ImgRegression(void)
{}

double ImgRegression::getRMSE(ImgRasterGdal& imgReader1, ImgRasterGdal& imgReader2, double& c0, double& c1, unsigned int band1, unsigned int band2, short verbose) const{
  c0=0;
  c1=1;
  unsigned int icol1=0,irow1=0;
  std::vector<double> rowBuffer1(imgReader1.nrOfCol());
  std::vector<double> rowBuffer2(imgReader2.nrOfCol());
  std::vector<double> buffer1;
  std::vector<double> buffer2;

  srand(time(NULL));
  for(irow1=0;irow1<imgReader1.nrOfRow();++irow1){
    if(irow1%m_down)
      continue;
    icol1=0;
    double icol2=0,irow2=0;
    double geox=0,geoy=0;
    imgReader1.readData(rowBuffer1,GDT_Float64,irow1,band1);
    imgReader1.image2geo(icol1,irow1,geox,geoy);
    imgReader2.geo2image(geox,geoy,icol2,irow2);
    icol2=static_cast<int>(icol2);
    irow2=static_cast<int>(irow2);
    if(irow2<0||irow2>=imgReader2.nrOfRow())
      continue;
    imgReader2.readData(rowBuffer2,GDT_Float64,irow2,band2);
    for(icol1=0;icol1<imgReader1.nrOfCol();++icol1){
      if(icol1%m_down)
	continue;
      if(m_threshold>0){//percentual value
	double p=static_cast<double>(rand())/(RAND_MAX);
	p*=100.0;
	if(p>m_threshold)
	  continue;//do not select for now, go to next column
      }
      imgReader1.image2geo(icol1,irow1,geox,geoy);
      imgReader2.geo2image(geox,geoy,icol2,irow2);
      if(icol2<0||icol2>=imgReader2.nrOfCol())
	continue;
      icol2=static_cast<int>(icol2);
      irow2=static_cast<int>(irow2);
      //check for nodata
      double value1=rowBuffer1[icol1];
      double value2=rowBuffer2[icol2];
      if(imgReader1.isNoData(value1)||imgReader2.isNoData(value2))
	continue;

      buffer1.push_back(value1);
      buffer2.push_back(value2);
      if(verbose>1)
	std::cout << geox << " " << geoy << " " << icol1 << " " << irow1 << " " << icol2 << " " << irow2 << " " << buffer1.back() << " " << buffer2.back() << std::endl;
    }
  }
  double err=0;
  if(buffer1.size()&&buffer2.size()){
    statfactory::StatFactory stat;
    err=stat.linear_regression_err(buffer1,buffer2,c0,c1);
  }
  if(verbose)
    std::cout << "linear regression based on " << buffer1.size() << " points: " << c0 << "+" << c1 << " * x " << " with rmse: " << err << std::endl;
  return err;
}

double ImgRegression::getR2(ImgRasterGdal& imgReader1, ImgRasterGdal& imgReader2, double& c0, double& c1, unsigned int band1, unsigned int band2, short verbose) const{
  c0=0;
  c1=1;
  int icol1=0,irow1=0;
  std::vector<double> rowBuffer1(imgReader1.nrOfCol());
  std::vector<double> rowBuffer2(imgReader2.nrOfCol());
  std::vector<double> buffer1;
  std::vector<double> buffer2;

  srand(time(NULL));
  for(irow1=0;irow1<imgReader1.nrOfRow();++irow1){
    if(irow1%m_down)
      continue;
    icol1=0;
    double icol2=0,irow2=0;
    double geox=0,geoy=0;
    imgReader1.readData(rowBuffer1,GDT_Float64,irow1,band1);
    imgReader1.image2geo(icol1,irow1,geox,geoy);
    imgReader2.geo2image(geox,geoy,icol2,irow2);
    icol2=static_cast<int>(icol2);
    irow2=static_cast<int>(irow2);
    if(irow2<0||irow2>=imgReader2.nrOfRow())
      continue;
    imgReader2.readData(rowBuffer2,GDT_Float64,irow2,band2);
    for(icol1=0;icol1<imgReader1.nrOfCol();++icol1){
      if(icol1%m_down)
	continue;
      if(m_threshold>0){//percentual value
	double p=static_cast<double>(rand())/(RAND_MAX);
	p*=100.0;
	if(p>m_threshold)
	  continue;//do not select for now, go to next column
      }
      imgReader1.image2geo(icol1,irow1,geox,geoy);
      imgReader2.geo2image(geox,geoy,icol2,irow2);
      if(icol2<0||icol2>=imgReader2.nrOfCol())
	continue;
      icol2=static_cast<int>(icol2);
      irow2=static_cast<int>(irow2);
      //check for nodata
      double value1=rowBuffer1[icol1];
      double value2=rowBuffer2[icol2];
      if(imgReader1.isNoData(value1)||imgReader2.isNoData(value2))
	continue;

      buffer1.push_back(value1);
      buffer2.push_back(value2);
      if(verbose>1)
	std::cout << geox << " " << geoy << " " << icol1 << " " << irow1 << " " << icol2 << " " << irow2 << " " << buffer1.back() << " " << buffer2.back() << std::endl;
    }
  }
  double r2=0;
  if(buffer1.size()&&buffer2.size()){
    statfactory::StatFactory stat;
    r2=stat.linear_regression(buffer1,buffer2,c0,c1);
  }
  if(verbose)
    std::cout << "linear regression based on " << buffer1.size() << " points: " << c0 << "+" << c1 << " * x " << " with r^2: " << r2 << std::endl;
  return r2;
}

double ImgRegression::pgetR2(ImgRasterGdal& imgReader1, ImgRasterGdal& imgReader2, double& c0, double& c1, unsigned int band1, unsigned int band2, short verbose) const{
  c0=0;
  c1=1;
  int icol1=0,irow1=0;
  std::vector<double> rowBuffer1(imgReader1.nrOfCol());
  std::vector<double> rowBuffer2(imgReader2.nrOfCol());
  std::vector<double> buffer1;
  std::vector<double> buffer2;

  srand(time(NULL));
  for(irow1=0;irow1<imgReader1.nrOfRow();++irow1){
    if(irow1%m_down)
      continue;
    icol1=0;
    double icol2=0,irow2=0;
    double geox=0,geoy=0;
    imgReader1.readData(rowBuffer1,GDT_Float64,irow1,band1);
    imgReader1.image2geo(icol1,irow1,geox,geoy);
    imgReader2.geo2image(geox,geoy,icol2,irow2);
    icol2=static_cast<int>(icol2);
    irow2=static_cast<int>(irow2);
    if(irow2<0||irow2>=imgReader2.nrOfRow())
      continue;
    imgReader2.readData(rowBuffer2,GDT_Float64,irow2,band2);
    for(icol1=0;icol1<imgReader1.nrOfCol();++icol1){
      if(icol1%m_down)
	continue;
      if(m_threshold>0){//percentual value
	double p=static_cast<double>(rand())/(RAND_MAX);
	p*=100.0;
	if(p>m_threshold)
	  continue;//do not select for now, go to next column
      }
      imgReader1.image2geo(icol1,irow1,geox,geoy);
      imgReader2.geo2image(geox,geoy,icol2,irow2);
      if(icol2<0||icol2>=imgReader2.nrOfCol())
	continue;
      icol2=static_cast<int>(icol2);
      irow2=static_cast<int>(irow2);
      //check for nodata
      double value1=rowBuffer1[icol1];
      double value2=rowBuffer2[icol2];
      if(imgReader1.isNoData(value1)||imgReader2.isNoData(value2))
	continue;

      buffer1.push_back(value1);
      buffer2.push_back(value2);
      if(verbose>1)
	std::cout << geox << " " << geoy << " " << icol1 << " " << irow1 << " " << icol2 << " " << irow2 << " " << buffer1.back() << " " << buffer2.back() << std::endl;
    }
  }
  double r=0;
  if(buffer1.size()&&buffer2.size()){
    statfactory::StatFactory stat;
    r=stat.correlation(buffer1,buffer2);
    // r=stat.gsl_correlation(buffer1,buffer2);
    double m1=0;
    double v1=0;
    double m2=0;
    double v2=0;
    stat.meanVar(buffer1,m1,v1);
    stat.meanVar(buffer2,m2,v2);
    if(v1>0){
      if(r>=0)
	c1=v2/v1;
      else
	c1=-v2/v1;
    }
    c0=m2-c1*m1;
  }
  if(verbose)
    std::cout << "orthogonal regression based on " << buffer1.size() << " points: " << c0 << "+" << c1 << " * x " << " with r^2: " << r*r << std::endl;
  return r*r;
}

double ImgRegression::getRMSE(ImgRasterGdal& imgReader, unsigned int band1, unsigned int band2, double& c0, double& c1, short verbose) const{
  c0=0;
  c1=1;
  int icol=0,irow=0;
  std::vector<double> rowBuffer1(imgReader.nrOfCol());
  std::vector<double> rowBuffer2(imgReader.nrOfCol());
  std::vector<double> buffer1;
  std::vector<double> buffer2;

  srand(time(NULL));
  assert(band1>=0);
  assert(band1<imgReader.nrOfBand());
  assert(band2>=0);
  assert(band2<imgReader.nrOfBand());
  for(irow=0;irow<imgReader.nrOfRow();++irow){
    if(irow%m_down)
      continue;
    icol=0;
    imgReader.readData(rowBuffer1,GDT_Float64,irow,band1);
    imgReader.readData(rowBuffer2,GDT_Float64,irow,band2);
    for(icol=0;icol<imgReader.nrOfCol();++icol){
      if(icol%m_down)
	continue;
      if(m_threshold>0){//percentual value
	double p=static_cast<double>(rand())/(RAND_MAX);
	p*=100.0;
	if(p>m_threshold)
	  continue;//do not select for now, go to next column
      }
      //check for nodata
      double value1=rowBuffer1[icol];
      double value2=rowBuffer2[icol];
      if(imgReader.isNoData(value1)||imgReader.isNoData(value2))
	continue;

      buffer1.push_back(value1);
      buffer2.push_back(value2);
      if(verbose>1)
	std::cout << icol << " " << irow << " " << buffer1.back() << " " << buffer2.back() << std::endl;
    }
  }
  double err=0;
  if(buffer1.size()&&buffer2.size()){
    statfactory::StatFactory stat;
    err=stat.linear_regression_err(buffer1,buffer2,c0,c1);
  }
  if(verbose)
    std::cout << "linear regression based on " << buffer1.size() << " points: " << c0 << "+" << c1 << " * x " << " with rmse: " << err << std::endl;
  return err;
}

double ImgRegression::getR2(ImgRasterGdal& imgReader, unsigned int band1, unsigned int band2, double& c0, double& c1, short verbose) const{
  c0=0;
  c1=1;
  int icol=0,irow=0;
  std::vector<double> rowBuffer1(imgReader.nrOfCol());
  std::vector<double> rowBuffer2(imgReader.nrOfCol());
  std::vector<double> buffer1;
  std::vector<double> buffer2;

  srand(time(NULL));
  assert(band1>=0);
  assert(band1<imgReader.nrOfBand());
  assert(band2>=0);
  assert(band2<imgReader.nrOfBand());
  for(irow=0;irow<imgReader.nrOfRow();++irow){
    if(irow%m_down)
      continue;
    icol=0;
    imgReader.readData(rowBuffer1,GDT_Float64,irow,band1);
    imgReader.readData(rowBuffer2,GDT_Float64,irow,band2);
    for(icol=0;icol<imgReader.nrOfCol();++icol){
      if(icol%m_down)
	continue;
      if(m_threshold>0){//percentual value
	double p=static_cast<double>(rand())/(RAND_MAX);
	p*=100.0;
	if(p>m_threshold)
	  continue;//do not select for now, go to next column
      }
      //check for nodata
      double value1=rowBuffer1[icol];
      double value2=rowBuffer2[icol];
      if(imgReader.isNoData(value1)||imgReader.isNoData(value2))
	continue;

      buffer1.push_back(value1);
      buffer2.push_back(value2);
      if(verbose>1)
	std::cout << icol << " " << irow << " " << buffer1.back() << " " << buffer2.back() << std::endl;
    }
  }
  double r2=0;
  if(buffer1.size()&&buffer2.size()){
    statfactory::StatFactory stat;
    r2=stat.linear_regression(buffer1,buffer2,c0,c1);
  }
  if(verbose)
    std::cout << "linear regression based on " << buffer1.size() << " points: " << c0 << "+" << c1 << " * x " << " with r^2: " << r2 << std::endl;
  return r2;
}

double ImgRegression::pgetR2(ImgRasterGdal& imgReader, unsigned int band1, unsigned int band2, double& c0, double& c1, short verbose) const{
  c0=0;
  c1=1;
  int icol=0,irow=0;
  std::vector<double> rowBuffer1(imgReader.nrOfCol());
  std::vector<double> rowBuffer2(imgReader.nrOfCol());
  std::vector<double> buffer1;
  std::vector<double> buffer2;

  srand(time(NULL));
  assert(band1>=0);
  assert(band1<imgReader.nrOfBand());
  assert(band2>=0);
  assert(band2<imgReader.nrOfBand());
  for(irow=0;irow<imgReader.nrOfRow();++irow){
    if(irow%m_down)
      continue;
    icol=0;
    imgReader.readData(rowBuffer1,GDT_Float64,irow,band1);
    imgReader.readData(rowBuffer2,GDT_Float64,irow,band2);
    for(icol=0;icol<imgReader.nrOfCol();++icol){
      if(icol%m_down)
	continue;
      if(m_threshold>0){//percentual value
	double p=static_cast<double>(rand())/(RAND_MAX);
	p*=100.0;
	if(p>m_threshold)
	  continue;//do not select for now, go to next column
      }
      //check for nodata
      double value1=rowBuffer1[icol];
      double value2=rowBuffer2[icol];
      if(imgReader.isNoData(value1)||imgReader.isNoData(value2))
	continue;

      buffer1.push_back(value1);
      buffer2.push_back(value2);
      if(verbose>1)
	std::cout << icol << " " << irow << " " << buffer1.back() << " " << buffer2.back() << std::endl;
    }
  }
  double r=0;
  if(buffer1.size()&&buffer2.size()){
    statfactory::StatFactory stat;
    r=stat.correlation(buffer1,buffer2);
    // r=stat.gsl_correlation(buffer1,buffer2);
    double m1=0;
    double v1=0;
    double m2=0;
    double v2=0;
    stat.meanVar(buffer1,m1,v1);
    stat.meanVar(buffer2,m2,v2);
    if(v1>0){
      if(r>=0)
	c1=v2/v1;
      else
	c1=-v2/v1;
    }
    c0=m2-c1*m1;
  }
  if(verbose)
    std::cout << "orthogonal regression based on " << buffer1.size() << " points: " << c0 << "+" << c1 << " * x " << " with r^2: " << r*r << std::endl;
  return r*r;
}
