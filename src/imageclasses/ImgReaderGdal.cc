/**********************************************************************
ImgReaderGdal.cc: class to read raster files using GDAL API library
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
#include <sstream>
#include <iostream>
#include <gsl/gsl_cdf.h>
#include "cpl_conv.h" // for CPLMalloc()
#include "ImgReaderGdal.h"

ImgReaderGdal::ImgReaderGdal(void){};

ImgReaderGdal::~ImgReaderGdal(void){};
// {
//   if(m_data.size()&&m_filename.size()){
//     for(int iband=0;iband<m_nband;++iband)
//       free(m_data[iband]);
//   }
// }

/**
 * @param filename Open a raster dataset with this filename
 * @param readMode Open dataset in ReadOnly or Update mode
 * @param memory Available memory to cache image raster data (in MB)
 **/
void ImgReaderGdal::open(const std::string& filename, const GDALAccess& readMode, unsigned long int memory)
{
  m_filename = filename;
  setCodec(readMode);
  initMem(memory);
  for(int iband=0;iband<m_nband;++iband){
    m_begin[iband]=0;
    m_end[iband]=0;
  }
}

void ImgReaderGdal::close(void)
{
  ImgRasterGdal::close();
}

/**
 * @param readMode Open dataset in ReadOnly or Update mode
 **/
void ImgReaderGdal::setCodec(const GDALAccess& readMode)
{
  GDALAllRegister();
  // m_gds = (GDALDataset *) GDALOpen(m_filename.c_str(), readMode );
#if GDAL_VERSION_MAJOR < 2
  GDALAllRegister();
  m_gds = (GDALDataset *) GDALOpen(m_filename.c_str(), readMode );
#else
  GDALAllRegister();
  if(readMode==GA_ReadOnly)
    m_gds = (GDALDataset*) GDALOpenEx(m_filename.c_str(), GDAL_OF_READONLY|GDAL_OF_RASTER, NULL, NULL, NULL);
  else if(readMode==GA_Update)
    m_gds = (GDALDataset*) GDALOpenEx(m_filename.c_str(), GDAL_OF_UPDATE|GDAL_OF_RASTER, NULL, NULL, NULL);
#endif

  if(m_gds == NULL){
    std::string errorString="FileOpenError";
    throw(errorString);
  }
  m_ncol= m_gds->GetRasterXSize();
  m_nrow= m_gds->GetRasterYSize();
  m_nband= m_gds->GetRasterCount();
  double adfGeoTransform[6];
  m_gds->GetGeoTransform( adfGeoTransform );
  m_gt[0]=adfGeoTransform[0];
  m_gt[1]=adfGeoTransform[1];
  m_gt[2]=adfGeoTransform[2];
  m_gt[3]=adfGeoTransform[3];
  m_gt[4]=adfGeoTransform[4];
  m_gt[5]=adfGeoTransform[5];
  m_projection=m_gds->GetProjectionRef();
}

/**
 * @param row Read a new block for caching this row (if needed)
 * @param band Band that must be read to cache
 * @return true if block was read
 **/
bool ImgReaderGdal::readNewBlock(int row, int band)
{
  if(m_end[band]<m_blockSize)//first time
    m_end[band]=m_blockSize;
  while(row>=m_end[band]&&m_begin[band]<nrOfRow()){
    m_begin[band]+=m_blockSize;
    m_end[band]=m_begin[band]+m_blockSize;
  }
  if(m_end[band]>nrOfRow())
    m_end[band]=nrOfRow();
  for(int iband=0;iband<m_nband;++iband){
    //fetch raster band
    GDALRasterBand  *poBand;
    assert(iband<nrOfBand()+1);
    poBand = m_gds->GetRasterBand(iband+1);//GDAL uses 1 based index
    poBand->RasterIO(GF_Read,0,m_begin[iband],nrOfCol(),m_end[iband]-m_begin[iband],m_data[iband],nrOfCol(),m_end[iband]-m_begin[iband],getDataType(),0,0);
  }
  return true;//new block was read
}

/**
 * @param x Reported column where minimum value in image was found (start counting from 0)
 * @param y Reported row where minimum value in image was found (start counting from 0)
 * @param band Search mininum value in image for this band
 * @return minimum value in image for the selected band
 **/
double ImgReaderGdal::getMin(int& x, int& y, int band){
  double minValue=0;
  std::vector<double> lineBuffer(nrOfCol());
  bool isValid=false;
  for(int irow=0;irow<nrOfRow();++irow){
    readData(lineBuffer,irow,band);
    for(int icol=0;icol<nrOfCol();++icol){
      if(isNoData(lineBuffer[icol]))
	continue;
      if(isValid){
	if(lineBuffer[icol]<minValue){
          y=irow;
          x=icol;
          minValue=lineBuffer[icol];
        }
      }
      else{
	y=irow;
	x=icol;
	minValue=lineBuffer[icol];
	isValid=true;
      }
    }
  }
  if(isValid)
    return minValue;
  else
    throw(static_cast<std::string>("Warning: not initialized"));
}

/**
 * @param x Reported column where maximum value in image was found (start counting from 0)
 * @param y Reported row where maximum value in image was found (start counting from 0)
 * @param band Search mininum value in image for this band
 * @return maximum value in image for the selected band
 **/
double ImgReaderGdal::getMax(int& x, int& y, int band){
  double maxValue=0;
  std::vector<double> lineBuffer(nrOfCol());
  bool isValid=false;
  for(int irow=0;irow<nrOfRow();++irow){
    readData(lineBuffer,irow,band);
    for(int icol=0;icol<nrOfCol();++icol){
      if(isNoData(lineBuffer[icol]))
	continue;
      if(isValid){
	if(lineBuffer[icol]>maxValue){
          y=irow;
          x=icol;
          maxValue=lineBuffer[icol];
        }
      }
      else{
	y=irow;
	x=icol;
	maxValue=lineBuffer[icol];
	isValid=true;
      }
    }
  }
  if(isValid)
    return maxValue;
  else
    throw(static_cast<std::string>("Warning: not initialized"));
}

/**
 * @param startCol, endCol, startRow, endRow Search extreme value in this region of interest (all indices start counting from 0)
 * @param band Search extreme value in image for this band
 * @param minValue Reported minimum value within searched region
 * @param maxValue Reported maximum value within searched region
 **/
void ImgReaderGdal::getMinMax(int startCol, int endCol, int startRow, int endRow, int band, double& minValue, double& maxValue)
{
  bool isConstraint=(maxValue>minValue);
  double minConstraint=minValue;
  double maxConstraint=maxValue;
  std::vector<double> lineBuffer(endCol-startCol+1);
  bool isValid=false;
  assert(endRow<nrOfRow());
  for(int irow=startCol;irow<endRow+1;++irow){
    readData(lineBuffer,startCol,endCol,irow,band);
    for(int icol=0;icol<lineBuffer.size();++icol){
      if(isNoData(lineBuffer[icol]))
	continue;
      if(isValid){
	if(isConstraint){
	  if(lineBuffer[icol]<minConstraint)
	    continue;
	  if(lineBuffer[icol]>maxConstraint)
	    continue;
	}
	if(lineBuffer[icol]<minValue)
	  minValue=lineBuffer[icol];
	if(lineBuffer[icol]>maxValue)
	  maxValue=lineBuffer[icol];
      }
      else{
	if(isConstraint){
	  if(lineBuffer[icol]<minConstraint)
	    continue;
	  if(lineBuffer[icol]>maxConstraint)
	    continue;
	}
	minValue=lineBuffer[icol];
	maxValue=lineBuffer[icol];
	isValid=true;
      }
    }
  }
  if(!isValid)
    throw(static_cast<std::string>("Warning: not initialized"));
}

/**
 * @param minValue Reported minimum value in image
 * @param maxValue Reported maximum value in image
 * @param band Search extreme value in image for this band
 **/
void ImgReaderGdal::getMinMax(double& minValue, double& maxValue, int band)
{
  bool isConstraint=(maxValue>minValue);
  double minConstraint=minValue;
  double maxConstraint=maxValue;
  std::vector<double> lineBuffer(nrOfCol());
  bool isValid=false;
  for(int irow=0;irow<nrOfRow();++irow){
    readData(lineBuffer,irow,band);
    for(int icol=0;icol<nrOfCol();++icol){
      if(isNoData(lineBuffer[icol]))
	continue;
      if(isValid){
	if(isConstraint){
	  if(lineBuffer[icol]<minConstraint)
	    continue;
	  if(lineBuffer[icol]>maxConstraint)
	    continue;
	}
	if(lineBuffer[icol]<minValue)
	  minValue=lineBuffer[icol];
	if(lineBuffer[icol]>maxValue)
	  maxValue=lineBuffer[icol];
      }
      else{
	if(isConstraint){
	  if(lineBuffer[icol]<minConstraint)
	    continue;
	  if(lineBuffer[icol]>maxConstraint)
	    continue;
	}
	minValue=lineBuffer[icol];
	maxValue=lineBuffer[icol];
	isValid=true;
      }
    }
  }
  if(!isValid)
    throw(static_cast<std::string>("Warning: not initialized"));
}

/**
 * @param histvector The reported histogram with counts per bin
 * @param min, max Only calculate histogram for values between min and max. If min>=max, calculate min and max from the image
 * @param nbin Number of bins used for calculating the histogram. If nbin is 0, the number of bins is  automatically calculated from min and max
 * @param theBand The band for which to calculate the histogram (start counting from 0)
 * @param kde Apply kernel density function for a Gaussian basis function
 * @return number of valid pixels in this dataset for the the selected band
 **/
double ImgReaderGdal::getHistogram(std::vector<double>& histvector, double& min, double& max, unsigned int& nbin, int theBand, bool kde){
  double minValue=0;
  double maxValue=0;
      
  if(min>=max)
    getMinMax(minValue,maxValue,theBand);
  else{
    minValue=min;
    maxValue=max;
  }
  if(min<max&&min>minValue)
    minValue=min;
  if(min<max&&max<maxValue)
    maxValue=max;
  min=minValue;
  max=maxValue;

  double sigma=0;
  if(kde){
    double meanValue=0;
    double stdDev=0;
    GDALProgressFunc pfnProgress;
    void* pProgressData;
    GDALRasterBand* rasterBand;
    rasterBand=getRasterBand(theBand);
    rasterBand->ComputeStatistics(0,&minValue,&maxValue,&meanValue,&stdDev,pfnProgress,pProgressData);
    //rest minvalue and MaxValue as ComputeStatistics does not account for nodata, scale and offset
    minValue=min;
    maxValue=max;

    if(m_scale.size()>theBand){
      stdDev*=m_scale[theBand];
    }
    sigma=1.06*stdDev*pow(getNvalid(theBand),-0.2);
  }

  double scale=0;
  if(maxValue>minValue){
    if(nbin==0)
      nbin=maxValue-minValue+1;
    scale=static_cast<double>(nbin-1)/(maxValue-minValue);
  }
  else
    nbin=1;
  assert(nbin>0);
  if(histvector.size()!=nbin){
    histvector.resize(nbin);
    for(int i=0;i<nbin;histvector[i++]=0);
  }
  double nvalid=0;
  unsigned long int nsample=0;
  unsigned long int ninvalid=0;
  std::vector<double> lineBuffer(nrOfCol());
  for(int irow=0;irow<nrOfRow();++irow){
    readData(lineBuffer,irow,theBand);
    for(int icol=0;icol<nrOfCol();++icol){
      if(isNoData(lineBuffer[icol]))
        ++ninvalid;
      else if(lineBuffer[icol]>maxValue)
        ++ninvalid;
      else if(lineBuffer[icol]<minValue)
        ++ninvalid;
      else if(nbin==1)
	++histvector[0];
      else{//scale to [0:nbin]
	if(sigma>0){
	  //create kde for Gaussian basis function
	  //todo: speed up by calculating first and last bin with non-zero contriubtion...
	  //todo: calculate real surface below pdf by using gsl_cdf_gaussian_P(x-mean+binsize,sigma)-gsl_cdf_gaussian_P(x-mean,sigma)
	  for(int ibin=0;ibin<nbin;++ibin){
	    double icenter=minValue+static_cast<double>(maxValue-minValue)*(ibin+0.5)/nbin;
	    double thePdf=gsl_ran_gaussian_pdf(lineBuffer[icol]-icenter, sigma);
	    histvector[ibin]+=thePdf;
	    nvalid+=thePdf;
	  }
	}
	else{
	  int theBin=static_cast<unsigned long int>(scale*(lineBuffer[icol]-minValue));
	  assert(theBin>=0);
	  assert(theBin<nbin);
	  ++histvector[theBin];
	  ++nvalid;
	}
      // else if(lineBuffer[icol]==maxValue)
      //   ++histvector[nbin-1];
      // else
      //   ++histvector[static_cast<int>(static_cast<double>(lineBuffer[icol]-minValue)/(maxValue-minValue)*(nbin-1))];
      }
    }
  }
  // unsigned long int nvalid=nrOfCol()*nrOfRow()-ninvalid;
  return nvalid;
}

/**
 * @param range Sorted vector containing the range of image values
 * @param band The band for which to calculate the range
 **/
void ImgReaderGdal::getRange(std::vector<short>& range, int band)
{
  std::vector<short> lineBuffer(nrOfCol());
  range.clear();
  for(int irow=0;irow<nrOfRow();++irow){
    readData(lineBuffer,irow,band);
    for(int icol=0;icol<nrOfCol();++icol){
      if(find(range.begin(),range.end(),lineBuffer[icol])==range.end())
        range.push_back(lineBuffer[icol]);
    }
  }
  sort(range.begin(),range.end());
}

/**
 * @param band The band for which to calculate the number of valid pixels
 * @return number of valid pixels in this dataset for the the selected band
 **/
unsigned long int ImgReaderGdal::getNvalid(int band)
{
  unsigned long int nvalid=0;
  if(m_noDataValues.size()){
    std::vector<double> lineBuffer(nrOfCol());
    for(int irow=0;irow<nrOfRow();++irow){
      readData(lineBuffer,irow,band);
      for(int icol=0;icol<nrOfCol();++icol){
	if(isNoData(lineBuffer[icol]))
	  continue;
	else
	  ++nvalid;
      }
    }
    return nvalid;
  }
  else
    return(nrOfCol()*nrOfRow());
}

/**
 * @param refX, refY Calculated reference pixel position in geo-refererenced coordinates
 * @param band The band for which to calculate the number of valid pixels
 **/

void ImgReaderGdal::getRefPix(double& refX, double &refY, int band)
{
  std::vector<double> lineBuffer(nrOfCol());
  double validCol=0;
  double validRow=0;
  int nvalidCol=0;
  int nvalidRow=0;
  for(int irow=0;irow<nrOfRow();++irow){
    readData(lineBuffer,irow,band);
    for(int icol=0;icol<nrOfCol();++icol){
      // bool valid=(find(m_noDataValues.begin(),m_noDataValues.end(),lineBuffer[icol])==m_noDataValues.end());
      // if(valid){
      if(!isNoData(lineBuffer[icol])){
        validCol+=icol+1;
        ++nvalidCol;
        validRow+=irow+1;
        ++nvalidRow;
      }
    }
  }
  if(isGeoRef()){
    //reference coordinate is lower left corner of pixel in center of gravity
    //we need geo coordinates for exactly this location: validCol(Row)/nvalidCol(Row)-0.5
    double cgravi=validCol/nvalidCol-0.5;
    double cgravj=validRow/nvalidRow-0.5;
    double refpixeli=floor(cgravi);
    double refpixelj=ceil(cgravj-1);
    //but image2geo provides location at center of pixel (shifted half pixel right down)
    image2geo(refpixeli,refpixelj,refX,refY);
    //refX and refY now refer to center of gravity pixel
    refX-=0.5*getDeltaX();//shift to left corner
    refY-=0.5*getDeltaY();//shift to lower left corner
  }
  else{
    refX=floor(validCol/nvalidCol-0.5);//left corner
    refY=floor(validRow/nvalidRow-0.5);//upper corner
    //shift to lower left corner of pixel
    refY+=1;
  }
}
