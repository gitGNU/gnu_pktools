/**********************************************************************
ImgReaderGdal.h: class to read raster files using GDAL API library
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
#ifndef _IMGREADERGDAL_H_
#define _IMGREADERGDAL_H_

#include "ImgRasterGdal.h"
#include <assert.h>
#include <fstream>
#include <string>
#include <sstream>
#include "gdal_priv.h"
#include "base/Vector2d.h"

//--------------------------------------------------------------------------
class ImgReaderGdal : public virtual ImgRasterGdal
{
public:
  ImgReaderGdal(void);
  ImgReaderGdal(const std::string& filename, const GDALAccess& readMode=GA_ReadOnly){open(filename, readMode);};
  ~ImgReaderGdal(void);
  void open(const std::string& filename, const GDALAccess& readMode=GA_ReadOnly);
  void close(void);

  void setScale(double theScale, int band=0){
    /* if(getRasterBand(band)->SetScale(theScale)==CE_Failure){ */
    if(m_scale.size()!=nrOfBand()){//initialize
      m_scale.resize(nrOfBand());
      for(int iband=0;iband<nrOfBand();++iband)
	m_scale[iband]=1.0;
    }
    m_scale[band]=theScale;
    /* }; */
  }
  void setOffset(double theOffset, int band=0){
    /* if(getRasterBand(band)->SetOffset(theOffset)==CE_Failure){ */
    if(m_offset.size()!=nrOfBand()){
      m_offset.resize(nrOfBand());
      for(int iband=0;iband<nrOfBand();++iband)
	m_offset[iband]=0.0;
    }
      m_offset[band]=theOffset;
    /* }; */
  }
  template<typename T> void readData(T& value, const GDALDataType& dataType, int col, int row, int band=0) const;
  template<typename T> void readData(std::vector<T>& buffer, const GDALDataType& dataType , int minCol, int maxCol, int row, int band=0) const;
  template<typename T> void readData(std::vector<T>& buffer, const GDALDataType& dataType , int minCol, int maxCol, double row, int band=0, RESAMPLE resample=NEAR) const;
  template<typename T> void readDataBlock(Vector2d<T>& buffer2d, const GDALDataType& dataType , int minCol, int maxCol, int minRow, int maxRow, int band=0) const;
  template<typename T> void readDataBlock(std::vector<T>& buffer, const GDALDataType& dataType , int minCol, int maxCol, int minRow, int maxRow, int band=0) const;
  template<typename T> void readData(std::vector<T>& buffer, const GDALDataType& dataType, int row, int band=0) const;
  template<typename T> void readData(std::vector<T>& buffer, const GDALDataType& dataType, double row, int band=0, RESAMPLE resample=NEAR) const;
  void getMinMax(int startCol, int endCol, int startRow, int endRow, int band, double& minValue, double& maxValue) const;
  void getMinMax(double& minValue, double& maxValue, int band=0) const;
  double getMin(int& col, int& row, int band=0) const;
  double getHistogram(std::vector<double>& histvector, double& min, double& max,unsigned int& nbin, int theBand=0, bool kde=false);
  double getMax(int& col, int& row, int band=0) const;
  void getRefPix(double& refX, double &refY, int band=0) const;
  void getRange(std::vector<short>& range, int Band=0) const;
  unsigned long int getNvalid(int band) const;

protected:
  void setCodec(const GDALAccess& readMode=GA_ReadOnly);

  std::vector<double> m_scale;
  std::vector<double> m_offset;
};

//     adfGeoTransform[0] /* top left x */
//     adfGeoTransform[1] /* w-e pixel resolution */
//     adfGeoTransform[2] /* rotation, 0 if image is "north up" */
//     adfGeoTransform[3] /* top left y */
//     adfGeoTransform[4] /* rotation, 0 if image is "north up" */
//     adfGeoTransform[5] /* n-s pixel resolution */

template<typename T> void ImgReaderGdal::readData(T& value, const GDALDataType& dataType, int col, int row, int band) const
{
  //fetch raster band
  GDALRasterBand  *poBand;
  assert(band<nrOfBand()+1);
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  assert(col<nrOfCol());
  assert(col>=0);
  assert(row<nrOfRow());
  assert(row>=0);
  poBand->RasterIO(GF_Read,col,row,1,1,&value,1,1,dataType,0,0);
  if(m_scale.size()>band)
    value=static_cast<double>(value)*m_scale[band];
  if(m_offset.size()>band)
    value=static_cast<double>(value)+m_offset[band];
}

template<typename T> void ImgReaderGdal::readData(std::vector<T>& buffer, const GDALDataType& dataType, int minCol, int maxCol, int row, int band) const
{
  //fetch raster band
  GDALRasterBand  *poBand;
  assert(band<nrOfBand()+1);
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  assert(minCol<nrOfCol());
  assert(minCol>=0);
  assert(maxCol<nrOfCol());
  assert(minCol<=maxCol);
  assert(row<nrOfRow());
  assert(row>=0);
  if(buffer.size()!=maxCol-minCol+1)
    buffer.resize(maxCol-minCol+1);
  poBand->RasterIO(GF_Read,minCol,row,buffer.size(),1,&(buffer[0]),buffer.size(),1,dataType,0,0);
  if(m_scale.size()>band||m_offset.size()>band){
    double theScale=1;
    double theOffset=0;
    if(m_scale.size()>band)
      theScale=m_scale[band];
    if(m_offset.size()>band)
      theOffset=m_offset[band];
    for(int index=0;index<buffer.size();++index)
      buffer[index]=theScale*static_cast<double>(buffer[index])+theOffset;
  }
}

template<typename T> void ImgReaderGdal::readData(std::vector<T>& buffer, const GDALDataType& dataType , int minCol, int maxCol, double row, int band, RESAMPLE resample) const
{
  //todo: make upper and lower row depend on isGeo...
  std::vector<T> readBuffer_upper;
  std::vector<T> readBuffer_lower;
  if(buffer.size()!=maxCol-minCol+1)
    buffer.resize(maxCol-minCol+1);
  double upperRow=row-0.5;
  upperRow=static_cast<int>(upperRow);
  double lowerRow=row+0.5;
  lowerRow=static_cast<int>(lowerRow);
  switch(resample){
  case(BILINEAR):
    if(lowerRow>=nrOfRow())
      lowerRow=nrOfRow()-1;
    if(upperRow<0)
      upperRow=0;
    readData(readBuffer_upper,GDT_Float64,minCol,maxCol,static_cast<int>(upperRow),band);
    readData(readBuffer_lower,GDT_Float64,minCol,maxCol,static_cast<int>(lowerRow),band);
    //do interpolation in y
    for(int icol=0;icol<maxCol-minCol+1;++icol){
      buffer[icol]=(lowerRow-row+0.5)*readBuffer_upper[icol]+(1-lowerRow+row-0.5)*readBuffer_lower[icol];
    }
    break;
  default:
    readData(buffer,GDT_Float64,minCol,maxCol,static_cast<int>(row),band);
    break;
  }
}

template<typename T> void ImgReaderGdal::readDataBlock(Vector2d<T>& buffer2d, const GDALDataType& dataType , int minCol, int maxCol, int minRow, int maxRow, int band) const
{
  buffer2d.resize(maxRow-minRow+1);
  typename std::vector<T> buffer;
  readDataBlock(buffer,dataType,minCol,maxCol,minRow,maxRow,band);
  typename std::vector<T>::const_iterator startit=buffer.begin();
  typename std::vector<T>::const_iterator endit=startit;
  for(int irow=minRow;irow<=maxRow;++irow){
    buffer2d[irow-minRow].resize(maxCol-minCol+1);
    endit+=maxCol-minCol+1;
    buffer2d[irow-minRow].assign(startit,endit);
    startit+=maxCol-minCol+1;
  }
}
  
template<typename T> void ImgReaderGdal::readDataBlock(std::vector<T>& buffer, const GDALDataType& dataType , int minCol, int maxCol, int minRow, int maxRow, int band) const
{
  double theScale=1;
  double theOffset=0;
  if(m_scale.size()>band)
    theScale=m_scale[band];
  if(m_offset.size()>band)
    theOffset=m_offset[band];
  //fetch raster band
  GDALRasterBand  *poBand;
  assert(band<nrOfBand()+1);
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  if(minCol>=nrOfCol() ||
     (minCol<0) ||
     (maxCol>=nrOfCol()) ||
     (minCol>maxCol) ||
     (minRow>=nrOfRow()) ||
     (minRow<0) ||
     (maxRow>=nrOfRow()) ||
     (minRow>maxRow)){
    std::string errorString="block not within image boundaries";
    throw(errorString);
  }
  /* assert(minCol<nrOfCol()); */
  /* assert(minCol>=0); */
  /* assert(maxCol<nrOfCol()); */
  /* assert(minCol<=maxCol); */
  /* assert(minRow<nrOfRow()); */
  /* assert(minRow>=0); */
  /* assert(maxRow<nrOfRow()); */
  /* assert(minRow<=maxRow); */
  if(buffer.size()!=(maxRow-minRow+1)*(maxCol-minCol+1))
    buffer.resize((maxRow-minRow+1)*(maxCol-minCol+1));
  poBand->RasterIO(GF_Read,minCol,minRow,maxCol-minCol+1,maxRow-minRow+1,&(buffer[0]),(maxCol-minCol+1),(maxRow-minRow+1),dataType,0,0);
  if(m_scale.size()>band||m_offset.size()>band){
    for(int index=0;index<buffer.size();++index)
      buffer[index]=theScale*buffer[index]+theOffset;
  }
}

// template<typename T> void ImgReaderGdal::readDataBlock(vector<T>& buffer, const GDALDataType& dataType , int minCol, int maxCol, int minRow, int maxRow, int band) const
// {
//   assert(band<nrOfBand()+1);
//   assert(minCol<nrOfCol());
//   assert(minCol>=0);
//   assert(maxCol<nrOfCol());
//   assert(minCol<=maxCol);
//   assert(minRow<nrOfRow());
//   assert(minRow>=0);
//   assert(maxRow<nrOfRow());
//   assert(minRow<=maxRow);
//   if(buffer.size()!=(maxRow-minRow+1)*(maxCol-minCol+1))
//     buffer.resize((maxRow-minRow+1)*(maxCol-minCol+1));
//   //fetch raster band
//   GDALRasterBand  *poBand;
//   assert(band<nrOfBand()+1);
//   poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
//   for(int irow=0;irow<maxRow-minRow+1;++irow)
//     poBand->RasterIO(GF_Read,minCol,minRow+irow,maxCol-minCol+1,1,&(buffer[irow*(maxCol-minCol+1)]),maxCol-minCol+1,1,dataType,0,0);
// }
  
template<typename T> void ImgReaderGdal::readData(std::vector<T>& buffer, const GDALDataType& dataType, int row, int band) const
{
  readData(buffer,dataType,0,nrOfCol()-1,row,band);
}

template<typename T> void ImgReaderGdal::readData(std::vector<T>& buffer, const GDALDataType& dataType, double row, int band, RESAMPLE resample) const
{
  readData(buffer,dataType,0,nrOfCol()-1,row,band,resample);
}


#endif // _IMGREADERGDAL_H_

//       //fetch raster band
//   GDALRasterBand  *poBand;
//   assert(band<nrOfBand()+1);
//   poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
//   buffer.resize(maxCol-minCol+1);
//   assert(minCol<nrOfCol());
//   assert(row<nrOfRow());
//   poBand->RasterIO(GF_Read,minCol,row,buffer.size(),1,&(buffer[0]),buffer.size(),1,GDT_Int16,0,0);
