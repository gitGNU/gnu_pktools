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

#include <assert.h>
#include <fstream>
#include <string>
#include <sstream>
#include "gdal_priv.h"
#include "base/Vector2d.h"

using namespace std;

enum RESAMPLE { NEAR = 0, BILINEAR = 1, BICUBIC = 2 };

//--------------------------------------------------------------------------
class ImgReaderGdal
{
public:
  ImgReaderGdal(void);
  ImgReaderGdal(const string& filename){open(filename);};
  ~ImgReaderGdal(void);
  void open(const string& filename);//, double magicX=1, double magicY=1);
  void close(void);
  int nrOfCol(void) const { return m_ncol;};
  int nrOfRow(void) const { return m_nrow;};
  int nrOfBand(void) const { return m_nband;};
  bool isGeoRef() const {return m_isGeoRef;};
  string getProjection(void) const;
  string getProjectionRef(void) const;
  string getGeoTransform() const;
  void getGeoTransform(double& ulx, double& uly, double& deltaX, double& deltaY, double& rot1, double& rot2) const;
  string getDescription() const;
  string getMetadataItem() const;
  string getImageDescription() const;
  bool getBoundingBox (double& ulx, double& uly, double& lrx, double& lry) const;
  bool getCentrePos(double& x, double& y) const;
  double getUlx() const {return (m_isGeoRef)? m_ulx : 0;};
  double getUly() const {return (m_isGeoRef)? m_uly : nrOfRow()-1;};
  double getLrx() const {return (m_isGeoRef)? m_ulx+nrOfCol()*m_delta_x : nrOfCol()-1;};
  double getLry() const {return (m_isGeoRef)? m_uly-nrOfRow()*m_delta_y : 0;};
  // bool getMagicPixel(double& magicX, double& magicY) const {magicX=m_magic_x;magicY=m_magic_y;};
  int getNoDataValues(vector<double>& noDataValues) const;
  bool isNoData(double value) const{return find(m_noDataValues.begin(),m_noDataValues.end(),value)!=m_noDataValues.end();};
  int pushNoDataValue(double noDataValue);
  CPLErr GDALSetNoDataValue(double noDataValue, int band=0) {getRasterBand(band)->SetNoDataValue(noDataValue);};
  bool covers(double x, double y) const;
  bool covers(double ulx, double  uly, double lrx, double lry) const;
  bool geo2image(double x, double y, double& i, double& j) const;
  bool image2geo(double i, double j, double& x, double& y) const;
  double getDeltaX(void) const {return m_delta_x;};
  double getDeltaY(void) const {return m_delta_y;};
  template<typename T> void readData(T& value, const GDALDataType& dataType, int col, int row, int band=0) const;
  template<typename T> void readData(vector<T>& buffer, const GDALDataType& dataType , int minCol, int maxCol, int row, int band=0) const;
  template<typename T> void readData(vector<T>& buffer, const GDALDataType& dataType , int minCol, int maxCol, double row, int band=0, RESAMPLE resample=0) const;
  template<typename T> void readDataBlock(Vector2d<T>& buffer, const GDALDataType& dataType , int minCol, int maxCol, int minRow, int maxRow, int band=0) const;
  template<typename T> void readDataBlock(vector<T>& buffer, const GDALDataType& dataType , int minCol, int maxCol, int minRow, int maxRow, int band=0) const;
  template<typename T> void readData(vector<T>& buffer, const GDALDataType& dataType, int row, int band=0) const;
  template<typename T> void readData(vector<T>& buffer, const GDALDataType& dataType, double row, int band=0, RESAMPLE resample=0) const;
  void getMinMax(double& minValue, double& maxValue, int band=0, bool exhaustiveSearch=false) const;
  double getMin(int& col, int& row, int band=0) const;
  double getMax(int& col, int& row, int band=0) const;
  void getRefPix(double& refX, double &refY, int band=0) const;
  void getRange(vector<short>& range, int Band=0) const;
  GDALDataType getDataType(int band=0) const;
  GDALRasterBand* getRasterBand(int band=0);
  GDALColorTable* getColorTable(int band=0) const;
  string getDriverDescription() const;
  string getImageType() const{return getDriverDescription();};
//   string getImageType() const{return "GTiff";};
  string getInterleave() const;
  string getCompression() const;
  GDALDataset* getDataset(){return m_gds;};
  char** getMetadata();
  char** getMetadata() const;
  void getMetadata(list<string>& metadata) const;

protected:
  void setCodec();//double magicX, double magicY);

  string m_filename;
  GDALDataset *m_gds;
  int m_ncol;
  int m_nrow;
  int m_nband;
  double m_ulx;
  double m_uly;
  double m_delta_x;
  double m_delta_y;
  // double m_magic_x;
  // double m_magic_y;
  bool m_isGeoRef;
  vector<double> m_noDataValues;
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
}

template<typename T> void ImgReaderGdal::readData(vector<T>& buffer, const GDALDataType& dataType, int minCol, int maxCol, int row, int band) const
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
}

template<typename T> void ImgReaderGdal::readData(vector<T>& buffer, const GDALDataType& dataType , int minCol, int maxCol, double row, int band, RESAMPLE resample) const
{
  vector<T> readBuffer_upper;
  vector<T> readBuffer_lower;
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

template<typename T> void ImgReaderGdal::readDataBlock(Vector2d<T>& buffer, const GDALDataType& dataType , int minCol, int maxCol, int minRow, int maxRow, int band) const
{
  buffer.resize(maxRow-minRow+1);
  for(int irow=minRow;irow<=maxRow;++irow){
    buffer[irow-minRow].resize(maxCol-minCol+1);
    readData(buffer[irow-minRow],dataType,minCol,maxCol,irow,band);
  }
}
  
template<typename T> void ImgReaderGdal::readDataBlock(vector<T>& buffer, const GDALDataType& dataType , int minCol, int maxCol, int minRow, int maxRow, int band) const
{
  //fetch raster band
  GDALRasterBand  *poBand;
  assert(band<nrOfBand()+1);
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  assert(minCol<nrOfCol());
  assert(minCol>=0);
  assert(maxCol<nrOfCol());
  assert(minCol<=maxCol);
  assert(minRow<nrOfRow());
  assert(minRow>=0);
  assert(maxRow<nrOfRow());
  assert(minRow<=maxRow);
  if(buffer.size()!=(maxRow-minRow+1)*(maxCol-minCol+1))
    buffer.resize((maxRow-minRow+1)*(maxCol-minCol+1));
  poBand->RasterIO(GF_Read,minCol,minRow,maxCol-minCol+1,maxRow-minRow+1,&(buffer[0]),(maxCol-minCol+1),(maxRow-minRow+1),dataType,0,0);
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
  
template<typename T> void ImgReaderGdal::readData(vector<T>& buffer, const GDALDataType& dataType, int row, int band) const
{
  readData(buffer,dataType,0,nrOfCol()-1,row,band);
}

template<typename T> void ImgReaderGdal::readData(vector<T>& buffer, const GDALDataType& dataType, double row, int band, RESAMPLE resample) const
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
