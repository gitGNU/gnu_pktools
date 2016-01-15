/**********************************************************************
ImgWriterGdal.h: class to write raster files using GDAL API library
Copyright (C) 2008-2012 Pieter Kempeneers

This file is part of pktools
n
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
#ifndef _IMGWRITERGDAL_H_
#define _IMGWRITERGDAL_H_

#include <assert.h>
#include <fstream>
#include <string>
#include <sstream>
#include "gdal_priv.h"
#include "ImgRasterGdal.h"
#include "ImgReaderGdal.h"
#include "ImgReaderOgr.h"

//--------------------------------------------------------------------------
class ImgWriterGdal : public virtual ImgRasterGdal
{
public:
  ImgWriterGdal(void);
  ~ImgWriterGdal(void);
  void open(const std::string& filename);
  void open(const std::string& filename, const ImgReaderGdal& imgSrc, const std::vector<std::string>& options=std::vector<std::string>());
  void open(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, const std::vector<std::string>& options=std::vector<std::string>());
  void close(void);

  void copyGeoTransform(const ImgReaderGdal& imgSrc);
  void setProjection(const std::string& projection);
  std::string setProjectionProj4(const std::string& projection);

  void setImageDescription(const std::string& imageDescription){m_gds->SetMetadataItem( "TIFFTAG_IMAGEDESCRIPTION",imageDescription.c_str());};
  void setGeoTransform(double* gt);

  template<typename T> bool writeData(T& value, const GDALDataType& dataType, int col, int row, int band=0) const;
  template<typename T> bool writeData(std::vector<T>& buffer, const GDALDataType& dataType , int minCol, int maxCol, int row, int band=0) const;
  template<typename T> bool writeData(std::vector<T>& buffer, const GDALDataType& dataType, int row, int band=0) const;
  bool writeData(void* pdata, const GDALDataType& dataType, int band=0) const;
  template<typename T> bool writeDataBlock(Vector2d<T>& buffer, const GDALDataType& dataType , int minCol, int maxCol, int minRow, int maxRow, int band=0) const;
  // std::string getInterleave(){return m_interleave;};
  // std::string getCompression(){return m_compression;};
  void setColorTable(const std::string& filename, int band=0);
  void setColorTable(GDALColorTable* colorTable, int band=0);
  void setMetadata(char** metadata);
  void rasterizeOgr(ImgReaderOgr& ogrReader, const std::vector<double>& burnValues=std::vector<double>(), const std::vector<std::string>& layernames=std::vector<std::string>());

protected:
  void setCodec(const GDALDataType& dataType, const std::string& imageType);
  void setCodec(const ImgReaderGdal& ImgSrc);

  /* double m_ulx; */
  /* double m_uly; */
  /* double m_delta_x; */
  /* double m_delta_y; */
  /* bool m_isGeoRef; */
  // std::string m_interleave;
  // std::string m_compression;
  std::vector<std::string> m_options;
};

template<typename T> bool ImgWriterGdal::writeData(T& value, const GDALDataType& dataType, int col, int row, int band) const
{
  //fetch raster band
  GDALRasterBand  *poBand;
  if(band>=nrOfBand()+1){
    std::ostringstream s;
    s << "band (" << band << ") exceeds nrOfBand (" << nrOfBand() << ")";
    throw(s.str());
  }
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  if(col>=nrOfCol()){
    std::ostringstream s;
    s << "col (" << col << ") exceeds nrOfCol (" << nrOfCol() << ")";
    throw(s.str());
  }
  if(col<0){
    std::ostringstream s;
    s << "col (" << col << ") is negative";
    throw(s.str());
  }
  if(row>=nrOfRow()){
    std::ostringstream s;
    s << "row (" << row << ") exceeds nrOfRow (" << nrOfRow() << ")";
    throw(s.str());
  }
  if(row<0){
    std::ostringstream s;
    s << "row (" << row << ") is negative";
    throw(s.str());
  }
  poBand->RasterIO(GF_Write,col,row,1,1,&value,1,1,dataType,0,0);
  return true;
}

template<typename T> bool ImgWriterGdal::writeData(std::vector<T>& buffer, const GDALDataType& dataType, int minCol, int maxCol, int row, int band) const
{
  //fetch raster band
  GDALRasterBand  *poBand;
  if(band>=nrOfBand()+1){
    std::ostringstream s;
    s << "band (" << band << ") exceeds nrOfBand (" << nrOfBand() << ")";
    throw(s.str());
  }
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  if(buffer.size()!=maxCol-minCol+1){
    std::string errorstring="invalid buffer size";
    throw(errorstring);
  }
  if(minCol>=nrOfCol()){
    std::ostringstream s;
    s << "minCol (" << minCol << ") exceeds nrOfCol (" << nrOfCol() << ")";
    throw(s.str());
  }
  if(minCol<0){
    std::ostringstream s;
    s << "mincol (" << minCol << ") is negative";
    throw(s.str());
  }
  if(maxCol>=nrOfCol()){
    std::ostringstream s;
    s << "maxCol (" << maxCol << ") exceeds nrOfCol (" << nrOfCol() << ")";
    throw(s.str());
  }
  if(maxCol<minCol){
    std::ostringstream s;
    s << "maxCol (" << maxCol << ") is less than minCol (" << minCol << ")";
    throw(s.str());
  }

  if(row>=nrOfRow()){
    std::ostringstream s;
    s << "row (" << row << ") exceeds nrOfRow (" << nrOfRow() << ")";
    throw(s.str());
  }
  if(row<0){
    std::ostringstream s;
    s << "row (" << row << ") is negative";
    throw(s.str());
  }
  poBand->RasterIO(GF_Write,minCol,row,buffer.size(),1,&(buffer[0]),buffer.size(),1,dataType,0,0);
  return true;
}

template<typename T> bool ImgWriterGdal::writeData(std::vector<T>& buffer, const GDALDataType& dataType, int row, int band) const
{
  //fetch raster band
  GDALRasterBand  *poBand;
  if(band>=nrOfBand()+1){
    std::ostringstream s;
    s << "band (" << band << ") exceeds nrOfBand (" << nrOfBand() << ")";
    throw(s.str());
  }
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  if(buffer.size()!=nrOfCol()){
    std::string errorstring="invalid buffer size";
    throw(errorstring);
  }
  if(row>=nrOfRow()){
    std::ostringstream s;
    s << "row (" << row << ") exceeds nrOfRow (" << nrOfRow() << ")";
    throw(s.str());
  }
  poBand->RasterIO(GF_Write,0,row,buffer.size(),1,&(buffer[0]),buffer.size(),1,dataType,0,0);
  return true;
}

template<typename T> bool ImgWriterGdal::writeDataBlock(Vector2d<T>& buffer, const GDALDataType& dataType , int minCol, int maxCol, int minRow, int maxRow, int band) const
{
  //fetch raster band
  GDALRasterBand  *poBand;
  if(band>=nrOfBand()+1){
    std::ostringstream s;
    s << "band (" << band << ") exceeds nrOfBand (" << nrOfBand() << ")";
    throw(s.str());
  }
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  assert(buffer.size()==maxRow-minRow+1);
  for(int irow=minRow;irow<=maxRow;++irow)
    writeData(buffer[irow-minRow], dataType, minCol, maxCol, irow, band);
  return true;
}

#endif // _IMGWRITERGDAL_H_


//     adfGeoTransform[0] /* top left x */
//     adfGeoTransform[1] /* w-e pixel resolution */
//     adfGeoTransform[2] /* rotation, 0 if image is "north up" */
//     adfGeoTransform[3] /* top left y */
//     adfGeoTransform[4] /* rotation, 0 if image is "north up" */
//     adfGeoTransform[5] /* n-s pixel resolution */
