/**********************************************************************
ImgWriterGdal.h: class to write raster files using GDAL API library
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
#ifndef _IMGWRITERGDAL_H_
#define _IMGWRITERGDAL_H_

#include <assert.h>
#include <fstream>
#include <string>
#include <sstream>
#include "gdal_priv.h"
#include "ImgReaderGdal.h"


using namespace std;

//--------------------------------------------------------------------------
class ImgWriterGdal
{
public:
  ImgWriterGdal(void);
  ~ImgWriterGdal(void);
  void open(const string& filename);
  void open(const string& filename, const ImgReaderGdal& imgSrc, const vector<string>& options=vector<string>());
  // void open(const string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const string& imageType="GTiff", const string& interleave="BAND", const string& compression="LZW", int magicX=1, int magicY=1);
  void open(const string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const string& imageType, const vector<string>& options=vector<string>());
  void close(void);
  int nrOfCol(void) const { return m_ncol;};
  int nrOfRow(void) const { return m_nrow;};
  int nrOfBand(void) const { return m_nband;};
  void setGeoTransform(double ulx, double uly, double deltaX, double deltaY, double rot1=0, double rot2=0);
  void copyGeoTransform(const ImgReaderGdal& imgSrc);
  string setProjection(void);//set (and return) default projection ETSR-LAEA
  void setProjection(const string& projection);
  string setProjectionProj4(const string& projection);
  void setImageDescription(const string& imageDescription){m_gds->SetMetadataItem( "TIFFTAG_IMAGEDESCRIPTION",imageDescription.c_str());};
  string getProjection(void) const;
  string getGeoTransform() const;
  void getGeoTransform(double& ulx, double& uly, double& deltaX, double& deltaY, double& rot1, double& rot2) const;
  bool getBoundingBox(double& ulx, double& uly, double& lrx, double& lry) const;
  bool getCentrePos(double& x, double& y) const;
  bool covers(double x, double y) const;
  bool covers(double ulx, double  uly, double lrx, double lry) const;
  bool geo2image(double x, double y, double& i, double& j) const;
  bool image2geo(double i, double j, double& x, double& y) const;
  bool isGeoRef() const {return m_isGeoRef;};
  // void getMagicPixel(double& x, double& y) const {x=m_magic_x;y=m_magic_y;};
  double getDeltaX(void) const {return m_delta_x;};
  double getDeltaY(void) const{return m_delta_y;};
  template<typename T> bool writeData(T& value, const GDALDataType& dataType, int col, int row, int band=0) const;
  template<typename T> bool writeData(vector<T>& buffer, const GDALDataType& dataType , int minCol, int maxCol, int row, int band=0) const;
  template<typename T> bool writeData(vector<T>& buffer, const GDALDataType& dataType, int row, int band=0) const;
  bool writeData(void* pdata, const GDALDataType& dataType, int band=0) const;
  // string getInterleave(){return m_interleave;};
  // string getCompression(){return m_compression;};
  GDALDataType getDataType(int band=0) const;
  GDALRasterBand* getRasterBand(int band);
  void setColorTable(const string& filename, int band=0);
  void setColorTable(GDALColorTable* colorTable, int band=0);
  void setMetadata(char** metadata);

protected:
  void setCodec(const string& imageType);
  void setCodec(const ImgReaderGdal& ImgSrc);

  string m_filename;
  GDALDataset *m_gds;
  int m_ncol;
  int m_nrow;
  int m_nband;
  GDALDataType m_type;
  double m_ulx;
  double m_uly;
  double m_delta_x;
  double m_delta_y;
  // double m_magic_x;
  // double m_magic_y;
  bool m_isGeoRef;
  // string m_interleave;
  // string m_compression;
  vector<string> m_options;
};

template<typename T> bool ImgWriterGdal::writeData(T& value, const GDALDataType& dataType, int col, int row, int band) const
{
  //fetch raster band
  GDALRasterBand  *poBand;
  if(band>=nrOfBand()+1){
    ostringstream s;
    s << "band (" << band << ") exceeds nrOfBand (" << nrOfBand() << ")";
    throw(s.str());
  }
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  if(col>=nrOfCol()){
    ostringstream s;
    s << "col (" << col << ") exceeds nrOfCol (" << nrOfCol() << ")";
    throw(s.str());
  }
  if(col<0){
    ostringstream s;
    s << "col (" << col << ") is negative";
    throw(s.str());
  }
  if(row>=nrOfRow()){
    ostringstream s;
    s << "row (" << row << ") exceeds nrOfRow (" << nrOfRow() << ")";
    throw(s.str());
  }
  if(row<0){
    ostringstream s;
    s << "row (" << row << ") is negative";
    throw(s.str());
  }
  poBand->RasterIO(GF_Write,col,row,1,1,&value,1,1,dataType,0,0);
}

template<typename T> bool ImgWriterGdal::writeData(vector<T>& buffer, const GDALDataType& dataType, int minCol, int maxCol, int row, int band) const
{
  //fetch raster band
  GDALRasterBand  *poBand;
  if(band>=nrOfBand()+1){
    ostringstream s;
    s << "band (" << band << ") exceeds nrOfBand (" << nrOfBand() << ")";
    throw(s.str());
  }
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  if(buffer.size()!=maxCol-minCol+1){
    string errorstring="invalid buffer size";
    throw(errorstring);
  }
  if(minCol>=nrOfCol()){
    ostringstream s;
    s << "minCol (" << minCol << ") exceeds nrOfCol (" << nrOfCol() << ")";
    throw(s.str());
  }
  if(minCol<0){
    ostringstream s;
    s << "mincol (" << minCol << ") is negative";
    throw(s.str());
  }
  if(maxCol>=nrOfCol()){
    ostringstream s;
    s << "maxCol (" << maxCol << ") exceeds nrOfCol (" << nrOfCol() << ")";
    throw(s.str());
  }
  if(maxCol<minCol){
    ostringstream s;
    s << "maxCol (" << maxCol << ") is less than minCol (" << minCol << ")";
    throw(s.str());
  }

  if(row>=nrOfRow()){
    ostringstream s;
    s << "row (" << row << ") exceeds nrOfRow (" << nrOfRow() << ")";
    throw(s.str());
  }
  if(row<0){
    ostringstream s;
    s << "row (" << row << ") is negative";
    throw(s.str());
  }
  poBand->RasterIO(GF_Write,minCol,row,buffer.size(),1,&(buffer[0]),buffer.size(),1,dataType,0,0);
}

template<typename T> bool ImgWriterGdal::writeData(vector<T>& buffer, const GDALDataType& dataType, int row, int band) const
{
  //fetch raster band
  GDALRasterBand  *poBand;
  if(band>=nrOfBand()+1){
    ostringstream s;
    s << "band (" << band << ") exceeds nrOfBand (" << nrOfBand() << ")";
    throw(s.str());
  }
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  if(buffer.size()!=nrOfCol()){
    string errorstring="invalid buffer size";
    throw(errorstring);
  }
  if(row>=nrOfRow()){
    ostringstream s;
    s << "row (" << row << ") exceeds nrOfRow (" << nrOfRow() << ")";
    throw(s.str());
  }
  poBand->RasterIO(GF_Write,0,row,buffer.size(),1,&(buffer[0]),buffer.size(),1,dataType,0,0);
}

#endif // _IMGWRITERGDAL_H_


//     adfGeoTransform[0] /* top left x */
//     adfGeoTransform[1] /* w-e pixel resolution */
//     adfGeoTransform[2] /* rotation, 0 if image is "north up" */
//     adfGeoTransform[3] /* top left y */
//     adfGeoTransform[4] /* rotation, 0 if image is "north up" */
//     adfGeoTransform[5] /* n-s pixel resolution */
