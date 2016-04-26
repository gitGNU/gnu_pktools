/**********************************************************************
ImgRasterGdal.h: class to read raster files using GDAL API library
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
#ifndef _IMGRASTERGDAL_H_
#define _IMGRASTERGDAL_H_

#include <typeinfo>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <assert.h>
#include "gdal_priv.h"

enum RESAMPLE { NEAR = 0, BILINEAR = 1, BICUBIC = 2 };

template<typename T1> GDALDataType getGDALDataType(){
  if (typeid(T1) == typeid(char))
    return GDT_Byte;
  else if (typeid(T1) == typeid(unsigned char))
    return GDT_Byte;
  else if (typeid(T1) == typeid(unsigned short))
    return GDT_UInt16;
  else if (typeid(T1) == typeid(short))
    return GDT_Int16;
  else if (typeid(T1) == typeid(int))
    return GDT_Int32;
  else if (typeid(T1) == typeid(unsigned int))
    return GDT_UInt32;
  else if (typeid(T1) == typeid(long))
    return GDT_Int32;
  else if (typeid(T1) == typeid(unsigned long))
    return GDT_UInt32;
  else if (typeid(T1) == typeid(float))
    return GDT_Float32;
  else if (typeid(T1) == typeid(double))
    return GDT_Float64;
  else
    return GDT_Byte;
};
//--------------------------------------------------------------------------
class ImgRasterGdal
{
public:
  ImgRasterGdal(void);
  virtual ~ImgRasterGdal(void){};
  virtual void close(void);
  std::string getFileName() const {return m_filename;};
  int nrOfCol(void) const { return m_ncol;};
  int nrOfRow(void) const { return m_nrow;};
  int nrOfBand(void) const { return m_nband;};
  bool isGeoRef() const {double gt[6];getGeoTransform(gt);if(gt[5]<0) return true;else return false;};
  std::string getProjection(void) const;
  std::string getProjectionRef(void) const;
  std::string getGeoTransform() const;
  void getGeoTransform(double* gt) const;
  void setGeoTransform(double* gt);
  void setProjection(const std::string& projection);
  std::string setProjectionProj4(const std::string& projection);

  bool getBoundingBox (double& ulx, double& uly, double& lrx, double& lry) const;
  bool getCenterPos(double& x, double& y) const;
  double getUlx() const {double ulx, uly, lrx,lry;getBoundingBox(ulx,uly,lrx,lry);return(ulx);};
  double getUly() const {double ulx, uly, lrx,lry;getBoundingBox(ulx,uly,lrx,lry);return(uly);};
  double getLrx() const {double ulx, uly, lrx,lry;getBoundingBox(ulx,uly,lrx,lry);return(lrx);};
  double getLry() const {double ulx, uly, lrx,lry;getBoundingBox(ulx,uly,lrx,lry);return(lry);};

  int getNoDataValues(std::vector<double>& noDataValues) const;
  bool isNoData(double value) const{if(m_noDataValues.empty()) return false;else return find(m_noDataValues.begin(),m_noDataValues.end(),value)!=m_noDataValues.end();};
  int pushNoDataValue(double noDataValue);
  int setNoData(const std::vector<double> nodata){m_noDataValues=nodata; return(m_noDataValues.size());};
  CPLErr GDALSetNoDataValue(double noDataValue, int band=0) {return getRasterBand(band)->SetNoDataValue(noDataValue);};
  bool covers(double x, double y) const;
  bool covers(double ulx, double  uly, double lrx, double lry) const;
  bool geo2image(double x, double y, double& i, double& j) const;
  bool image2geo(double i, double j, double& x, double& y) const;
  double getDeltaX(void) const {double gt[6];getGeoTransform(gt);return gt[1];};
  double getDeltaY(void) const {double gt[6];getGeoTransform(gt);return -gt[5];};

  GDALDataType getDataType(int band=0) const;
  GDALRasterBand* getRasterBand(int band=0);
  GDALColorTable* getColorTable(int band=0) const;
  std::string getDriverDescription() const;
  std::string getImageType() const{return getDriverDescription();};
//   std::string getImageType() const{return "GTiff";};
  std::string getInterleave() const;
  std::string getCompression() const;
  GDALDataset* getDataset(){return m_gds;};
  char** getMetadata();
  char** getMetadata() const;
  void getMetadata(std::list<std::string>& metadata) const;

  std::string getDescription() const;
  std::string getMetadataItem() const;
  std::string getImageDescription() const;

  friend class ImgReaderGdal;
  friend class ImgWriterGdal;

protected:
  std::string m_filename;
  GDALDataset *m_gds;
  int m_ncol;
  int m_nrow;
  int m_nband;
  double m_gt[6];
  std::vector<double> m_noDataValues;

private:
};

#endif // _IMGRASTERGDAL_H_
