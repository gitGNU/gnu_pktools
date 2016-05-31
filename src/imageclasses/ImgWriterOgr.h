/**********************************************************************
ImgWriterOgr.h: class to write vector files using OGR API library
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
#ifndef _IMGWRITEROGR_H_
#define _IMGWRITEROGR_H_

#include <assert.h>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "ogrsf_frmts.h"
#include "ImgRaster.h"
#include "ImgReaderOgr.h"

//--------------------------------------------------------------------------
class ImgWriterOgr
{
public:
  ImgWriterOgr(void);
  ImgWriterOgr(const std::string& filename, const std::string& imageType="ESRI Shapefile");
  ImgWriterOgr(const std::string& filename, ImgReaderOgr& imgReaderOgr);
  ImgWriterOgr(const std::string& filename, ImgReaderOgr& imgReaderOgr, bool copyData);
  ~ImgWriterOgr(void);
  void open(const std::string& filename, ImgReaderOgr& imgReaderOgr);
  void open(const std::string& filename, const std::string& imageType="ESRI Shapefile");
  void close(void);
  int ascii2ogr(const std::string& filename, const std::string &layername, const std::vector<std::string>& fieldName, const std::vector<OGRFieldType>& fieldType, short colX=1, short colY=2, const std::string& theProjection="", const OGRwkbGeometryType& eGType=wkbPoint, const char fs=' ');
  OGRLayer* createLayer(const std::string& layername="New layer", const std::string& theProjection="", const OGRwkbGeometryType& eGType=wkbUnknown, char** papszOptions=NULL);
  OGRLayer* copyLayer(OGRLayer* poSrcLayer, const std::string& layername, char** papszOptions=NULL);
  void createField(const std::string& fieldname, const OGRFieldType& fieldType, int theLayer=0);//default: get back layer
  OGRLayer* getLayer(int layer=0) const {return m_datasource->GetLayer(layer);};
  std::string getLayerName(int layer=0){return m_datasource->GetLayer(layer)->GetLayerDefn()->GetName();};
  int getFields(std::vector<std::string>& fields, int layer=0) const;
  int getFields(std::vector<OGRFieldDefn*>& fields, int layer=0) const;
  void copyFields(const ImgReaderOgr& imgReaderOgr, int srcLayer=0, int targetLayer=0);//default: get back layer
  void addLineString(std::vector<OGRPoint*>& points, const std::string& fieldName, const std::string& theId, int layer=0);
  void addRing(std::vector<OGRPoint*>& points, const std::string& fieldName, int theId, int layer=0);
  void addLineString(std::vector<OGRPoint*>& points, const std::string& fieldName, int theId, int layer=0);
  void addPoint(double x, double y, const std::map<std::string,double>& pointAttributes, std::string fieldName, const std::string& theId, int layer=0);
  void addPoint(double x, double y, const std::map<std::string,double>& pointAttributes, std::string fieldName, int theId, int layer=0);
  int addData(ImgRaster& imgReader, int layer=0, bool verbose=false);
  OGRFeature* createFeature(int layer=0);
  OGRErr createFeature(OGRFeature* theFeature, int layer=0);
  int getFieldCount(int layer=0) const;
  int getFeatureCount(int layer=0) const;
#if GDAL_VERSION_MAJOR < 2
  OGRDataSource* getDataSource(void) {return m_datasource;};
  OGRSFDriver* getDriver(void) const {return m_datasource->GetDriver();};
#else
  GDALDataset* getDataSource(void) {return m_datasource;};
  GDALDriver* getDriver(void) const {return m_datasource->GetDriver();};
#endif

protected:
  void setCodec(const std::string& imageType);
#if GDAL_VERSION_MAJOR < 2
  void setCodec(OGRSFDriver *poDriver);
#else
  void setCodec(GDALDriver *poDriver);
#endif
  
//   OGRLayer* getLayer(int layer=0);
    
  std::string m_filename;
#if GDAL_VERSION_MAJOR < 2
  OGRDataSource *m_datasource;
#else
  GDALDataset *m_datasource;
#endif
//   vector<OGRLayer*> m_layers;
};

#endif // _IMGWRITEROGR_H_
