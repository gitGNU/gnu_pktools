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
#include "ImgReaderGdal.h"
#include "ImgWriterGdal.h"
#include "ImgReaderOgr.h"

using namespace std;

//--------------------------------------------------------------------------
class ImgWriterOgr
{
public:
  ImgWriterOgr(void);
  ImgWriterOgr(const string& filename);
  ImgWriterOgr(const string& filename, ImgReaderOgr& imgReaderOgr);
  ImgWriterOgr(const string& filename, ImgReaderOgr& imgReaderOgr, bool copyData);
  ~ImgWriterOgr(void);
  void open(const string& filename, const string& imageType="ESRI Shapefile");
  void close(void);
  int ascii2shape(const string& filename, const string &layername, const vector<string>& fieldName, const vector<OGRFieldType>& fieldType, short colX=1, short colY=2, const string& theProjection="", const OGRwkbGeometryType& eGType=wkbPoint);
  OGRLayer* createLayer(const string& layername="New layer", const string& theProjection="", const OGRwkbGeometryType& eGType=wkbUnknown, char** papszOptions=NULL);
  OGRLayer* copyLayer(OGRLayer* poSrcLayer, const string& layername, char** papszOptions=NULL);
  void createField(const string& fieldname, const OGRFieldType& fieldType, int theLayer=-1);//default: get back layer
  OGRLayer* getLayer(int layer=0) const {return m_datasource->GetLayer(layer);};
  string getLayerName(int layer=0){return m_datasource->GetLayer(layer)->GetLayerDefn()->GetName();};
  int getFields(vector<string>& fields, int layer=0) const;
  int getFields(vector<OGRFieldDefn*>& fields, int layer=0) const;
  void copyFields(const ImgReaderOgr& imgReaderOgr, int theLayer=-1);//default: get back layer
  void addLineString(vector<OGRPoint*>& points, const string& fieldName, const string& theId);
  void addRing(vector<OGRPoint*>& points, const string& fieldName, int theId);
  void addLineString(vector<OGRPoint*>& points, const string& fieldName, int theId);
  void addPoint(double x, double y, const map<string,double>& pointAttributes, string fieldName, const string& theId);
  void addPoint(double x, double y, const map<string,double>& pointAttributes, string fieldName, int theId);
  int addData(const ImgReaderGdal& imgReader, int layer=-1, bool verbose=false);
  OGRFeature* createFeature();
  OGRErr createFeature(OGRFeature* theFeature);
  int getFieldCount(int layer=0) const;
  int getFeatureCount() const;
  OGRDataSource* getDataSource(void) {return m_datasource;};
  OGRSFDriver* getDriver(void) const {return m_datasource->GetDriver();};

protected:
  void setCodec(const string& imageType);
  void setCodec(OGRSFDriver *poDriver);
  
//   OGRLayer* getLayer(int layer=-1);
    
  string m_filename;
  OGRDataSource *m_datasource;
//   vector<OGRLayer*> m_layers;
};

#endif // _IMGWRITEROGR_H_
