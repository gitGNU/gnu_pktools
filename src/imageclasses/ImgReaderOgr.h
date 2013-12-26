/**********************************************************************
ImgReaderOgr.h: class to read vector files using OGR API library
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
#ifndef _IMGREADEROGR_H_
#define _IMGREADEROGR_H_

#include <assert.h>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>
#include "ogrsf_frmts.h"
#include "base/Vector2d.h"
#include "ImgReaderGdal.h"

//--------------------------------------------------------------------------
class ImgReaderOgr
{
public:
  ImgReaderOgr(void);
  ImgReaderOgr(const std::string& filename);
  ~ImgReaderOgr(void);
  void open(const std::string& filename);
  void close(void);
  //   int readData(Vector2d<double>& data, const OGRFieldType& fieldType, std::vector<std::string>& fields, int layer=0, bool pos=false);//default layer 0 and no pos information in data
  template <typename T> int readXY(std::vector<T>& xVector, std::vector<T>& yVector, int layer=0, bool verbose=false);
  template <typename T> int readY(std::vector<T>& yVector, int layer=0, bool verbose=false);
  template <typename T> int readData(std::vector<T>& data, const OGRFieldType& fieldType, std::vector<std::string>& fields, OGRFeature *poFeature, int layer=0, bool pos=false, bool verbose=false);
  template <typename T> int readData(std::vector<T>& data, const OGRFieldType& fieldType, const std::string& theField, int layer=0, bool verbose=false);
  template <typename T> int readData(Vector2d<T>& data, const OGRFieldType& fieldType, std::vector<std::string>& fields, int layer=0, bool pos=false, bool verbose=false);//default layer 0 and no pos information in data
  template <typename T> int readData(std::map<int,Vector2d<T> >& data, const OGRFieldType& fieldType, std::vector<std::string>& fields, const std::string& label, int layer=0, bool pos=false, bool verbose=false);//default layer 0 and no pos information in data
  template <typename T> int readData(std::map<std::string,Vector2d<T> >& data, const OGRFieldType& fieldType, std::vector<std::string>& fields, const std::string& label, int layer=0, bool pos=false, bool verbose=false);//default layer 0 and no pos information in data
  unsigned int readDataImageShape(std::map<std::string,Vector2d<float> > &mapPixels, //[classNr][pixelNr][bandNr],
						       std::vector<std::string>& fields,
						       const std::vector<short>& bands,
						       const std::string& label,
						       int verbose=false);

  unsigned int readDataImageShape(std::map<std::string,Vector2d<float> > &mapPixels, //[classNr][pixelNr][bandNr],
						       std::vector<std::string>& fields,
						       double start,
						       double end,
						       const std::string& label,
						       int verbose=false);

  void shape2ascii(std::ostream& theOstream, const std::string& pointname, int layer=0, bool verbose=false);
  unsigned long int getFeatureCount(int layer=0) const;
  int getFieldCount(int layer=0) const;
  OGRLayer* getLayer(int layer=0){return m_datasource->GetLayer(layer);};
  std::string getProjection(int layer=0) const;
  OGRwkbGeometryType getGeometryType(int layer=0) const;
  std::string getLayerName(int layer=0){return m_datasource->GetLayer(layer)->GetLayerDefn()->GetName();};
  //  int getLayer(int layer=0) const;
  int getFields(std::vector<std::string>& fields, int layer=0) const;
  int getFields(std::vector<OGRFieldDefn*>& fields, int layer=0) const;
  OGRDataSource* getDataSource(void) {return m_datasource;};
  OGRSFDriver* getDriver(void) const {return m_datasource->GetDriver();};
//   OGRLayer *executeSql(const std::string& output,const std::string& sqlStatement, OGRGeometry* spatialFilter=NULL);
  template<typename T> int readSql(Vector2d<T>& data, const OGRFieldType& fieldType, std::vector<std::string>& fields, const std::string& sqlStatement, OGRGeometry* spatialFilter=NULL, int layer=0, bool pos=false, bool verbose=false);
  template<typename T> int readSql(std::map<int,Vector2d<T> >& data, const OGRFieldType& fieldType, std::vector<std::string>& fields, const std::string& label, const std::string& sqlStatement, OGRGeometry* spatialFilter, int layer=0, bool pos=false, bool verbose=false);
  bool getExtent(double& ulx, double& uly, double& lrx, double& lry, int layer=0);

  friend std::ostream& operator<<(std::ostream& theOstream, ImgReaderOgr& theImageReader);
  
protected:
  void setCodec(void);

  std::string m_filename;
  OGRDataSource *m_datasource;
};

//read data from all features in a map, organized by classes
template <typename T> int ImgReaderOgr::readData(std::map<int,Vector2d<T> >& data, const OGRFieldType& fieldType, std::vector<std::string>& fields, const std::string& label, int layer, bool pos, bool verbose)
{
  if(layer<0)
    layer=m_datasource->GetLayerCount()-1;
  assert(m_datasource->GetLayerCount()>layer);
  OGRLayer  *poLayer;
  if(verbose)
    std::cout << "number of layers: " << m_datasource->GetLayerCount() << std::endl;
  poLayer = m_datasource->GetLayer(layer);
  if(poLayer!=NULL){
    OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
    if(fields.empty()){
      fields.resize(poFDefn->GetFieldCount());
      if(verbose)
        std::cout << "resized fields to " << fields.size() << std::endl;
    }
    //start reading features from the layer
    OGRFeature *poFeature;
    if(verbose)
      std::cout << "reset reading" << std::endl;
    poLayer->ResetReading();
    unsigned long int ifeature=0;
    int posOffset=(pos)?2:0;
    if(verbose)
      std::cout << "going through features" << std::endl << std::flush;
    int theClass=0;
    while( (poFeature = poLayer->GetNextFeature()) != NULL ){
      std::vector<T> theFeature;//(fields.size()+posOffset);//x,y+selectedfields
      if(verbose)
        std::cout << "reading feature " << ifeature << std::endl << std::flush;
      OGRGeometry *poGeometry;
      poGeometry = poFeature->GetGeometryRef();
      if(verbose){
        if(poGeometry == NULL)
          std::cerr << "no geometry defined" << std::endl << std::flush;
        else if(wkbFlatten(poGeometry->getGeometryType()) != wkbPoint)
          std::cerr << "Warning: poGeometry type: " << wkbFlatten(poGeometry->getGeometryType()) << std::endl << std::flush;
      }
      assert(poGeometry != NULL );
             // && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint);
      if(pos){
        if(wkbFlatten(poGeometry->getGeometryType()) == wkbPoint){
          OGRPoint *poPoint;
          poPoint = (OGRPoint *) poGeometry;
          theFeature.push_back(poPoint->getX());
          theFeature.push_back(poPoint->getY());
        }
        else if(wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon){
          OGRPoint thePoint;
          OGRPolygon * poPolygon = (OGRPolygon *) poGeometry;
          poPolygon->Centroid(&thePoint);
          theFeature.push_back(thePoint.getX());
          theFeature.push_back(thePoint.getY());
        }        
        else{
          //Centroid for non polygon geometry not supported until OGR 1.8.0, comment out if version < 1.8.0 is installed...";
          OGRPoint thePoint;
          poGeometry->Centroid(&thePoint);
          theFeature.push_back(thePoint.getX());
          theFeature.push_back(thePoint.getY());
        }       
      }
      // OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
      std::string featurename;
      for(int iField=0;iField<poFDefn->GetFieldCount();++iField){
        OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);
        std::string fieldname=poFieldDefn->GetNameRef();
        if(fieldname==label)
          theClass=poFeature->GetFieldAsInteger(iField);
        else{
          switch(fieldType){
          case(OFTReal):
            if(fields.size()<poFDefn->GetFieldCount()){
              if(find(fields.begin(),fields.end(),fieldname)!=fields.end())
                theFeature.push_back(poFeature->GetFieldAsDouble(iField));
            }
            else{
              fields[iField]=fieldname;
              theFeature.push_back(poFeature->GetFieldAsDouble(iField));
            }
            break;
          case(OFTInteger):
            if(fields.size()<poFDefn->GetFieldCount()){
              if(find(fields.begin(),fields.end(),fieldname)!=fields.end())
                theFeature.push_back(poFeature->GetFieldAsDouble(iField));
            }
            else{
              fields[iField]=fieldname;
              theFeature.push_back(poFeature->GetFieldAsDouble(iField));
            }
            break;
          default:
            {
              std::string errorstring="field type not supported in ImgReaderOgr::ReadData";
              throw(errorstring);
            }
            break;
          }
        }
      }
      data[theClass].push_back(theFeature);
      ++ifeature;
      ++ifeature;
    }
    if(verbose)
      std::cout << "number of features read: " << ifeature << std::endl << std::flush;
    typename std::map<int,Vector2d<T> >::const_iterator mit=data.begin();
    int nband=0;
    if(verbose)
      std::cout << "read classes: " << std::flush;
    while(mit!=data.end()){
      if(verbose)
        std::cout << mit->first << " " << std::flush;
      if(!nband)
        nband=fields.size();
      if(pos)
        assert((mit->second)[0].size()==nband+2);
      else
        assert((mit->second)[0].size()==nband);
      ++mit;
    }
    if(verbose)
      std::cout << std::endl << std::flush;
    return(nband);
  }
  else{
    std::ostringstream ess;
    ess << "no layer in " << m_filename;
    throw(ess.str());
  }
}

//read data from all features in a map, organized by class names
template <typename T> int ImgReaderOgr::readData(std::map<std::string,Vector2d<T> >& data, const OGRFieldType& fieldType, std::vector<std::string>& fields, const std::string& label, int layer, bool pos, bool verbose)
{
  if(layer<0)
    layer=m_datasource->GetLayerCount()-1;
  assert(m_datasource->GetLayerCount()>layer);
  OGRLayer  *poLayer;
  if(verbose)
    std::cout << "number of layers: " << m_datasource->GetLayerCount() << std::endl;
  poLayer = m_datasource->GetLayer(layer);
  if(poLayer!=NULL){
    OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
    assert(poFDefn!=NULL);

    if(fields.empty()){
      fields.resize(poFDefn->GetFieldCount());
      if(verbose)
        std::cout << "resized fields to " << fields.size() << std::endl;
    }

    //start reading features from the layer
    OGRFeature *poFeature;
    if(verbose)
      std::cout << "reset reading" << std::endl;
    poLayer->ResetReading();
    unsigned long int ifeature=0;
    int posOffset=(pos)?2:0;
    if(verbose)
      std::cout << "going through features to fill in string map" << std::endl << std::flush;
    std::string theClass;
    while( (poFeature = poLayer->GetNextFeature()) != NULL ){
      std::vector<T> theFeature;//(fields.size()+posOffset);//x,y+selectedfields
      if(verbose)
        std::cout << "reading feature " << ifeature << std::endl << std::flush;
      OGRGeometry *poGeometry;
      poGeometry = poFeature->GetGeometryRef();
      if(verbose){
        if(poGeometry == NULL)
          std::cerr << "no geometry defined" << std::endl << std::flush;
        else if(wkbFlatten(poGeometry->getGeometryType()) != wkbPoint)
          std::cerr << "Warning: poGeometry type: " << wkbFlatten(poGeometry->getGeometryType()) << std::endl << std::flush;
      }
      assert(poGeometry != NULL );
             // && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint);
      if(pos){
        if(wkbFlatten(poGeometry->getGeometryType()) == wkbPoint){
          OGRPoint *poPoint;
          poPoint = (OGRPoint *) poGeometry;
          theFeature.push_back(poPoint->getX());
          theFeature.push_back(poPoint->getY());
        }
        else if(wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon){
          OGRPoint thePoint;
          poGeometry->Centroid(&thePoint);
          theFeature.push_back(thePoint.getX());
          theFeature.push_back(thePoint.getY());
        }        
        else{
          //Centroid for non polygon geometry not supported until OGR 1.8.0, comment out if version < 1.8.0 is installed...";
          OGRPoint thePoint;
          poGeometry->Centroid(&thePoint);
          theFeature.push_back(thePoint.getX());
          theFeature.push_back(thePoint.getY());
        }       
      }
      // OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();//got LayerDefn already...
      std::string featurename;
      for(int iField=0;iField<poFDefn->GetFieldCount();++iField){
        OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);
        std::string fieldname=poFieldDefn->GetNameRef();
        if(fieldname==label){
          theClass=poFeature->GetFieldAsString(iField);
          if(verbose)
            std::cout << "read feature for " << theClass << std::endl;
        }
        else{
          switch(fieldType){
          case(OFTReal):
            if(fields.size()<poFDefn->GetFieldCount()){
              if(find(fields.begin(),fields.end(),fieldname)!=fields.end())
                theFeature.push_back(poFeature->GetFieldAsDouble(iField));
            }
            else{
              fields[iField]=fieldname;
              theFeature.push_back(poFeature->GetFieldAsDouble(iField));
            }
            break;
          case(OFTInteger):
            if(fields.size()<poFDefn->GetFieldCount()){
              if(find(fields.begin(),fields.end(),fieldname)!=fields.end())
                theFeature.push_back(poFeature->GetFieldAsDouble(iField));
            }
            else{
              fields[iField]=fieldname;
              theFeature.push_back(poFeature->GetFieldAsDouble(iField));
            }
            break;
          default:
            {
              std::string errorstring="field type not supported in ImgReaderOgr::ReadData";
              throw(errorstring);
            }
            break;
          }
        }
        assert(poFDefn!=NULL);
      }
      data[theClass].push_back(theFeature);
      ++ifeature;
    }
    if(verbose)
      std::cout << "number of features read: " << ifeature << std::endl << std::flush;
    typename std::map<std::string,Vector2d<T> >::const_iterator mit=data.begin();
    int nband=0;
    if(verbose)
      std::cout << "read classes: " << std::flush;
    while(mit!=data.end()){
      if(verbose)
        std::cout << mit->first << " " << std::flush;
      if(!nband)
        nband=fields.size();
      if(pos)
        assert((mit->second)[0].size()==nband+2);
      else
        assert((mit->second)[0].size()==nband);
      ++mit;
    }
    if(verbose)
      std::cout << std::endl << std::flush;
    return(nband);
  }
  else{
    std::ostringstream ess;
    ess << "no layer in " << m_filename;
    throw(ess.str());
  }
}

//read x positions
template <typename T> int ImgReaderOgr::readXY(std::vector<T>& xVector, std::vector<T>& yVector, int layer, bool verbose){
  if(layer<0)
    layer=m_datasource->GetLayerCount()-1;
  assert(m_datasource->GetLayerCount()>layer);
  OGRLayer  *poLayer;
  if(verbose)
    std::cout << "number of layers: " << m_datasource->GetLayerCount() << std::endl;
  poLayer = m_datasource->GetLayer(layer);
  OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
  //start reading features from the layer
  OGRFeature *poFeature;
  if(verbose)
    std::cout << "reset reading" << std::endl;
  poLayer->ResetReading();
  unsigned long int ifeature=0;
  if(verbose)
    std::cout << "going through features" << std::endl << std::flush;
  while( (poFeature = poLayer->GetNextFeature()) != NULL ){
    if(verbose)
      std::cout << "reading feature " << ifeature << std::endl << std::flush;
    OGRGeometry *poGeometry;
    poGeometry = poFeature->GetGeometryRef();
    if(verbose){
      if(poGeometry == NULL)
        std::cerr << "no geometry defined" << std::endl << std::flush;
      else// if(wkbFlatten(poGeometry->getGeometryType()) != wkbPoint)
        std::cout << "poGeometry type: " << wkbFlatten(poGeometry->getGeometryType()) << std::endl;
    }
    // assert(poGeometry != NULL 
    //        && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint);
    OGRPoint *poPoint = (OGRPoint *) poGeometry;
    xVector.push_back(poPoint->getX());
    yVector.push_back(poPoint->getY());
    ++ifeature;
  }
  assert(xVector.size()==yVector.size());
  if(xVector.size()){
    return xVector.size();
  }
  else{
    std::ostringstream ess;
    ess << "no layer in " << m_filename;
    throw(ess.str());
  }
}

//read data from a single feature
template <typename T> int ImgReaderOgr::readData(std::vector<T>& data, const OGRFieldType& fieldType, std::vector<std::string>& fields, OGRFeature *poFeature, int layer, bool pos, bool verbose)
{
  if(layer<0)
    layer=m_datasource->GetLayerCount()-1;
  assert(m_datasource->GetLayerCount()>layer);
  OGRLayer  *poLayer;
  if(verbose)
    std::cout << "number of layers: " << m_datasource->GetLayerCount() << std::endl;
  poLayer = m_datasource->GetLayer(layer);
  OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
  if(fields.empty()){
    fields.resize(poFDefn->GetFieldCount());
    if(verbose)
      std::cout << "resized fields to " << fields.size() << std::endl;
  }
  OGRGeometry *poGeometry;
  poGeometry = poFeature->GetGeometryRef();
  if(verbose){
    if(poGeometry == NULL)
      std::cerr << "no geometry defined" << std::endl << std::flush;
    else if(wkbFlatten(poGeometry->getGeometryType()) != wkbPoint)
      std::cerr << "poGeometry type: " << wkbFlatten(poGeometry->getGeometryType()) << std::endl << std::flush;
  }
  assert(poGeometry != NULL);
         // && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint);
  OGRPoint *poPoint = (OGRPoint *) poGeometry;
  if(pos){
    if(wkbFlatten(poGeometry->getGeometryType()) == wkbPoint){
      OGRPoint *poPoint;
      poPoint = (OGRPoint *) poGeometry;
      data.push_back(poPoint->getX());
      data.push_back(poPoint->getY());
    }
    else if(wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon){
      OGRPoint thePoint;
      poGeometry->Centroid(&thePoint);
      data.push_back(thePoint.getX());
      data.push_back(thePoint.getY());
    }        
    else{
      //Centroid for non polygon geometry not supported until OGR 1.8.0, comment out if version < 1.8.0 is installed...";
      OGRPoint thePoint;
      poGeometry->Centroid(&thePoint);
      data.push_back(thePoint.getX());
      data.push_back(thePoint.getY());
    }       

  }
  std::string featurename;
  for(int iField=0;iField<poFDefn->GetFieldCount();++iField){
    OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);
    std::string fieldname=poFieldDefn->GetNameRef();
    switch(fieldType){
    case(OFTReal):
      if(fields.size()<poFDefn->GetFieldCount()){
        if(find(fields.begin(),fields.end(),fieldname)!=fields.end())
          data.push_back(poFeature->GetFieldAsDouble(iField));
      }
      else{
        fields[iField]=fieldname;
        data.push_back(poFeature->GetFieldAsDouble(iField));
      }
      break;
    case(OFTInteger):
      if(fields.size()<poFDefn->GetFieldCount()){
        if(find(fields.begin(),fields.end(),fieldname)!=fields.end())
          data.push_back(poFeature->GetFieldAsDouble(iField));
      }
      else{
        fields[iField]=fieldname;
        data.push_back(poFeature->GetFieldAsDouble(iField));
      }
      break;
    default:
      {
        std::string errorstring="field type not supported in ImgReaderOgr::ReadData";
        throw(errorstring);
      }
      break;
    }
  }
  //   assert(data.size()==ifeature);
  if(data.size()){
    if(pos)
      assert(data.size()==fields.size()+2);
    else
      assert(data.size()==fields.size());
    return fields.size();
  }
  else{
    std::ostringstream ess;
    ess << "no layer in " << m_filename;
    throw(ess.str());
  }
}

//read one field from all features
template <typename T> int ImgReaderOgr::readData(std::vector<T>& data, const OGRFieldType& fieldType, const std::string& theField, int layer, bool verbose)
{
  if(layer<0)
    layer=m_datasource->GetLayerCount()-1;
  assert(m_datasource->GetLayerCount()>layer);
  OGRLayer  *poLayer;
  if(verbose)
    std::cout << "number of layers: " << m_datasource->GetLayerCount() << std::endl;
  poLayer = m_datasource->GetLayer(layer);
  OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
  int nfield=(theField!="")? poFDefn->GetFieldCount() : 1;
  if(theField==""){
    //read first field available 
    if(verbose)
      std::cout << "read first field from total of " << nfield << std::endl;
  }

  //start reading features from the layer
  OGRFeature *poFeature;
  if(verbose)
    std::cout << "reset reading" << std::endl;
  poLayer->ResetReading();
  unsigned long int ifeature=0;
  if(verbose)
    std::cout << "going through features" << std::endl << std::flush;
  while( (poFeature = poLayer->GetNextFeature()) != NULL ){
    // std::vector<T> theFeature;//(fields.size()+posOffset);//x,y+selectedfields
    T theFeature;
    if(verbose)
      std::cout << "reading feature " << ifeature << std::endl << std::flush;
    OGRGeometry *poGeometry;
    poGeometry = poFeature->GetGeometryRef();
    if(verbose){
      if(poGeometry == NULL)
        std::cerr << "no geometry defined" << std::endl << std::flush;
      else// if(wkbFlatten(poGeometry->getGeometryType()) != wkbPoint)
        std::cout << "poGeometry type: " << wkbFlatten(poGeometry->getGeometryType()) << std::endl;
    }
    // assert(poGeometry != NULL 
    //        && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint);
    OGRPoint *poPoint = (OGRPoint *) poGeometry;

    for(int iField=0;iField<nfield;++iField){
      OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);
      std::string fieldname=poFieldDefn->GetNameRef();
      if(fieldname!=theField)
        continue;
      switch(fieldType){
      case(OFTInteger):
      case(OFTReal):
        theFeature=poFeature->GetFieldAsDouble(iField);
      break;
      default:
        {
          std::string errorstring="field type not supported in ImgReaderOgr::ReadData";
          throw(errorstring);
        }
        break;
      }
    }
    data.push_back(theFeature);
    if(verbose)
      std::cout << "feature is: " << theFeature << std::endl;
    ++ifeature;
  }
  if(data.size()){
    return data.size();
  }
  else{
    std::ostringstream ess;
    ess << "no layer in " << m_filename;
    throw(ess.str());
  }
}

//read data from all features  
template <typename T> int ImgReaderOgr::readData(Vector2d<T>& data, const OGRFieldType& fieldType, std::vector<std::string>& fields, int layer, bool pos, bool verbose)
{
  if(layer<0)
    layer=m_datasource->GetLayerCount()-1;
  assert(m_datasource->GetLayerCount()>layer);
  OGRLayer  *poLayer;
  if(verbose)
    std::cout << "number of layers: " << m_datasource->GetLayerCount() << std::endl;
  poLayer = m_datasource->GetLayer(layer);
  OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
  if(fields.empty()){
    fields.resize(poFDefn->GetFieldCount());
    if(verbose)
      std::cout << "resized fields to " << fields.size() << std::endl;
  }
  //start reading features from the layer
  OGRFeature *poFeature;
  if(verbose)
    std::cout << "reset reading" << std::endl;
  poLayer->ResetReading();
  unsigned long int ifeature=0;
  int posOffset=(pos)?2:0;
  if(verbose)
    std::cout << "going through features" << std::endl << std::flush;
  while( (poFeature = poLayer->GetNextFeature()) != NULL ){
    std::vector<T> theFeature;//(fields.size()+posOffset);//x,y+selectedfields
    if(verbose)
      std::cout << "reading feature " << ifeature << std::endl << std::flush;
    OGRGeometry *poGeometry;
    poGeometry = poFeature->GetGeometryRef();
    if(verbose){
      if(poGeometry == NULL)
        std::cerr << "no geometry defined" << std::endl << std::flush;
      else if(wkbFlatten(poGeometry->getGeometryType()) != wkbPoint)
        std::cerr << "poGeometry type: " << wkbFlatten(poGeometry->getGeometryType()) << std::endl << std::flush;
    }
    assert(poGeometry != NULL 
           && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint);
    OGRPoint *poPoint = (OGRPoint *) poGeometry;
    if(pos){
      theFeature.push_back(poPoint->getX());
      theFeature.push_back(poPoint->getY());
    }
    OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
    std::string featurename;
    for(int iField=0;iField<poFDefn->GetFieldCount();++iField){
      OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);
      std::string fieldname=poFieldDefn->GetNameRef();
      switch(fieldType){
      case(OFTReal):
        if(fields.size()<poFDefn->GetFieldCount()){
          if(find(fields.begin(),fields.end(),fieldname)!=fields.end())
            theFeature.push_back(poFeature->GetFieldAsDouble(iField));
        }
        else{
          fields[iField]=fieldname;
          theFeature.push_back(poFeature->GetFieldAsDouble(iField));
        }
        break;
      case(OFTInteger):
        if(fields.size()<poFDefn->GetFieldCount()){
          if(find(fields.begin(),fields.end(),fieldname)!=fields.end())
            theFeature.push_back(poFeature->GetFieldAsDouble(iField));
        }
        else{
          fields[iField]=fieldname;
          theFeature.push_back(poFeature->GetFieldAsDouble(iField));
        }
        break;
      default:
        {
          std::string errorstring="field type not supported in ImgReaderOgr::ReadData";
          throw(errorstring);
        }
        break;
      }
    }
    data.push_back(theFeature);
    ++ifeature;
  }
//   assert(data.size()==ifeature);
  if(data.size()){
    if(pos)
      assert(data[0].size()==fields.size()+2);
    else
      assert(data[0].size()==fields.size());
    return fields.size();
  }
  else{
    std::ostringstream ess;
    ess << "no layer in " << m_filename;
    throw(ess.str());
  }
}

template<typename T> int ImgReaderOgr::readSql(std::map<int, Vector2d<T> >& data, const OGRFieldType& fieldType, std::vector<std::string>& fields, const std::string& label, const std::string& sqlStatement, OGRGeometry* spatialFilter, int layer, bool pos, bool verbose)
{
  if(layer<0)
    layer=m_datasource->GetLayerCount()-1;
  assert(m_datasource->GetLayerCount()>layer);
  OGRLayer *poLayer;
  poLayer = m_datasource->ExecuteSQL(sqlStatement.c_str(), spatialFilter,NULL );
  if(poLayer!=NULL){
    OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
    if(fields.empty()){
      fields.resize(poFDefn->GetFieldCount());
      if(verbose)
        std::cout << "resized fields to " << fields.size() << std::endl;
    }
    //start reading features from the layer
    OGRFeature *poFeature;
    if(verbose)
      std::cout << "reset reading" << std::endl;
    poLayer->ResetReading();
    unsigned long int ifeature=0;
    int posOffset=(pos)?2:0;
    if(verbose)
      std::cout << "going through features" << std::endl << std::flush;
    int theClass=0;
    while( (poFeature = poLayer->GetNextFeature()) != NULL ){
      std::vector<T> theFeature;//(fields.size()+posOffset);//x,y+selectedfields
      if(verbose)
        std::cout << "reading feature " << ifeature << std::endl << std::flush;
      OGRGeometry *poGeometry;
      poGeometry = poFeature->GetGeometryRef();
      if(verbose){
        if(poGeometry == NULL)
          std::cerr << "no geometry defined" << std::endl << std::flush;
        else if(wkbFlatten(poGeometry->getGeometryType()) != wkbPoint)
          std::cerr << "poGeometry type: " << wkbFlatten(poGeometry->getGeometryType()) << std::endl << std::flush;
      }
      assert(poGeometry != NULL 
             && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint);
      OGRPoint *poPoint = (OGRPoint *) poGeometry;
      if(pos){
        theFeature.push_back(poPoint->getX());
        theFeature.push_back(poPoint->getY());
      }
      OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
      std::string featurename;
      for(int iField=0;iField<poFDefn->GetFieldCount();++iField){
        OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);
        std::string fieldname=poFieldDefn->GetNameRef();
        if(fieldname==label)
          theClass=poFeature->GetFieldAsInteger(iField);
        else{
          switch(fieldType){
          case(OFTReal):
            if(fields.size()<poFDefn->GetFieldCount()){
              if(find(fields.begin(),fields.end(),fieldname)!=fields.end())
                theFeature.push_back(poFeature->GetFieldAsDouble(iField));
            }
            else{
              fields[iField]=fieldname;
              theFeature.push_back(poFeature->GetFieldAsDouble(iField));
            }
            break;
          case(OFTInteger):
            if(fields.size()<poFDefn->GetFieldCount()){
              if(find(fields.begin(),fields.end(),fieldname)!=fields.end())
                theFeature.push_back(poFeature->GetFieldAsDouble(iField));
            }
            else{
              fields[iField]=fieldname;
              theFeature.push_back(poFeature->GetFieldAsDouble(iField));
            }
            break;
          default:
            {
              std::string errorstring="field type not supported in ImgReaderOgr::ReadData";
              throw(errorstring);
            }
            break;
          }
        }
      }
      data[theClass].push_back(theFeature);
      ++ifeature;
    }
    if(verbose)
      std::cout << "number of features read: " << ifeature << std::endl << std::flush;
    typename std::map<int,Vector2d<T> >::const_iterator mit=data.begin();
    int nband=0;
    if(verbose)
      std::cout << "read classes: " << std::flush;
    while(mit!=data.end()){
      if(verbose)
        std::cout << mit->first << " " << std::flush;
      if(!nband)
        nband=fields.size();
      if(pos)
        assert((mit->second)[0].size()==nband+2);
      else
        assert((mit->second)[0].size()==nband);
      ++mit;
    }
    if(verbose)
      std::cout << std::endl << std::flush;
    return(nband);
  }
  else{
    std::ostringstream ess;
    ess << "no layer in " << m_filename;
    throw(ess.str());
  }
}

template<typename T> int ImgReaderOgr::readSql(Vector2d<T>& data, const OGRFieldType& fieldType, std::vector<std::string>& fields, const std::string& sqlStatement, OGRGeometry* spatialFilter, int layer, bool pos, bool verbose)
{
  if(layer<0)
    layer=m_datasource->GetLayerCount()-1;
  assert(m_datasource->GetLayerCount()>layer);
  OGRLayer *poLayer;
  poLayer = m_datasource->ExecuteSQL(sqlStatement.c_str(), spatialFilter,NULL );
  if(poLayer!=NULL){
    OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
    if(fields.empty()){
      fields.resize(poFDefn->GetFieldCount());
      if(verbose)
        std::cout << "resized fields to " << fields.size() << std::endl;
    }
    //start reading features from the layer
    OGRFeature *poFeature;
    if(verbose)
      std::cout << "reset reading" << std::endl;
    poLayer->ResetReading();
    unsigned long int ifeature=0;
    int posOffset=(pos)?2:0;
    if(verbose)
      std::cout << "going through features" << std::endl << std::flush;
    while( (poFeature = poLayer->GetNextFeature()) != NULL ){
      std::vector<T> theFeature;//(fields.size()+posOffset);//x,y+selectedfields
      if(verbose)
        std::cout << "reading feature " << ifeature << std::endl << std::flush;
      OGRGeometry *poGeometry;
      poGeometry = poFeature->GetGeometryRef();
      if(verbose){
        if(poGeometry == NULL)
          std::cerr << "no geometry defined" << std::endl << std::flush;
        else if(wkbFlatten(poGeometry->getGeometryType()) != wkbPoint)
          std::cerr << "poGeometry type: " << wkbFlatten(poGeometry->getGeometryType()) << std::endl << std::flush;
      }
      assert(poGeometry != NULL 
             && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint);
      OGRPoint *poPoint = (OGRPoint *) poGeometry;
      if(pos){
        theFeature.push_back(poPoint->getX());
        theFeature.push_back(poPoint->getY());
      }
      OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
      std::string featurename;
      for(int iField=0;iField<poFDefn->GetFieldCount();++iField){
        OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);
        std::string fieldname=poFieldDefn->GetNameRef();
        switch(fieldType){
        case(OFTReal):
          if(fields.size()<poFDefn->GetFieldCount()){
            if(find(fields.begin(),fields.end(),fieldname)!=fields.end())
              theFeature.push_back(poFeature->GetFieldAsDouble(iField));
          }
          else{
            fields[iField]=fieldname;
            theFeature.push_back(poFeature->GetFieldAsDouble(iField));
          }
          break;
        case(OFTInteger):
          if(fields.size()<poFDefn->GetFieldCount()){
            if(find(fields.begin(),fields.end(),fieldname)!=fields.end())
              theFeature.push_back(poFeature->GetFieldAsDouble(iField));
          }
          else{
            fields[iField]=fieldname;
            theFeature.push_back(poFeature->GetFieldAsDouble(iField));
          }
          break;
        default:
          {
            std::string errorstring="field type not supported in ImgReaderOgr::ReadData";
            throw(errorstring);
          }
          break;
        }
      }
      data.push_back(theFeature);
      ++ifeature;
    }
    m_datasource->ReleaseResultSet( poLayer );
    //   assert(data.size()==ifeature);
    if(data.size()){
      if(pos)
        assert(data[0].size()==fields.size()+2);
      else
        assert(data[0].size()==fields.size());
      return fields.size();
    }
    else
      return(0);
  }
  else{
    std::ostringstream ess;
    ess << "no layer in " << m_filename;
    throw(ess.str());
  }
}

#endif // _IMGREADEROGR_H_
