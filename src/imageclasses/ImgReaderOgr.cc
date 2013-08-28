/**********************************************************************
ImgReaderOgr.cc: class to read vector files using OGR API library
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
#include <iostream>
#include <fstream>
#include "ImgReaderOgr.h"
#include "ImgWriterOgr.h"
#include "cpl_string.h"
//---------------------------------------------------------------------------
ImgReaderOgr::ImgReaderOgr(void)
{}

ImgReaderOgr::ImgReaderOgr(const std::string& filename)
{
  open(filename);
}

ImgReaderOgr::~ImgReaderOgr(void)
{
}

//---------------------------------------------------------------------------

void ImgReaderOgr::open(const std::string& filename)
{
  m_filename = filename;
  setCodec();
}

//---------------------------------------------------------------------------
void ImgReaderOgr::close(void)
{
  OGRDataSource::DestroyDataSource(m_datasource);
}

//---------------------------------------------------------------------------
void ImgReaderOgr::setCodec(void){
  //register the drivers
  OGRRegisterAll();
  //open the input OGR datasource. Datasources can be files, RDBMSes, directories full of files, or even remote web services depending on the driver being used. However, the datasource name is always a single string.
  m_datasource = OGRSFDriverRegistrar::Open(m_filename.c_str(), FALSE);//FAlSE: do not update
  if( m_datasource == NULL ){
    std::string errorString="Open failed";
    throw(errorString);
  }
}

bool ImgReaderOgr::getExtent(double& ulx, double& uly, double& lrx, double& lry, int layer)
{
  OGREnvelope oExt;
  if(getLayer(layer)->GetExtent(&oExt,TRUE)==OGRERR_NONE){
    ulx=oExt.MinX;
    uly=oExt.MaxY;
    lrx=oExt.MaxX;
    lry=oExt.MinY;
    return true;
  }
  else
    return false;
}

unsigned long int ImgReaderOgr::getFeatureCount(int layer) const
{
  return(m_datasource->GetLayer(layer)->GetFeatureCount());
}

int ImgReaderOgr::getFieldCount(int layer) const
{
  if(layer<0)
    layer=m_datasource->GetLayerCount()-1;
  assert(m_datasource->GetLayerCount()>layer);
  OGRLayer  *poLayer;
  if((poLayer = m_datasource->GetLayer(layer))==NULL){
    std::string errorstring="Could not get layer";
    throw(errorstring);
  }
  OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
  return(poFDefn->GetFieldCount());
}

int ImgReaderOgr::getFields(std::vector<std::string>& fields, int layer) const
{
  if(layer<0)
    layer=m_datasource->GetLayerCount()-1;
  assert(m_datasource->GetLayerCount()>layer);
  OGRLayer  *poLayer;
  if((poLayer = m_datasource->GetLayer(layer))==NULL){
    std::string errorstring="Could not get layer";
    throw(errorstring);
  }
  OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
  fields.clear();
  fields.resize(poFDefn->GetFieldCount());
  for(int iField=0;iField<poFDefn->GetFieldCount();++iField){
    OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);
    fields[iField]=poFieldDefn->GetNameRef();
  }
  return(fields.size());
}

int ImgReaderOgr::getFields(std::vector<OGRFieldDefn*>& fields, int layer) const
{
  if(layer<0)
    layer=m_datasource->GetLayerCount()-1;
  assert(m_datasource->GetLayerCount()>layer);
  OGRLayer  *poLayer;
  if((poLayer = m_datasource->GetLayer(layer))==NULL){
    std::string errorstring="Could not get layer";
    throw(errorstring);
  }
  OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
  fields.clear();
  fields.resize(poFDefn->GetFieldCount());
  for(int iField=0;iField<poFDefn->GetFieldCount();++iField){
    OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);
    fields[iField]=poFDefn->GetFieldDefn(iField);
  }
  assert(fields.size()==getFieldCount());
  return(fields.size());
}

std::string ImgReaderOgr::getProjection(int layer) const
{
  if(m_datasource->GetLayer(layer)->GetSpatialRef()){
    char* ppszResult;
    m_datasource->GetLayer(layer)->GetSpatialRef()->exportToWkt(&ppszResult);
    return(ppszResult);
  }
  else
    return "";
}

OGRwkbGeometryType ImgReaderOgr::getGeometryType(int layer) const
{
  return m_datasource->GetLayer(layer)->GetLayerDefn()->GetGeomType();
}

std::ostream& operator<<(std::ostream& theOstream, ImgReaderOgr& theImageReader){
  //An OGRDataSource can potentially have many layers associated with it. The number of layers available can be queried with OGRDataSource::GetLayerCount() and individual layers fetched by index using OGRDataSource::GetLayer(). However, we wil just fetch the layer by name.
  //todo: try to open and catch if failure...
  // ofstream fpoints(filename.c_str(),ios::out);
  OGRLayer  *poLayer;
  poLayer = theImageReader.getDataSource()->GetLayer(0);
  OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();

  poLayer->ResetReading();

  theOstream << "#";
  int iField=0;
  // theOstream << "X" << " " << "Y" << " ";
  for(int iField=0;iField<poFDefn->GetFieldCount();++iField){
      OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);
      std::string fieldname=poFieldDefn->GetNameRef();
      theOstream << fieldname << " ";
  }
  theOstream << std::endl;

  poLayer->ResetReading();
  
  //start reading features from the layer
  OGRFeature *poFeature;
  unsigned long int ifeature=0;
  while( (poFeature = poLayer->GetNextFeature()) != NULL ){
    OGRGeometry *poGeometry;
    poGeometry = poFeature->GetGeometryRef();
    assert(poGeometry != NULL);
    double x,y;
    if(wkbFlatten(poGeometry->getGeometryType()) == wkbPoint){
      OGRPoint *poPoint = (OGRPoint *) poGeometry;
      x=poPoint->getX();
      y=poPoint->getY();
    }
    std::vector<std::string> vfields(poFDefn->GetFieldCount());
    std::string featurename;
    std::vector<std::string>::iterator fit=vfields.begin();
    for(int iField=0;iField<poFDefn->GetFieldCount();++iField){
      *(fit++)=poFeature->GetFieldAsString(iField);
    }
    theOstream.precision(12);
    if(wkbFlatten(poGeometry->getGeometryType()) == wkbPoint)
      theOstream << x << " " << y;
    for(fit=vfields.begin();fit!=vfields.end();++fit)
      theOstream << " " << *fit;
    theOstream << endl;
    ++ifeature;
  }
  return(theOstream);
}

// OGRLayer * ImgReaderOgr::executeSql(const std::string& output, const std::string& sqlStatement, OGRGeometry* spatialFilter)
// {
//   OGRLayer *poResultSet;
//   poResultSet = m_datasource->ExecuteSQL(sqlStatement.c_str(), spatialFilter,NULL );

//   if( poResultSet != NULL ){
//     ImgWriterOgr imgWriter;
//     imgWriter.open(output);
//     imgWriter.copyLayer(poResultSet,output);
//     m_datasource->ReleaseResultSet( poResultSet );
//     imgWriter.close();
//   }
// }
