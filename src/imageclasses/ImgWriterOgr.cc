/**********************************************************************
ImgWriterOgr.cc: class to write vector files using OGR API library
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
#include <sstream>
#include "ImgReaderOgr.h"
#include "ImgWriterOgr.h"
#include "ImgRasterGdal.h"
#include "cpl_string.h"
//---------------------------------------------------------------------------
ImgWriterOgr::ImgWriterOgr(void)
{}

ImgWriterOgr::~ImgWriterOgr(void)
{
}

ImgWriterOgr::ImgWriterOgr(const std::string& filename, const std::string& imageType)
{
  open(filename,imageType);
}

ImgWriterOgr::ImgWriterOgr(const std::string& filename, ImgReaderOgr& imgReaderOgr)
{
  m_filename=filename;
  setCodec(imgReaderOgr.getDriver());
  int nlayer=imgReaderOgr.getDataSource()->GetLayerCount();
  for(int ilayer=0;ilayer<nlayer;++ilayer){
    std::string layername = imgReaderOgr.getLayer(ilayer)->GetName();
    createLayer(layername,imgReaderOgr.getProjection(),imgReaderOgr.getGeometryType(),NULL);
    copyFields(imgReaderOgr,ilayer,ilayer);
  }
}

ImgWriterOgr::ImgWriterOgr(const std::string& filename, ImgReaderOgr& imgReaderOgr, bool copyData)
{
  CPLErrorReset();
  m_filename=filename;
  setCodec(imgReaderOgr.getDriver());
  int nlayer=imgReaderOgr.getDataSource()->GetLayerCount();
  for(int ilayer=0;ilayer<nlayer;++ilayer){
    std::string layername = imgReaderOgr.getLayer(ilayer)->GetName();
    createLayer(layername,imgReaderOgr.getProjection(),imgReaderOgr.getGeometryType(),NULL);
    copyFields(imgReaderOgr,ilayer,ilayer);
    if(copyData){
      OGRFeature *poFeature;
      while(true){// (poFeature = imgReaderOgr.getLayer()->GetNextFeature()) != NULL ){
	poFeature = imgReaderOgr.getLayer(ilayer)->GetNextFeature();
	if( poFeature == NULL )
	  break;
	OGRFeature *poDstFeature = NULL;

	poDstFeature=createFeature(ilayer);
	//todo: check here if SetFrom works (experienced segmentation fault)
	if( poDstFeature->SetFrom( poFeature, TRUE ) != OGRERR_NONE ){
	  const char* fmt;
	  std::string errorString="Unable to translate feature %d from layer %s.\n";
	  fmt=errorString.c_str();
	  CPLError( CE_Failure, CPLE_AppDefined,
		    fmt,
		    poFeature->GetFID(), getLayerName().c_str() );
	  // CPLError( CE_Failure, CPLE_AppDefined,
	  //           "Unable to translate feature %d from layer %s.\n",
	  //           poFeature->GetFID(), getLayerName().c_str() );
            
	  OGRFeature::DestroyFeature( poFeature );
	  OGRFeature::DestroyFeature( poDstFeature );
	}
	poDstFeature->SetFID( poFeature->GetFID() );
	OGRFeature::DestroyFeature( poFeature );

	CPLErrorReset();
	if(createFeature( poDstFeature,ilayer ) != OGRERR_NONE){
	  const char* fmt;
	  std::string errorString="Unable to translate feature %d from layer %s.\n";
	  fmt=errorString.c_str();
	  CPLError( CE_Failure, CPLE_AppDefined,
		    fmt,
		    poFeature->GetFID(), getLayerName().c_str() );
	  OGRFeature::DestroyFeature( poDstFeature );
	}
	OGRFeature::DestroyFeature( poDstFeature );
      }
    }
  }
}

//---------------------------------------------------------------------------

void ImgWriterOgr::open(const std::string& filename, ImgReaderOgr& imgReaderOgr)
{
  m_filename=filename;
  setCodec(imgReaderOgr.getDriver());
  int nlayer=imgReaderOgr.getDataSource()->GetLayerCount();
  for(int ilayer=0;ilayer<nlayer;++ilayer){
    std::string layername = imgReaderOgr.getLayer(ilayer)->GetName();
    createLayer(layername,imgReaderOgr.getProjection(),imgReaderOgr.getGeometryType(),NULL);
    copyFields(imgReaderOgr,ilayer,ilayer);
  }
}

void ImgWriterOgr::open(const std::string& filename, const std::string& imageType)
{
  m_filename = filename;
  setCodec(imageType);
}

//---------------------------------------------------------------------------
void ImgWriterOgr::close(void)
{
#if GDAL_VERSION_MAJOR < 2
  OGRDataSource::DestroyDataSource(m_datasource);
#else
  GDALClose(m_datasource);
#endif
}

//---------------------------------------------------------------------------
void ImgWriterOgr::setCodec(const std::string& imageType){
#if GDAL_VERSION_MAJOR < 2
  //register the drivers
  OGRRegisterAll();
  //fetch the OGR file driver
  OGRSFDriver *poDriver;
  poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(imageType.c_str());
#else
  //register the drivers
  GDALAllRegister();
  GDALDriver *poDriver;
  poDriver = GetGDALDriverManager()->GetDriverByName(imageType.c_str());
#endif
  if( poDriver == NULL ){
    std::string errorString="FileOpenError";
    throw(errorString);
  }
#if GDAL_VERSION_MAJOR < 2
  m_datasource = OGRSFDriverRegistrar::Open( m_filename.c_str(), TRUE );
#else
  m_datasource = (GDALDataset*) GDALOpenEx(m_filename.c_str(), GDAL_OF_UPDATE||GDAL_OF_VECTOR, NULL, NULL, NULL);
#endif
  if( m_datasource == NULL ){
#if GDAL_VERSION_MAJOR < 2
    m_datasource = OGRSFDriverRegistrar::Open( m_filename.c_str(), FALSE );
#else
    m_datasource = (GDALDataset*) GDALOpenEx(m_filename.c_str(), GDAL_OF_READONLY||GDAL_OF_VECTOR, NULL, NULL, NULL);
#endif
    if ( m_datasource != NULL){// we can only open in not update mode
      std::string errorString="Update mode not supported, delete output dataset first";
      throw(errorString);
#if GDAL_VERSION_MAJOR < 2
      OGRDataSource::DestroyDataSource(m_datasource);
#else
      GDALClose(m_datasource);
#endif
      m_datasource = NULL;
    }
    else //create the data source
#if GDAL_VERSION_MAJOR < 2
      m_datasource=poDriver->CreateDataSource(m_filename.c_str(),NULL);
#else
      m_datasource=poDriver->Create(m_filename.c_str(),0,0,0,GDT_Unknown,NULL);
#endif
  }
  else{//datasets exists, always overwrite all layers (no update append for now)
    int nLayerCount = m_datasource->GetLayerCount();
    for(int iLayer = 0; iLayer < nLayerCount; ++iLayer){
      if(m_datasource->DeleteLayer(iLayer)!=OGRERR_NONE){
	std::string errorstring="DeleteLayer() failed when overwrite requested";
	throw(errorstring);
      }
    }
  }
  if(m_datasource==NULL){
    std::string errorString="Creation of output file failed";
    throw(errorString);
  }
}

#if GDAL_VERSION_MAJOR < 2
void ImgWriterOgr::setCodec(OGRSFDriver *poDriver){
  OGRRegisterAll();
#else
void ImgWriterOgr::setCodec(GDALDriver *poDriver){
  GDALAllRegister();
#endif
  if( poDriver == NULL ){
    std::string errorString="FileOpenError";
    throw(errorString);
  }
#if GDAL_VERSION_MAJOR < 2
  m_datasource = OGRSFDriverRegistrar::Open( m_filename.c_str(), TRUE );
#else
  m_datasource = (GDALDataset*) GDALOpenEx(m_filename.c_str(), GDAL_OF_UPDATE||GDAL_OF_VECTOR, NULL, NULL, NULL);
#endif
  if( m_datasource == NULL ){
#if GDAL_VERSION_MAJOR < 2
    m_datasource = OGRSFDriverRegistrar::Open( m_filename.c_str(), FALSE );
#else
    m_datasource = (GDALDataset*) GDALOpenEx(m_filename.c_str(), GDAL_OF_READONLY||GDAL_OF_VECTOR, NULL, NULL, NULL);
#endif
    if ( m_datasource != NULL){// we can only open in not update mode
      std::string errorString="Update mode not supported, delete output dataset first";
      throw(errorString);
#if GDAL_VERSION_MAJOR < 2
      OGRDataSource::DestroyDataSource(m_datasource);
#else
      GDALClose(m_datasource);
#endif
      m_datasource = NULL;
    }
    else //create the data source
#if GDAL_VERSION_MAJOR < 2
      m_datasource=poDriver->CreateDataSource(m_filename.c_str(),NULL);
#else
      m_datasource=poDriver->Create(m_filename.c_str(),0,0,0,GDT_Unknown,NULL);
#endif
  }
  else{//datasets exists, always overwrite all layers (no update append for now)
    int nLayerCount = m_datasource->GetLayerCount();
    for(int iLayer = 0; iLayer < nLayerCount; ++iLayer){
      if(m_datasource->DeleteLayer(iLayer)!=OGRERR_NONE){
	std::string errorstring="DeleteLayer() failed when overwrite requested";
	throw(errorstring);
      }
    }
  }
  if(m_datasource==NULL){
    std::string errorString="Creation of output file failed";
    throw(errorString);
  }
}

// OGRLayer* ImgWriterOgr::copyLayer(OGRLayer* poSrcLayer, const std::string& layername, char** papszOptions)
// {
//   return(m_datasource->CopyLayer(poSrcLayer, layername.c_str(),papszOptions));
// }

OGRLayer* ImgWriterOgr::createLayer(const std::string& layername, const std::string& theProjection, const OGRwkbGeometryType& eGType, char** papszOptions)
{
  if( !m_datasource->TestCapability( ODsCCreateLayer ) ){
    std::string errorString="Test capability to create layer failed";
    throw(errorString);
  }
  //papszOptions = CSLSetNameValue( papszOptions, "DIM", "1" );
  //if points: use wkbPoint
  //if no constraints on the types geometry to be written: use wkbUnknown 
  OGRLayer* poLayer;

  OGRSpatialReference oSRS;

  if(theProjection!=""){
    oSRS.SetFromUserInput(theProjection.c_str());
    poLayer=m_datasource->CreateLayer( layername.c_str(), &oSRS, eGType,papszOptions );
    //   if(theProjection.find("EPSPG:")!=std::string::npos){
    //     int epsg_code=atoi(theProjection.substr(theProjection.find_first_not_of("EPSG:")).c_str());
    //     OGRSpatialReference oSRS;
    //     oSRS.importFromEPSG(epsg_code);
    //     poLayer=m_datasource->CreateLayer( layername.c_str(), &oSRS, eGType,papszOptions );
    //   }
    //   else{
    //     OGRSpatialReference oSRS(theProjection.c_str());
    //     poLayer=m_datasource->CreateLayer( layername.c_str(), &oSRS, eGType,papszOptions );
    //   }
    // }
    // oSRS.importFromProj4(theProjection);
  }
  else
    poLayer=m_datasource->CreateLayer( layername.c_str(), NULL, eGType,papszOptions );
  //check if destroy is needed?!
  CSLDestroy( papszOptions );
  if( poLayer == NULL ){
    std::string errorstring="Layer creation failed";
    throw(errorstring);
  }
  // OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
  return poLayer;
}

void ImgWriterOgr::createField(const std::string& fieldname, const OGRFieldType& fieldType, int theLayer)
{
  OGRFieldDefn oField( fieldname.c_str(), fieldType );
  if(fieldType==OFTString)
    oField.SetWidth(32);
  if(theLayer<0)
    theLayer=m_datasource->GetLayerCount()-1;//get back layer
  if(m_datasource->GetLayer(theLayer)->CreateField( &oField ) != OGRERR_NONE ){
      std::ostringstream es;
      es << "Creating field " << fieldname << " failed";
      std::string errorString=es.str();
      throw(errorString);
  }
}

int ImgWriterOgr::getFields(std::vector<std::string>& fields, int layer) const
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

int ImgWriterOgr::getFields(std::vector<OGRFieldDefn*>& fields, int layer) const
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
    // OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);
    fields[iField]=poFDefn->GetFieldDefn(iField);
  }
  assert(fields.size()==getFieldCount(layer));
  return(fields.size());
}

void ImgWriterOgr::copyFields(const ImgReaderOgr& imgReaderOgr, int srcLayer, int targetLayer){
  if(targetLayer<0)
    targetLayer=m_datasource->GetLayerCount()-1;//get back layer
  //get fields from imgReaderOgr
  std::vector<OGRFieldDefn*> fields;
  
  imgReaderOgr.getFields(fields,srcLayer);
  for(unsigned int iField=0;iField<fields.size();++iField){
    if(m_datasource->GetLayer(targetLayer)->CreateField(fields[iField]) != OGRERR_NONE ){
      std::ostringstream es;
      es << "Creating field " << fields[iField]->GetNameRef() << " failed";
      std::string errorString=es.str();
      throw(errorString);
    }
  }
}

void ImgWriterOgr::addPoint(double x, double y, const std::map<std::string,double>& pointAttributes, std::string fieldName, const std::string& theId, int layer){
  OGRFeature *poFeature;
  poFeature=createFeature(layer);
  OGRPoint pt;
  poFeature->SetField( fieldName.c_str(), theId.c_str());
  for(std::map<std::string,double>::const_iterator mit=pointAttributes.begin();mit!=pointAttributes.end();++mit){
    poFeature->SetField((mit->first).c_str(), mit->second);
  }
  pt.setX(x);
  pt.setY(y);
  poFeature->SetGeometry( &pt );
  if(createFeature(poFeature,layer)!=OGRERR_NONE){
    std::string errorString="Failed to create feature";
    throw(errorString);
  }
  OGRFeature::DestroyFeature( poFeature );
}

void ImgWriterOgr::addPoint(double x, double y, const std::map<std::string,double>& pointAttributes, std::string fieldName, int theId, int layer){
  OGRFeature *poFeature;
  poFeature = createFeature(layer);
  OGRPoint pt;
  if(pointAttributes.size()+1!=poFeature->GetFieldCount()){
    std::ostringstream ess;
    ess << "Failed to add feature: " << pointAttributes.size() << " !=" << poFeature->GetFieldCount() << std::endl;
    throw(ess.str());
  }
  assert(pointAttributes.size()+1==poFeature->GetFieldCount());
  poFeature->SetField( fieldName.c_str(), theId);
  for(std::map<std::string,double>::const_iterator mit=pointAttributes.begin();mit!=pointAttributes.end();++mit){
    poFeature->SetField((mit->first).c_str(),mit->second);
  }
  pt.setX(x);
  pt.setY(y);
  poFeature->SetGeometry( &pt );
  if(createFeature(poFeature,layer)!=OGRERR_NONE){
    std::string errorString="Failed to create feature";
    throw(errorString);
  }
  OGRFeature::DestroyFeature( poFeature );
}

//add a line std::string (polygon), caller is responsible to close the line (end point=start point)
void ImgWriterOgr::addLineString(std::vector<OGRPoint*>& points, const std::string& fieldName, int theId, int layer){
  OGRFeature *poFeature;
  poFeature = createFeature(layer);
  poFeature->SetStyleString("PEN(c:#FF0000,w:5px)");//see also http://www.gdal.org/ogr/ogr_feature_style.html
  poFeature->SetField( fieldName.c_str(), theId);
  OGRLineString theLineString;
  theLineString.setNumPoints(points.size());
  for(unsigned int ip=0;ip<points.size();++ip)
    theLineString.setPoint(ip,points[ip]);
  if(poFeature->SetGeometry( &theLineString )!=OGRERR_NONE){
    std::string errorString="Failed to set line OGRLineString as feature geometry";
    throw(errorString);
  }
  if(createFeature(poFeature,layer)!=OGRERR_NONE){
    std::string errorString="Failed to create feature";
    throw(errorString);
  }
  OGRFeature::DestroyFeature( poFeature );
}

//add a ring (polygon), caller is responsible to close the line (end point=start point)?
void ImgWriterOgr::addRing(std::vector<OGRPoint*>& points, const std::string& fieldName, int theId, int layer){
  OGRFeature *poFeature;
  poFeature = createFeature(layer);
  poFeature->SetStyleString("PEN(c:#FF0000,w:5px)");//see also http://www.gdal.org/ogr/ogr_feature_style.html
  poFeature->SetField( fieldName.c_str(), theId);
  // OGRLineString theLineString;
  // theLineString.setNumPoints(points.size());
  OGRPolygon thePolygon;
  OGRLinearRing theRing;
  for(unsigned int ip=0;ip<points.size();++ip)
    theRing.addPoint(points[ip]);
  //  theRing.addPoint(points[0]);//close the ring
  theRing.closeRings();//redundent with previous line?
  thePolygon.addRing(&theRing);
  // SetSpatialFilter(&thePolygon)
  poFeature->SetGeometry( &thePolygon );
  if(createFeature(poFeature,layer)!=OGRERR_NONE){
    std::string errorString="Failed to create feature";
    throw(errorString);
    OGRFeature::DestroyFeature( poFeature );
  }
  if(poFeature->SetGeometry( &thePolygon )!=OGRERR_NONE){
    std::string errorString="Failed to set polygon as feature geometry";
    throw(errorString);
  }
  OGRFeature::DestroyFeature( poFeature );
}

//add a line string (polygon), caller is responsible to close the line (end point=start point)
void ImgWriterOgr::addLineString(std::vector<OGRPoint*>& points, const std::string& fieldName, const std::string& theId, int layer){
  OGRFeature *poFeature;
  poFeature = createFeature(layer);
  poFeature->SetField( fieldName.c_str(), theId.c_str());
  OGRLineString theLineString;
  theLineString.setNumPoints(points.size());
  for(unsigned int ip=0;ip<points.size();++ip)
    theLineString.setPoint(ip,points[ip]);
  if(poFeature->SetGeometry( &theLineString )!=OGRERR_NONE){
    std::string errorString="Failed to set line OGRLineString as feature geometry";
    throw(errorString);
  }
  if(createFeature(poFeature,layer)!=OGRERR_NONE){
    std::string errorString="Failed to create feature";
    throw(errorString);
  }
  OGRFeature::DestroyFeature( poFeature );
}

OGRFeature* ImgWriterOgr::createFeature(int layer){
  return(OGRFeature::CreateFeature(m_datasource->GetLayer(layer)->GetLayerDefn()));
}

OGRErr ImgWriterOgr::createFeature(OGRFeature *theFeature,int layer){
  return m_datasource->GetLayer(layer)->CreateFeature(theFeature);
}

unsigned int ImgWriterOgr::getFieldCount(int layer) const
{
  assert(m_datasource->GetLayerCount()>layer);
  OGRLayer  *poLayer;
  if((poLayer = m_datasource->GetLayer(layer))==NULL){
    std::string errorstring="Could not get layer";
    throw(errorstring);
  }
  OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
  return(poFDefn->GetFieldCount());
}

unsigned int ImgWriterOgr::getFeatureCount(int layer) const
{
  return(getLayer(layer)->GetFeatureCount());
}

int ImgWriterOgr::ascii2ogr(const std::string& filename, const std::string &layername, const std::vector<std::string>& fieldName, const std::vector<OGRFieldType>& fieldType, short colX, short colY, const std::string& theProjection, const OGRwkbGeometryType& eGType, const char fs)
{
  char     **papszOptions=NULL;
  createLayer(layername, theProjection, eGType, papszOptions);
  //create attribute fields that should appear on the layer. Fields must be added to the layer before any features are written. To create a field we initialize an OGRField object with the information about the field. In the case of Shapefiles, the field width and precision is significant in the creation of the output .dbf file, so we set it specifically, though generally the defaults are OK
  int ncol=fieldName.size();
  assert(colX>=0);
  assert(colY>=0);
  assert(colX<ncol+2);
  assert(colY<ncol+2);
  for(int ifield=0;ifield<ncol;++ifield)
    createField(fieldName[ifield],fieldType[ifield]);
  //create a local OGRFeature, set attributes and attach geometry before trying to write it to the layer. It is imperative that this feature be instantiated from the OGRFeatureDefn associated with the layer it will be written to.
  //todo: try to open and catch if failure...
  std::ifstream fpoints(filename.c_str(),std::ios::in);
  std::string line;
  OGRPolygon thePolygon;
  OGRLinearRing theRing;
  OGRPoint firstPoint;
  OGRFeature *polyFeature=0;
  if(eGType!=wkbPoint)
    polyFeature=createFeature();


  if(fs>' '&&fs<='~'){//field separator is a regular character (minimum ASCII code is space, maximum ASCII code is tilde)
    std::string csvRecord;
    while(getline(fpoints,csvRecord)){//read a line
      OGRFeature *pointFeature=0;
      if(eGType==wkbPoint)
        pointFeature=createFeature();
      OGRPoint thePoint;
      std::istringstream csvstream(csvRecord);
      std::string value;
      int colId=0;
      int fieldId=0;
      while(getline(csvstream,value,fs)){//read a column
        if(colId==colX)
          thePoint.setX(atof(value.c_str()));
        else if(colId==colY)
          thePoint.setY(atof(value.c_str()));
        else{
          switch(fieldType[fieldId]){
          case(OFTReal):
            if(eGType==wkbPoint)
              pointFeature->SetField(fieldId,atof(value.c_str()));
            else if(firstPoint.IsEmpty())
              polyFeature->SetField(fieldId,atof(value.c_str()));
            break;
          case(OFTInteger):
            if(eGType==wkbPoint)
              pointFeature->SetField(fieldId,atoi(value.c_str()));
            else if(firstPoint.IsEmpty())
              polyFeature->SetField(fieldId,atoi(value.c_str()));
            break;
          case(OFTString):
            if(eGType==wkbPoint)
              pointFeature->SetField(fieldId,value.c_str());
            else if(firstPoint.IsEmpty())
              polyFeature->SetField(fieldId,value.c_str());
            break;
          default:
            break;
          }
          ++fieldId;
        }
        ++colId;
      }
      if(colId!=fieldId+2){
        std::ostringstream ess;
        ess << "Error: colId = " << colId << " is different from fieldId+2 = " << fieldId;
        throw(ess.str());
      }
      if(eGType==wkbPoint){
        pointFeature->SetGeometry( &thePoint );
        if(createFeature(pointFeature)!=OGRERR_NONE){
          std::string errorString="Failed to create feature";
          throw(errorString);
          OGRFeature::DestroyFeature( pointFeature );
        }
      }
      else{
        if(firstPoint.IsEmpty()){
          firstPoint=thePoint;
        }
        theRing.addPoint(&thePoint);
      }
    }
  }
  else{//space or tab delimited fields
    while(getline(fpoints,line)){
      OGRFeature *pointFeature=0;
      if(eGType==wkbPoint)
        pointFeature=createFeature();
      OGRPoint thePoint;
      std::istringstream ist(line);
      std::string value;
      int colId=0;
      int fieldId=0;
      while(ist >> value){
        if(colId==colX)
          thePoint.setX(atof(value.c_str()));
        else if(colId==colY)
          thePoint.setY(atof(value.c_str()));
        else{
          switch(fieldType[fieldId]){
          case(OFTReal):
            if(eGType==wkbPoint)
              pointFeature->SetField(fieldId,atof(value.c_str()));
            else if(firstPoint.IsEmpty())
              polyFeature->SetField(fieldId,atof(value.c_str()));
            break;
          case(OFTInteger):
            if(eGType==wkbPoint)
              pointFeature->SetField(fieldId,atoi(value.c_str()));
            else if(firstPoint.IsEmpty())
              polyFeature->SetField(fieldId,atoi(value.c_str()));
            break;
          case(OFTString):
            if(eGType==wkbPoint)
              pointFeature->SetField(fieldId,value.c_str());
            else if(firstPoint.IsEmpty())
              polyFeature->SetField(fieldId,value.c_str());
            break;
          default:
            break;
          }
          ++fieldId;
        }
        ++colId;
      }
      if(colId!=fieldId+2){
        std::ostringstream ess;
        ess << "Error: colId = " << colId << " is different from fieldId+2 = " << fieldId;
        throw(ess.str());
      }
      if(eGType==wkbPoint){
        pointFeature->SetGeometry( &thePoint );
        if(createFeature(pointFeature)!=OGRERR_NONE){
          std::string errorString="Failed to create feature";
          throw(errorString);
          OGRFeature::DestroyFeature( pointFeature );
        }
      }
      else{
        if(firstPoint.IsEmpty()){
          firstPoint=thePoint;
        }
        theRing.addPoint(&thePoint);
      }
    }
  }
  if(eGType!=wkbPoint){
    theRing.addPoint(&firstPoint);//close the ring
    thePolygon.addRing(&theRing);
    // SetSpatialFilter(&thePolygon)
    polyFeature->SetGeometry( &thePolygon );
    if(createFeature(polyFeature)!=OGRERR_NONE){
      std::string errorString="Failed to create feature";
      throw(errorString);
      OGRFeature::DestroyFeature( polyFeature );
    }
  }
  return getFeatureCount();
}

int ImgWriterOgr::addData(ImgRasterGdal& imgReader, int theLayer, bool verbose)
{
  OGRLayer  *poLayer;
  assert(m_datasource->GetLayerCount()>theLayer);
  if(verbose)
    std::cout << "number of layers: " << m_datasource->GetLayerCount() << std::endl;
  if(verbose)
    std::cout << "get layer " << theLayer << std::endl;
  poLayer = m_datasource->GetLayer(theLayer);
  //start reading features from the layer
  OGRFeature *poFeature;
  if(verbose)
    std::cout << "reset reading" << std::endl;
  poLayer->ResetReading();
  for(unsigned int iband=0;iband<imgReader.nrOfBand();++iband){
    std::ostringstream fs;
    fs << "band" << iband;
    createField(fs.str(),OFTReal,theLayer);
  }
  OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
  int nfield=poFDefn->GetFieldCount();
  if(verbose)
    std::cout << "new number of fields: " << nfield << std::endl;
  while( (poFeature = poLayer->GetNextFeature()) != NULL ){
    OGRGeometry *poGeometry;
    poGeometry = poFeature->GetGeometryRef();
    assert(poGeometry != NULL 
           && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint);
    OGRPoint *poPoint = (OGRPoint *) poGeometry;
    double x=poPoint->getX();
    double y=poPoint->getY();
    for(unsigned int iband=0;iband<imgReader.nrOfBand();++iband){
      double imgData;
      imgReader.readData(imgData,x,y,iband);
      //todo: put imgdata in field
      std::ostringstream fs;
      fs << "band" << iband;
      poFeature->SetField(fs.str().c_str(),imgData);
    }
  }
  return(nfield);
}

