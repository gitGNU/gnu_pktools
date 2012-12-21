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
#include "ImgReaderGdal.h"
#include "cpl_string.h"
//---------------------------------------------------------------------------
ImgWriterOgr::ImgWriterOgr(void)
{}

ImgWriterOgr::~ImgWriterOgr(void)
{
}

ImgWriterOgr::ImgWriterOgr(const string& filename)
{
  open(filename);
}

ImgWriterOgr::ImgWriterOgr(const string& filename, ImgReaderOgr& imgReaderOgr)
{
  m_filename=filename;
  setCodec(imgReaderOgr.getDriver());
  createLayer(filename,imgReaderOgr.getProjection(),imgReaderOgr.getGeometryType(),NULL);
  copyFields(imgReaderOgr);
}

ImgWriterOgr::ImgWriterOgr(const string& filename, ImgReaderOgr& imgReaderOgr, bool copyData)
{
  CPLErrorReset();
  m_filename=filename;
  setCodec(imgReaderOgr.getDriver());
  createLayer(filename,imgReaderOgr.getProjection(),imgReaderOgr.getGeometryType(),NULL);
  copyFields(imgReaderOgr);
  if(copyData){
    OGRFeature *poFeature;
    while( (poFeature = imgReaderOgr.getLayer()->GetNextFeature()) != NULL ){
      poFeature = imgReaderOgr.getLayer()->GetNextFeature();
      if( poFeature == NULL )
        break;
      OGRFeature *poDstFeature = NULL;

      poDstFeature=createFeature();
//       poDstFeature = OGRFeature::CreateFeature(m_datasource->GetLayer(m_datasource->GetLayerCount()-1)->GetLayerDefn());
      if( poDstFeature->SetFrom( poFeature, TRUE ) != OGRERR_NONE ){
	const char* fmt;
	string errorString="Unable to translate feature %d from layer %s.\n";
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
      if(createFeature( poDstFeature ) != OGRERR_NONE){
	const char* fmt;
	string errorString="Unable to translate feature %d from layer %s.\n";
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

//---------------------------------------------------------------------------

void ImgWriterOgr::open(const string& filename, const string& imageType)
{
  m_filename = filename;
  setCodec(imageType);
}

//---------------------------------------------------------------------------
void ImgWriterOgr::close(void)
{
  OGRDataSource::DestroyDataSource(m_datasource);
}

//---------------------------------------------------------------------------
void ImgWriterOgr::setCodec(const string& imageType){
  //register the drivers
  OGRRegisterAll();
  //fetch the shape file driver
  OGRSFDriver *poDriver;
  poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(imageType.c_str());
  if( poDriver == NULL ){
    string errorString="FileOpenError";
    throw(errorString);
  }
  //create the data source
  m_datasource=poDriver->CreateDataSource(m_filename.c_str(),NULL);
  if(m_datasource==NULL){
    string errorString="Creation of output file failed";
    throw(errorString);
  }
}

void ImgWriterOgr::setCodec(OGRSFDriver *poDriver){
  OGRRegisterAll();
  if( poDriver == NULL ){
    string errorString="FileOpenError";
    throw(errorString);
  }
  //create the data source
  m_datasource=poDriver->CreateDataSource(m_filename.c_str(),NULL);
  if(m_datasource==NULL){
    string errorString="Creation of output file failed";
    throw(errorString);
  }
}

// OGRLayer* ImgWriterOgr::copyLayer(OGRLayer* poSrcLayer, const string& layername, char** papszOptions)
// {
//   return(m_datasource->CopyLayer(poSrcLayer, layername.c_str(),papszOptions));
// }

OGRLayer* ImgWriterOgr::createLayer(const string& layername, const string& theProjection, const OGRwkbGeometryType& eGType, char** papszOptions)
{
  if( !m_datasource->TestCapability( ODsCCreateLayer ) ){
    string errorString="Test capability to create layer failed";
    throw(errorString);
  }
  //papszOptions = CSLSetNameValue( papszOptions, "DIM", "1" );
  //if points: use wkbPoint
  //if no constraints on the types geometry to be written: use wkbUnknown 
  OGRLayer* poLayer;
  if(theProjection!=""){
    if(theProjection.find("EPSPG:")!=string::npos){
      int epsg_code=atoi(theProjection.substr(theProjection.find_first_not_of("EPSG:")).c_str());
      OGRSpatialReference oSRS;
      oSRS.importFromEPSG(epsg_code);
      poLayer=m_datasource->CreateLayer( layername.c_str(), &oSRS, eGType,papszOptions );
    }
    else{
      OGRSpatialReference oSRS(theProjection.c_str());
      poLayer=m_datasource->CreateLayer( layername.c_str(), &oSRS, eGType,papszOptions );
    }
  }
  else
    poLayer=m_datasource->CreateLayer( layername.c_str(), NULL, eGType,papszOptions );
  //check if destroy is needed?!
  CSLDestroy( papszOptions );
  if( poLayer == NULL ){
    string errorstring="Layer creation failed";
    throw(errorstring);
  }
  OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
  return poLayer;
}

void ImgWriterOgr::createField(const string& fieldname, const OGRFieldType& fieldType, int theLayer)
{
  OGRFieldDefn oField( fieldname.c_str(), fieldType );
  if(fieldType==OFTString)
    oField.SetWidth(32);
  if(theLayer<0)
    theLayer=m_datasource->GetLayerCount()-1;//get back layer
  if(m_datasource->GetLayer(theLayer)->CreateField( &oField ) != OGRERR_NONE ){
      ostringstream es;
      es << "Creating field " << fieldname << " failed";
      string errorString=es.str();
      throw(errorString);
  }
}

int ImgWriterOgr::getFields(vector<string>& fields, int layer) const
{
  if(layer<0)
    layer=m_datasource->GetLayerCount()-1;
  assert(m_datasource->GetLayerCount()>layer);
  OGRLayer  *poLayer;
  if((poLayer = m_datasource->GetLayer(layer))==NULL){
    string errorstring="Could not get layer";
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

int ImgWriterOgr::getFields(vector<OGRFieldDefn*>& fields, int layer) const
{
  if(layer<0)
    layer=m_datasource->GetLayerCount()-1;
  assert(m_datasource->GetLayerCount()>layer);
  OGRLayer  *poLayer;
  if((poLayer = m_datasource->GetLayer(layer))==NULL){
    string errorstring="Could not get layer";
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

void ImgWriterOgr::copyFields(const ImgReaderOgr& imgReaderOgr, int theLayer){
  if(theLayer<0)
    theLayer=m_datasource->GetLayerCount()-1;//get back layer
  //get fields from imgReaderOgr
  vector<OGRFieldDefn*> fields;
  
  imgReaderOgr.getFields(fields);
//   OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
  for(int iField=0;iField<fields.size();++iField){
    if(m_datasource->GetLayer(theLayer)->CreateField(fields[iField]) != OGRERR_NONE ){
      ostringstream es;
      es << "Creating field " << fields[iField]->GetNameRef() << " failed";
      string errorString=es.str();
      throw(errorString);
    }
  }
}

void ImgWriterOgr::addPoint(double x, double y, const map<string,double>& pointAttributes, string fieldName, const string& theId){
  OGRFeature *poFeature;
  poFeature=createFeature();
  OGRPoint pt;
  poFeature->SetField( fieldName.c_str(), theId.c_str());
  for(map<string,double>::const_iterator mit=pointAttributes.begin();mit!=pointAttributes.end();++mit){
    poFeature->SetField((mit->first).c_str(), mit->second);
  }
  pt.setX(x);
  pt.setY(y);
  poFeature->SetGeometry( &pt );
  if(createFeature(poFeature)!=OGRERR_NONE){
    string errorString="Failed to create feature in shapefile";
    throw(errorString);
  }
  OGRFeature::DestroyFeature( poFeature );
}

void ImgWriterOgr::addPoint(double x, double y, const map<string,double>& pointAttributes, string fieldName, int theId){
  OGRFeature *poFeature;
  poFeature = createFeature();
  OGRPoint pt;
  if(pointAttributes.size()+1!=poFeature->GetFieldCount()){
    ostringstream ess;
    ess << "Failed to add feature: " << pointAttributes.size() << " !=" << poFeature->GetFieldCount() << endl;
    throw(ess.str());
  }
  assert(pointAttributes.size()+1==poFeature->GetFieldCount());
  poFeature->SetField( fieldName.c_str(), theId);
  int fid=0;
  for(map<string,double>::const_iterator mit=pointAttributes.begin();mit!=pointAttributes.end();++mit){
    poFeature->SetField((mit->first).c_str(),mit->second);
  }
  pt.setX(x);
  pt.setY(y);
  poFeature->SetGeometry( &pt );
  if(createFeature(poFeature)!=OGRERR_NONE){
    string errorString="Failed to create feature in shapefile";
    throw(errorString);
  }
  OGRFeature::DestroyFeature( poFeature );
}

//add a line string (polygon), caller is responsible to close the line (end point=start point)
void ImgWriterOgr::addLineString(vector<OGRPoint*>& points, const string& fieldName, int theId){
  OGRFeature *poFeature;
  poFeature = createFeature();
  poFeature->SetStyleString("PEN(c:#FF0000,w:5px)");//see also http://www.gdal.org/ogr/ogr_feature_style.html
  poFeature->SetField( fieldName.c_str(), theId);
  OGRLineString theLineString;
  theLineString.setNumPoints(points.size());
  for(int ip=0;ip<points.size();++ip)
    theLineString.setPoint(ip,points[ip]);
  if(poFeature->SetGeometry( &theLineString )!=OGRERR_NONE){
    string errorString="Failed to set line OGRLineString as feature geometry";
    throw(errorString);
  }
  if(createFeature(poFeature)!=OGRERR_NONE){
    string errorString="Failed to create feature in shapefile";
    throw(errorString);
  }
  OGRFeature::DestroyFeature( poFeature );
}

//add a ring (polygon), caller is responsible to close the line (end point=start point)?
void ImgWriterOgr::addRing(vector<OGRPoint*>& points, const string& fieldName, int theId){
  OGRFeature *poFeature;
  poFeature = createFeature();
  poFeature->SetStyleString("PEN(c:#FF0000,w:5px)");//see also http://www.gdal.org/ogr/ogr_feature_style.html
  poFeature->SetField( fieldName.c_str(), theId);
  // OGRLineString theLineString;
  // theLineString.setNumPoints(points.size());
  OGRPolygon thePolygon;
  OGRLinearRing theRing;
  for(int ip=0;ip<points.size();++ip)
    theRing.addPoint(points[ip]);
  //  theRing.addPoint(points[0]);//close the ring
  theRing.closeRings();//redundent with previous line?
  thePolygon.addRing(&theRing);
  // SetSpatialFilter(&thePolygon)
  poFeature->SetGeometry( &thePolygon );
  if(createFeature(poFeature)!=OGRERR_NONE){
    string errorString="Failed to create feature in shapefile";
    throw(errorString);
    OGRFeature::DestroyFeature( poFeature );
  }
  if(poFeature->SetGeometry( &thePolygon )!=OGRERR_NONE){
    string errorString="Failed to set polygon as feature geometry";
    throw(errorString);
  }
  OGRFeature::DestroyFeature( poFeature );
}

//add a line string (polygon), caller is responsible to close the line (end point=start point)
void ImgWriterOgr::addLineString(vector<OGRPoint*>& points, const string& fieldName, const string& theId){
  OGRFeature *poFeature;
  poFeature = createFeature();
  poFeature->SetField( fieldName.c_str(), theId.c_str());
  OGRLineString theLineString;
  theLineString.setNumPoints(points.size());
  for(int ip=0;ip<points.size();++ip)
    theLineString.setPoint(ip,points[ip]);
  if(poFeature->SetGeometry( &theLineString )!=OGRERR_NONE){
    string errorString="Failed to set line OGRLineString as feature geometry";
    throw(errorString);
  }
  if(createFeature(poFeature)!=OGRERR_NONE){
    string errorString="Failed to create feature in shapefile";
    throw(errorString);
  }
  OGRFeature::DestroyFeature( poFeature );
}

OGRFeature* ImgWriterOgr::createFeature(){
  return(OGRFeature::CreateFeature(m_datasource->GetLayer(m_datasource->GetLayerCount()-1)->GetLayerDefn()));
}

OGRErr ImgWriterOgr::createFeature(OGRFeature *theFeature){
  return m_datasource->GetLayer(m_datasource->GetLayerCount()-1)->CreateFeature(theFeature);
}

int ImgWriterOgr::getFieldCount(int layer) const
{
  if(layer<0)
    layer=m_datasource->GetLayerCount()-1;
  assert(m_datasource->GetLayerCount()>layer);
  OGRLayer  *poLayer;
  if((poLayer = m_datasource->GetLayer(layer))==NULL){
    string errorstring="Could not get layer";
    throw(errorstring);
  }
  OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
  return(poFDefn->GetFieldCount());
}

int ImgWriterOgr::getFeatureCount() const
{
  return(getLayer()->GetFeatureCount());
}

int ImgWriterOgr::ascii2ogr(const string& filename, const string &layername, const vector<string>& fieldName, const vector<OGRFieldType>& fieldType, short colX, short colY, const string& theProjection, const OGRwkbGeometryType& eGType, const char fs)
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
  ifstream fpoints(filename.c_str(),ios::in);
  string line;
  OGRPolygon thePolygon;
  OGRLinearRing theRing;
  OGRPoint firstPoint;
  OGRFeature *polyFeature;
  if(eGType!=wkbPoint)
    polyFeature=createFeature();


  if(fs>' '&&fs<='~'){//field separator is a regular character (minimum ASCII code is space, maximum ASCII code is tilde)
    string csvRecord;
    while(getline(fpoints,csvRecord)){//read a line
      OGRFeature *pointFeature;
      if(eGType==wkbPoint)
        pointFeature=createFeature();
      OGRPoint thePoint;
      bool skip=false;
      istringstream csvstream(csvRecord);
      string value;
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
              polyFeature->SetField(fieldId,atof(value.c_str()));
            break;
          case(OFTString):
            if(eGType==wkbPoint)
              pointFeature->SetField(fieldId,value.c_str());
            else if(firstPoint.IsEmpty())
              polyFeature->SetField(fieldId,atof(value.c_str()));
            break;
          default:
            break;
          }
          ++fieldId;
        }
        ++colId;
      }
      if(colId!=fieldId+2){
        ostringstream ess;
        ess << "Error: colId = " << colId << " is different from fieldId+2 = " << fieldId;
        throw(ess.str());
      }
      if(eGType==wkbPoint){
        pointFeature->SetGeometry( &thePoint );
        if(createFeature(pointFeature)!=OGRERR_NONE){
          string errorString="Failed to create feature in shapefile";
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
      OGRFeature *pointFeature;
      if(eGType==wkbPoint)
        pointFeature=createFeature();
      OGRPoint thePoint;
      bool skip=false;
      istringstream ist(line);
      string value;
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
              polyFeature->SetField(fieldId,atof(value.c_str()));
            break;
          case(OFTString):
            if(eGType==wkbPoint)
              pointFeature->SetField(fieldId,value.c_str());
            else if(firstPoint.IsEmpty())
              polyFeature->SetField(fieldId,atof(value.c_str()));
            break;
          default:
            break;
          }
          ++fieldId;
        }
        ++colId;
      }
      if(colId!=fieldId+2){
        ostringstream ess;
        ess << "Error: colId = " << colId << " is different from fieldId+2 = " << fieldId;
        throw(ess.str());
      }
      if(eGType==wkbPoint){
        pointFeature->SetGeometry( &thePoint );
        if(createFeature(pointFeature)!=OGRERR_NONE){
          string errorString="Failed to create feature in shapefile";
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
      string errorString="Failed to create feature in shapefile";
      throw(errorString);
      OGRFeature::DestroyFeature( polyFeature );
    }
  }
  return getFeatureCount();
}

int ImgWriterOgr::addData(const ImgReaderGdal& imgReader, int theLayer, bool verbose)
{
  OGRLayer  *poLayer;
  if(theLayer<0)
    theLayer=m_datasource->GetLayerCount()-1;//get back layer
  assert(m_datasource->GetLayerCount()>theLayer);
  if(verbose)
    cout << "number of layers: " << m_datasource->GetLayerCount() << endl;
  if(verbose)
    cout << "get layer " << theLayer << endl;
  poLayer = m_datasource->GetLayer(theLayer);
  //start reading features from the layer
  OGRFeature *poFeature;
  if(verbose)
    cout << "reset reading" << endl;
  poLayer->ResetReading();
  for(int iband=0;iband<imgReader.nrOfBand();++iband){
    ostringstream fs;
    fs << "band" << iband;
    createField(fs.str(),OFTReal);
  }
  OGRFeatureDefn *poFDefn = poLayer->GetLayerDefn();
  int nfield=poFDefn->GetFieldCount();
  if(verbose)
    cout << "new number of fields: " << nfield << endl;
  while( (poFeature = poLayer->GetNextFeature()) != NULL ){
    OGRGeometry *poGeometry;
    poGeometry = poFeature->GetGeometryRef();
    assert(poGeometry != NULL 
           && wkbFlatten(poGeometry->getGeometryType()) == wkbPoint);
    OGRPoint *poPoint = (OGRPoint *) poGeometry;
    double x=poPoint->getX();
    double y=poPoint->getY();
    for(int iband=0;iband<imgReader.nrOfBand();++iband){
      double imgData;
      imgReader.readData(imgData,GDT_Float64,x,y,iband);
      //todo: put imgdata in field
      ostringstream fs;
      fs << "band" << iband;
      poFeature->SetField(fs.str().c_str(),imgData);
    }
  }
  return(nfield);
}

