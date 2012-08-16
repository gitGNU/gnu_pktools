/**********************************************************************
ImgWriterGdal.cc: class to write raster files using GDAL API library
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
#include <iomanip>
#include <time.h>
#include "ogr_spatialref.h"
#include "ImgWriterGdal.h"
//---------------------------------------------------------------------------
ImgWriterGdal::ImgWriterGdal(void)
  : m_gds(NULL), m_isGeoRef(false), m_ncol(0), m_nrow(0), m_nband(0)
{}

// ImgWriterGdal::ImgWriterGdal(void)
//   : m_gds(NULL), m_magic_x(1), m_magic_y(1), m_isGeoRef(false), m_ncol(0), m_nrow(0), m_nband(0), m_interleave("BAND"), m_compression("LZW")
// {}

ImgWriterGdal::~ImgWriterGdal(void)
{
  // delete m_gds;
//   GDALDumpOpenDatasets(stderr);
//   GDALDestroyDriverManager();//could still be be used by other objects
}

//---------------------------------------------------------------------------
void ImgWriterGdal::open(const string& filename, const ImgReaderGdal& imgSrc, const vector<string>& options)
{
  m_isGeoRef=imgSrc.isGeoRef();
  m_filename=filename;
  m_ncol=imgSrc.nrOfCol();
  m_nrow=imgSrc.nrOfRow();
  m_nband=imgSrc.nrOfBand();
  m_type=imgSrc.getDataType();
  m_options=options;
  // m_interleave=imgSrc.getInterleave();
  // m_compression=imgSrc.getCompression();
  // imgSrc.getMagicPixel(m_magic_x,m_magic_y);
  setCodec(imgSrc);
}

// void ImgWriterGdal::open(const string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const string& imageType, const string& interleave, const string& compression, int magicX, int magicY)
// {
//   m_isGeoRef=false;
//   m_filename = filename;
//   m_ncol = ncol;
//   m_nrow = nrow;
//   m_nband = nband;
//   m_type=dataType;
//   m_interleave = interleave;
//   m_compression=compression;
//   m_magic_x=magicX;
//   m_magic_y=magicY;
//   setCodec(imageType);
// }

void ImgWriterGdal::open(const string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const string& imageType, const vector<string>& options)
{
  m_isGeoRef=false;
  m_filename = filename;
  m_ncol = ncol;
  m_nrow = nrow;
  m_nband = nband;
  m_type=dataType;
  // m_interleave = interleave;
  // m_compression=compression;
  m_options=options;
  // m_magic_x=magicX;
  // m_magic_y=magicY;
  setCodec(imageType);
}

//---------------------------------------------------------------------------
void ImgWriterGdal::close(void)
{
  GDALClose(m_gds);
  char **papszOptions=NULL;
  for(vector<string>::const_iterator optionIt=m_options.begin();optionIt!=m_options.end();++optionIt)
    papszOptions=CSLAddString(papszOptions,optionIt->c_str());
  CSLDestroy(papszOptions);
}

//---------------------------------------------------------------------------
void ImgWriterGdal::setCodec(const ImgReaderGdal& imgSrc){
  GDALAllRegister();
  GDALDriver *poDriver;
  poDriver = GetGDALDriverManager()->GetDriverByName(imgSrc.getDriverDescription().c_str());
  if( poDriver == NULL ){
    string errorString="FileOpenError";
    throw(errorString);
  }
  char **papszMetadata;
  papszMetadata = poDriver->GetMetadata();
  assert( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE ));
  char **papszOptions=NULL;
  for(vector<string>::const_iterator optionIt=m_options.begin();optionIt!=m_options.end();++optionIt)
    papszOptions=CSLAddString(papszOptions,optionIt->c_str());
  // char **papszOptions=NULL;
  // ostringstream compressList;
  // compressList << "COMPRESS=" << m_compression;
  // papszOptions = CSLAddString(papszOptions,(compressList.str()).c_str());
  // ostringstream interleaveList;
  // interleaveList << "INTERLEAVE=" << m_interleave;
  // papszOptions = CSLAddString(papszOptions,(interleaveList.str()).c_str());
  m_gds=poDriver->Create(m_filename.c_str(),m_ncol,m_nrow,m_nband,m_type,papszOptions);
  if(imgSrc.isGeoRef()){
    setProjection(imgSrc.getProjection());
    double ulx,uly,deltaX,deltaY,rot1,rot2;
    imgSrc.getGeoTransform(ulx,uly,deltaX,deltaY,rot1,rot2);
    setGeoTransform(ulx,uly,deltaX,deltaY,rot1,rot2);
  }
  m_gds->SetMetadata(imgSrc.getMetadata() ); 
  m_gds->SetMetadataItem( "TIFFTAG_DOCUMENTNAME", m_filename.c_str());
  m_gds->SetMetadataItem( "TIFFTAG_SOFTWARE", "pktools v2.0 by Pieter Kempeneers");
  time_t rawtime;
  time ( &rawtime );

  time_t tim=time(NULL);
  tm *now=localtime(&tim);
  ostringstream datestream;
  //date string must be 20 characters long...
  datestream << now->tm_year+1900;
  if(now->tm_mon+1<10)
    datestream << ":0" << now->tm_mon+1;
  else
    datestream << ":" << now->tm_mon+1;
  if(now->tm_mday<10)
    datestream << ":0" << now->tm_mday;
  else
    datestream << ":" << now->tm_mday;
  if(now->tm_hour<10)
    datestream << " 0" << now->tm_hour;
  else
    datestream << " " << now->tm_hour;
  if(now->tm_min<10)
    datestream << ":0" << now->tm_min;
  else
    datestream << ":" << now->tm_min;
  if(now->tm_sec<10)
    datestream << ":0" << now->tm_sec;
  else
    datestream << ":" << now->tm_sec;
  m_gds->SetMetadataItem( "TIFFTAG_DATETIME", datestream.str().c_str());
//   list<string> lmeta;
//   imgReader.getMetadata(lmeta);
//   list<string>::const_iterator lit=lmeta.begin();
//   while(lit!=lmeta.end()){
//     cout << *lit << endl;
//     ++lit;
//   }
  // m_gds->SetMetadataItem( "INTERLEAVE", m_interleave.c_str(), "IMAGE_STRUCTURE" );
  // m_gds->SetMetadataItem( "COMPRESS", m_compression.c_str(), "IMAGE_STRUCTURE" );
  if(imgSrc.getColorTable()!=NULL)
    setColorTable(imgSrc.getColorTable());
}

void ImgWriterGdal::setCodec(const string& imageType)
{
  GDALAllRegister();
  GDALDriver *poDriver;
  poDriver = GetGDALDriverManager()->GetDriverByName(imageType.c_str());
  if( poDriver == NULL ){
    ostringstream s;
    s << "FileOpenError (" << imageType << ")";
    throw(s.str());
  }
  char **papszMetadata;
  papszMetadata = poDriver->GetMetadata();
  assert( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE ));
  char **papszOptions=NULL;
  for(vector<string>::const_iterator optionIt=m_options.begin();optionIt!=m_options.end();++optionIt)
    papszOptions=CSLAddString(papszOptions,optionIt->c_str());
  // ostringstream compressList;
  // compressList << "COMPRESS=" << m_compression;
  // papszOptions = CSLAddString(papszOptions,(compressList.str()).c_str());
  // ostringstream interleaveList;
  // interleaveList << "INTERLEAVE=" << m_interleave;
  // papszOptions = CSLAddString(papszOptions,(interleaveList.str()).c_str());
  m_gds=poDriver->Create(m_filename.c_str(),m_ncol,m_nrow,m_nband,m_type,papszOptions);

  // m_gds->SetMetadataItem( "INTERLEAVE", m_interleave.c_str(), "IMAGE_STRUCTURE" );
  // m_gds->SetMetadataItem( "COMPRESSION", m_compression.c_str(), "IMAGE_STRUCTURE" );
  m_gds->SetMetadataItem( "TIFFTAG_DOCUMENTNAME", m_filename.c_str());
  m_gds->SetMetadataItem( "TIFFTAG_SOFTWARE", "pktools v2.0 by Pieter Kempeneers");
  time_t rawtime;
  time ( &rawtime );

  time_t tim=time(NULL);
  tm *now=localtime(&tim);
  ostringstream datestream;
  //date string must be 20 characters long...
  datestream << now->tm_year+1900;
  if(now->tm_mon+1<10)
    datestream << ":0" << now->tm_mon+1;
  else
    datestream << ":" << now->tm_mon+1;
  if(now->tm_mday<10)
    datestream << ":0" << now->tm_mday;
  else
    datestream << ":" << now->tm_mday;
  if(now->tm_hour<10)
    datestream << " 0" << now->tm_hour;
  else
    datestream << " " << now->tm_hour;
  if(now->tm_min<10)
    datestream << ":0" << now->tm_min;
  else
    datestream << ":" << now->tm_min;
  if(now->tm_sec<10)
    datestream << ":0" << now->tm_sec;
  else
    datestream << ":" << now->tm_sec;
  m_gds->SetMetadataItem( "TIFFTAG_DATETIME", datestream.str().c_str());
//   m_gds->SetMetadataItem( "TIFFTAG_DATETIME", ctime(&rawtime));
}

void ImgWriterGdal::setMetadata(char** metadata)
{
  assert(m_gds);
  m_gds->SetMetadata(metadata); 
}

string ImgWriterGdal::getProjection(void) const 
{
  assert(m_gds);
  string theProjection=m_gds->GetProjectionRef();
  //due to error in Gdal? AUTHORITY fields do not seem to work!
  // size_t startpos,endpos;
  // while((startpos=theProjection.find(",AUTHORITY"))!=string::npos){
  //   endpos=theProjection.find("]",startpos+1,1)+1;
  //   theProjection.erase(startpos,endpos-startpos);
  // }
  return theProjection;
}

//---------------------------------------------------------------------------
void ImgWriterGdal::setGeoTransform(double ulx, double uly, double deltaX, double deltaY, double rot1, double rot2)
{
  m_isGeoRef=true;
  m_ulx=ulx;
  m_uly=uly;
  m_delta_x=deltaX;
  m_delta_y=deltaY;
  double adfGeoTransform[6];// { 444720, 30, 0, 3751320, 0, -30 };
  adfGeoTransform[0]=ulx;
  adfGeoTransform[1]=deltaX;
  adfGeoTransform[2]=rot1;
  adfGeoTransform[3]=uly;
  adfGeoTransform[4]=rot2;
  adfGeoTransform[5]=-deltaY;//convention of GDAL!
  if(m_gds)
    m_gds->SetGeoTransform(adfGeoTransform);
}

void ImgWriterGdal::copyGeoTransform(const ImgReaderGdal& imgSrc)
{
  setProjection(imgSrc.getProjection());
  double ulx,uly,deltaX,deltaY,rot1,rot2;
  imgSrc.getGeoTransform(ulx,uly,deltaX,deltaY,rot1,rot2);
  setGeoTransform(ulx,uly,deltaX,deltaY,rot1,rot2);
}

string ImgWriterGdal::setProjectionProj4(const string& projection)
{
  if(!m_isGeoRef)
    m_isGeoRef=true;

    OGRSpatialReferenceH hSRS;  
    char *pszResult = NULL;  
  
    CPLErrorReset();  
      
    hSRS = OSRNewSpatialReference( NULL );  
    if( OSRSetFromUserInput( hSRS, projection.c_str() ) == OGRERR_NONE )  
        OSRExportToWkt( hSRS, &pszResult );  
    else  
    {  
        // CPLError( CE_Failure, CPLE_AppDefined,  
        //           "Translating source or target SRS failed:\n%s",  
        //           projection.c_str() );  
        ostringstream s;
        s << "Error in set projection " << projection;
        throw(s.str());
        // exit( 1 );
    }  
    // m_gds->SetProjection(pszResult);

    //due to error in Gdal? AUTHORITY fields do not seem to work!
    string theProjection=pszResult;
    // size_t startpos,endpos;
    // while((startpos=theProjection.find(",AUTHORITY"))!=string::npos){
    //   endpos=theProjection.find("]",startpos+1,1)+1;
    //   theProjection.erase(startpos,endpos-startpos);
    // }
    assert(m_gds);
    m_gds->SetProjection(theProjection.c_str());
    OSRDestroySpatialReference( hSRS );  
  
    // return pszResult;  
    return theProjection;  

//   OGRSpatialReference oSRS;
//   char *pszSRS_WKT = NULL;
//   oSRS.importFromProj4(projection.c_str());
// //   oSRS.SetLAEA(52,10,4321000,3210000);//redundant but according to JRC standard
//   oSRS.exportToWkt(&pszSRS_WKT);
//   m_gds->SetProjection(pszSRS_WKT);
//   CPLFree(pszSRS_WKT);
}

void ImgWriterGdal::setProjection(const string& projection)
{
  if(!m_isGeoRef)
    m_isGeoRef=true;
  OGRSpatialReference oSRS;
  char *pszSRS_WKT = NULL;
  assert(m_gds);
  m_gds->SetProjection(projection.c_str());
  CPLFree(pszSRS_WKT);
}

//default projection: ETSR-LAEA
string ImgWriterGdal::setProjection(void)
{
  string theProjection;
  OGRSpatialReference oSRS;
  char *pszSRS_WKT = NULL;
  //// oSRS.importFromEPSG(3035);
  oSRS.SetGeogCS("ETRS89","European_Terrestrial_Reference_System_1989","GRS 1980",6378137,298.2572221010042,"Greenwich",0,"degree",0.0174532925199433);
  // cout << setprecision(16) << "major axis: " << oSRS.GetSemiMajor(NULL) << endl;//notice that major axis can be set to a different value than the default to the well known standard corresponding to the name (European_Terrestrial_Reference_System_1989), but that new value, while recognized by GetSemiMajor, will not be written in the geotiff tag!
  oSRS.SetProjCS( "ETRS89 / ETRS-LAEA" );
  oSRS.SetLAEA(52,10,4321000,3210000);
  oSRS.exportToWkt( &pszSRS_WKT );
  theProjection=pszSRS_WKT;
  CPLFree( pszSRS_WKT );
  assert(m_gds);
  m_gds->SetProjection(theProjection.c_str());
  return(theProjection);
}

bool ImgWriterGdal::getBoundingBox(double& ulx, double& uly, double& lrx, double& lry) const
{
  if(m_isGeoRef){
//     ulx=m_ulx-(m_magic_x-1.0)*m_delta_x;
//     uly=m_uly+(m_magic_y-1.0)*m_delta_y;
//     lrx=ulx+(nrOfCol()+1.0-m_magic_x)*m_delta_x;
//     lry=uly-(nrOfRow()+1.0-m_magic_y)*m_delta_y;
    ulx=m_ulx;
    uly=m_uly;
    lrx=ulx+nrOfCol()*m_delta_x;
    lry=uly-nrOfRow()*m_delta_y;
    return true;
  }
  else{
    ulx=0;
    uly=nrOfRow()-1;
    lrx=nrOfCol()-1;
    lry=0;
    return false;
  }
}

bool ImgWriterGdal::getCentrePos(double& x, double& y) const
{
  if(m_isGeoRef){
//     x=m_ulx+(nrOfCol()/2.0-(m_magic_x-1.0))*m_delta_x;
//     y=m_uly-(nrOfRow()/2.0+(m_magic_y-1.0))*m_delta_y;
    x=m_ulx+nrOfCol()/2.0*m_delta_x;
    y=m_uly-nrOfRow()/2.0*m_delta_y;
    return true;
  }
  else
    return false;
}

bool ImgWriterGdal::geo2image(double x, double y, double& i, double& j) const
{
  //double values are returned, caller is responsible for interpolation step
  if(m_isGeoRef){
//     double ulx=m_ulx-(m_magic_x-1.0)*m_delta_x;
//     double uly=m_uly+(m_magic_y-1.0)*m_delta_y;
    double ulx=m_ulx;
    double uly=m_uly;
    i=(x-ulx)/m_delta_x;
    j=(uly-y)/m_delta_y;
    return true;
  }
  else{
    i=x;
    j=nrOfRow()-y;
    return false;
  }
}

//centre of pixel is always returned (regardless of magic pixel reference)!
bool ImgWriterGdal::image2geo(double i, double j, double& x, double& y) const
{
  if(m_isGeoRef){
//     x=m_ulx+(1.5-m_magic_x+i)*m_delta_x;
//     y=m_uly-(1.5-m_magic_y+j)*m_delta_y;
    x=m_ulx+(0.5+i)*m_delta_x;
    y=m_uly-(0.5+j)*m_delta_y;
    return true;
  }
  else
    return false;
}

bool ImgWriterGdal::covers(double x, double  y) const
{
  double theULX, theULY, theLRX, theLRY;
  getBoundingBox(theULX,theULY,theLRX,theLRY);
  return((x > theULX)&&
	 (x < theLRX)&&
	 (y < theULY)&&
	 (y >theLRY));
}

bool ImgWriterGdal::covers(double ulx, double  uly, double lrx, double lry) const
{
  double theULX, theULY, theLRX, theLRY;
  getBoundingBox(theULX,theULY,theLRX,theLRY);
  return((ulx < theLRX)&&(lrx > theULX)&&(lry < theULY)&&(uly > theLRY));
}

string ImgWriterGdal::getGeoTransform() const
{
  double adfGeoTransform[6];// { 444720, 30, 0, 3751320, 0, -30 };
  double ulx;
  double deltaX;
  double rot1;
  double uly;
  double rot2;
  double deltaY;
  if(m_gds){
    m_gds->GetGeoTransform(adfGeoTransform);
    ulx=adfGeoTransform[0];
    deltaX=adfGeoTransform[1];
    rot1=adfGeoTransform[2];
    uly=adfGeoTransform[3];
    rot2=adfGeoTransform[4];
    deltaY=-adfGeoTransform[5];//convention of GDAL!
  }
  else{//virtual writer
    ulx=m_ulx;
    uly=m_uly;
    deltaX=m_delta_x;
    deltaY=m_delta_y;
    rot1=0;
    rot2=0;
  }
  ostringstream s;
  s << "[" << ulx << "," << deltaX << "," << rot1 << "," << uly << "," << rot2 << "," << -deltaY << "]";
  return(s.str());
}

void ImgWriterGdal::getGeoTransform(double& ulx, double& uly, double& deltaX, double& deltaY, double& rot1, double& rot2) const
{
  if(m_gds){
    double adfGeoTransform[6];// { 444720, 30, 0, 3751320, 0, -30 };
    m_gds->GetGeoTransform(adfGeoTransform);
    ulx=adfGeoTransform[0];
    deltaX=adfGeoTransform[1];
    rot1=adfGeoTransform[2];
    uly=adfGeoTransform[3];
    rot2=adfGeoTransform[4];
    deltaY=-adfGeoTransform[5];//convention of GDAL!
  }
  else{//virtual writer
    ulx=m_ulx;
    uly=m_uly;
    deltaX=m_delta_x;
    deltaY=m_delta_y;
    rot1=0;
    rot2=0;
  }
}

GDALDataType ImgWriterGdal::getDataType(int band) const
{
  assert(band<m_nband+1);
  return (m_gds->GetRasterBand(band+1))->GetRasterDataType();
}

GDALRasterBand* ImgWriterGdal::getRasterBand(int band)
{
  assert(band<m_nband+1);
  return (m_gds->GetRasterBand(band+1));
}

//filename is ascii file containing 5 columns: index R G B ALFA (0:transparent, 255:solid)
void ImgWriterGdal::setColorTable(const string& filename, int band)
{
  //todo: fool proof table in file (no checking currently done...)
  ifstream ftable(filename.c_str(),ios::in);
  string line;
//   poCT=new GDALColorTable();
  GDALColorTable colorTable;
  short nline=0;
  while(getline(ftable,line)){
    ++nline;
    istringstream ist(line);
    GDALColorEntry sEntry;
    short id;
    ist >> id >> sEntry.c1 >> sEntry.c2 >> sEntry.c3 >> sEntry.c4;
//     poCT->SetColorEntry(id,&sEntry);
    colorTable.SetColorEntry(id,&sEntry);
  }
  // assert(nline==colorTable.GetColorEntryCount());
//   (m_gds->GetRasterBand(band+1))->SetColorTable(poCT);
  (m_gds->GetRasterBand(band+1))->SetColorTable(&colorTable);
}

void ImgWriterGdal::setColorTable(GDALColorTable* colorTable, int band)
{
  (m_gds->GetRasterBand(band+1))->SetColorTable(colorTable);
}

//write an entire image from memory to file
bool ImgWriterGdal::writeData(void* pdata, const GDALDataType& dataType, int band) const{
  //fetch raster band
  GDALRasterBand  *poBand;
  if(band>=nrOfBand()+1){
    ostringstream s;
    s << "band (" << band << ") exceeds nrOfBand (" << nrOfBand() << ")";
    throw(s.str());
  }
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  poBand->RasterIO(GF_Write,0,0,nrOfCol(),nrOfRow(),pdata,nrOfCol(),nrOfRow(),dataType,0,0);
  return true;
}  
