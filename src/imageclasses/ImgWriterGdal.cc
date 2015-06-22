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
#include <algorithm>
#include "ogr_spatialref.h"
extern "C" {
#include "gdal_alg.h"
}
#include "ImgWriterGdal.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

//---------------------------------------------------------------------------
ImgWriterGdal::ImgWriterGdal(void)
  : m_gds(NULL), m_ncol(0), m_nrow(0), m_nband(0)
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
void ImgWriterGdal::open(const std::string& filename, const ImgReaderGdal& imgSrc, const std::vector<std::string>& options)
{
  // m_isGeoRef=imgSrc.isGeoRef();
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

// void ImgWriterGdal::open(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, const std::string& interleave, const std::string& compression, int magicX, int magicY)
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

void ImgWriterGdal::open(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, const std::vector<std::string>& options)
{
  // m_isGeoRef=false;
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
  for(std::vector<std::string>::const_iterator optionIt=m_options.begin();optionIt!=m_options.end();++optionIt)
    papszOptions=CSLAddString(papszOptions,optionIt->c_str());
  CSLDestroy(papszOptions);
}

//---------------------------------------------------------------------------
void ImgWriterGdal::setCodec(const ImgReaderGdal& imgSrc){
  GDALAllRegister();
  GDALDriver *poDriver;
  poDriver = GetGDALDriverManager()->GetDriverByName(imgSrc.getDriverDescription().c_str());
  if( poDriver == NULL ){
    std::string errorString="FileOpenError";
    throw(errorString);
  }
  char **papszMetadata;
  papszMetadata = poDriver->GetMetadata();
  //todo: try and catch if CREATE is not supported (as in PNG)
  assert( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE ));
  char **papszOptions=NULL;
  for(std::vector<std::string>::const_iterator optionIt=m_options.begin();optionIt!=m_options.end();++optionIt)
    papszOptions=CSLAddString(papszOptions,optionIt->c_str());
  // char **papszOptions=NULL;
  // std::ostringstream compressList;
  // compressList << "COMPRESS=" << m_compression;
  // papszOptions = CSLAddString(papszOptions,(compressList.str()).c_str());
  // std::ostringstream interleaveList;
  // interleaveList << "INTERLEAVE=" << m_interleave;
  // papszOptions = CSLAddString(papszOptions,(interleaveList.str()).c_str());
  m_gds=poDriver->Create(m_filename.c_str(),m_ncol,m_nrow,m_nband,m_type,papszOptions);
  // if(imgSrc.isGeoRef()){
    setProjection(imgSrc.getProjection());
    double gt[6];
    imgSrc.getGeoTransform(gt);
    setGeoTransform(gt);
  // }
  m_gds->SetMetadata(imgSrc.getMetadata() ); 
  m_gds->SetMetadataItem( "TIFFTAG_DOCUMENTNAME", m_filename.c_str());
  std::string versionString="pktools ";
  versionString+=VERSION;
  versionString+=" by Pieter Kempeneers";
  m_gds->SetMetadataItem( "TIFFTAG_SOFTWARE", versionString.c_str());
  time_t rawtime;
  time ( &rawtime );

  time_t tim=time(NULL);
  tm *now=localtime(&tim);
  std::ostringstream datestream;
  //date std::string must be 20 characters long...
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
//   list<std::string> lmeta;
//   imgReader.getMetadata(lmeta);
//   list<std::string>::const_iterator lit=lmeta.begin();
//   while(lit!=lmeta.end()){
//     cout << *lit << endl;
//     ++lit;
//   }
  // m_gds->SetMetadataItem( "INTERLEAVE", m_interleave.c_str(), "IMAGE_STRUCTURE" );
  // m_gds->SetMetadataItem( "COMPRESS", m_compression.c_str(), "IMAGE_STRUCTURE" );
  if(imgSrc.getColorTable()!=NULL)
    setColorTable(imgSrc.getColorTable());
}

void ImgWriterGdal::setCodec(const std::string& imageType)
{
  GDALAllRegister();
  GDALDriver *poDriver;
  poDriver = GetGDALDriverManager()->GetDriverByName(imageType.c_str());
  if( poDriver == NULL ){
    std::ostringstream s;
    s << "FileOpenError (" << imageType << ")";
    throw(s.str());
  }
  char **papszMetadata;
  papszMetadata = poDriver->GetMetadata();
  //todo: try and catch if CREATE is not supported (as in PNG)
  assert( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE ));
  char **papszOptions=NULL;
  for(std::vector<std::string>::const_iterator optionIt=m_options.begin();optionIt!=m_options.end();++optionIt)
    papszOptions=CSLAddString(papszOptions,optionIt->c_str());
  // std::ostringstream compressList;
  // compressList << "COMPRESS=" << m_compression;
  // papszOptions = CSLAddString(papszOptions,(compressList.str()).c_str());
  // std::ostringstream interleaveList;
  // interleaveList << "INTERLEAVE=" << m_interleave;
  // papszOptions = CSLAddString(papszOptions,(interleaveList.str()).c_str());
  m_gds=poDriver->Create(m_filename.c_str(),m_ncol,m_nrow,m_nband,m_type,papszOptions);

  // m_gds->SetMetadataItem( "INTERLEAVE", m_interleave.c_str(), "IMAGE_STRUCTURE" );
  // m_gds->SetMetadataItem( "COMPRESSION", m_compression.c_str(), "IMAGE_STRUCTURE" );
  m_gds->SetMetadataItem( "TIFFTAG_DOCUMENTNAME", m_filename.c_str());
  std::string versionString="pktools ";
  versionString+=VERSION;
  versionString+=" by Pieter Kempeneers";
  m_gds->SetMetadataItem( "TIFFTAG_SOFTWARE", versionString.c_str());
  time_t rawtime;
  time ( &rawtime );

  time_t tim=time(NULL);
  tm *now=localtime(&tim);
  std::ostringstream datestream;
  //date std::string must be 20 characters long...
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

std::string ImgWriterGdal::getProjection(void) const 
{
  assert(m_gds);
  std::string theProjection=m_gds->GetProjectionRef();
  //due to error in Gdal? AUTHORITY fields do not seem to work!
  // size_t startpos,endpos;
  // while((startpos=theProjection.find(",AUTHORITY"))!=string::npos){
  //   endpos=theProjection.find("]",startpos+1,1)+1;
  //   theProjection.erase(startpos,endpos-startpos);
  // }
  return theProjection;
}

//---------------------------------------------------------------------------
void ImgWriterGdal::setGeoTransform(double* gt){
  // m_isGeoRef=true;
  m_gt[0]=gt[0];
  m_gt[1]=gt[1];
  m_gt[2]=gt[2];
  m_gt[3]=gt[3];
  m_gt[4]=gt[4];
  m_gt[5]=gt[5];
  if(m_gds)
    m_gds->SetGeoTransform(m_gt);
}

// void ImgWriterGdal::setGeoTransform(double ulx, double uly, double deltaX, double deltaY, double rot1, double rot2)
// {
//   m_isGeoRef=true;
//   m_ulx=ulx;
//   m_uly=uly;
//   m_delta_x=deltaX;
//   m_delta_y=deltaY;
//   double adfGeoTransform[6];// { 444720, 30, 0, 3751320, 0, -30 };
//   adfGeoTransform[0]=ulx;
//   adfGeoTransform[1]=deltaX;
//   adfGeoTransform[2]=rot1;
//   adfGeoTransform[3]=uly;
//   adfGeoTransform[4]=rot2;
//   adfGeoTransform[5]=-deltaY;//convention of GDAL!
//   if(m_gds)
//     m_gds->SetGeoTransform(adfGeoTransform);
// }

void ImgWriterGdal::copyGeoTransform(const ImgReaderGdal& imgSrc)
{
  setProjection(imgSrc.getProjection());
  double gt[6];
  imgSrc.getGeoTransform(gt);
  setGeoTransform(gt);
  // imgSrc.getGeoTransform(ulx,uly,deltaX,deltaY,rot1,rot2);
  // setGeoTransform(ulx,uly,deltaX,deltaY,rot1,rot2);
}

std::string ImgWriterGdal::setProjectionProj4(const std::string& projection)
{
  // if(!m_isGeoRef)
  //   m_isGeoRef=true;

  OGRSpatialReference theRef;
  theRef.SetFromUserInput(projection.c_str());
  char *wktString;
  theRef.exportToWkt(&wktString);
  assert(m_gds);
  m_gds->SetProjection(wktString);
  return(wktString);

    // OGRSpatialReferenceH hSRS;  
    // char *pszResult = NULL;  
  
    // CPLErrorReset();  
      
    // hSRS = OSRNewSpatialReference( NULL );  
    // if( OSRSetFromUserInput( hSRS, projection.c_str() ) == OGRERR_NONE )  
    //     OSRExportToWkt( hSRS, &pszResult );  
    // else  
    // {  
    //     std::ostringstream s;
    //     s << "Error in set projection " << projection;
    //     throw(s.str());
    // }  
    // std::string theProjection=pszResult;
    // assert(m_gds);
    // m_gds->SetProjection(theProjection.c_str());
    // OSRDestroySpatialReference( hSRS );  
  
    // return theProjection;  
}

void ImgWriterGdal::setProjection(const std::string& projection)
{
  // if(!m_isGeoRef)
  //   m_isGeoRef=true;
  OGRSpatialReference oSRS;
  char *pszSRS_WKT = NULL;
  assert(m_gds);
  m_gds->SetProjection(projection.c_str());
  CPLFree(pszSRS_WKT);
}

//default projection: ETSR-LAEA
std::string ImgWriterGdal::setProjection(void)
{
  std::string theProjection;
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
  double gt[6];// { 444720, 30, 0, 3751320, 0, -30 };
  m_gds->GetGeoTransform(gt);

  //assuming
  //adfGeotransform[0]: ULX (upper left X coordinate)
  //adfGeotransform[1]: $cos(\alpha)\cdot\textrm{Xres}$
  //adfGeotransform[2]: $-sin(\alpha)\cdot\textrm{Xres}$
  //adfGeotransform[3]: ULY (upper left Y coordinate)
  //adfGeotransform[4]: $-sin(\alpha)\cdot\textrm{Yres}$
  //adfGeotransform[5]: $-cos(\alpha)\cdot\textrm{Yres}$

  ulx=gt[0];
  uly=gt[3];
  lrx=gt[0]+nrOfCol()*gt[1]+nrOfRow()*gt[2];
  lry=gt[3]+nrOfCol()*gt[4]+nrOfRow()*gt[5];
  if(isGeoRef()){
    // ulx=m_ulx;
    // uly=m_uly;
    // lrx=ulx+nrOfCol()*m_delta_x;
    // lry=uly-nrOfRow()*m_delta_y;
    return true;
  }
  else{
    // ulx=0;
    // uly=nrOfRow()-1;
    // lrx=nrOfCol()-1;
    // lry=0;
    return false;
  }
}

bool ImgWriterGdal::getCentrePos(double& x, double& y) const
{
  double gt[6];// { 444720, 30, 0, 3751320, 0, -30 };
  m_gds->GetGeoTransform(gt);

  //assuming
  //adfGeotransform[0]: ULX (upper left X coordinate)
  //adfGeotransform[1]: $cos(\alpha)\cdot\textrm{Xres}$
  //adfGeotransform[2]: $-sin(\alpha)\cdot\textrm{Xres}$
  //adfGeotransform[3]: ULY (upper left Y coordinate)
  //adfGeotransform[4]: $-sin(\alpha)\cdot\textrm{Yres}$
  //adfGeotransform[5]: $-cos(\alpha)\cdot\textrm{Yres}$
  x=gt[0]+(nrOfCol()/2.0)*gt[1]+(nrOfRow()/2.0)*gt[2];
  y=gt[3]+(nrOfCol()/2.0)*gt[4]+(nrOfRow()/2.0)*gt[5];
  if(isGeoRef()){
    // x=m_ulx+nrOfCol()/2.0*m_delta_x;
    // y=m_uly-nrOfRow()/2.0*m_delta_y;
    return true;
  }
  else
    return false;
}

bool ImgWriterGdal::geo2image(double x, double y, double& i, double& j) const
{
  //double values are returned, caller is responsible for interpolation step
  //double values are returned, caller is responsible for interpolation step
  double gt[6];// { 444720, 30, 0, 3751320, 0, -30 };
  m_gds->GetGeoTransform(gt);
  //assuming
  //adfGeotransform[0]: ULX (upper left X coordinate)
  //adfGeotransform[1]: $cos(\alpha)\cdot\textrm{Xres}$
  //adfGeotransform[2]: $-sin(\alpha)\cdot\textrm{Xres}$
  //adfGeotransform[3]: ULY (upper left Y coordinate)
  //adfGeotransform[4]: $-sin(\alpha)\cdot\textrm{Yres}$
  //adfGeotransform[5]: $-cos(\alpha)\cdot\textrm{Yres}$
  double denom=(gt[1]-gt[2]*gt[4]/gt[5]);
  double eps=0.00001;
  if(fabs(denom)>eps){
    i=(x-gt[0]-gt[2]/gt[5]*(y-gt[3]))/denom;
    j=(y-gt[3]-gt[4]*(x-gt[0]-gt[2]/gt[5]*(y-gt[3]))/denom)/gt[5];
  }
  if(isGeoRef()){
    // double ulx=m_ulx;
    // double uly=m_uly;
    // i=(x-ulx)/m_delta_x;
    // j=(uly-y)/m_delta_y;
    return true;
  }
  else{
    // i=x;
    // j=nrOfRow()-y;
    return false;
  }
}

//centre of pixel is always returned (regardless of magic pixel reference)!
bool ImgWriterGdal::image2geo(double i, double j, double& x, double& y) const
{
  double gt[6];// { 444720, 30, 0, 3751320, 0, -30 };
  m_gds->GetGeoTransform(gt);

  //assuming
  //adfGeotransform[0]: ULX (upper left X coordinate)
  //adfGeotransform[1]: $cos(\alpha)\cdot\textrm{Xres}$
  //adfGeotransform[2]: $-sin(\alpha)\cdot\textrm{Xres}$
  //adfGeotransform[3]: ULY (upper left Y coordinate)
  //adfGeotransform[4]: $-sin(\alpha)\cdot\textrm{Yres}$
  //adfGeotransform[5]: $-cos(\alpha)\cdot\textrm{Yres}$

  x=gt[0]+(0.5+i)*gt[1]+(0.5+j)*gt[2];
  y=gt[3]+(0.5+i)*gt[4]+(0.5+j)*gt[5];

  if(isGeoRef()){
    // x=m_ulx+(0.5+i)*m_delta_x;
    // y=m_uly-(0.5+j)*m_delta_y;
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

std::string ImgWriterGdal::getGeoTransform() const
{
  double gt[6];// { 444720, 30, 0, 3751320, 0, -30 };
  if(m_gds)
    m_gds->GetGeoTransform(gt);
  else{//virtual writer
    gt[0]=m_gt[0];
    gt[1]=m_gt[1];
    gt[2]=m_gt[2];
    gt[3]=m_gt[3];
    gt[4]=m_gt[4];
    gt[5]=m_gt[5];
  }
  std::ostringstream s;
  s << "[" << gt[0] << "," << gt[1] << "," << gt[2] << "," << gt[3] << "," << gt[4] << "," << gt[5] << "]";
  return(s.str());
}

void ImgWriterGdal::getGeoTransform(double* gt) const{
  if(m_gds)
    m_gds->GetGeoTransform(gt);
  else{//virtual writer
    gt[0]=m_gt[0];
    gt[1]=m_gt[1];
    gt[2]=m_gt[2];
    gt[3]=m_gt[3];
    gt[4]=m_gt[4];
    gt[5]=m_gt[5];
  }
}

// void ImgWriterGdal::getGeoTransform(double& ulx, double& uly, double& deltaX, double& deltaY, double& rot1, double& rot2) const
// {
//   if(m_gds){
//     double adfGeoTransform[6];// { 444720, 30, 0, 3751320, 0, -30 };
//     m_gds->GetGeoTransform(adfGeoTransform);
//     ulx=adfGeoTransform[0];
//     deltaX=adfGeoTransform[1];
//     rot1=adfGeoTransform[2];
//     uly=adfGeoTransform[3];
//     rot2=adfGeoTransform[4];
//     deltaY=-adfGeoTransform[5];//convention of GDAL!
//   }
//   else{//virtual writer
//     ulx=m_ulx;
//     uly=m_uly;
//     deltaX=m_delta_x;
//     deltaY=m_delta_y;
//     rot1=0;
//     rot2=0;
//   }
// }

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
void ImgWriterGdal::setColorTable(const std::string& filename, int band)
{
  //todo: fool proof table in file (no checking currently done...)
  std::ifstream ftable(filename.c_str(),std::ios::in);
  std::string line;
//   poCT=new GDALColorTable();
  GDALColorTable colorTable;
  short nline=0;
  while(getline(ftable,line)){
    ++nline;
    std::istringstream ist(line);
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
    std::ostringstream s;
    s << "band (" << band << ") exceeds nrOfBand (" << nrOfBand() << ")";
    throw(s.str());
  }
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  poBand->RasterIO(GF_Write,0,0,nrOfCol(),nrOfRow(),pdata,nrOfCol(),nrOfRow(),dataType,0,0);
  return true;
}  

void ImgWriterGdal::rasterizeOgr(ImgReaderOgr& ogrReader,const std::vector<std::string> layernames ){
  std::vector<int> bands;
  for(int iband=0;iband<nrOfBand();++iband)
    bands.push_back(iband+1);
  std::vector<OGRLayerH> layers;
  int nlayer=0;
  std::vector<double> burnValues;//todo: replace hard coded 1 with arg
  for(int ilayer=0;ilayer<ogrReader.getLayerCount();++ilayer){
    std::string currentLayername=ogrReader.getLayer(ilayer)->GetName();
    if(layernames.size())
      if(find(layernames.begin(),layernames.end(),currentLayername)==layernames.end())
	continue;
    std::cout << "processing layer " << currentLayername << std::endl;
    layers.push_back((OGRLayerH)ogrReader.getLayer(ilayer));
    ++nlayer;
    for(int iband=0;iband<nrOfBand();++iband)
      burnValues.push_back(1);
  }
  void *pTransformArg;
  char **papszOptions;
  double dfComplete=0.0;
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  pfnProgress(dfComplete,pszMessage,pProgressArg);
  if(GDALRasterizeLayers( (GDALDatasetH)m_gds,nrOfBand(),&(bands[0]),layers.size(),&(layers[0]),NULL,pTransformArg,&(burnValues[0]),papszOptions,pfnProgress,pProgressArg)!=CE_None){
    std::cerr << CPLGetLastErrorMsg() << std::endl;
    exit(1);
  }
  else{
    dfComplete=1.0;
    pfnProgress(dfComplete,pszMessage,pProgressArg);
  }
}
