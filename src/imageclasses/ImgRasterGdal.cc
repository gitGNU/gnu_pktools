/**********************************************************************
ImgRasterGdal.cc: class to read raster files using GDAL API library
Copyright (C) 2008-2016 Pieter Kempeneers

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
#include "ogr_spatialref.h"
extern "C" {
#include "gdal_alg.h"
}
#include <config.h>
#include "ImgRasterGdal.h"
#include "algorithms/StatFactory.h"

ImgRasterGdal::ImgRasterGdal(){
  reset();
}

void ImgRasterGdal::reset(void)
{
  m_gds=0;
  m_ncol=0;
  m_nrow=0;
  m_nband=0;
  m_dataType=GDT_Unknown;
  m_projection="";
  m_noDataValues.clear();
  m_scale.clear();
  m_offset.clear();
  m_options.clear();
  // m_writeMode=false;
  m_access=READ_ONLY;
#if GDAL_VERSION_MAJOR >= 2
  m_resample=GRIORA_NearestNeighbour;
#endif
  m_filename.clear();
  m_data.clear();
}

/**
 * @param memory Available memory to cache image raster data (in MB)
 **/
CPLErr ImgRasterGdal::initMem()
{
  freeMem();
  m_data.resize(nrOfBand());
  for(int iband=0;iband<m_nband;++iband){
    m_data[iband]=(void *) malloc(((GDALGetDataTypeSize(getDataType()))>>3)*nrOfCol()*m_nrow);
    if(!(m_data[iband])){
      std::string errorString="Error: could not allocate memory in initMem";
      throw(errorString);
    }
  }
  return(CE_None);
}

/**
   /**
   * @param memory Available memory to cache image raster data (in MB)
   **/
void ImgRasterGdal::freeMem()
{
  for(int iband=0;iband<m_data.size();++iband){
    free(m_data[iband]);
  }
  m_data.clear();
}

/**
 * @param imgSrc Use this source image as a template to copy image attributes
 **/
ImgRasterGdal& ImgRasterGdal::operator=(ImgRasterGdal& imgSrc)
{
  //check for assignment to self (of the form v=v)
  if(this==&imgSrc)
     return *this;
  else{
    open(imgSrc);
    return *this;
  }
}

void ImgRasterGdal::close(void)
{
  // if(writeMode()){
  if(writeMode()||updateMode()){
    char **papszOptions=NULL;
    for(std::vector<std::string>::const_iterator optionIt=m_options.begin();optionIt!=m_options.end();++optionIt)
      papszOptions=CSLAddString(papszOptions,optionIt->c_str());
    if(papszOptions)
      CSLDestroy(papszOptions);
  }
  if(m_gds)
    GDALClose(m_gds);
  reset();
}

/**
 * @return the projection of this data set in string format
 **/
std::string ImgRasterGdal::getProjection(void) const
{
  // if(m_gds)
  //   return(m_gds->GetProjectionRef());
  // else
  return(m_projection);
}

/**
 * @return the projection of this data set in string format
 **/
std::string ImgRasterGdal::getProjectionRef(void) const
{
  // if(m_gds)
  //   return(m_gds->GetProjectionRef());
  // else
  return(m_projection);
}

/**
 * @param projection projection string to be used for this dataset
 * @return the projection of this data set in string format
 **/
CPLErr ImgRasterGdal::setProjectionProj4(const std::string& projection)
{
  OGRSpatialReference theRef;
  theRef.SetFromUserInput(projection.c_str());
  char *wktString;
  theRef.exportToWkt(&wktString);
  m_projection=wktString;
  if(m_gds)
    return(m_gds->SetProjection(wktString));
  else
    return(CE_Failure);
}

/**
 * @param projection projection string to be used for this dataset
 **/
CPLErr ImgRasterGdal::setProjection(const std::string& projection)
{
  m_projection=projection;
  if(m_gds){
    return(m_gds->SetProjection(projection.c_str()));
  }
  else{
    return(CE_Failure);
  }
}

/**
 * @param band get data type for this band (start counting from 0)
 * @return the GDAL data type of this data set for the selected band
 **/
GDALDataType ImgRasterGdal::getDataType(int band) const
{
  if(nrOfBand()<=band){
    std::string errorString="Error: band number exceeds available bands in getDataType";
    throw(errorString);
  }
  if(getRasterBand(band))
    return((getRasterBand(band)->GetRasterDataType()));
  else
    return(m_dataType);
}

/**
 * @param band get GDAL raster band for this band (start counting from 0)
 * @return the GDAL raster band of this data set for the selected band
 **/
GDALRasterBand* ImgRasterGdal::getRasterBand(int band) const
{
  if(nrOfBand()<=band){
    std::string errorString="Error: band number exceeds available bands in getRasterBand";
    throw(errorString);
  }
  if(m_gds)
    return((m_gds->GetRasterBand(band+1)));
  else
    return(0);
}

/**
 * @param band get GDAL color table for this band (start counting from 0)
 * @return the GDAL color table of this data set for the selected band
 **/
GDALColorTable* ImgRasterGdal::getColorTable(int band) const
{
  if(nrOfBand()<=band){
    std::string errorString="Error: band number exceeds available bands in getColorTable";
    throw(errorString);
  }
  GDALRasterBand* theRasterBand=getRasterBand(band);
  if(theRasterBand)
    return(theRasterBand->GetColorTable());
  else
    return(0);
}

/**
 * @return the driver description of this data set in string format
 **/
std::string ImgRasterGdal::getDriverDescription() const
{
  std::string driverDescription;
  if(m_gds)
    driverDescription=m_gds->GetDriver()->GetDescription();
  return(driverDescription);
}

/**
 * @param gt pointer to the six geotransform parameters:
 * @param adfGeoTransform[0] top left x
 * @param GeoTransform[1] w-e pixel resolution
 * @param GeoTransform[2] rotation, 0 if image is "north up"
 * @param GeoTransform[3] top left y
 * @param GeoTransform[4] rotation, 0 if image is "north up"
 * @param GeoTransform[5] n-s pixel resolution
 **/
CPLErr ImgRasterGdal::setGeoTransform(double* gt){
  // m_isGeoRef=true;
  m_gt[0]=gt[0];
  m_gt[1]=gt[1];
  m_gt[2]=gt[2];
  m_gt[3]=gt[3];
  m_gt[4]=gt[4];
  m_gt[5]=gt[5];
  if(m_gds&&m_access==WRITE)
    return(m_gds->SetGeoTransform(m_gt));
  else
    return(CE_Failure);
}

/**
 * @param imgSrc Use this source image as a template to copy geotranform information
 **/
void ImgRasterGdal::copyGeoTransform(const ImgRasterGdal& imgSrc)
{
  double gt[6];
  imgSrc.getGeoTransform(gt);
  setGeoTransform(gt);
}

/**
 * @param imgSrc Use this pointer to source image as a template to copy geotranform information
 **/
// void ImgRasterGdal::copyGeoTransform(const std::shared_ptr<ImgRasterGdal>& imgSrc)
// {
//   double gt[6];
//   imgSrc->getGeoTransform(gt);
//   setGeoTransform(gt);
// }

/**
 * @param gt pointer to the six geotransform parameters:
 * @param adfGeoTransform[0] top left x
 * @param GeoTransform[1] w-e pixel resolution
 * @param GeoTransform[2] rotation, 0 if image is "north up"
 * @param GeoTransform[3] top left y
 * @param GeoTransform[4] rotation, 0 if image is "north up"
 * @param GeoTransform[5] n-s pixel resolution
 **/
void ImgRasterGdal::getGeoTransform(double* gt) const{
  // if(m_gds){
  //   m_gds->GetGeoTransform(gt);
  // }
  // else{
    gt[0]=m_gt[0];
    gt[1]=m_gt[1];
    gt[2]=m_gt[2];
    gt[3]=m_gt[3];
    gt[4]=m_gt[4];
    gt[5]=m_gt[5];
  // }
}

/**
 * @return the geotransform of this data set in string format
 **/
std::string ImgRasterGdal::getGeoTransform() const
{
  std::string gtString;
  double gt[6];// { 444720, 30, 0, 3751320, 0, -30 };
  getGeoTransform(gt);
  std::ostringstream s;
  s << "[" << gt[0] << "," << gt[1] << "," << gt[2] << "," << gt[3] << "," << gt[4] << "," << gt[5] << "]";
  gtString=s.str();
  return(s.str());
}

/**
 * @return the metadata of this data set in C style string format (const version)
 **/
char** ImgRasterGdal::getMetadata() const
{
  if(m_gds){
    if(m_gds->GetMetadata()!=NULL)
      return(m_gds->GetMetadata());
    else
      return(0);
  }
  else
    return(0);
    // return (char**)"";
}

/**
 * @return the metadata of this data set in standard template library (stl) string format
 **/
void ImgRasterGdal::getMetadata(std::list<std::string>& metadata) const
{
  if(m_gds){
    char** cmetadata=m_gds->GetMetadata();
    while(*cmetadata!=NULL){
      metadata.push_back(*(cmetadata));
      ++cmetadata;
    }
  }
}

/**
 * @return the description of this data set in string format
 **/
std::string ImgRasterGdal::getDescription() const
{
  if(m_gds){
    if(m_gds->GetDriver()->GetDescription()!=NULL)
      return m_gds->GetDriver()->GetDescription();
    else
      return("");
  }
  else
    return("");
}

/**
 * @return the meta data item of this data set in string format
 **/
std::string ImgRasterGdal::getMetadataItem() const
{
  if(m_gds){
    if(m_gds->GetDriver()->GetMetadataItem( GDAL_DMD_LONGNAME )!=NULL)
      return m_gds->GetDriver()->GetMetadataItem( GDAL_DMD_LONGNAME );
    return("");
  }
  else
    return("");
}

/**
 * @return the image description (TIFFTAG) of this data set in string format
 **/
std::string ImgRasterGdal::getImageDescription() const
{
  if(m_gds){
    if(m_gds->GetDriver()->GetMetadataItem("TIFFTAG_IMAGEDESCRIPTION")!=NULL)
      return m_gds->GetDriver()->GetMetadataItem("TIFFTAG_IMAGEDESCRIPTION");
    return("");
  }
  else
    return("");
}

/**
 * @return the band coding interleave of this data set in string format
 **/
std::string ImgRasterGdal::getInterleave() const
{
  if(m_gds){
    if(m_gds->GetMetadataItem( "INTERLEAVE", "IMAGE_STRUCTURE"))
      return m_gds->GetMetadataItem( "INTERLEAVE", "IMAGE_STRUCTURE");
    else
      return("BAND");
  }
  else
    return("");
}

/**
 * @return the compression meta data of this data set in string format
 **/
std::string ImgRasterGdal::getCompression() const
{
  if(m_gds){
    if(m_gds->GetMetadataItem( "COMPRESSION", "IMAGE_STRUCTURE"))
      return m_gds->GetMetadataItem( "COMPRESSION", "IMAGE_STRUCTURE");
    return("NONE");
  }
  else
    return("NONE");
}

/**
 * assuming
 * adfGeotransform[0]: ULX (upper left X coordinate)
 * adfGeotransform[1]: $cos(\alpha)\cdot\textrm{Xres}$
 * adfGeotransform[2]: $-sin(\alpha)\cdot\textrm{Xres}$
 * adfGeotransform[3]: ULY (upper left Y coordinate)
 * adfGeotransform[4]: $-sin(\alpha)\cdot\textrm{Yres}$
 * adfGeotransform[5]: $-cos(\alpha)\cdot\textrm{Yres}$

 * @param ulx upper left coordinate in x
 * @param uly upper left coordinate in y
 * @param lrx lower left coordinate in x
 * @param lry lower left coordinate in y
 * @return true if image is georeferenced
 **/
bool ImgRasterGdal::getBoundingBox(double& ulx, double& uly, double& lrx, double& lry) const
{
  double gt[6];// { 444720, 30, 0, 3751320, 0, -30 };
  getGeoTransform(gt);

  ulx=gt[0];
  uly=gt[3];
  lrx=gt[0]+nrOfCol()*gt[1]+nrOfRow()*gt[2];
  lry=gt[3]+nrOfCol()*gt[4]+nrOfRow()*gt[5];
  if(isGeoRef())
    return true;
  else
    return false;
}

/**
 * assuming
 * adfGeotransform[0]: ULX (upper left X coordinate)
 * adfGeotransform[1]: $cos(\alpha)\cdot\textrm{Xres}$
 * adfGeotransform[2]: $-sin(\alpha)\cdot\textrm{Xres}$
 * adfGeotransform[3]: ULY (upper left Y coordinate)
 * adfGeotransform[4]: $-sin(\alpha)\cdot\textrm{Yres}$
 * adfGeotransform[5]: $-cos(\alpha)\cdot\textrm{Yres}$
 * @param x, y centre coordinates in x and y
 * @return true if image is georeferenced
 **/
bool ImgRasterGdal::getCenterPos(double& x, double& y) const
{
  double gt[6];// { 444720, 30, 0, 3751320, 0, -30 };
  getGeoTransform(gt);

  x=gt[0]+(nrOfCol()/2.0)*gt[1]+(nrOfRow()/2.0)*gt[2];
  y=gt[3]+(nrOfCol()/2.0)*gt[4]+(nrOfRow()/2.0)*gt[5];
  if(isGeoRef()){
    // x=m_ulx+(nrOfCol()/2.0)*m_delta_x;
    // y=m_uly-(nrOfRow()/2.0)*m_delta_y;
    return true;
  }
  else
    return false;
}

/**
 * assuming
 * adfGeotransform[0]: ULX (upper left X coordinate)
 * adfGeotransform[1]: $cos(\alpha)\cdot\textrm{Xres}$
 * adfGeotransform[2]: $-sin(\alpha)\cdot\textrm{Xres}$
 * adfGeotransform[3]: ULY (upper left Y coordinate)
 * adfGeotransform[4]: $-sin(\alpha)\cdot\textrm{Yres}$
 * adfGeotransform[5]: $-cos(\alpha)\cdot\textrm{Yres}$
 * @param x,y georeferenced coordinates in x and y
 * @param i,j image coordinates (can be fraction of pixels)
 * @return true if image is georeferenced
 **/
bool ImgRasterGdal::geo2image(double x, double y, double& i, double& j) const
{
  double gt[6];// { 444720, 30, 0, 3751320, 0, -30 };
  getGeoTransform(gt);

  double denom=(gt[1]-gt[2]*gt[4]/gt[5]);
  double eps=0.00001;
  if(fabs(denom)>eps){
    i=(x-gt[0]-gt[2]/gt[5]*(y-gt[3]))/denom;
    j=(y-gt[3]-gt[4]*(x-gt[0]-gt[2]/gt[5]*(y-gt[3]))/denom)/gt[5];
  }
  if(isGeoRef())
    return true;
  else
    return false;
}

/**
 * assuming
 * adfGeotransform[0]: ULX (upper left X coordinate)
 * adfGeotransform[1]: $cos(\alpha)\cdot\textrm{Xres}$
 * adfGeotransform[2]: $-sin(\alpha)\cdot\textrm{Xres}$
 * adfGeotransform[3]: ULY (upper left Y coordinate)
 * adfGeotransform[4]: $-sin(\alpha)\cdot\textrm{Yres}$
 * adfGeotransform[5]: $-cos(\alpha)\cdot\textrm{Yres}$
 * @param i,j image coordinates (can be fraction of pixels)
 * @param x,y georeferenced coordinates in x and y (can be fraction of pixels)
 * @return true if image is georeferenced
 **/
bool ImgRasterGdal::image2geo(double i, double j, double& x, double& y) const
{
  double gt[6];// { 444720, 30, 0, 3751320, 0, -30 };
  getGeoTransform(gt);

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

/**
 * assuming
 * adfGeotransform[0]: ULX (upper left X coordinate)
 * adfGeotransform[1]: $cos(\alpha)\cdot\textrm{Xres}$
 * adfGeotransform[2]: $-sin(\alpha)\cdot\textrm{Xres}$
 * adfGeotransform[3]: ULY (upper left Y coordinate)
 * adfGeotransform[4]: $-sin(\alpha)\cdot\textrm{Yres}$
 * adfGeotransform[5]: $-cos(\alpha)\cdot\textrm{Yres}$
 * @param x,y georeferenced coordinates in x and y
 * @return true if image covers the georeferenced location
 **/
bool ImgRasterGdal::covers(double x, double  y) const
{
  double theULX, theULY, theLRX, theLRY;
  getBoundingBox(theULX,theULY,theLRX,theLRY);
  return((x > theULX)&&
         (x < theLRX)&&
         (y < theULY)&&
         (y >theLRY));
}

/**
 * assuming
 * adfGeotransform[0]: ULX (upper left X coordinate)
 * adfGeotransform[1]: $cos(\alpha)\cdot\textrm{Xres}$
 * adfGeotransform[2]: $-sin(\alpha)\cdot\textrm{Xres}$
 * adfGeotransform[3]: ULY (upper left Y coordinate)
 * adfGeotransform[4]: $-sin(\alpha)\cdot\textrm{Yres}$
 * adfGeotransform[5]: $-cos(\alpha)\cdot\textrm{Yres}$
 * @param ulx upper left coordinate in x
 * @param uly upper left coordinate in y
 * @param lrx lower left coordinate in x
 * @param lry lower left coordinate in y
 * @return true if image (partially) covers the bounding box
 **/
bool ImgRasterGdal::covers(double ulx, double  uly, double lrx, double lry) const
{
  double theULX, theULY, theLRX, theLRY;
  getBoundingBox(theULX,theULY,theLRX,theLRY);
  return((ulx < theLRX)&&(lrx > theULX)&&(lry < theULY)&&(uly > theLRY));
}

/**
 * @param noDataValues standard template library (stl) vector containing no data values
 * @return number of no data values in this dataset
 **/
int ImgRasterGdal::getNoDataValues(std::vector<double>& noDataValues) const
{
  if(m_noDataValues.size()){
    noDataValues=m_noDataValues;
    return m_noDataValues.size();
  }
  else
    return 0;
}

/**
 * @param noDataValue no data value to be pushed for this dataset
 * @return number of no data values in this dataset
 **/
int ImgRasterGdal::pushNoDataValue(double noDataValue)
{
  if(find(m_noDataValues.begin(),m_noDataValues.end(),noDataValue)==m_noDataValues.end())
    m_noDataValues.push_back(noDataValue);
  return(m_noDataValues.size());
}

//From Reader
/**
 * @param filename Open a raster dataset with this filename
 **/
CPLErr ImgRasterGdal::open(const std::string& filename)
// void ImgRasterGdal::open(const std::string& filename, const GDALAccess& readMode)
{
  // m_writeMode=false;
  m_access=READ_ONLY;
  m_filename = filename;
  registerDriver();
  return(CE_None);
}

/**
 **/
void ImgRasterGdal::registerDriver()
{
  if(writeMode()){
    GDALAllRegister();
    GDALDriver *poDriver;
    poDriver = GetGDALDriverManager()->GetDriverByName(m_imageType.c_str());
    if( poDriver == NULL ){
      std::ostringstream s;
      s << "FileOpenError (" << m_imageType << ")";
      throw(s.str());
    }
    char **papszMetadata;
    papszMetadata = poDriver->GetMetadata();
    //todo: try and catch if CREATE is not supported (as in PNG)
    if( ! CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE )){
      std::ostringstream s;
      s << "Error: image type " << m_imageType << " not supported";
      throw(s.str());
    }
    char **papszOptions=NULL;
    for(std::vector<std::string>::const_iterator optionIt=m_options.begin();optionIt!=m_options.end();++optionIt)
      papszOptions=CSLAddString(papszOptions,optionIt->c_str());

    m_gds=poDriver->Create(m_filename.c_str(),m_ncol,m_nrow,m_nband,m_dataType,papszOptions);
    double gt[6];
    getGeoTransform(gt);
    if(setGeoTransform(gt)!=CE_None)
      std::cerr << "Warning: could not write geotransform information in " << m_filename << std::endl;
    if(setProjection(m_projection)!=CE_None)
      std::cerr << "Warning: could not write projection information in " << m_filename << std::endl;


    if(m_noDataValues.size()){
      for(int iband=0;iband<nrOfBand();++iband)
        GDALSetNoDataValue(m_noDataValues[0],iband);
    }

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
  }
  else{
    GDALAllRegister();
    // m_gds = (GDALDataset *) GDALOpen(m_filename.c_str(), readMode );
#if GDAL_VERSION_MAJOR < 2
    GDALAllRegister();
    if(m_access==UPDATE)
      m_gds = (GDALDataset *) GDALOpen(m_filename.c_str(), GA_Update);
    else
      m_gds = (GDALDataset *) GDALOpen(m_filename.c_str(), GA_ReadOnly );
    // m_gds = (GDALDataset *) GDALOpen(m_filename.c_str(), readMode );
#else
    GDALAllRegister();
    // if(readMode==GA_ReadOnly)
    if(m_access==UPDATE)
      m_gds = (GDALDataset*) GDALOpenEx(m_filename.c_str(), GDAL_OF_UPDATE|GDAL_OF_RASTER, NULL, NULL, NULL);
    else
      m_gds = (GDALDataset*) GDALOpenEx(m_filename.c_str(), GDAL_OF_READONLY|GDAL_OF_RASTER, NULL, NULL, NULL);
    // m_gds = (GDALDataset*) GDALOpenEx(m_filename.c_str(), GDAL_OF_READONLY|GDAL_OF_RASTER, NULL, NULL, NULL);
    // else if(readMode==GA_Update)
    //   m_gds = (GDALDataset*) GDALOpenEx(m_filename.c_str(), GDAL_OF_UPDATE|GDAL_OF_RASTER, NULL, NULL, NULL);
#endif

    if(m_gds == NULL){
      std::string errorString="FileOpenError";
      throw(errorString);
    }
    m_ncol= m_gds->GetRasterXSize();
    m_nrow= m_gds->GetRasterYSize();
    m_nband= m_gds->GetRasterCount();
    m_dataType=getDataType();
    m_imageType=getImageType();
    double adfGeoTransform[6];
    m_gds->GetGeoTransform( adfGeoTransform );
    m_gt[0]=adfGeoTransform[0];
    m_gt[1]=adfGeoTransform[1];
    m_gt[2]=adfGeoTransform[2];
    m_gt[3]=adfGeoTransform[3];
    m_gt[4]=adfGeoTransform[4];
    m_gt[5]=adfGeoTransform[5];
    m_projection=m_gds->GetProjectionRef();
  }
}

/**
 * @param app application options
 **/
ImgRasterGdal::ImgRasterGdal(const app::AppFactory &app) {
  ImgRasterGdal();
  //input
  Optionpk<std::string> input_opt("i", "input", "input filename");
  Optionpk<std::string> resample_opt("r", "r", "resample: GRIORA_NearestNeighbour|GRIORA_Bilinear|GRIORA_Cubic|GRIORA_CubicSpline|GRIORA_Lanczos|GRIORA_Average|GRIORA_Average|GRIORA_Gauss (check http://www.gdal.org/gdal_8h.html#a640ada511cbddeefac67c548e009d5a)","GRIORA_NearestNeighbour");
  Optionpk<std::string> extra_opt("extra", "extra", "RGDALRasterIOExtraArg (check http://www.gdal.org/structGDALRasterIOExtraArg.html)");
  // Optionpk<std::string> targetSRS_opt("t_srs", "t_srs", "Target spatial reference system in EPSG format (e.g., epsg:3035)");//todo
  //output
  Optionpk<int> nsample_opt("ns", "nsample", "Number of samples");
  Optionpk<int> nline_opt("nl", "nline", "Number of lines");
  Optionpk<int> band_opt("b", "band", "Number of bands");
  Optionpk<std::string> otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64})","Byte");
  Optionpk<double> nodata_opt("nodata", "nodata", "Nodata value to put in image if out of bounds.");
  Optionpk<unsigned long int> seed_opt("seed", "seed", "seed value for random generator",0);
  Optionpk<double> mean_opt("mean", "mean", "Mean value for random generator",0);
  Optionpk<double> sigma_opt("sigma", "sigma", "Sigma value for random generator",0);
  Optionpk<std::string> assignSRS_opt("a_srs", "a_srs", "Assign the spatial reference for the output file, e.g., psg:3035 to use European projection and force to European grid");
  Optionpk<std::string> description_opt("d", "description", "Set image description");
  //input and output
  Optionpk<double> ulx_opt("ulx", "ulx", "Upper left x value bounding box");
  Optionpk<double> uly_opt("uly", "uly", "Upper left y value bounding box");
  Optionpk<double> lrx_opt("lrx", "lrx", "Lower right x value bounding box");
  Optionpk<double> lry_opt("lry", "lry", "Lower right y value bounding box");
  Optionpk<double> dx_opt("dx", "dx", "Resolution in x");
  Optionpk<double> dy_opt("dy", "dy", "Resolution in y");
  Optionpk<std::string> access_opt("access", "access", "access (READ_ONLY, UPDATE)","READ_ONLY",2);//todo

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  doProcess=input_opt.retrieveOption(app.getArgc(),app.getArgv());
  nsample_opt.retrieveOption(app.getArgc(),app.getArgv());
  nline_opt.retrieveOption(app.getArgc(),app.getArgv());
  band_opt.retrieveOption(app.getArgc(),app.getArgv());
  ulx_opt.retrieveOption(app.getArgc(),app.getArgv());
  uly_opt.retrieveOption(app.getArgc(),app.getArgv());
  lrx_opt.retrieveOption(app.getArgc(),app.getArgv());
  lry_opt.retrieveOption(app.getArgc(),app.getArgv());
  dx_opt.retrieveOption(app.getArgc(),app.getArgv());
  dy_opt.retrieveOption(app.getArgc(),app.getArgv());
  resample_opt.retrieveOption(app.getArgc(),app.getArgv());
  extra_opt.retrieveOption(app.getArgc(),app.getArgv());
  otype_opt.retrieveOption(app.getArgc(),app.getArgv());
  nodata_opt.retrieveOption(app.getArgc(),app.getArgv());
  seed_opt.retrieveOption(app.getArgc(),app.getArgv());
  mean_opt.retrieveOption(app.getArgc(),app.getArgv());
  sigma_opt.retrieveOption(app.getArgc(),app.getArgv());
  description_opt.retrieveOption(app.getArgc(),app.getArgv());
  assignSRS_opt.retrieveOption(app.getArgc(),app.getArgv());
  // targetSRS_opt.retrieveOption(app.getArgc(),app.getArgv());
  access_opt.retrieveOption(app.getArgc(),app.getArgv());
  if(!doProcess){
    std::cout << std::endl;
    std::ostringstream helpStream;
    helpStream << "help info: ";
    throw(helpStream.str());//help was invoked, stop processing
  }
  if(app.empty()){
    ImgRasterGdal();
  }
  else if(input_opt.empty()){
    if(nsample_opt.empty()){
      std::ostringstream errorStream;
      errorStream << "Warning: no number of samples (use option -ns). Returning empty image" << std::endl;
      ImgRasterGdal();
      throw(errorStream.str());
    }
    if(nline_opt.empty()){
      std::ostringstream errorStream;
      errorStream << "Warning: no number of lines (use option -nl). Returning empty image" << std::endl;
      ImgRasterGdal();
      throw(errorStream.str());
    }
    GDALDataType theType=getGDALDataType(otype_opt[0]);
    open(nsample_opt[0],nline_opt[0],band_opt[0],theType);
    setNoData(nodata_opt);
    if(description_opt.size())
      setImageDescription(description_opt[0]);
    double gt[6];
    if(ulx_opt[0]<lrx_opt[0])
      gt[0]=ulx_opt[0];
    else
      gt[0]=0;
    if(dx_opt.size())
      gt[1]=dx_opt[0];
    else if(lrx_opt[0]>0){
      gt[1]=lrx_opt[0]-ulx_opt[0];
      gt[1]/=nrOfCol();
    }
    else
      gt[1]=1;
    gt[2]=0;
    if(uly_opt[0]>lry_opt[0])
      gt[3]=uly_opt[0];
    else
      gt[3]=0;
    gt[4]=0;
    if(dy_opt.size())
      gt[5]=-dy_opt[0];
    else if(lry_opt[0]>0){
      gt[5]=lry_opt[0]-uly_opt[0];
      gt[5]/=nrOfRow();
    }
    else
      gt[5]=1;
    setGeoTransform(gt);
    if(assignSRS_opt.size())
      setProjectionProj4(assignSRS_opt[0]);
    statfactory::StatFactory stat;
    gsl_rng* rndgen=stat.getRandomGenerator(seed_opt[0]);
    std::vector<double> lineBuffer(nrOfCol());
    double value=stat.getRandomValue(rndgen,"gaussian",mean_opt[0],sigma_opt[0]);
    for(unsigned int iband=0;iband<nrOfBand();++iband){
      for(unsigned int irow=0;irow<nrOfRow();++irow){
        for(unsigned int icol=0;icol<nrOfCol();++icol){
          if(sigma_opt[0]>0||(!irow&&!iband)){
            value=stat.getRandomValue(rndgen,"gaussian",mean_opt[0],sigma_opt[0]);
            lineBuffer[icol]=value;
          }
        }
        writeData(lineBuffer,irow,iband);
      }
    }
  }
  else{
    setAccess(access_opt[0]);
    m_filename=input_opt[0];
    registerDriver();
    if(band_opt.empty()){
      while(band_opt.size()<nrOfBand())
        band_opt.push_back(band_opt.size());
    }
    if(ulx_opt.empty())
      ulx_opt.push_back(getUlx());
    if(uly_opt.empty())
      uly_opt.push_back(getUly());
    if(lrx_opt.empty())
      lrx_opt.push_back(getLrx());
    if(lry_opt.empty())
      lry_opt.push_back(getLry());
    if(dx_opt.empty())
      dx_opt.push_back(getDeltaX());
    if(dy_opt.empty())
      dy_opt.push_back(getDeltaY());

    //force bounding box to be within dataset
    if(ulx_opt[0]<getUlx())
      ulx_opt[0]=getUlx();
    if(uly_opt[0]>getUly())
      uly_opt[0]=getUly();
    if(lrx_opt[0]>getLrx())
      lrx_opt[0]=getLrx();
    if(lry_opt[0]<getLry())
      lry_opt[0]=getLry();

    //todo: reproject on the fly using
    // OGRSpatialReference::SetFromUserInput

#if GDAL_VERSION_MAJOR >= 2
    m_resample=getGDALResample(resample_opt[0]);
#endif

    double gt[6];
    gt[0]=ulx_opt[0];
    gt[3]=uly_opt[0];
    gt[1]=dx_opt[0];//todo: adfGeotransform[1]: $cos(\alpha)\cdot\textrm{Xres}$
    gt[2]=0;//todo: $-sin(\alpha)\cdot\textrm{Xres}$
    gt[4]=0;//todo: $-sin(\alpha)\cdot\textrm{Yres}$
    gt[5]=-dy_opt[0];//todo: a$-cos(\alpha)\cdot\textrm{Yres}
    setGeoTransform(gt);

    int nBufXSize=abs(static_cast<unsigned int>(ceil((lrx_opt[0]-ulx_opt[0])/dx_opt[0])));
    int nBufYSize=abs(static_cast<unsigned int>(ceil((uly_opt[0]-lry_opt[0])/dy_opt[0])));
    m_ncol=nBufXSize;
    m_nrow=nBufYSize;
    //todo: support user defined selection of bands
    // m_nband=band_opt.size();
  }
}

/**
 * @param x Reported column where minimum value in image was found (start counting from 0)
 * @param y Reported row where minimum value in image was found (start counting from 0)
 * @param band Search mininum value in image for this band
 * @return minimum value in image for the selected band
 **/
double ImgRasterGdal::getMin(int& x, int& y, int band){
  double minValue=0;
  std::vector<double> lineBuffer(nrOfCol());
  bool isValid=false;
  for(int irow=0;irow<nrOfRow();++irow){
    readData(lineBuffer,irow,band);
    for(int icol=0;icol<nrOfCol();++icol){
      if(isNoData(lineBuffer[icol]))
  continue;
      if(isValid){
  if(lineBuffer[icol]<minValue){
          y=irow;
          x=icol;
          minValue=lineBuffer[icol];
        }
      }
      else{
  y=irow;
  x=icol;
  minValue=lineBuffer[icol];
  isValid=true;
      }
    }
  }
  if(isValid)
    return minValue;
  else
    throw(static_cast<std::string>("Warning: not initialized"));
}

/**
 * @param x Reported column where maximum value in image was found (start counting from 0)
 * @param y Reported row where maximum value in image was found (start counting from 0)
 * @param band Search mininum value in image for this band
 * @return maximum value in image for the selected band
 **/
double ImgRasterGdal::getMax(int& x, int& y, int band){
  double maxValue=0;
  std::vector<double> lineBuffer(nrOfCol());
  bool isValid=false;
  for(int irow=0;irow<nrOfRow();++irow){
    readData(lineBuffer,irow,band);
    for(int icol=0;icol<nrOfCol();++icol){
      if(isNoData(lineBuffer[icol]))
  continue;
      if(isValid){
  if(lineBuffer[icol]>maxValue){
          y=irow;
          x=icol;
          maxValue=lineBuffer[icol];
        }
      }
      else{
  y=irow;
  x=icol;
  maxValue=lineBuffer[icol];
  isValid=true;
      }
    }
  }
  if(isValid)
    return maxValue;
  else
    throw(static_cast<std::string>("Warning: not initialized"));
}

/**
 * @param startCol, endCol, startRow, endRow Search extreme value in this region of interest (all indices start counting from 0)
 * @param band Search extreme value in image for this band
 * @param minValue Reported minimum value within searched region
 * @param maxValue Reported maximum value within searched region
 **/
void ImgRasterGdal::getMinMax(int startCol, int endCol, int startRow, int endRow, int band, double& minValue, double& maxValue)
{
  bool isConstraint=(maxValue>minValue);
  double minConstraint=minValue;
  double maxConstraint=maxValue;
  std::vector<double> lineBuffer(endCol-startCol+1);
  bool isValid=false;
  //todo: replace assert with exception
  assert(endRow<nrOfRow());
  for(int irow=startCol;irow<endRow+1;++irow){
    readData(lineBuffer,startCol,endCol,irow,band);
    for(int icol=0;icol<lineBuffer.size();++icol){
      if(isNoData(lineBuffer[icol]))
  continue;
      if(isValid){
  if(isConstraint){
    if(lineBuffer[icol]<minConstraint)
      continue;
    if(lineBuffer[icol]>maxConstraint)
      continue;
  }
  if(lineBuffer[icol]<minValue)
    minValue=lineBuffer[icol];
  if(lineBuffer[icol]>maxValue)
    maxValue=lineBuffer[icol];
      }
      else{
  if(isConstraint){
    if(lineBuffer[icol]<minConstraint)
      continue;
    if(lineBuffer[icol]>maxConstraint)
      continue;
  }
  minValue=lineBuffer[icol];
  maxValue=lineBuffer[icol];
  isValid=true;
      }
    }
  }
  if(!isValid)
    throw(static_cast<std::string>("Warning: not initialized"));
}

/**
 * @param minValue Reported minimum value in image
 * @param maxValue Reported maximum value in image
 * @param band Search extreme value in image for this band
 **/
void ImgRasterGdal::getMinMax(double& minValue, double& maxValue, int band)
{
  bool isConstraint=(maxValue>minValue);
  double minConstraint=minValue;
  double maxConstraint=maxValue;
  std::vector<double> lineBuffer(nrOfCol());
  bool isValid=false;
  for(int irow=0;irow<nrOfRow();++irow){
    readData(lineBuffer,irow,band);
    for(int icol=0;icol<nrOfCol();++icol){
      if(isNoData(lineBuffer[icol]))
  continue;
      if(isValid){
  if(isConstraint){
    if(lineBuffer[icol]<minConstraint)
      continue;
    if(lineBuffer[icol]>maxConstraint)
      continue;
  }
  if(lineBuffer[icol]<minValue)
    minValue=lineBuffer[icol];
  if(lineBuffer[icol]>maxValue)
    maxValue=lineBuffer[icol];
      }
      else{
  if(isConstraint){
    if(lineBuffer[icol]<minConstraint)
      continue;
    if(lineBuffer[icol]>maxConstraint)
      continue;
  }
  minValue=lineBuffer[icol];
  maxValue=lineBuffer[icol];
  isValid=true;
      }
    }
  }
  if(!isValid)
    throw(static_cast<std::string>("Warning: not initialized"));
}


/**
 * @param histvector The reported histogram with counts per bin
 * @param min, max Only calculate histogram for values between min and max. If min>=max, calculate min and max from the image
 * @param nbin Number of bins used for calculating the histogram. If nbin is 0, the number of bins is  automatically calculated from min and max
 * @param theBand The band for which to calculate the histogram (start counting from 0)
 * @param kde Apply kernel density function for a Gaussian basis function
 * @return number of valid pixels in this dataset for the the selected band
 **/
double ImgRasterGdal::getHistogram(std::vector<double>& histvector, double& min, double& max, int& nbin, int theBand, bool kde){
  double minValue=0;
  double maxValue=0;

  if(min>=max)
    getMinMax(minValue,maxValue,theBand);
  else{
    minValue=min;
    maxValue=max;
  }
  if(min<max&&min>minValue)
    minValue=min;
  if(min<max&&max<maxValue)
    maxValue=max;
  min=minValue;
  max=maxValue;

  double sigma=0;
  if(kde){
    double meanValue=0;
    double stdDev=0;
    GDALProgressFunc pfnProgress;
    void* pProgressData;
    GDALRasterBand* rasterBand;
    rasterBand=getRasterBand(theBand);
    rasterBand->ComputeStatistics(0,&minValue,&maxValue,&meanValue,&stdDev,pfnProgress,pProgressData);
    //rest minvalue and MaxValue as ComputeStatistics does not account for nodata, scale and offset
    minValue=min;
    maxValue=max;

    if(m_scale.size()>theBand){
      stdDev*=m_scale[theBand];
    }
    sigma=1.06*stdDev*pow(getNvalid(theBand),-0.2);
  }

  double scale=0;
  if(maxValue>minValue){
    if(nbin==0)
      nbin=maxValue-minValue+1;
    scale=static_cast<double>(nbin-1)/(maxValue-minValue);
  }
  else
    nbin=1;
  //todo: replace assert with exception
  assert(nbin>0);
  if(histvector.size()!=nbin){
    histvector.resize(nbin);
    for(int i=0;i<nbin;histvector[i++]=0);
  }
  double nvalid=0;
  unsigned long int ninvalid=0;
  std::vector<double> lineBuffer(nrOfCol());
  for(int irow=0;irow<nrOfRow();++irow){
    readData(lineBuffer,irow,theBand);
    for(int icol=0;icol<nrOfCol();++icol){
      if(isNoData(lineBuffer[icol]))
        ++ninvalid;
      else if(lineBuffer[icol]>maxValue)
        ++ninvalid;
      else if(lineBuffer[icol]<minValue)
        ++ninvalid;
      else if(nbin==1)
        ++histvector[0];
      else{//scale to [0:nbin]
        if(sigma>0){
          //create kde for Gaussian basis function
          //todo: speed up by calculating first and last bin with non-zero contriubtion...
          //todo: calculate real surface below pdf by using gsl_cdf_gaussian_P(x-mean+binsize,sigma)-gsl_cdf_gaussian_P(x-mean,sigma)
          for(int ibin=0;ibin<nbin;++ibin){
            double icenter=minValue+static_cast<double>(maxValue-minValue)*(ibin+0.5)/nbin;
            double thePdf=gsl_ran_gaussian_pdf(lineBuffer[icol]-icenter, sigma);
            histvector[ibin]+=thePdf;
            nvalid+=thePdf;
          }
        }
        else{
          int theBin=static_cast<unsigned long int>(scale*(lineBuffer[icol]-minValue));
          //todo: replace assert with exception
          assert(theBin>=0);
          assert(theBin<nbin);
          ++histvector[theBin];
          ++nvalid;
        }
        // else if(lineBuffer[icol]==maxValue)
        //   ++histvector[nbin-1];
        // else
        //   ++histvector[static_cast<int>(static_cast<double>(lineBuffer[icol]-minValue)/(maxValue-minValue)*(nbin-1))];
      }
    }
  }
  // unsigned long int nvalid=nrOfCol()*nrOfRow()-ninvalid;
  return nvalid;
}

/**
 * @param range Sorted vector containing the range of image values
 * @param band The band for which to calculate the range
 **/
void ImgRasterGdal::getRange(std::vector<short>& range, int band)
{
  std::vector<short> lineBuffer(nrOfCol());
  range.clear();
  for(int irow=0;irow<nrOfRow();++irow){
    readData(lineBuffer,irow,band);
    for(int icol=0;icol<nrOfCol();++icol){
      if(find(range.begin(),range.end(),lineBuffer[icol])==range.end())
        range.push_back(lineBuffer[icol]);
    }
  }
  sort(range.begin(),range.end());
}

/**
 * @param band The band for which to calculate the number of valid pixels
 * @return number of valid pixels in this dataset for the the selected band
 **/
unsigned long int ImgRasterGdal::getNvalid(int band)
{
  unsigned long int nvalid=0;
  if(m_noDataValues.size()){
    std::vector<double> lineBuffer(nrOfCol());
    for(int irow=0;irow<nrOfRow();++irow){
      readData(lineBuffer,irow,band);
      for(int icol=0;icol<nrOfCol();++icol){
  if(isNoData(lineBuffer[icol]))
    continue;
  else
    ++nvalid;
      }
    }
    return nvalid;
  }
  else
    return(nrOfCol()*nrOfRow());
}

/**
 * @param band The band for which to calculate the number of valid pixels
 * @return number of invalid pixels in this dataset for the the selected band
 **/
unsigned long int ImgRasterGdal::getNinvalid(int band)
{
  unsigned long int nvalid=0;
  if(m_noDataValues.size()){
    std::vector<double> lineBuffer(nrOfCol());
    for(int irow=0;irow<nrOfRow();++irow){
      readData(lineBuffer,irow,band);
      for(int icol=0;icol<nrOfCol();++icol){
  if(isNoData(lineBuffer[icol]))
    continue;
  else
    ++nvalid;
      }
    }
    return (nrOfCol()*nrOfRow())-nvalid;
  }
  else
    return(0);
}

/**
 * @param refX, refY Calculated reference pixel position in geo-refererenced coordinates
 * @param band The band for which to calculate the number of valid pixels
 **/

void ImgRasterGdal::getRefPix(double& refX, double &refY, int band)
{
  std::vector<double> lineBuffer(nrOfCol());
  double validCol=0;
  double validRow=0;
  int nvalidCol=0;
  int nvalidRow=0;
  for(int irow=0;irow<nrOfRow();++irow){
    readData(lineBuffer,irow,band);
    for(int icol=0;icol<nrOfCol();++icol){
      // bool valid=(find(m_noDataValues.begin(),m_noDataValues.end(),lineBuffer[icol])==m_noDataValues.end());
      // if(valid){
      if(!isNoData(lineBuffer[icol])){
        validCol+=icol+1;
        ++nvalidCol;
        validRow+=irow+1;
        ++nvalidRow;
      }
    }
  }
  if(isGeoRef()){
    //reference coordinate is lower left corner of pixel in center of gravity
    //we need geo coordinates for exactly this location: validCol(Row)/nvalidCol(Row)-0.5
    double cgravi=validCol/nvalidCol-0.5;
    double cgravj=validRow/nvalidRow-0.5;
    double refpixeli=floor(cgravi);
    double refpixelj=ceil(cgravj-1);
    //but image2geo provides location at center of pixel (shifted half pixel right down)
    image2geo(refpixeli,refpixelj,refX,refY);
    //refX and refY now refer to center of gravity pixel
    refX-=0.5*getDeltaX();//shift to left corner
    refY-=0.5*getDeltaY();//shift to lower left corner
  }
  else{
    refX=floor(validCol/nvalidCol-0.5);//left corner
    refY=floor(validRow/nvalidRow-0.5);//upper corner
    //shift to lower left corner of pixel
    refY+=1;
  }
}

// /**
//  * @param filename Open a raster dataset with this filename
//  * @param imgSrc Use this source image as a template to copy image attributes
//  * @param options Creation options
//  **/
// void ImgRasterGdal::open(const std::string& filename, const ImgReaderGdal& imgSrc, const std::vector<std::string>& options)
// {
//   m_ncol=imgSrc.nrOfCol();
//   m_nrow=imgSrc.nrOfRow();
//   m_nband=imgSrc.nrOfBand();
//   m_dataType=imgSrc.getDataType();
//   setFile(filename,imgSrc,options);
//   // m_filename=filename;
//   // m_options=options;
//   // setDriver(imgSrc);
// }

/**
 * @param filename Open a raster dataset with this filename
 * @param imgSrc Use this source image as a template to copy image attributes
 * @param options Creation options
 **/
CPLErr ImgRasterGdal::open(const std::string& filename, const ImgRasterGdal& imgSrc, const std::vector<std::string>& options)
{
  m_ncol=imgSrc.nrOfCol();
  m_nrow=imgSrc.nrOfRow();
  m_nband=imgSrc.nrOfBand();
  m_dataType=imgSrc.getDataType();
  setProjection(imgSrc.getProjection());
  copyGeoTransform(imgSrc);
  if(setFile(filename,imgSrc.getImageType(),options)!=CE_None)
    return(CE_Failure);
  m_gds->SetMetadata(imgSrc.getMetadata());
  if(imgSrc.getColorTable()!=NULL)
    setColorTable(imgSrc.getColorTable());
  return(CE_None);
}

/**
 * @param imgSrc Use this source image as a template to copy image attributes
 **/
CPLErr ImgRasterGdal::open(ImgRasterGdal& imgSrc)
{
  m_ncol=imgSrc.nrOfCol();
  m_nrow=imgSrc.nrOfRow();
  m_nband=imgSrc.nrOfBand();
  m_dataType=imgSrc.getDataType();
  setProjection(imgSrc.getProjection());
  copyGeoTransform(imgSrc);
  imgSrc.getNoDataValues(m_noDataValues);
  imgSrc.getScale(m_scale);
  imgSrc.getOffset(m_offset);
  if(m_filename!=""){
    // m_writeMode=true;
    m_access=WRITE;
    registerDriver();
  }
  else
    m_access=READ_ONLY;
    // m_writeMode=false;
  //todo: check if filename needs to be set, but as is it is used for writing, I don't think so.
  // if(imgSrc.getFileName()!=""){
  //   m_filename=imgSrc.getFileName();
    // std::cerr << "Warning: filename not set, dataset not defined yet" << std::endl;
  // }
  return(CE_None);
}

/**
 * @param filename Open a raster dataset with this filename
 * @param ncol Number of columns in image
 * @param nrow Number of rows in image
 * @param nband Number of bands in image
 * @param dataType The data type of the image (one of the GDAL supported datatypes: GDT_Byte, GDT_[U]Int[16|32], GDT_Float[32|64])
 * @param imageType Image type. Currently only those formats where the drivers support the Create method can be written
 * @param options Creation options
 **/
CPLErr ImgRasterGdal::open(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, const std::vector<std::string>& options)
{
  m_ncol = ncol;
  m_nrow = nrow;
  m_nband = nband;
  m_dataType = dataType;
  return(setFile(filename,imageType,options));
}

/**
 * @param ncol Number of columns in image
 * @param nrow Number of rows in image
 * @param nband Number of bands in image
 * @param dataType The data type of the image (one of the GDAL supported datatypes: GDT_Byte, GDT_[U]Int[16|32], GDT_Float[32|64])
 **/
CPLErr ImgRasterGdal::open(int ncol, int nrow, int nband, const GDALDataType& dataType)
{
  m_ncol = ncol;
  m_nrow = nrow;
  m_nband = nband;
  m_dataType = dataType;
  if(m_filename!=""){
    // m_writeMode=true;
    m_access=WRITE;
    registerDriver();
  }
  else
    initMem();
  return(CE_None);
}

/**
 * @param imgSrc Use this source image as a template to copy image attributes
 **/
// void ImgRasterGdal::setDriver(const ImgRasterGdal& imgSrc){
//   GDALAllRegister();
//   GDALDriver *poDriver;
//   poDriver = GetGDALDriverManager()->GetDriverByName(imgSrc.getDriverDescription().c_str());
//   if( poDriver == NULL ){
//     std::string errorString="FileOpenError";
//     throw(errorString);
//   }

//   char **papszMetadata = poDriver->GetMetadata();
//   //todo: try and catch if CREATE is not supported (as in PNG)
//   if( ! CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE )){
//     std::ostringstream s;
//     s << "Error: image type " << imgSrc.getImageType() << " not supported";
//     throw(s.str());
//   }
//   char **papszOptions=NULL;
//   for(std::vector<std::string>::const_iterator optionIt=m_options.begin();optionIt!=m_options.end();++optionIt)
//     papszOptions=CSLAddString(papszOptions,optionIt->c_str());

//   m_gds=poDriver->Create(m_filename.c_str(),m_ncol,m_nrow,m_nband,imgSrc.getDataType(),papszOptions);
//   double gt[6];
//   imgSrc.getGeoTransform(gt);
//   if(setGeoTransform(gt)!=CE_None)
//     std::cerr << "Warning: could not write geotransform information in " << m_filename << std::endl;
//   setProjection(imgSrc.getProjection());
//   if(setProjection(imgSrc.getProjection())!=CE_None)
//     std::cerr << "Warning: could not write projection information in " << m_filename << std::endl;

//   if(m_noDataValues.size()){
//     for(int iband=0;iband<nrOfBand();++iband)
//       GDALSetNoDataValue(m_noDataValues[0],iband);
//   }

//   m_gds->SetMetadata(imgSrc.getMetadata());
//   m_gds->SetMetadataItem( "TIFFTAG_DOCUMENTNAME", m_filename.c_str());
//   std::string versionString="pktools ";
//   versionString+=VERSION;
//   versionString+=" by Pieter Kempeneers";
//   m_gds->SetMetadataItem( "TIFFTAG_SOFTWARE", versionString.c_str());
//   time_t rawtime;
//   time ( &rawtime );

//   time_t tim=time(NULL);
//   tm *now=localtime(&tim);
//   std::ostringstream datestream;
//   //date std::string must be 20 characters long...
//   datestream << now->tm_year+1900;
//   if(now->tm_mon+1<10)
//     datestream << ":0" << now->tm_mon+1;
//   else
//     datestream << ":" << now->tm_mon+1;
//   if(now->tm_mday<10)
//     datestream << ":0" << now->tm_mday;
//   else
//     datestream << ":" << now->tm_mday;
//   if(now->tm_hour<10)
//     datestream << " 0" << now->tm_hour;
//   else
//     datestream << " " << now->tm_hour;
//   if(now->tm_min<10)
//     datestream << ":0" << now->tm_min;
//   else
//     datestream << ":" << now->tm_min;
//   if(now->tm_sec<10)
//     datestream << ":0" << now->tm_sec;
//   else
//     datestream << ":" << now->tm_sec;
//   m_gds->SetMetadataItem( "TIFFTAG_DATETIME", datestream.str().c_str());
//   if(imgSrc.getColorTable()!=NULL)
//     setColorTable(imgSrc.getColorTable());
// }

/**
 * @param dataType The data type of the image (one of the GDAL supported datatypes: GDT_Byte, GDT_[U]Int[16|32], GDT_Float[32|64])
 * @param imageType Image type. Currently only those formats where the drivers support the Create method can be written
 **/
// void ImgRasterGdal::setDriver(const std::string& imageType){
//   GDALAllRegister();
//   GDALDriver *poDriver;
//   poDriver = GetGDALDriverManager()->GetDriverByName(imageType.c_str());
//   if( poDriver == NULL ){
//     std::ostringstream s;
//     s << "FileOpenError (" << imageType << ")";
//     throw(s.str());
//   }
//   char **papszMetadata;
//   papszMetadata = poDriver->GetMetadata();
//   //todo: try and catch if CREATE is not supported (as in PNG)
//   if( ! CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE )){
//     std::ostringstream s;
//     s << "Error: image type " << imageType << " not supported";
//     throw(s.str());
//   }
//   char **papszOptions=NULL;
//   for(std::vector<std::string>::const_iterator optionIt=m_options.begin();optionIt!=m_options.end();++optionIt)
//     papszOptions=CSLAddString(papszOptions,optionIt->c_str());

//   m_gds=poDriver->Create(m_filename.c_str(),m_ncol,m_nrow,m_nband,m_dataType,papszOptions);
//   double gt[6];
//   getGeoTransform(gt);
//   if(setGeoTransform(gt)!=CE_None)
//     std::cerr << "Warning: could not write geotransform information in " << m_filename << std::endl;
//   if(setProjection(m_projection)!=CE_None)
//     std::cerr << "Warning: could not write projection information in " << m_filename << std::endl;


//   if(m_noDataValues.size()){
//     for(int iband=0;iband<nrOfBand();++iband)
//       GDALSetNoDataValue(m_noDataValues[0],iband);
//   }

//   m_gds->SetMetadataItem( "TIFFTAG_DOCUMENTNAME", m_filename.c_str());
//   std::string versionString="pktools ";
//   versionString+=VERSION;
//   versionString+=" by Pieter Kempeneers";
//   m_gds->SetMetadataItem( "TIFFTAG_SOFTWARE", versionString.c_str());
//   time_t rawtime;
//   time ( &rawtime );

//   time_t tim=time(NULL);
//   tm *now=localtime(&tim);
//   std::ostringstream datestream;
//   //date std::string must be 20 characters long...
//   datestream << now->tm_year+1900;
//   if(now->tm_mon+1<10)
//     datestream << ":0" << now->tm_mon+1;
//   else
//     datestream << ":" << now->tm_mon+1;
//   if(now->tm_mday<10)
//     datestream << ":0" << now->tm_mday;
//   else
//     datestream << ":" << now->tm_mday;
//   if(now->tm_hour<10)
//     datestream << " 0" << now->tm_hour;
//   else
//     datestream << " " << now->tm_hour;
//   if(now->tm_min<10)
//     datestream << ":0" << now->tm_min;
//   else
//     datestream << ":" << now->tm_min;
//   if(now->tm_sec<10)
//     datestream << ":0" << now->tm_sec;
//   else
//     datestream << ":" << now->tm_sec;
//   m_gds->SetMetadataItem( "TIFFTAG_DATETIME", datestream.str().c_str());
// }

/**
 * @param filename Open a raster dataset with this filename
 * @param imageType Image type. Currently only those formats where the drivers support the Create method can be written
 **/
CPLErr ImgRasterGdal::setFile(const std::string& filename, const std::string& imageType, const std::vector<std::string>& options)
{
  m_access=WRITE;
  // m_writeMode=true;
  m_filename=filename;
  m_options=options;
  m_imageType=imageType;
  if(nrOfCol()&&nrOfRow()&&nrOfBand()){
    registerDriver();
  }
  return(CE_None);
}

CPLErr ImgRasterGdal::setFile(const app::AppFactory &app){
  Optionpk<std::string> input_opt("fn", "fn", "filename");
  Optionpk<std::string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate).","GTiff");
  Optionpk<std::string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");

  option_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  doProcess=input_opt.retrieveOption(app.getArgc(),app.getArgv());
  oformat_opt.retrieveOption(app.getArgc(),app.getArgv());
  option_opt.retrieveOption(app.getArgc(),app.getArgv());
  if(!doProcess){
    std::cout << std::endl;
    std::ostringstream helpStream;
    helpStream << "help info: ";
    throw(helpStream.str());//help was invoked, stop processing
  }
  return(setFile(input_opt[0],oformat_opt[0],option_opt));
}

/**
 * @param metadata Set this metadata when writing the image (if supported byt the driver)
 **/
void ImgRasterGdal::setMetadata(char** metadata)
{
  if(m_gds)
    m_gds->SetMetadata(metadata);
}

//default projection: ETSR-LAEA
// std::string ImgRasterGdal::setProjection(void)
// {
//   std::string theProjection;
//   OGRSpatialReference oSRS;
//   char *pszSRS_WKT = NULL;
//   //// oSRS.importFromEPSG(3035);
//   oSRS.SetGeogCS("ETRS89","European_Terrestrial_Reference_System_1989","GRS 1980",6378137,298.2572221010042,"Greenwich",0,"degree",0.0174532925199433);
//   // cout << setprecision(16) << "major axis: " << oSRS.GetSemiMajor(NULL) << endl;//notice that major axis can be set to a different value than the default to the well known standard corresponding to the name (European_Terrestrial_Reference_System_1989), but that new value, while recognized by GetSemiMajor, will not be written in the geotiff tag!
//   oSRS.SetProjCS( "ETRS89 / ETRS-LAEA" );
//   oSRS.SetLAEA(52,10,4321000,3210000);
//   oSRS.exportToWkt( &pszSRS_WKT );
//   theProjection=pszSRS_WKT;
//   CPLFree( pszSRS_WKT );
//   assert(m_gds);
//   m_gds->SetProjection(theProjection.c_str());
//   return(theProjection);
// }


/**
 * @param filename ASCII file containing 5 columns: index R G B ALFA (0:transparent, 255:solid)
 * @param band band number to set color table (starting counting from 0)
 **/
void ImgRasterGdal::setColorTable(const std::string& filename, int band)
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
    colorTable.SetColorEntry(id,&sEntry);
  }
  if(m_gds)
    (m_gds->GetRasterBand(band+1))->SetColorTable(&colorTable);
}

/**
 * @param colorTable Instance of the GDAL class GDALColorTable
 * @param band band number to set color table (starting counting from 0)
 **/
void ImgRasterGdal::setColorTable(GDALColorTable* colorTable, int band)
{
  if(m_gds)
    (m_gds->GetRasterBand(band+1))->SetColorTable(colorTable);
}

/**
 * @param ogrReader Vector dataset as an instance of the ImgReaderOgr that must be rasterized
 * @param burnValues Values to burn into raster cells (one value for each band)
 * @param controlOptions special options controlling rasterization (ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG)
 * "ATTRIBUTE":
 * Identifies an attribute field on the features to be used for a burn in value. The value will be burned into all output bands. If specified, padfLayerBurnValues will not be used and can be a NULL pointer.
 * "CHUNKYSIZE":
 * The height in lines of the chunk to operate on. The larger the chunk size the less times we need to make a pass through all the shapes. If it is not set or set to zero the default chunk size will be used. Default size will be estimated based on the GDAL cache buffer size using formula: cache_size_bytes/scanline_size_bytes, so the chunk will not exceed the cache.
 * "ALL_TOUCHED":
 * May be set to TRUE to set all pixels touched by the line or polygons, not just those whose center is within the polygon or that are selected by brezenhams line algorithm. Defaults to FALSE.
 "BURN_VALUE_
 * May be set to "Z" to use the Z values of the geometries. The value from padfLayerBurnValues or the attribute field value is added to this before burning. In default case dfBurnValue is burned as it is. This is implemented properly only for points and lines for now. Polygons will be burned using the Z value from the first point. The M value may be supported in the future.
 * "MERGE_ALG":
 * May be REPLACE (the default) or ADD. REPLACE results in overwriting of value, while ADD adds the new value to the existing raster, suitable for heatmaps for instance.
 * @param layernames Names of the vector dataset layers to process. Leave empty to process all layers
 **/
void ImgRasterGdal::rasterizeOgr(ImgReaderOgr& ogrReader, const std::vector<double>& burnValues, const std::vector<std::string>& controlOptions, const std::vector<std::string>& layernames ){
  std::vector<int> bands;
  if(burnValues.empty()&&controlOptions.empty()){
    std::string errorString="Error: either burn values or control options must be provided";
    throw(errorString);
  }
  for(int iband=0;iband<nrOfBand();++iband)
    bands.push_back(iband+1);
  std::vector<OGRLayerH> layers;
  int nlayer=0;

  std::vector<double> burnBands;//burn values for all bands in a single layer
  std::vector<double> burnLayers;//burn values for all bands and all layers
  if(burnValues.size()){
    burnBands=burnValues;
    while(burnBands.size()<nrOfBand())
      burnBands.push_back(burnValues[0]);
  }
  for(int ilayer=0;ilayer<ogrReader.getLayerCount();++ilayer){
    std::string currentLayername=ogrReader.getLayer(ilayer)->GetName();
    if(layernames.size())
      if(find(layernames.begin(),layernames.end(),currentLayername)==layernames.end())
        continue;
    std::cout << "processing layer " << currentLayername << std::endl;
    layers.push_back((OGRLayerH)ogrReader.getLayer(ilayer));
    ++nlayer;
    if(burnValues.size()){
      for(int iband=0;iband<nrOfBand();++iband)
        burnLayers.insert(burnLayers.end(),burnBands.begin(),burnBands.end());
    }
  }
  void* pTransformArg=NULL;
  GDALProgressFunc pfnProgress=NULL;
  void* pProgressArg=NULL;

  char **coptions=NULL;
  for(std::vector<std::string>::const_iterator optionIt=controlOptions.begin();optionIt!=controlOptions.end();++optionIt)
    coptions=CSLAddString(coptions,optionIt->c_str());

  if(controlOptions.size()){
    if(GDALRasterizeLayers( (GDALDatasetH)m_gds,nrOfBand(),&(bands[0]),layers.size(),&(layers[0]),NULL,pTransformArg,NULL,coptions,pfnProgress,pProgressArg)!=CE_None){
      std::string errorString(CPLGetLastErrorMsg());
      throw(errorString);
    }
  }
  else if(burnValues.size()){
    if(GDALRasterizeLayers( (GDALDatasetH)m_gds,nrOfBand(),&(bands[0]),layers.size(),&(layers[0]),NULL,pTransformArg,&(burnLayers[0]),NULL,pfnProgress,pProgressArg)!=CE_None){
      std::string errorString(CPLGetLastErrorMsg());
      throw(errorString);
    }
  }
  else{
    std::string errorString="Error: either attribute fieldname or burn values must be set to rasterize vector dataset";
    throw(errorString);
  }
  m_access=READ_ONLY;
  // m_writeMode=false;
}

/**
 * @param ogrReader Vector dataset as an instance of the ImgReaderOgr that must be rasterized
 * @param burnValues Values to burn into raster cells (one value for each band)
 * @param layernames Names of the vector dataset layers to process. Leave empty to process all layers
 **/
void ImgRasterGdal::rasterizeBuf(ImgReaderOgr& ogrReader, double burnValue, const std::vector<std::string>& layernames ){
  std::vector<OGRLayerH> layers;
  int nlayer=0;

  std::vector<double> burnBands;//burn values for all bands in a single layer
  while(burnBands.size()<nrOfBand())
    burnBands.push_back(burnValue);
  for(int ilayer=0;ilayer<ogrReader.getLayerCount();++ilayer){
    std::string currentLayername=ogrReader.getLayer(ilayer)->GetName();
    if(layernames.size())
      if(find(layernames.begin(),layernames.end(),currentLayername)==layernames.end())
        continue;
    std::cout << "processing layer " << currentLayername << std::endl;
    layers.push_back((OGRLayerH)ogrReader.getLayer(ilayer));
    ++nlayer;
  }
  void* pTransformArg=NULL;
  GDALProgressFunc pfnProgress=NULL;
  void* pProgressArg=NULL;

  if(m_data.size()!=nrOfBand()){
    std::string errorString="Error: m_data not initialized";
    throw(errorString);
  }
  for(int iband=0;iband<nrOfBand();++iband){
    if(!(m_data[iband])){
      std::string errorString="Error: m_data not initialized";
      throw(errorString);
    }
    Vector2d<double> initBlock(nrOfRow(),nrOfCol());
    writeDataBlock(initBlock,0,nrOfCol()-1,0,nrOfRow()-1,iband);
    double gt[6];
    getGeoTransform(gt);
    if(GDALRasterizeLayersBuf(m_data[iband],nrOfCol(),nrOfRow(),getDataType(),GDALGetDataTypeSize(getDataType())>>3,0,layers.size(),&(layers[0]), getProjectionRef().c_str(),gt,NULL, pTransformArg, burnBands[iband],NULL,pfnProgress,pProgressArg)!=CE_None){
      std::string errorString(CPLGetLastErrorMsg());
      throw(errorString);
    }
  }
}

/**
 * @param ogrReader Vector dataset as an instance of the ImgReaderOgr that must be rasterized
 * @param controlOptions special options controlling rasterization (ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG)
 * "ATTRIBUTE":
 * Identifies an attribute field on the features to be used for a burn in value. The value will be burned into all output bands. If specified, padfLayerBurnValues will not be used and can be a NULL pointer.
 * "ALL_TOUCHED":
 * May be set to TRUE to set all pixels touched by the line or polygons, not just those whose center is within the polygon or that are selected by brezenhams line algorithm. Defaults to FALSE.
 "BURN_VALUE_FROM":
 * May be set to "Z" to use the Z values of the geometries. The value from padfLayerBurnValues or the attribute field value is added to this before burning. In default case dfBurnValue is burned as it is. This is implemented properly only for points and lines for now. Polygons will be burned using the Z value from the first point. The M value may be supported in the future.
 * "MERGE_ALG":
 * May be REPLACE (the default) or ADD. REPLACE results in overwriting of value, while ADD adds the new value to the existing raster, suitable for heatmaps for instance.
 * @param layernames Names of the vector dataset layers to process. Leave empty to process all layers
 **/
void ImgRasterGdal::rasterizeBuf(ImgReaderOgr& ogrReader, const std::vector<std::string>& controlOptions, const std::vector<std::string>& layernames ){
  std::vector<OGRLayerH> layers;
  int nlayer=0;

  for(int ilayer=0;ilayer<ogrReader.getLayerCount();++ilayer){
    std::string currentLayername=ogrReader.getLayer(ilayer)->GetName();
    if(layernames.size())
      if(find(layernames.begin(),layernames.end(),currentLayername)==layernames.end())
        continue;
    std::cout << "processing layer " << currentLayername << std::endl;
    layers.push_back((OGRLayerH)ogrReader.getLayer(ilayer));
    ++nlayer;
  }
  void* pTransformArg=NULL;
  GDALProgressFunc pfnProgress=NULL;
  void* pProgressArg=NULL;

  char **coptions=NULL;
  for(std::vector<std::string>::const_iterator optionIt=controlOptions.begin();optionIt!=controlOptions.end();++optionIt)
    coptions=CSLAddString(coptions,optionIt->c_str());

  if(m_data.size()!=nrOfBand()){
    std::string errorString="Error: m_data not initialized";
    throw(errorString);
  }
  for(int iband=0;iband<nrOfBand();++iband){
    if(!(m_data[iband])){
      std::string errorString="Error: m_data not initialized";
      throw(errorString);
    }
    Vector2d<double> initBlock(nrOfRow(),nrOfCol());
    writeDataBlock(initBlock,0,nrOfCol()-1,0,nrOfRow()-1,iband);
    double gt[6];
    getGeoTransform(gt);
    if(GDALRasterizeLayersBuf(m_data[iband],nrOfCol(),nrOfRow(),getDataType(),GDALGetDataTypeSize(getDataType())>>3,0,layers.size(),&(layers[0]), getProjectionRef().c_str(),gt,NULL, pTransformArg, 0,coptions,pfnProgress,pProgressArg)!=CE_None){
      std::string errorString(CPLGetLastErrorMsg());
      throw(errorString);
    }
  }
}

/**
 *
 *
 * @param t1 minimum threshold
 * @param t2 maximum threshold
 * @param bg value if outside thresholds
 *
 * @return CE_None if success, CE_Failure if failed
 */CPLErr ImgRasterGdal::setThreshold(double t1, double t2){
  try{
    if(m_noDataValues.empty()){
      std::string errorString="Error: no data value not set";
      throw(errorString);
    }
    std::vector<double> lineInput(nrOfCol());
    for(int iband=0;iband<nrOfBand();++iband){
      for(int irow=0;irow<nrOfRow();++irow){
        readData(lineInput,irow,iband);
        for(int icol=0;icol<nrOfCol();++icol){
          if(lineInput[icol]>=t1&&lineInput[icol]<=t2)
            continue;
          else
            lineInput[icol]=m_noDataValues[0];
        }
        writeData(lineInput,irow,iband);
      }
    }
  }
  catch(std::string errorstring){
    std::cerr << errorstring << std::endl;
    return(CE_Failure);
  }
  catch(...){
    return(CE_Failure);
  }
  return(CE_None);
}

/**
 *
 *
 * @param t1 minimum threshold
 * @param t2 maximum threshold
 * @param fg value if within thresholds
 * @param bg value if outside thresholds
 *
 * @return CE_None if success, CE_Failure if failed
 */CPLErr ImgRasterGdal::setThreshold(double t1, double t2, double value){
  try{
    if(m_noDataValues.empty()){
      std::string errorString="Error: no data value not set";
      throw(errorString);
    }
    std::vector<double> lineInput(nrOfCol());
    for(int iband=0;iband<nrOfBand();++iband){
      for(int irow=0;irow<nrOfRow();++irow){
        readData(lineInput,irow,iband);
        for(int icol=0;icol<nrOfCol();++icol){
          if(lineInput[icol]>=t1&&lineInput[icol]<=t2)
            lineInput[icol]=value;
          else
            lineInput[icol]=m_noDataValues[0];
        }
        writeData(lineInput,irow,iband);
      }
    }
  }
  catch(std::string errorstring){
    std::cerr << errorstring << std::endl;
    return(CE_Failure);
  }
  catch(...){
    return(CE_Failure);
  }
  return(CE_None);
}
