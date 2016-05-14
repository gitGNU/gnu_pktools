/**********************************************************************
ImgRasterGdal.cc: class to read raster files using GDAL API library
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
#include "ogr_spatialref.h"
#include "ImgRasterGdal.h"

ImgRasterGdal::ImgRasterGdal(void)
  : m_gds(NULL), m_ncol(0), m_nrow(0), m_nband(0), m_dataType(GDT_Unknown)
{}

ImgRasterGdal::~ImgRasterGdal(void)
{
}

void ImgRasterGdal::close(void)
{
  GDALClose(m_gds);
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
  if(m_gds)
    return(m_gds->SetProjection(projection.c_str()));
  else
    return(CE_Failure);
}

/**
 * @param band get data type for this band (start counting from 0)
 * @return the GDAL data type of this data set for the selected band
 **/
GDALDataType ImgRasterGdal::getDataType(int band) const
{
  assert(band<m_nband+1);
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
  assert(band<m_nband+1);
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
  assert(band<m_nband+1);
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
  if(m_gds)
    return(m_gds->SetGeoTransform(m_gt));
  else
    return(CE_Failure);
      
}

/**
 * @param imgSrc Use this source image as a template to copy geotranform information
 **/
void ImgRasterGdal::copyGeoTransform(const ImgRasterGdal& imgSrc)
{
  setProjection(imgSrc.getProjection());
  double gt[6];
  imgSrc.getGeoTransform(gt);
  setGeoTransform(gt);
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
 * @return the metadata of this data set in string format
 **/
char** ImgRasterGdal::getMetadata()
{
  if(m_gds){
    if(m_gds->GetMetadata()!=NULL)
      return(m_gds->GetMetadata());
  }
  else
    return (char**)"";
}

/**
 * @return the metadata of this data set in C style string format (const version)
 **/
char** ImgRasterGdal::getMetadata() const
{
  if(m_gds){
    if(m_gds->GetMetadata()!=NULL)
      return(m_gds->GetMetadata());
  }
  else 
    return (char**)"";
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
