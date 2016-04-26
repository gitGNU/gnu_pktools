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
  : m_gds(NULL), m_ncol(0), m_nrow(0), m_nband(0)
{}

void ImgRasterGdal::close(void)
{
  GDALClose(m_gds);
}

std::string ImgRasterGdal::getProjection(void) const 
{
  std::string theProjection=m_gds->GetProjectionRef();
  // size_t startpos,endpos;
  // while((startpos=theProjection.find(",AUTHORITY"))!=std::string::npos){
  //   endpos=theProjection.find("]",startpos+1,1)+1;
  //   theProjection.erase(startpos,endpos-startpos);
  // }
  return theProjection;
}

std::string ImgRasterGdal::getProjectionRef(void) const 
{
  std::string theProjection;
  if(m_gds->GetProjectionRef())
    return(m_gds->GetProjectionRef());
  else
    return "";
}

void ImgRasterGdal::setGeoTransform(double* gt){
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

std::string ImgRasterGdal::setProjectionProj4(const std::string& projection)
{
  OGRSpatialReference theRef;
  theRef.SetFromUserInput(projection.c_str());
  char *wktString;
  theRef.exportToWkt(&wktString);
  assert(m_gds);
  m_gds->SetProjection(wktString);
  return(wktString);
}

void ImgRasterGdal::setProjection(const std::string& projection)
{
  OGRSpatialReference oSRS;
  char *pszSRS_WKT = NULL;
  assert(m_gds);
  m_gds->SetProjection(projection.c_str());
  CPLFree(pszSRS_WKT);
}


GDALDataType ImgRasterGdal::getDataType(int band) const
{
  assert(band<m_nband+1);
  return (m_gds->GetRasterBand(band+1))->GetRasterDataType();
}

GDALRasterBand* ImgRasterGdal::getRasterBand(int band)
{
  assert(band<m_nband+1);
  return (m_gds->GetRasterBand(band+1));
}

GDALColorTable* ImgRasterGdal::getColorTable(int band) const
{
  assert(band<m_nband+1);
  return (m_gds->GetRasterBand(band+1))->GetColorTable();
}

std::string ImgRasterGdal::getDriverDescription() const
{
  return m_gds->GetDriver()->GetDescription();
}

void ImgRasterGdal::getGeoTransform(double* gt) const{
  m_gds->GetGeoTransform(gt);
}

// void ImgRasterGdal::getGeoTransform(double& ulx, double& uly, double& deltaX, double& deltaY, double& rot1, double& rot2) const
// {
//   double adfGeoTransform[6];// { 444720, 30, 0, 3751320, 0, -30 };
//   m_gds->GetGeoTransform(adfGeoTransform);
//   ulx=adfGeoTransform[0];
//   deltaX=adfGeoTransform[1];
//   rot1=adfGeoTransform[2];
//   uly=adfGeoTransform[3];
//   rot2=adfGeoTransform[4];
//   deltaY=-adfGeoTransform[5];//convention of GDAL!
// }

std::string ImgRasterGdal::getGeoTransform() const
{
  double gt[6];// { 444720, 30, 0, 3751320, 0, -30 };
  m_gds->GetGeoTransform(gt);
  std::ostringstream s;
  s << "[" << gt[0] << "," << gt[1] << "," << gt[2] << "," << gt[3] << "," << gt[4] << "," << gt[5] << "]";
  return(s.str());
  // if(!isGeoRef())
  //   return("");
  // else{
  //   double adfGeoTransform[6];// { 444720, 30, 0, 3751320, 0, -30 };
  //   m_gds->GetGeoTransform(adfGeoTransform);
  //   double ulx=adfGeoTransform[0];
  //   double deltaX=adfGeoTransform[1];
  //   double rot1=adfGeoTransform[2];
  //   double uly=adfGeoTransform[3];
  //   double rot2=adfGeoTransform[4];
  //   double deltaY=-adfGeoTransform[5];//convention of GDAL!
  //   std::ostringstream s;
  //   s << "[" << ulx << "," << deltaX << "," << rot1 << "," << uly << "," << rot2 << "," << -deltaY << "]";
  //   return(s.str());
  // }
}

char** ImgRasterGdal::getMetadata()
{
  if(m_gds->GetMetadata()!=NULL)
    return(m_gds->GetMetadata());
  else
    return (char**)"";
}

char** ImgRasterGdal::getMetadata() const
{
  if(m_gds->GetMetadata()!=NULL)
    return(m_gds->GetMetadata());
  else 
    return (char**)"";
}

void ImgRasterGdal::getMetadata(std::list<std::string>& metadata) const
{
  char** cmetadata=m_gds->GetMetadata();
  while(*cmetadata!=NULL){
    metadata.push_back(*(cmetadata));
    ++cmetadata;
  }
}

std::string ImgRasterGdal::getDescription() const
{
  if(m_gds->GetDriver()->GetDescription()!=NULL)
    return m_gds->GetDriver()->GetDescription();
  else
    return "";
}

std::string ImgRasterGdal::getMetadataItem() const 
{
  if(m_gds->GetDriver()->GetMetadataItem( GDAL_DMD_LONGNAME )!=NULL)
    return m_gds->GetDriver()->GetMetadataItem( GDAL_DMD_LONGNAME );
  else
    return "";
}
std::string ImgRasterGdal::getImageDescription() const 
{
  if(m_gds->GetDriver()->GetMetadataItem("TIFFTAG_IMAGEDESCRIPTION")!=NULL)
    return m_gds->GetDriver()->GetMetadataItem("TIFFTAG_IMAGEDESCRIPTION");
  else
    return "";
}

std::string ImgRasterGdal::getInterleave() const
{
  if(m_gds->GetMetadataItem( "INTERLEAVE", "IMAGE_STRUCTURE"))
    return m_gds->GetMetadataItem( "INTERLEAVE", "IMAGE_STRUCTURE");
  else
    return("BAND");
}

std::string ImgRasterGdal::getCompression() const
{
  if(m_gds->GetMetadataItem( "COMPRESSION", "IMAGE_STRUCTURE"))
    return m_gds->GetMetadataItem( "COMPRESSION", "IMAGE_STRUCTURE");
  else
    return("NONE");
}

bool ImgRasterGdal::getBoundingBox(double& ulx, double& uly, double& lrx, double& lry) const
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

bool ImgRasterGdal::getCenterPos(double& x, double& y) const
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
    // x=m_ulx+(nrOfCol()/2.0)*m_delta_x;
    // y=m_uly-(nrOfRow()/2.0)*m_delta_y;
    return true;
  }
  else{
    // x=nrOfCol()/2.0;
    // y=nrOfRow()/2.0;
    return false;
  }
}

//i and j represent fraction of pixels, return true if image is georeferenced
bool ImgRasterGdal::geo2image(double x, double y, double& i, double& j) const
{
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

//x and y represent center of pixel, return true if image is georeferenced
bool ImgRasterGdal::image2geo(double i, double j, double& x, double& y) const
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
  else{
    // x=0.5+i;
    // y=nrOfRow()-(0.5+j);
    return false;
  }
}

bool ImgRasterGdal::covers(double x, double  y) const
{
  double theULX, theULY, theLRX, theLRY;
  getBoundingBox(theULX,theULY,theLRX,theLRY);
  return((x > theULX)&&
         (x < theLRX)&&
         (y < theULY)&&
         (y >theLRY));
}

bool ImgRasterGdal::covers(double ulx, double  uly, double lrx, double lry) const
{
  double theULX, theULY, theLRX, theLRY;
  getBoundingBox(theULX,theULY,theLRX,theLRY);
  return((ulx < theLRX)&&(lrx > theULX)&&(lry < theULY)&&(uly > theLRY));
}

int ImgRasterGdal::getNoDataValues(std::vector<double>& noDataValues) const
{
  if(m_noDataValues.size()){
    noDataValues=m_noDataValues;
    return m_noDataValues.size();
  }
  else
    return 0;
}

int ImgRasterGdal::pushNoDataValue(double noDataValue)
{
  if(find(m_noDataValues.begin(),m_noDataValues.end(),noDataValue)==m_noDataValues.end())
    m_noDataValues.push_back(noDataValue);
  return(m_noDataValues.size());
}

// bool ImgRasterGdal::setNoDataValue(double noDataValue,int band)
// {
//   GDALRasterBand  *poBand;
//   poBand = m_gds->GetRasterBand(band+1);
//   if(poBand->SetNoDataValue(noDataValue)!=CE_None)
//     return false;
//   else
//     return true;
// }
